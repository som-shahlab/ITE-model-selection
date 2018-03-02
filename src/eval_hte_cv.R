library(magrittr)
library(tidyverse)

# group data into all test folds
# for each test fold, do the internal matching
# record the matched patient
find_matches = function(data) {
	treated_match = Match(Tr=data$treatment, 
	   					  X=data.matrix(data%>%select(starts_with("covariate"))),
	    				  replace=T, estimand="ATT") 
	control_match = Match(Tr=!data$treatment, # (separate so one treated doesn't get two controls ambiguously)
	    				  X=data.matrix(data%>%select(starts_with("covariate"))),
	    				  replace=T, estimand="ATT")
	data.frame(subject = data$subject[c(treated_match$index.treated, control_match$index.treated)], # all subjects
			   match = data$subject[c(treated_match$index.control, control_match$index.control)]) # their matches
}


setup_data = function(scenario, training_percent, n_folds) {
	setup = setup_simulation(scenario) 
	setup$params$p = 10 # for efficiency

	data = create_data(setup) %>% 
	    mutate(set = ifelse(subject > n()*training_percent, "test", "training"))
	true_effect = data %>% 
	    select(starts_with("covariate")) %>% 
	    (setup$functions$effect)
	true_mean = data %>% 
	    select(starts_with("covariate")) %>% 
	    (setup$functions$mean)
	true_propensity = data %>% 
	    select(starts_with("covariate")) %>% 
	    (setup$functions$propensity)
	# opt_value = true_hte_value(true_effect, true_effect, true_mean)

	cv_index = data %>% 
	    filter(set=="training") %>%
	    create_cv_index(n_folds=n_folds)
	test_index = list("training" = (data %>% filter(set=="training") %>% pull(subject)))

	cv_held_out = cv_index %>% 
    	map(~test_index$training[!(test_index$training %in% .)])
	held_out = c(cv_held_out, list("training" = (data %>% filter(set=="test") %>% pull(subject))))
	matches = held_out %>% 
	    map(~ data.frame(subject=.) %>%
	    	inner_join(data, by="subject") %>% 
	    	find_matches()) %>%
	    bind_rows(.id="fold")

	aux_data = data.frame(subject=data$subject, true_mean=true_mean, true_effect=true_effect, true_propensity=true_propensity) %>%
		inner_join(matches, by="subject") %>%
		inner_join(data %>% select(subject, treatment), by='subject') %>%
		mutate(ip_weights = 1/(1-treatment + 2*treatment*true_propensity - true_propensity)) %>% # evaluates to either 1/p1 (if w=1) or 1/p0 (if w=0)
		select(-treatment)
	return(list(data=data, aux_data=aux_data, cv_index=cv_index, test_index=test_index))
}

get_estimates = function(data, model, cv_index, test_index, aux_data) {
	cv_estimates = data %>% 
	    filter(set=="training") %>%
	    cross_estimate_hte(model$method, model$tune_grid, cv_index) %>% 
	    inner_join(aux_data, by=c("subject", "fold"))
	test_estimates = data %>% 
	    cross_estimate_hte(model$method, model$tune_grid, test_index) %>% 
	    inner_join(aux_data, by=c("subject", "fold"))
	return(list(cv_estimates=cv_estimates, test_estimates=test_estimates))
}

compute_cv_metrics = function(estimates) {
	estimates  %>%
	    # dplyr::group_by(!!!syms(c(param_names, "fold"))) %>% # I do this for each fold
	    dplyr::group_by(model, fold) %>% # I do this for each fold
	    dplyr::summarize(prediction_accuracy = -mse(est_outcome, outcome),
	                     transformed_accuracy = -transformed_mse(est_effect, treatment, outcome, weights=ip_weights),
	                     matching_accuracy = -matching_mse(est_effect, treatment, outcome, subject, match),
	                     # uplift = uplift(est_effect, treatment, outcome),
	                     # decile = decile(est_effect, outcome, treatment),
	                     # c_benefit = c_benefit(est_effect,treatment,outcome),
	                     est_gain = gain(est_effect, treatment, outcome, weights=ip_weights),
	                     est_value = value(est_effect, treatment, outcome, weights=ip_weights),
	                     # value_max = value_max(est_effect, treatment, outcome, weights=ip_weights),
	                     # c_benefit_k = c_benefit_k(est_effect, treatment, outcome, weights=ip_weights)
	                     ) %>%
	    # dplyr::ungroup() %>% dplyr::group_by(!!!syms(param_names)) %>% # Then average over the folds
	    dplyr::ungroup() %>% dplyr::group_by(model) %>% # Then average over the folds
	    dplyr::summarize(prediction_accuracy = mean(prediction_accuracy, na.rm=T),
	                     est_value = mean(est_value, na.rm=T),
	                     transformed_accuracy = mean(transformed_accuracy, na.rm=T),
	                     matching_accuracy = mean(matching_accuracy, na.rm=T),
	                     # uplift = mean(uplift, na.rm=T),
	                     est_gain = mean(est_gain, na.rm=T),
	                     # decile = mean(decile),
	                     # c_benefit = mean(c_benefit, na.rm=T),
	                     # value_max = mean(value_max, na.rm=T),
	                     # c_benefit_k = mean(c_benefit_k, na.rm=T)
	                     )
}

compute_test_metrics = function(estimates) {
	estimates  %>%
	    dplyr::group_by(model) %>% # there should only be one fold
	    dplyr::summarize(true_hte_error = mse(est_effect, true_effect),
	                     true_value = true_hte_value(est_effect, true_effect, true_mean)) 
}

get_errors = function(cv_estimates, test_estimates) {
	cv_accuracy = cv_estimates %>% 
    	compute_cv_metrics() %>%
        gather(selection_method, accuracy, -model)
	min_cv_accuracy = cv_accuracy %>%
	    group_by(selection_method) %>%
	    filter(accuracy == max(accuracy, na.rm=T)) %>%
	    sample_n(1) %>% # if there are ties for the highest accuracy, break at random
	    select(-accuracy)
	test_accuracy = test_estimates %>% 
	    compute_test_metrics() 
	true_selection_accuracy = min_cv_accuracy %>%
	    inner_join(test_accuracy, by="model") 
	    # mutate(optimal_deficiency = -true_hte_value(true_effect, true_effect, true_mean)) # %>%
	    # mutate(scenario=scenario, n_folds=n_folds, training_percent=training_percent, rep=rep)
	return(list(cv_accuracy=cv_accuracy, test_accuracy=test_accuracy, true_selection_accuracy=true_selection_accuracy))
}
