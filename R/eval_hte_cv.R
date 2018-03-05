#' @import dplyr
#' @import purrr
#' @import tidyr
#' @import magrittr
#' @import caret
#' @import Matching

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

#' Prepares simulated data for experiments
#'
#' @param DGP a list created by a call to dgp()
#' @param n the number of samples desired
#' @param training_percent how much data to use for training 
#' @param n_folds the number of folds to use for treatment effect cross-validation
#' @keywords
#' @export
#' @examples
setup_data = function(DGP, n, training_percent, n_folds) {
	simulation = create_data(DGP, n) 
	data = simulation$data %>% 
	    mutate(set = ifelse(subject > n*training_percent, "test", "training"))
	aux_data = simulation$aux_data
	# opt_value = true_value(true_effect, true_effect, true_mean)

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

	aux_data = aux_data %>%
		inner_join(matches, by="subject") %>%
		inner_join(data %>% select(subject, treatment), by='subject') %>%
		mutate(ip_weights = 1/(1-treatment + 2*treatment*true_propensity - true_propensity)) %>% # evaluates to either 1/p1 (if w=1) or 1/p0 (if w=0)
		select(-treatment)
	return(list(data=data, aux_data=aux_data, cv_index=cv_index, test_index=test_index))
}

#' Estimates conditonal mean regression models on the data via caret. 
#' Does cross-estimation on the training data and runs each model trained on the full training set on the test set.
#' Returns all out-of-sample estimates from each model
#'
#' @param data  data 
#' @param models  models 
#' @param cv_index cv_index
#' @param test_index  test_index 
#' @keywords
#' @export
#' @examples
get_estimates = function(data, models, cv_index, test_index) {
	training_data = data %>% 
	    filter(set=="training") 
	cv_estimates = models %>% # cross validate on CV data (ignore the test set)
		map(~cross_estimate_hte(training_data, .$method, .$tune_grid, cv_index)) %>%
	    bind_rows() 
	test_estimates = models %>% # now train on all CV data and test on the test set
		map(~cross_estimate_hte(data, .$method, .$tune_grid, test_index)) %>% 
		bind_rows() 
	return(list(cv_estimates=cv_estimates, test_estimates=test_estimates))
}

compute_cv_metrics = function(estimates) {
	estimates  %>%
	    # dplyr::group_by(!!!syms(c(param_names, "fold"))) %>% # I do this for each fold
	    dplyr::group_by(model, fold) %>% # I do this for each fold
	    # mutate(est_effect_test_match = est_effect_covariate_matching(treatment, outcome, subject, match),
	    	   # est_effect_test_trans = est_effect_transformed_outcome(treatment, outcome, ip_weights)) %>%
	    dplyr::summarize(#### framework: ###
	    				 match_mse = est_effect_covariate_matching(treatment, outcome, subject, match) %>% loss_squared_error(est_effect),
	    				 trans_mse = est_effect_transformed_outcome(treatment, outcome, ip_weights) %>% loss_squared_error(est_effect), 
	    				 match_decision = est_effect_covariate_matching(treatment, outcome, subject, match) %>% loss_decision(est_effect),
	    				 trans_decision = est_effect_transformed_outcome(treatment, outcome, ip_weights) %>% loss_decision(est_effect), # aka gain!
	    				 # value:
						 value = -value(est_effect, treatment, outcome, weights=ip_weights),
						 # gain = -gain(est_effect, treatment, outcome, weights=ip_weights),
						 #### broken: ####
	    				 prediction_error = loss_squared_error(est_outcome, outcome)
	                     # uplift = uplift(est_effect, treatment, outcome),
	                     # decile = decile(est_effect, outcome, treatment),
	                     # c_benefit = c_benefit(est_effect,treatment,outcome),
	                     # value_max = value_max(est_effect, treatment, outcome, weights=ip_weights),
	                     # c_benefit_k = c_benefit_k(est_effect, treatment, outcome, weights=ip_weights)
	                     ) %>%
	    # dplyr::ungroup() %>% dplyr::group_by(!!!syms(param_names)) %>% # Then average over the folds
	    dplyr::ungroup() %>% dplyr::select(-fold) %>% dplyr::group_by(model) %>% # Then average over the folds
	    dplyr::summarize_all(mean, na.rm=T)
}

compute_test_metrics = function(estimates) {
	estimates  %>%
	    dplyr::group_by(model) %>% # there should only be one fold
	    dplyr::summarize(true_hte_error = loss_squared_error(est_effect, true_effect),
	                     true_value = -true_value(est_effect, true_effect, true_mean)) 
}

#' Estimates estimated cv errors and true test errors (the latter via true values in aux_data). 
#'
#' @param cv_estimates  cv_estimates 
#' @param test_estimates test_estimates 
#' @param aux_data aux_data
#' @keywords
#' @export
#' @examples
get_errors = function(cv_estimates, test_estimates, aux_data) {
	cv_error = cv_estimates %>% 
		inner_join(aux_data, by=c("subject", "fold")) %>%
    	compute_cv_metrics() %>%
        gather(selection_method, error, -model)
	min_cv_error = cv_error %>%
	    group_by(selection_method) %>%
	    filter(error == min(error, na.rm=T)) %>%
	    sample_n(1) %>% # if there are ties for the lowest error, break at random
	    select(-error)
	test_error = test_estimates %>% 
		inner_join(aux_data, by=c("subject", "fold")) %>%
	    compute_test_metrics() 
	true_selection_error = min_cv_error %>%
	    inner_join(test_error, by="model") 
	    # mutate(optimal_deficiency = -true_hte_value(true_effect, true_effect, true_mean)) # %>%
	    # mutate(scenario=scenario, n_folds=n_folds, training_percent=training_percent, rep=rep)
	return(list(cv_error=cv_error, test_error=test_error, true_selection_error=true_selection_error))
}