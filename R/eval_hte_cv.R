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

create_cv_index = function(data, n_folds=5) {
		createFolds(data$treatment, k=n_folds) %>% # hold-out indices
    	map(~(data$subject)[-.]) #  complement of each index, in terms of the subject IDs
}


# see:
# Fewell, Z., Hernán, M. A., Wolfe, F., Tilling, K., Stata, H. C., 2004. (n.d.). Controlling for time-dependent confounding using marginal structural models. Pdfs.Semanticscholar.org
# Rodriguez, G. (2014). Survival Models. In Lecture Notes on Generalized Linear Models (pp. 1–34).
ipcw = function(data, n_time_intervals=10) {
	interval_data = data %>% 
	    mutate(event = !event) %>% # "event" is now censoring, not death
	    mutate(event_interval = ntile(time, n_time_intervals))
	pseudo_data = list(subject = interval_data %$% unique(subject), 
                   interval = interval_data %$% unique(event_interval)) %>%
    cross_df %>%
    inner_join(interval_data, by="subject") %>%
    filter(interval <= event_interval) %>% # keep only pseudo-observations before or at event time
    mutate(event = ifelse(interval==event_interval, event, FALSE)) # before the real time of censoring, set censoring=F
    pCt_XW = glm.fit( # no global intercept
	    x = pseudo_data %>% 
	            select(subject, interval, treatment, starts_with("covariate")) %>% 
	            mutate(interval_copy = interval) %>%
	            mutate(trash=1) %>%
	            spread(interval, trash, fill=0) %>% # give each interval its own intercept
	            select(-interval_copy, -subject) %>%
	            data.matrix(),
	    y = pseudo_data %>% pull(event),
	    family = binomial(link="cloglog")) %$% 
		fitted.values
    data.frame(pCt_XW = pCt_XW,
          	   subject = pseudo_data$subject) %>%
	    group_by(subject) %>%
	    summarize(cens_prob_at_event_time = prod(1-pCt_XW)) %>%
	    ungroup() %>%
	    mutate(est_ipcw = 1/cens_prob_at_event_time) %>%
	    inner_join(data, by="subject") %>%
	    select(subject, est_ipcw)
}

iptw = function(data) {
	pW_X = glm.fit(
		x = data %>% select(starts_with("covariate")) %>% as.matrix,
		y = data %>% pull("treatment"),
 		family = binomial(link="logit")) %$% 
	fitted.values
	data.frame(pW_X = pW_X, 
			   treatment=data$treatment,
			   subject=data$subject) %>%
		mutate(est_iptw = 1/(1-treatment + 2*treatment*pW_X - pW_X)) %>%
		select(subject, est_iptw)
}

#' Prepares simulated data for experiments
#'
#' @param DGP a list created by a call to dgp()
#' @param n_train number of samples to train and validate on
#' @param n_test number of samples to test on
#' @param n_folds the number of folds to use for treatment effect cross-validation
#' @keywords
#' @export
#' @examples
setup_data = function(DGP, n_train, n_test, n_folds) {
	simulation = create_data(DGP, n_train+n_test) 
	data = simulation$data %>% 
	    mutate(set = ifelse(subject > n_train, "test", "training"))
	aux_data = simulation$aux_data %>%
		mutate(set = ifelse(subject > n_train, "test", "training"))

	cv_index = data %>% 
	    filter(set=="training") %>%
	    create_cv_index(n_folds=n_folds)
	test_index = list("training" = (data %>% filter(set=="training") %>% pull(subject)))

	aux_data = aux_data %>%
		inner_join(data %>% select(subject, treatment), by='subject') %>%
		inner_join(iptw(data), by="subject") %>% 	# I'm cheating a little bit here because I use test data to fit these
		inner_join(ipcw(data), by="subject") %>% 	# but the models should be underfit anyways so that data shouldn't add much...
		select(-treatment) # can fix this later but won't really change it
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
		map(~cross_estimate_hte(training_data, .$method, .$method_name, .$tune_grid, cv_index)) %>%
	    bind_rows() 
	test_estimates = models %>% # now train on all CV data and test on the test set
		map(~cross_estimate_hte(data, .$method, .$method_name, .$tune_grid, test_index)) %>% 
		bind_rows() 
	return(list(cv_estimates=cv_estimates, test_estimates=test_estimates))
}

compute_cv_metrics = function(estimates) {
	estimates  %>%
	    # dplyr::group_by(!!!syms(c(param_names, "fold"))) %>% # I do this for each fold
	    dplyr::group_by(model, fold) %>% # I do this for each fold
	    # mutate(est_effect_test_match = est_effect_covariate_matching(treatment, outcome, subject, match),
	    	   # est_effect_test_trans = est_effect_transformed_outcome(treatment, outcome, ip_weights)) %>%
	    dplyr::summarize(value = value(est_effect, treatment, time, event, IPTW=est_iptw),
	    				 c_stat = c_stat(est_pseudo_outcome, time, event, treatment, IPCW=est_ipcw),
	                     random = random_metric()
	                     ) %>%
	   	dplyr::ungroup() %>% dplyr::select(-fold) %>% dplyr::group_by(model) %>% # Then average over the folds
	    dplyr::summarize_all(mean, na.rm=T)
}

compute_test_metrics = function(estimates) {
	estimates  %>%
	    dplyr::group_by(model) %>% # there should only be one fold
	    dplyr::summarize(true_value = true_value(est_effect, treated_mean, control_mean))
}

#' Estimates estimated cv metrics and true test metrics (the latter via true values in aux_data). 
#'
#' @param cv_estimates  cv_estimates 
#' @param test_estimates test_estimates 
#' @param aux_data aux_data
#' @keywords
#' @export
#' @examples
get_metrics = function(cv_estimates, test_estimates, aux_data) {
	cv_metrics = cv_estimates %>% 
		inner_join(aux_data, by="subject") %>%
    	compute_cv_metrics() %>%
        gather(selection_method, metrics, -model)
	best_cv_metrics = cv_metrics %>%
	    group_by(selection_method) %>%
	    filter(metrics == max(metrics, na.rm=T)) %>%
	    sample_n(1) %>% # if there are ties for the lowest metrics, break at random
	    select(-metrics) %>% ungroup()
	test_metrics = test_estimates %>% 
		inner_join(aux_data, by="subject") %>%
	    compute_test_metrics()
	oracle_metrics = test_metrics %>%
		gather(selection_method, metrics, -model) %>%
		group_by(selection_method) %>%
	    filter(metrics == max(metrics, na.rm=T)) %>%
	    sample_n(1) %>%
		select(-metrics) %>% ungroup() %>%
		mutate(selection_method = str_c("oracle_selector", selection_method, sep="_"))
	true_selection_metrics = best_cv_metrics %>%
		bind_rows(oracle_metrics) %>%
	    inner_join(test_metrics, by="model") %>%
	    bind_rows(data.frame(model="truth", selection_method="oracle",  # this is the true model
	    					 true_value=(aux_data %>% filter(set=="test") %$% true_value(-effect, treated_mean, control_mean)))) %>% 
	    bind_rows(data.frame(model="harm", selection_method="demon",  # this is the true "evil" model
	    					 true_value=(aux_data %>% filter(set=="test") %$% true_value(effect, treated_mean, control_mean)))) 
	    # mutate(optimal_deficiency = -true_hte_value(true_effect, true_effect, true_mean)) # %>%
	    # mutate(scenario=scenario, n_folds=n_folds, training_percent=training_percent, rep=rep)
	return(list(cv_metrics=cv_metrics, test_metrics=test_metrics, true_selection_metrics=true_selection_metrics))
}
