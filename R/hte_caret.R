#' @import dplyr
#' @import purrr
#' @import tidyr
#' @import magrittr
#' @import caret
#' @import Matching

# setup_fitting_1_model = function(method, tune_grid=NULL, train_index) {
# 	function(data) train(
# 		  x = data %>% dplyr::select(treatment, starts_with("covariate")) %>% as.matrix,
# 		  y = data$outcome,
# 		  method = method,
# 		  trControl = trainControl(method='cv', 
#                                  number=length(train_index),
# 		  						 index=train_index,
#                                  returnResamp="all",
#                                  savePredictions="all"),
# 		  tuneGrid = tune_grid)
# }

# test_estimate_hte_1_model = function(data, method, tune_grid, train_index) {
# 	fitting_function = setup_fitting_1_model(method, tune_grid, train_index)
# 	cf_data = data %>% # must be in this order for the train index to work
#     	mutate(treatment = !treatment)
# 	full_data = bind_rows(data, cf_data)
# 	models = full_data %>% fitting_function
# 	models$pred %>% # this carries with it columns with all the values of the hyperparameters
# 		mutate(method = method) %>%
# 		# unite(model, -pred, -obs, -rowIndex, -Resample, sep="^") %>%
# 		unite_("model", c("method", names(tune_grid)), sep="~") %>%
# 		# unite(model, !!!syms(c("method", names(tune_grid))), sep="~") %>%
# 	    mutate(cf = ifelse(rowIndex > nrow(data), "counterfactual", "factual")) %>%
# 	    mutate(subject = ifelse(cf=="factual", rowIndex, rowIndex-nrow(data))) %>%
# 	    select(-rowIndex, -obs) %>%
# 	    spread(cf, pred) %>% 
# 	    arrange(Resample, subject) %>%
# 	    filter(!is.na(factual)) %>%
# 	    inner_join(data %>% select(subject, treatment, outcome), by="subject") %>%
# 	    mutate(treated=ifelse(treatment, factual, counterfactual),
# 	           control=ifelse(treatment, counterfactual, factual),
# 	           effect=treated-control)
# }


gbm_ph_fit_predict = function(x_train, y_train, x_val, subject, interaction.depth, n.minobsinnode, shrinkage, n.trees) {
    model = gbm.fit(x_train, y_train, distribution="coxph", verbose=F,
                    n.trees=max(n.trees), interaction.depth=interaction.depth, 
                    shrinkage=shrinkage, n.minobsinnode=n.minobsinnode)
    predict(model, x_val, n.trees=n.trees) %>% 
        data.frame %>%
        mutate(subject=subject, interaction.depth=interaction.depth, 
               shrinkage=shrinkage, n.minobsinnode=n.minobsinnode) %>%
        gather(n.trees, pred, -subject, -interaction.depth, -shrinkage, -n.minobsinnode) %>%
        mutate(n.trees = str_replace(n.trees,"X","") %>% as.numeric)
}

# the thing that is estimated by coxph models is the log-relative (to the basline) risk 
fit_model = function(data, train_index, method, tune_grid) {
    x_train = data[train_index,] %>% dplyr::select(starts_with("covariate")) %>% as.matrix
    y_train = data[train_index,] %$% Surv(time, event)
    x_val = data[-train_index,] %>% dplyr::select(starts_with("covariate")) %>% as.matrix
    subject = data[-train_index,] %>% pull(subject)
    
    if(method=="gbm") {
        grouped_tune_grid = tune_grid %>% 
            group_by(interaction.depth, n.minobsinnode, shrinkage) %>%
            summarize(n.trees=list(n.trees))
        grouped_tune_grid %>%
            pmap(function(interaction.depth, n.minobsinnode, shrinkage, n.trees) {
                gbm_ph_fit_predict(x_train, y_train, x_val, subject,
                                   interaction.depth, n.minobsinnode, 
                                   shrinkage, n.trees)}) %>%
            bind_rows()
    }
}

prep_fold_data = function(training_data, validation_data) {
	data = bind_rows(training_data, validation_data) %>%
		mutate(rowIndex=row_number())
	index = data %>% 
		filter(sample_type=="training") %>% 
		pull(rowIndex)
	return(list(data=data, index=index))
}

# The two-model approach might not make sense for proportional hazards: the baseline hazard is different in both models,
# what is estimated is the difference over that baseline. 
test_estimate_hte = function(data, method, tune_grid, fold, fold_name) {
	training_data = data %>% 
		filter(subject %in% fold) %>% 
		mutate(sample_type="training")
	validation_data = data %>% 
		filter(!(subject %in% fold)) %>%
		mutate(sample_type="test")

	fold_data = training_data %>% 
		split(.$treatment) %>%
		map(~prep_fold_data(., validation_data)) #now have a list (treat => (data, fold))

	predictions = fold_data %>% # fit one model to each treatment group
		map(~fit_model(.$data, .$index, method, tune_grid)) # returns the big matrix with all test set predictions for each treatment
	test_estimates = fold_data %>%
	    map(~select(.$data, subject, treatment, time, event)) %>%
	    list(predictions) %>% 
	    pmap(function(data, predictions) inner_join(data, predictions, by="subject")) %>%
	    imap(function(df,name) rename_(df, .dots=setNames("pred", str_c("est_rel_risk", name, sep="_")))) %>% # df_TRUE$pred and df_FALSE$pred become (est_rel_risk_TRUE, ..._FALSE) in the same df
	    reduce(inner_join, by=c("subject", "treatment", "time", "event", names(tune_grid))) %>% # will join on all columns... if didn't want to join on params would also have to join across methods!
	    mutate(method=method) %>%
	    unite_("model", c("method", names(tune_grid)), sep="~") %>% 
	    mutate(fold=fold_name) %>%
	    mutate(est_effect=est_rel_risk_TRUE-est_rel_risk_FALSE, # this is the log-relative (to the control group) risk 
	    	   est_rel_risk=treatment*(est_rel_risk_TRUE) + (1-treatment)*est_rel_risk_FALSE) %>% # selects the estimated relative risk from the appropriate model
	    select(subject, model, treatment, time, event, est_effect, est_rel_risk, fold) 
}

cross_estimate_hte = function(data, method, tune_grid, train_index) {
	train_index %>%
	imap(function(fold, fold_name) test_estimate_hte(data, method, tune_grid, fold, fold_name)) %>%
	bind_rows()
}


# methods = list("xgbTree", "lm")
# tune_grids = list(expand.grid(nrounds = 1:150, max_depth = 2, eta = 0.1), NULL)
