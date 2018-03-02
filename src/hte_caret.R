library(caret)
library(magrittr)
library(tidyverse)
library(rlang)

create_cv_index = function(data, n_folds=5) {
		createFolds(data$treatment, k=n_folds) %>% # hold-out indices
    	map(~(data$subject)[-.]) #  complement of each index, in terms of the subject IDs
}

# fit_model = setup_fitting("gbm", expand.grid(...))
# model = fit_model(data)
setup_fitting_1_model = function(method, tune_grid=NULL, train_index) {
	function(data) train(
		  x = data %>% dplyr::select(treatment, starts_with("covariate")) %>% as.matrix,
		  y = data$outcome,
		  method = method,
		  trControl = trainControl(method='cv', 
                                 number=length(train_index),
		  						 index=train_index,
                                 returnResamp="all",
                                 savePredictions="all"),
		  tuneGrid = tune_grid)
}

test_estimate_hte_1_model = function(data, method, tune_grid, train_index) {
	fitting_function = setup_fitting_1_model(method, tune_grid, train_index)
	cf_data = data %>% # must be in this order for the train index to work
    	mutate(treatment = !treatment)
	full_data = bind_rows(data, cf_data)
	models = full_data %>% fitting_function
	models$pred %>% # this carries with it columns with all the values of the hyperparameters
		mutate(method = method) %>%
		# unite(model, -pred, -obs, -rowIndex, -Resample, sep="^") %>%
		unite_("model", c("method", names(tune_grid)), sep="~") %>%
		# unite(model, !!!syms(c("method", names(tune_grid))), sep="~") %>%
	    mutate(cf = ifelse(rowIndex > nrow(data), "counterfactual", "factual")) %>%
	    mutate(subject = ifelse(cf=="factual", rowIndex, rowIndex-nrow(data))) %>%
	    select(-rowIndex, -obs) %>%
	    spread(cf, pred) %>% 
	    arrange(Resample, subject) %>%
	    filter(!is.na(factual)) %>%
	    inner_join(data %>% select(subject, treatment, outcome), by="subject") %>%
	    mutate(treated=ifelse(treatment, factual, counterfactual),
	           control=ifelse(treatment, counterfactual, factual),
	           effect=treated-control)
}

# I'm going to have to hand-code the cross-validation here... i.e. split into train/test: 
# for each train, fit a bunch models on the treated and control people separately
# for each test, use (each combination?) of two models on all that data together to generate two columns: 
# pred treated outcome and pred control outcome. 
# Also a 3rd column pred outcome that is w*treated_outcome + (1-w)*control_outcome (for outcome CV)

fit_model = function(data, train_index, method, tune_grid=NULL) {
	train(
		  x = data %>% dplyr::select(starts_with("covariate")) %>% as.matrix,
		  y = data$outcome,
		  method = method,
		  trControl = trainControl(method='cv', # will feed in 1 fold at a time
                                 number=length(train_index),
		  						 index=train_index,
                                 returnResamp="all",
                                 savePredictions="all"),
		  tuneGrid = tune_grid)
}

prep_fold_data = function(training_data, test_data) {
	data = bind_rows(training_data, test_data) %>%
		mutate(rowIndex=row_number())
	index = list(c(data %>% 
					filter(sample_type=="training") %>% 
					pull(rowIndex)))
	return(list(data=data, index=index))
}

test_estimate_hte = function(data, method, tune_grid, fold, fold_name) {
	training_data = data[fold,] %>% mutate(sample_type="training")
	test_data = data[-fold,] %>% mutate(sample_type="test")

	fold_data = training_data %>% 
		split(.$treatment) %>%
		map(~prep_fold_data(., test_data)) #now have a list (treat => (data, fold))
	predictions = fold_data %>%
		map(~fit_model(.$data, .$index, method, tune_grid)$pred) # returns the big matrix with all test set predictions for each treatment
	test_estimates = fold_data %>%
	    map(~select(.$data, subject, treatment, outcome, rowIndex)) %>%
	    list(predictions) %>%
	    pmap(function(data, predictions) inner_join(data, predictions, by="rowIndex")) %>%
	    imap(function(df,name) rename_(df, .dots=setNames("pred", str_c("est_outcome", name, sep="_")))) %>%
	    reduce(inner_join, by=c("subject", "treatment", "outcome", names(tune_grid))) %>% # will join on all columns... if didn't want to join on params would also have to join across methods!
	    mutate(method=method) %>%
	    unite_("model", c("method", names(tune_grid)), sep="~") %>% 
	    mutate(fold=fold_name) %>%
	    mutate(est_effect=est_outcome_TRUE-est_outcome_FALSE, 
	    	   est_outcome=treatment*(est_outcome_TRUE) + (1-treatment)*est_outcome_FALSE) %>%
	    select(subject, model, treatment, outcome, est_effect, est_outcome, fold) 
}

cross_estimate_hte = function(data, method, tune_grid, train_index) {
	train_index %>%
	imap(function(fold, fold_name) test_estimate_hte(data, method, tune_grid, fold, fold_name)) %>%
	bind_rows()
}

# methods = list("xgbTree", "lm")
# tune_grids = list(expand.grid(nrounds = 1:150, max_depth = 2, eta = 0.1), NULL)
