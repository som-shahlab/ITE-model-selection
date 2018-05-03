#' @import dplyr
#' @import purrr
#' @import tidyr
#' @import magrittr
#' @import caret

# # returns a function that uses training data to return the expected expected values of y (real-valued or binary)
# # given x in a test set (under different hyperparameter values)
# # this can be used to generate different ML methods on-the-fly by passing in the method name in caret
# # will use this to estimate the final treatment effect function for the R-learner
# make_learner_caret = function(caret_method_name) {
# 	function(x, y, tune_grid, training_index, weights=NULL) {
# 		if(is.logical(y)) {y = factor(ifelse(y, "treated", "control"))}
# 		all_model_predictions = train(x = x, y = y, weights = weights, 
# 									  method = caret_method_name, tuneGrid = tune_grid,
# 			  						  trControl = trainControl(method='cv', number=1, # will feed in 1 fold at a time
# 									  						   index=list(index=training_index),
# 									                      	   returnResamp="all", savePredictions="all", classProbs=T))$pred
# 		if(is.factor(y)){
# 			all_model_predictions %>% select(prediction=treated, index=rowIndex, names(tune_grid))
# 		} else {
# 			all_model_predictions %>% select(prediction=pred, index=rowIndex, names(tune_grid))
# 		}
# 		# returns df with names (estimate, row, ... ) where ... are the names of the tune grid
# 	}
# }

# takes a list of fit caret models (hyperparams already optimized by caret), returns the one with lowest RMSE
pick_model = function(models, metric, opt=min) {
	if(length(models)==1) {
		return(models[[1]])
	} else {
		best_model_name = resamples(models)$values %>% # each row is a fold, columns are (model x metric)
		    gather(model_metric, value, -Resample) %>% 
		    separate(model_metric, c("model","metric"), sep="~") %>%
		    filter(metric==metric) %>%
		    group_by(model) %>%
		    summarize(mean_value = mean(as.numeric(value), na.rm=T)) %>% # as.numeric in case of weird things because of NAs
		    filter(mean_value==opt(mean_value)) %>%
		    pull(model) %>% first() # in case of ties
		return(models[[best_model_name]])
	}
}

# this function only returns the test-set predictions from the model determined best via internal k-fold cross-validation
# to be used to cross-estimate mu and p on the true validation set
cv_learner = function(x, y, model_specs, weights=NULL, k_folds=5) {
	if(is.logical(y)) {y = factor(ifelse(y, "treated", "control"))}
	best_models = model_specs %>% imap(function(tune_grid, method) {
		set.seed(1)
		train(x = x, y = y, weights = weights, 
		  	  method = method, tuneGrid = tune_grid,
			  trControl = trainControl(method='cv', number=k_folds,
                      	   returnResamp="final", savePredictions="final", classProbs=T,
                      	   verbose=F))
	})
	if(is.factor(y)){
		best_model = pick_model(best_models, "Accuracy", max)
	} else {
		best_model = pick_model(best_models, "RMSE", min)
	}
}

# model specs should be a list, each element of which is a list(method, tune_grid)
cross_validated_cross_estimation = function(x, y, model_specs, weights=NULL, k_folds_cv=3, k_folds_ce=5) {
	set.seed(1)
	test_indices = createFolds(y, k=k_folds_ce) 
	test_indices %>% map(function(test_index) {
		model = cv_learner(x[-test_index,], y[-test_index], model_specs, weights=weights, k_folds=k_folds_cv) 
		if(is.logical(y)) {
			predict(model, newdata=x[test_index,], type="prob") %>%
			data.frame(cross_estimate = .$treated, index=test_index)
		} else {
			predict(model, newdata=x[test_index,]) %>%
			data.frame(cross_estimate = ., index=test_index)
		}
	}) %>% bind_rows %>% arrange(index) %>% pull(cross_estimate)
}

# https://topepo.github.io/caret/model-training-and-tuning.html#alternate-performance-metrics
wRMSE = function(data, lev=NULL, model=NULL) {
	c(wRMSE = data%$%sqrt(sum(weights*(obs-pred)^2)/sum(weights)))
} # not sure why using this makes it weird

# returns a fit R-learner model object. mu_hat and p_hat are cross-validatedly cross-estimated, then fixed.
# The hyperparameters of the tau_hat model are selected with cross-validation and is refit to the full data
R_learner_cv = function(x, w, y, mu_model_specs, p_model_specs, tau_model_specs, k_folds_cv=3, k_folds_ce=5) {
	mu_hat = cross_validated_cross_estimation(x, y, mu_model_specs, 
											  k_folds_cv=k_folds_cv, k_folds_ce=k_folds_ce)
	p_hat = cross_validated_cross_estimation(x, w, p_model_specs, 
											 k_folds_cv=k_folds_cv, k_folds_ce=k_folds_ce)

	pseudo_outcome = (y - mu_hat)/(w - p_hat)
	weights = (w - p_hat)^2
	
	tau_model_specs %>% imap(function(tune_grid, method) {
		set.seed(1)
		train(x = x, y = pseudo_outcome, weights = weights, 
			  method = method, tuneGrid = tune_grid, metric="wRMSE", maximize=F,
			  trControl = trainControl(method='cv', number=k_folds_cv,
		                      	   returnResamp="final", savePredictions="final",
		                      	   summaryFunction=wRMSE
		                      	   ))
	}) %>% pick_model("wRMSE", min) 
}

# returns predictions from multiple R-learners. mu_hat and p_hat are cross-validatedly cross-estimated, then fixed.
# Test-set predictions from each model of tau_hat (derived from each set of hyperparameter values) are returned.
R_learners_pred_test = function(training_index, x, w, y, mu_model_specs, p_model_specs, tau_model_specs, k_folds_cv=3, k_folds_ce=5, default_weight=0) {
	mu_hat = cross_validated_cross_estimation(x[training_index,], y[training_index], mu_model_specs, 
											  k_folds_cv=k_folds_cv, k_folds_ce=k_folds_ce)
	p_hat = cross_validated_cross_estimation(x[training_index,], w[training_index], p_model_specs, 
											 k_folds_cv=k_folds_cv, k_folds_ce=k_folds_ce)
	pseudo_outcome = rep(default_weight,length(y)) # the default weight changes nothing. It is mutable for testing purposes.
	pseudo_outcome[training_index] = (y[training_index] - mu_hat)/(w[training_index] - p_hat)
	weights = rep(default_weight, length(w))
	weights[training_index] = (w[training_index] - p_hat)^2

	tau_model_specs %>% imap(function(tune_grid, method) {
		train(x = x, y = pseudo_outcome, weights = weights, 
			  method = method, tuneGrid = tune_grid,
			  trControl = trainControl(method='cv', number=1, index=list(index=training_index), # use the training data and predict on the test/validation data
		                      	   returnResamp="none", savePredictions="all")) %$%
		pred %>%
		unite(model, names(tune_grid), sep="~") %>%
		mutate(model = str_c(method, model, sep="@")) %>%
		select(est_effect=pred, index=rowIndex, model)
	}) %>% bind_rows()
}

couterfactual_test_obs = function(x, training_index) {
	if (is.vector(x)) {c(x[training_index], x[-training_index], x[-training_index])}
	else {rbind(x[training_index,], x[-training_index,], x[-training_index,])}
}

# returns predictions from multiple S-learners. 
# Test-set predictions from each model of tau_hat (derived from each set of hyperparameter values) are returned.
S_learners_pred_test = function(training_index, x, w, y, model_specs) {
	N = length(y)
	index = 1:N
	N_val = length(index[-training_index])
	list(x, y, index) %>% 
	    map(~couterfactual_test_obs(., training_index)) %->%
	    c(x_cf, y_cf, index_cf)
	w_cf = c(w[training_index], rep(0, N_val), rep(1, N_val))
	model_specs %>% imap(function(tune_grid, method) {
		all_model_predictions = train(x = cbind(x_cf,(w_cf-0.5)*x_cf), y = y_cf,  
									  method = method, tuneGrid = tune_grid,
			  						  trControl = trainControl(method='cv', number=1, index=list(index=1:length(training_index)),
									                      	   returnResamp="none", savePredictions="all"))$pred
		all_model_predictions %>% 
			mutate(couterfactual = ifelse(rowIndex <= N, "control", "treated")) %>%
			mutate(index = index_cf[rowIndex]) %>%
			select(-rowIndex, -obs) %>%
			spread(couterfactual, pred) %>%
			mutate(est_effect=treated-control, est_outcome=treated*w[index]+control*!w[index]) %>%
			unite(model, names(tune_grid), sep="~") %>%
			mutate(model = str_c(method, model, sep="@")) %>%
			select(est_effect, est_outcome, index, model)
	}) %>% bind_rows()
}

filter_treatment_index = function(w, training_index, condition) {
	list(index=training_index[w[training_index]==condition])
}

# returns predictions from multiple T-learners. 
# The two models in each pair are constrained to have the same hyperparameters for computational feasibility
# Test-set predictions from each model of tau_hat (derived from each set of hyperparameter values) are returned.
T_learners_pred_test	= function(training_index, x, w, y, model_specs) {
	model_specs %>% imap(function(tune_grid, method) {
		list(treated=TRUE, control=FALSE) %>% imap(function(condition, counterfactual) {
			all_model_predictions = train(x = x, y = y,  
										  method = method, tuneGrid = tune_grid,
				  						  trControl = trainControl(method='cv', number=1, 
				  						  						   index=filter_treatment_index(w, training_index, condition),
										                      	   returnResamp="none", savePredictions="all"))$pred
			all_model_predictions %>% 
				mutate(counterfactual = counterfactual)
		}) %>% 
		bind_rows %>%
		filter(!(rowIndex %in% training_index)) %>% # each model will have predicted for subjects in the training set of the other model
		select(-obs) %>%
		spread(counterfactual, pred) %>%
		mutate(est_effect=treated-control, est_outcome=treated*w[rowIndex]+control*!w[rowIndex]) %>%
		unite(model, names(tune_grid), sep="~") %>%
		mutate(model = str_c(method, model, sep="@")) %>%
		select(est_effect, est_outcome, index=rowIndex, model)
	}) %>% bind_rows()
}
