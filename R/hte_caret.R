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

# https://topepo.github.io/caret/model-training-and-tuning.html#alternate-performance-metrics
summary_metrics = function(data, lev=NULL, model=NULL) {
	if (is.null(data$weights)) {
		data %<>% mutate(weights=1)
	}
	data %>% select(obs, pred, weights) %->% c(obs, pred, weights) # makes sure the order is correct
	if (!is.factor(obs) && is.numeric(obs)) {
		c(wRMSE = sqrt(sum(weights*(obs-pred)^2)/sum(weights)))
	} else {
		pred = factor(pred, levels = levels(obs)) # not sure if needed, but this line is in caret::postResample
		c(wAccuracy = sum(weights*(obs == pred))/sum(weights))
	}
}

# takes a list of fit caret models (hyperparams already optimized by caret), returns the one with lowest wRMSE or 
# highest wAccuracy
pick_model = function(models) {
	if(length(models)==1) {
		return(models[[1]])
	} else {
		best_model_name = resamples(models)$values %>% # each row is a fold, columns are (model x metric)
		    gather(model_metric, value, -Resample) %>% 
		    separate(model_metric, c("model","metric"), sep="~") %>%
		    # filter(metric==metric) %>% # there should only be a single metric: wRMSE or wAccuracy
		    group_by(model) %>%
		    summarize(mean_value = mean(as.numeric(value), na.rm=T)) %>% # as.numeric in case of weird things because of NAs
		    filter(mean_value==ifelse(models[[1]]$maximize, max(mean_value), min(mean_value))) %>%
		    pull(model) %>% first() # in case of ties
		return(models[[best_model_name]])
	}
}

# this function only returns the test-set predictions from the model determined best via internal k-fold cross-validation
# to be used to cross-estimate mu and p on the true validation set
learner_cv = function(x, y, model_specs, weights=NULL, k_folds=4) {
	if(is.logical(y)) {y = factor(ifelse(y, "treated", "control"))}
	model_specs %>% imap(function(settings, method) {
		train_args = list(
			x = x, y = y, weights = weights, 
			metric = "wRMSE", maximize=F, # these will be changed if it is a classification problem
			method = method, tuneGrid = settings$tune_grid,
			trControl = trainControl(
				summaryFunction=summary_metrics,
				method='cv', number=k_folds,
			  	returnResamp="final", savePredictions="final"))
		if(is.factor(y)) {
			train_args$trControl$classProbs=T
			train_args$metric="wAccuracy"
			train_args$maximize=T}
		set.seed(1)
		do.call(train, c(train_args, settings$extra_args))
	}) %>% pick_model(
	)
}

learners_pred_test = function(training_index, x, y, model_specs, weights=NULL) {
	if(is.logical(y)) {y = factor(ifelse(y, "treated", "control"))}
	model_specs %>% imap(function(settings, method) {
		train_args = list(
			x = x, y = y, weights = weights, 
			metric = "wRMSE", maximize=F, # these will be changed if it is a classification problem
			method = method, tuneGrid = settings$tune_grid,
			trControl = trainControl(
				summaryFunction=summary_metrics,
				method='cv', number=1, index=list(index=training_index),
			  	returnResamp="none", savePredictions="all"))
		if(is.factor(y)) {
			train_args$trControl$classProbs=T
			train_args$metric="wAccuracy"
			train_args$maximize=T}
		set.seed(1)
		trained = do.call(train, c(train_args, settings$extra_args)) 

		trained$pred %>%
		unite(model, names(settings$tune_grid), sep="~") %>%
		mutate(model = str_c(method, model, sep="@")) %>%
		select(yhat=pred, index=rowIndex, model)
	}) %>% bind_rows()
}

# model specs should be a list, each element of which is a list(method=list(tune_grid, extra_args))
cross_validated_cross_estimation = function(x, y, model_specs, weights=NULL, k_folds_cv=4, k_folds_ce=4) {
	set.seed(1)
	test_indices = createFolds(y, k=k_folds_ce) 
	test_indices %>% map(function(test_index) {
		model = learner_cv(x[-test_index,], y[-test_index], model_specs, weights=weights, k_folds=k_folds_cv) 
		if(is.logical(y)) {
			predict(model, newdata=x[test_index,], type="prob") %>%
			data.frame(cross_estimate = .$treated, index=test_index)
		} else {
			predict(model, newdata=x[test_index,]) %>%
			data.frame(cross_estimate = ., index=test_index)
		}
	}) %>% bind_rows %>% arrange(index) %>% pull(cross_estimate)
}

# returns a fit R-learner model object. mu_hat and p_hat are cross-validatedly cross-estimated, then fixed.
# The hyperparameters of the tau_hat model are selected with cross-validation and is refit to the full data
R_learner_cv = function(x, w, y, mu_model_specs, p_model_specs, tau_model_specs, k_folds_cv=4, k_folds_ce=4) {
	mu_hat = cross_validated_cross_estimation(
		x, y, mu_model_specs, 
		k_folds_cv=k_folds_cv, 
		k_folds_ce=k_folds_ce)
	p_hat = cross_validated_cross_estimation(
		x, w, p_model_specs, 
		k_folds_cv=k_folds_cv, 
		k_folds_ce=k_folds_ce)

	pseudo_outcome = (y - mu_hat)/(w - p_hat)
	weights = (w - p_hat)^2

	learner_cv(x, pseudo_outcome, tau_model_specs, weights=weights, k_folds=k_folds_cv)
}

# returns predictions from multiple R-learners. mu_hat and p_hat are cross-validatedly cross-estimated, then fixed.
# Test-set predictions from each model of tau_hat (derived from each set of hyperparameter values) are returned.
R_learners_pred_test = function(training_index, x, w, y, mu_model_specs, p_model_specs, tau_model_specs, k_folds_cv=4, k_folds_ce=4) {
	mu_hat = cross_validated_cross_estimation(
		x[training_index,], y[training_index], 
		mu_model_specs, 
		k_folds_cv=k_folds_cv, 
		k_folds_ce=k_folds_ce)
	p_hat = cross_validated_cross_estimation(
		x[training_index,], w[training_index], 
		p_model_specs, 
		k_folds_cv=k_folds_cv, 
		k_folds_ce=k_folds_ce)

	c(pseudo_outcome, weights) %<-% rep(list(rep(1,length(y))),2) # filler values, no impact on result (but can't be 0 to avoid caret error)
	pseudo_outcome[training_index] = (y[training_index] - mu_hat)/(w[training_index] - p_hat)
	weights[training_index] = (w[training_index] - p_hat)^2

	learners_pred_test(training_index, x, pseudo_outcome, tau_model_specs, weights=weights)
}

# returns predictions from multiple S-learners. 
# Test-set predictions from each model of tau_hat (derived from each set of hyperparameter values) are returned.
S_learners_pred_test = function(training_index, x, w, y, model_specs) {
	N = length(y)
	all_index = 1:N
	N_val = length(all_index[-training_index])

	x_cf = rbind(x[training_index,], x[-training_index,], x[-training_index,])
		list(y, all_index) %>% 
	    map(~c(.[training_index], .[-training_index], .[-training_index])) %->%
	c(y_cf, index_cf)
	w_cf = c(w[training_index], rep(0, N_val), rep(1, N_val))

	learners_pred_test(1:length(training_index), cbind(x_cf,(w_cf-0.5)*x_cf), y_cf, model_specs) %>%
		mutate(couterfactual = ifelse(index <= N, "control", "treated")) %>%
		mutate(index = index_cf[index]) %>%
		spread(couterfactual, yhat) %>%
		mutate(est_effect=treated-control, est_outcome=treated*w[index]+control*!w[index])
}

# returns predictions from multiple T-learners. 
# The two models in each pair are constrained to have the same hyperparameters for computational feasibility
# Test-set predictions from each model of tau_hat (derived from each set of hyperparameter values) are returned.
T_learners_pred_test	= function(training_index, x, w, y, model_specs) {
	list(treated=TRUE, control=FALSE) %>% imap(function(condition, counterfactual) {
		learners_pred_test(
			training_index[w[training_index]==condition], # filters the training index to only include subjects treated under "condition"
			x, y, model_specs) %>%
			mutate(counterfactual = counterfactual)
	}) %>%
	bind_rows() %>%
	filter(!(index %in% training_index)) %>% # each model will have predicted for subjects in the training set of the other model
	spread(counterfactual, yhat) %>%
	mutate(est_effect=treated-control, est_outcome=treated*w[index]+control*!w[index]) 
}
