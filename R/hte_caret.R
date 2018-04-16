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
test_estimate_hte = function(data, method, method_name, tune_grid, fold, fold_name) {
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
		map(~fit_model(.$data, .$index, method, tune_grid)) %>% # returns the big matrix with all test set predictions for each treatment
	    imap(function(df,name) rename_(df, .dots=setNames("pred", str_c("est_pseudo_outcome", name, sep="_")))) %>% # df_TRUE$pred and df_FALSE$pred become (est_pseudo_outcome_TRUE, ..._FALSE) in the same df
	    reduce(inner_join, by=c("subject", names(tune_grid))) %>% # will join on all columns... if didn't want to join on params would also have to join across methods in some far outer loop!
	    mutate(method=method_name) %>%
	    unite_("model", c("method", names(tune_grid)), sep="~") %>% 
	    mutate(fold=fold_name)  # selects the estimated relative risk from the appropriate model
}

# rewrite a function like this for all necessary methods
gbm_ph = function(data, tune_grid, fold, fold_name) {
	training_data = data %>% 
		filter(subject %in% fold)
	validation_data = data %>% 
		filter(!(subject %in% fold)) 

	x_train = training_data %>% dplyr::select(treatment, starts_with("covariate")) %>% as.matrix
    y_train = training_data %$% Surv(time, event)
    x_val = validation_data %>% dplyr::select(treatment, starts_with("covariate"))
    subject = validation_data %>% pull(subject)

    x_val_cf = list("TRUE"=TRUE, "FALSE"=FALSE) %>% map(~mutate(x_val, treatment=.) %>% as.matrix) # makes list(x_val_0, x_val_1) where treatments have been replaced 

    tune_grid %>% 
        group_by(interaction.depth, n.minobsinnode, shrinkage) %>%
        summarize(n.trees=list(n.trees)) %>%
        pmap(function(interaction.depth, n.minobsinnode, shrinkage, n.trees) {
		    model = gbm.fit(x_train, y_train, distribution="coxph", verbose=F,
		                    n.trees=max(n.trees), interaction.depth=interaction.depth, 
		                    shrinkage=shrinkage, n.minobsinnode=n.minobsinnode)
		    x_val_cf %>% 
		    	map(~ predict(model, ., n.trees=n.trees) %>% 
		    		data.frame %>%
		    		mutate(subject=subject, interaction.depth=interaction.depth, 
	                   shrinkage=shrinkage, n.minobsinnode=n.minobsinnode) %>%
		            gather(n.trees, pred, -subject, -interaction.depth, -shrinkage, -n.minobsinnode) %>%
		            mutate(n.trees = str_replace(n.trees,"X","") %>% as.numeric)) %>% 
		    	imap(function(df,name) rename_(df, .dots=setNames("pred", str_c("relative_risk_to_basline", name, sep="_")))) %>%
		    	reduce(inner_join, by=c("subject", names(tune_grid)))
        }) %>%
    	bind_rows() %>% 
	    mutate(method="gbm_ph", fold=fold_name) %>%
	    unite_("model", c("method", names(tune_grid)), sep="~") %>%
	    inner_join(validation_data, by='subject') %>%
	    group_by(model) %>%
		mutate(est_treatment = relative_risk_to_basline_FALSE > relative_risk_to_basline_TRUE, # false > true so that positive favors treatment since this is log relative risk
	    	   est_pseudo_outcome=treatment*(relative_risk_to_basline_TRUE) + (1-treatment)*relative_risk_to_basline_FALSE,
	    	   est_ranking=row_number(-est_pseudo_outcome)) %>% #higher ranking means longer survival 
		select(subject, model, treatment, time, event, est_treatment, est_ranking, fold) 
}

# things to put in the tune_grid:
# num.trees, mtry, min.node.size,
# use treatment as always.split.variables? 

extract_mean_survival = function(ranger_pred) {
	dt = data.frame(t = ranger_pred$unique.death.times) %>% 
	    mutate(dt=lead(t)-t, dt=ifelse(is.na(dt), 0, dt)) %>%
	    pull(dt)
	data.frame(pred=(ranger_pred$survival %*% dt))
}

one_model_surv_rf = function(data, tune_grid, fold, fold_name) {
	training_data = data %>% 
		filter(subject %in% fold) %>%
		select(time, event, treatment, starts_with("covariate"))
	validation_data = data %>% 
		filter(!(subject %in% fold))
	validation_data_cf = list("TRUE"=TRUE, "FALSE"=FALSE) %>% 
		map(~mutate(validation_data, treatment=.))
	subject = validation_data %>% pull(subject)

	tune_grid %>%
		pmap(function(num.trees, mtry, min.node.size) {
			model = ranger(dependent.variable.name="time", status.variable.name="event", 
						   data=training_data, 
			   			   mtry=mtry, num.trees=num.trees, min.node.size=min.node.size)
			validation_data_cf %>% 
		    	map(~ predict(model, ., n.trees=n.trees) %>% 
		    		extract_mean_survival() %>%
		    		mutate(subject=subject, mtry=mtry, num.trees=num.trees, min.node.size=min.node.size)) %>%
		    	imap(function(df,name) rename_(df, .dots=setNames("pred", str_c("mean_survival", name, sep="_")))) %>%
		    	reduce(inner_join, by=c("subject", names(tune_grid))) 
		}) %>%
    	bind_rows() %>% 
	    mutate(method="one_model_surv_rf", fold=fold_name) %>%
	    unite_("model", c("method", names(tune_grid)), sep="~") %>%
	    inner_join(validation_data, by='subject') %>%
	    group_by(model) %>%
		mutate(est_effect=mean_survival_TRUE - mean_survival_FALSE, 
	    	   est_survival=treatment*(mean_survival_TRUE) + (1-treatment)*mean_survival_FALSE,
	    	   est_ranking=row_number(est_survival)) %>% #higher ranking means longer survival 
		select(subject, model, treatment, time, event, est_effect, est_survival, est_ranking, fold) 
}

two_model_surv_rf = function(data, tune_grid, fold, fold_name) {
	training_data = data %>% 
		filter(subject %in% fold) %>%
		select(time, event, treatment, starts_with("covariate"))
	validation_data = data %>% 
		filter(!(subject %in% fold))
	subject = validation_data %>% pull(subject)

	tune_grid %>%
		pmap(function(num.trees, mtry, min.node.size) {
			models = list("TRUE"=TRUE, "FALSE"=FALSE) %>%
				map(function(treatment) {
					ranger(dependent.variable.name="time", status.variable.name="event", 
						   data=(training_data %>% filter(treatment==treatment) %>% select(-treatment)), 
			   			   mtry=mtry, num.trees=num.trees, min.node.size=min.node.size) %>%
					predict(validation_data) %>%
					extract_mean_survival() %>%
					mutate(subject=subject, mtry=mtry, num.trees=num.trees, min.node.size=min.node.size)
				}) %>%
		    	imap(function(df,name) rename_(df, .dots=setNames("pred", str_c("mean_survival", name, sep="_")))) %>%
		    	reduce(inner_join, by=c("subject", names(tune_grid))) 
		}) %>%
    	bind_rows() %>% 
	    mutate(method="two_model_surv_rf", fold=fold_name) %>%
	    unite_("model", c("method", names(tune_grid)), sep="~") %>%
		inner_join(validation_data, by='subject') %>%
	    group_by(model) %>%
		mutate(est_effect=mean_survival_TRUE - mean_survival_FALSE, 
	    	   est_survival=treatment*(mean_survival_TRUE) + (1-treatment)*mean_survival_FALSE,
	    	   est_ranking=row_number(est_survival)) %>% #higher ranking means longer survival 
		select(subject, model, treatment, time, event, est_effect, est_survival, est_ranking, fold) 
}

# change this so that "method" is the function that the data and options get passed to and everything else happens internally
cross_estimate_hte = function(data, method, tune_grid, train_index) {
	train_index %>%
	imap(function(fold, fold_name) method(data, tune_grid, fold, fold_name)) %>%
	bind_rows() 
}


# methods = list("xgbTree", "lm")
# tune_grids = list(expand.grid(nrounds = 1:150, max_depth = 2, eta = 0.1), NULL)
