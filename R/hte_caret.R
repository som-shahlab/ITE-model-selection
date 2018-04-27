#' @import BART
#' @import ranger

# rewrite a function like this for all necessary methods
#' @export
gbm_ph = function(data, tune_grid, training_fold, training_fold_name) {
	training_data = data %>% 
		filter(subject %in% training_fold)
	validation_data = data %>% 
		filter(!(subject %in% training_fold)) 

	x_train = training_data %>% dplyr::select(treatment, starts_with("covariate")) %>% make_matrix
    y_train = training_data %$% Surv(time, event)
    x_val = validation_data %>% dplyr::select(treatment, starts_with("covariate"))
    subject = validation_data %>% pull(subject)

    x_val_cf = list("TRUE"=TRUE, "FALSE"=FALSE) %>% map(~mutate(x_val, treatment=.) %>% make_matrix) # makes list(x_val_0, x_val_1) where treatments have been replaced 

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
	    mutate(method="gbm_ph", training_fold=training_fold_name) %>%
	    unite_("model", c("method", names(tune_grid)), sep="~") %>%
	    inner_join(validation_data, by='subject') %>%
	    group_by(model) %>%
		mutate(log_relative_risk = relative_risk_to_basline_TRUE - relative_risk_to_basline_FALSE,
			   est_treatment = relative_risk_to_basline_FALSE > relative_risk_to_basline_TRUE, # false > true so that positive favors treatment since this is log relative risk
	    	   est_pseudo_outcome=treatment*(relative_risk_to_basline_TRUE) + (1-treatment)*relative_risk_to_basline_FALSE,
	    	   est_ranking=row_number(-est_pseudo_outcome)) %>% #higher ranking means longer survival 
		select(subject, model, treatment, time, event, est_treatment, est_ranking, training_fold) 
}

#' @export
treat_all = function(data, tune_grid, training_fold, training_fold_name) {
	data %>% 
		mutate(est_treatment = TRUE, est_ranking=1, training_fold=training_fold_name, model="treat_all") %>%
		select(subject, model, treatment, time, event, est_treatment, est_ranking, training_fold)
}

#' @export
treat_none = function(data, tune_grid, training_fold, training_fold_name) {
	data %>% 
		mutate(est_treatment = FALSE, est_ranking=1, training_fold=training_fold_name, model="treat_none") %>%
		select(subject, model, treatment, time, event, est_treatment, est_ranking, training_fold)
}

#' @export
extract_mean_survival_ranger = function(ranger_pred) {
	dt = data.frame(t = ranger_pred$unique.death.times) %>% 
		arrange(t) %>%
	    mutate(dt=lead(t)-t, dt=ifelse(is.na(dt), 0, dt)) %>%
	    pull(dt)
	data.frame(pred=(ranger_pred$survival %*% dt))
}

# use treatment as always.split.variables? 
# tune_grid should include: n_times, power, base, ntree, ndpost, nskip
#' @export
one_model_surv_bart = function(data, tune_grid, training_fold, training_fold_name) {
	training_data = data %>% 
		filter(subject %in% training_fold)
	validation_data = data %>% 
		filter(!(subject %in% training_fold)) 

	x_train = training_data %>% dplyr::select(treatment, starts_with("covariate")) %>% make_matrix
    d_train = training_data %>% pull(event)
    x_val = validation_data %>% dplyr::select(subject, treatment, starts_with("covariate"))
    x_val_cf_df = list("TRUE"=TRUE, "FALSE"=FALSE) %>% 
		map(~mutate(x_val, treatment=.)) %>%
		bind_rows() 
    subject = x_val_cf_df %>% pull(subject)
	n_subject = length(subject)
    x_val_cf = x_val_cf_df %>% select(-subject) %>% make_matrix

    tune_grid %>% 
    	pmap(function(n_times, power, base, ntree, ndpost, nskip) {
    		t_train = training_data %>% 
		    	mutate(event_interval = ntile(time, n_times)) %>%
		    	group_by(event_interval) %>%
		    	mutate(time = min(time)) %>% # group times together to prevent data explosion
		    	pull(time)

			subject_index = subject %>% map(~rep(.,n_times)) %>% unlist
			treatment_index = x_val_cf_df %>% pull(treatment) %>% map(~rep(.,n_times)) %>% unlist
			t_index = rep(t_train %>% unique %>% sort, n_subject)

			print(data.frame(t_train, d_train))

		    mc.surv.bart(times=t_train, delta=d_train, x.train=x_train, x.test=x_val_cf,
		    		  power=power, base=base, ntree=ntree, ndpost=ndpost, nskip=nskip) %$%
		    	surv.test %>%
		    	t %>% data.frame %>%
		    	mutate(subject=subject_index, time=t_index, treatment=treatment_index) %>%
		    	gather(mcmc_sample, survival, -subject, -time, -treatment) %>%
		    	group_by(subject, mcmc_sample, treatment) %>%
		    	arrange(time) %>%
		    	mutate(dt=lead(time)-time, dt=ifelse(is.na(dt), 0, dt)) %>%
		    	summarize(est_survival = sum(dt*survival)) %>%
		    	ungroup() %>% group_by(subject, treatment) %>%
		    	summarize(est_survival = mean(est_survival)) %>%
		    	spread(treatment, est_survival) %>%
		    	rename("est_survival_TRUE"="TRUE", "est_survival_FALSE"="FALSE") %>%
		    	mutate(n_times=n_times, power=power, base=base, ntree=ntree, ndpost=ndpost, nskip=nskip)
    	}) %>%
    	bind_rows() %>%
    	mutate(method="one_model_surv_bart", training_fold=training_fold_name) %>%
		unite_("model", c("method", names(tune_grid)), sep="~") %>%
	    inner_join(validation_data, by='subject') %>%
		group_by(model) %>%
		mutate(est_effect=est_survival_TRUE - est_survival_FALSE, 
	    	   est_survival=treatment*(est_survival_TRUE) + (1-treatment)*est_survival_FALSE,
	    	   est_treatment = est_effect > 0,
	    	   est_ranking=row_number(est_survival)) %>% #higher ranking means longer survival 
		select(subject, model, treatment, time, event, est_effect, est_survival, est_treatment, est_ranking, training_fold) 
}

#' @export
one_model_surv_rf = function(data, tune_grid, training_fold, training_fold_name) {
	training_data = data %>% 
		filter(subject %in% training_fold) %>%
		select(time, event, treatment, starts_with("covariate"))
	validation_data = data %>% 
		filter(!(subject %in% training_fold))
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
		    		extract_mean_survival_ranger() %>%
		    		mutate(subject=subject, mtry=mtry, num.trees=num.trees, min.node.size=min.node.size)) %>%
		    	imap(function(df,name) rename_(df, .dots=setNames("pred", str_c("mean_survival", name, sep="_")))) %>%
		    	reduce(inner_join, by=c("subject", names(tune_grid))) 
		}) %>%
    	bind_rows() %>% 
	    mutate(method="one_model_surv_rf", training_fold=training_fold_name) %>%
	    unite_("model", c("method", names(tune_grid)), sep="~") %>%
	    inner_join(validation_data, by='subject') %>%
	    group_by(model) %>%
		mutate(est_effect=mean_survival_TRUE - mean_survival_FALSE, 
	    	   est_survival=treatment*(mean_survival_TRUE) + (1-treatment)*mean_survival_FALSE,
	    	   est_treatment = est_effect > 0,
	    	   est_ranking=row_number(est_survival)) %>% #higher ranking means longer survival 
		select(subject, model, treatment, time, event, est_effect, est_survival, est_treatment, est_ranking, training_fold) 
}

#' @export
two_model_surv_rf = function(data, tune_grid, training_fold, training_fold_name) {
	training_data = data %>% 
		filter(subject %in% training_fold) %>%
		select(time, event, treatment, starts_with("covariate"))
	validation_data = data %>% 
		filter(!(subject %in% training_fold))
	subject = validation_data %>% pull(subject)

	tune_grid %>%
		pmap(function(num.trees, mtry, min.node.size) {
			models = list("TRUE"=TRUE, "FALSE"=FALSE) %>%
				map(function(treatment) {
					ranger(dependent.variable.name="time", status.variable.name="event", 
						   data=(training_data %>% filter(treatment==treatment) %>% select(-treatment)), 
			   			   mtry=mtry, num.trees=num.trees, min.node.size=min.node.size) %>%
					predict(validation_data) %>%
					extract_mean_survival_ranger() %>%
					mutate(subject=subject, mtry=mtry, num.trees=num.trees, min.node.size=min.node.size)
				}) %>%
		    	imap(function(df,name) rename_(df, .dots=setNames("pred", str_c("mean_survival", name, sep="_")))) %>%
		    	reduce(inner_join, by=c("subject", names(tune_grid))) 
		}) %>%
    	bind_rows() %>% 
	    mutate(method="two_model_surv_rf", training_fold=training_fold_name) %>%
	    unite_("model", c("method", names(tune_grid)), sep="~") %>%
		inner_join(validation_data, by='subject') %>%
	    group_by(model) %>%
		mutate(est_effect=mean_survival_TRUE - mean_survival_FALSE, 
	    	   est_survival=treatment*(mean_survival_TRUE) + (1-treatment)*mean_survival_FALSE,
	    	   est_treatment = est_effect > 0,
	    	   est_ranking=row_number(est_survival)) %>% #higher ranking means longer survival 
		select(subject, model, treatment, time, event, est_effect, est_survival, est_treatment, est_ranking, training_fold) 
}

# change this so that "method" is the function that the data and options get passed to and everything else happens internally
cross_estimate_hte = function(data, method, tune_grid, train_index) {
	data %<>% select(subject, time, event, treatment, starts_with("covariate"))
	train_index %>%
	imap(function(training_fold, training_fold_name) method(data, tune_grid, training_fold, fold_name)) %>%
	bind_rows() 
}


# methods = list("xgbTree", "lm")
# tune_grids = list(expand.grid(nrounds = 1:150, max_depth = 2, eta = 0.1), NULL)
