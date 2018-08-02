#' @import magrittr
#' @import caret
#' @import zeallot
#' @import Matching
#' @import tidyverse

#' @export
make_indices = function(n_train, n_val, n_test) {
    itrain = 1:n_train
    ival = (n_train+1):(n_train+n_val)
    itest = (n_train+n_val+n_train+1):(n_train+n_val+n_train)
    list(itrain, ival, itest)
}

#' @export
cv_indices = function(itrain, ival, n_folds) {
    train_total = c(itrain, ival)
    n_train = length(train_total)
    splits = sample(train_total) %>%
        split(rep(1:n_folds, ceiling(n_train/n_folds)))
    icv = 1:n_folds %>% map(function(fold){
        list(
            train = splits[-fold] %>% unlist(use.names=F),
            val = splits[[fold]])
    })
}

# returns the index of the observation most similar to the observation in that position:
# i.e. (1->5), (2->3), ... 
find_matches = function(x, w) {
    index = 1:length(w)
    treated_match = Match(Tr=w, X=x, replace=T, estimand="ATT") 
    control_match = Match(Tr=!w, X=x, replace=T, estimand="ATT")
    subject = index[c(treated_match$index.treated, control_match$index.treated)] # all subjects
    match = index[c(treated_match$index.control, control_match$index.control)] # their matches
 
    match_sorted = rep(NA,length(index))
    match_sorted[subject] = match
    return(match_sorted)
}

make_matrix = function(x) stats::model.matrix(~.-1, x)

subset_rows = function(items, index) {
    items %>% map(function(item) {
        if (is.vector(item)) item[index]
        else if (is.array(item)) item[index,]
        else stop("Don't know what to do with item of this type")
    })
}

iptw = function(p,w) 1/(1-w + 2*w*p - p)

# estimate many treatment effects models on the training set and 
# get validation + test set estimates from each
#' @export
estimate_val_test = function(data, itrain, model_specs) {
    c(x, w, y, ....) %<-% data
    list(
        S=S_learners_pred_test, 
        T=T_learners_pred_test, 
        R=R_learners_pred_test) %>%
    imap(function(learner, learner_type){
        learner(itrain, x, w, y, model_specs) %>% # uses train to predict on val/test
        mutate(model = str_c(learner_type, model, sep="$"))
    }) %>% bind_rows
}

# cross-estimate mean outcome and propensity score using cross-validated models 
# on the validation set and match on validation set (VALIDATION)
#' @export
learn_validation_auxiliaries = function(data, ival, model_specs, randomized=F) {
    data %>% 
    subset_rows(ival) %->%
    c(x_val, w_val, y_val, p_val, ....)

    mu_hat_val = cross_validated_cross_estimation(x_val, y_val, model_specs)
    c(mu0_hat_val, mu1_hat_val) %<-% S_learner_cv_ce(x_val, w_val, y_val, model_specs)
    if (randomized) {
        p_hat_val = p_val
    } else {
        p_hat_val = cross_validated_cross_estimation(x_val, w_val, model_specs)
    }

    list(
        matches = find_matches(x_val, w_val),
        w = w_val,
        y = y_val,
        p = p_hat_val,
        mu0 = mu0_hat_val,
        mu1 = mu1_hat_val,
        pseudo_outcome = (y_val - mu_hat_val)/(w_val - p_hat_val),
        weight = (w_val - p_hat_val)^2,
        iptw = iptw(p_hat_val, w_val))
}

# use those auxiliary data and estimates on the validation set to estimate all validation metrics
#' @export
estimate_val_metrics = function(estimates, val_bundle, metrics, ival) {
    val_estimates_grp = estimates %>% # estimates on test and validation sets
        filter(index %in% ival) %>%
        group_by(model)
    metrics %>% imap(function(metric, metric_name) {
        val_estimates_grp %>% 
            summarize(!!metric_name:=metric(val_bundle, y_hat, tau_hat))
    }) %>% reduce(inner_join, by="model")
}

#' @export
estimate_learner_test = function(data, learners, model_specs, itrain, ival, itest) {
    learners %>%
    imap(function(learner, learner_name) {
        data %>% 
            subset_rows(c(itrain, ival)) %->%
            c(xtrain, wtrain, ytrain, ...) 
        tau_hat = learner(xtrain, wtrain, ytrain, model_specs) %>%
            predict(data$x[itest,])
        data.frame(tau_hat=tau_hat, index=itest, model=learner_name, yhat=NA)
    }) %>% bind_rows()
}

# apply tMSE and value to test-set estimates from all models
#' @export
calc_test_metrics = function(data, estimates, itest) {
    data %>% 
    subset_rows(itest) %->%
    c(..., mu0_test, mu1_test, tau_test)
    
    estimates %>% 
        filter(index %in% itest) %>% # index here should be a column inside of estimates
        group_by(model) %>%
        summarize(
            tmse = mse(tau_test, tau_hat),
            value = mean(mu1_test*(tau_hat>0) + mu0_test*(tau_hat<=0)))
}
