library(tidyverse)

f1 = function(x) rep(0, nrow(x))

f2 = function(x) 6 * (x[, 1] > 1) - 1 - 
    (6*pnorm(-1) -1) # take out the expectation

f3 = function(x) 5 * x[, 1]

f4 = function(x) {
    1 * x[,2] * x[,4] * x[,6] + 
    2 * x[,2] * x[,4] * (1-x[,6]) +
    3 * x[,2] * (1-x[,4]) * x[,6] + 
    4 * x[,2] * (1-x[,4]) * (1-x[,6]) +
    5 * (1-x[,2]) * x[,4] * x[,6] + 
    6 * (1-x[,2]) * x[,4] * (1-x[,6]) +
    7 * (1-x[,2]) * (1-x[,4]) * x[,6] + 
    8 * (1-x[,2]) * (1-x[,4]) * (1-x[,6]) - 
    5 -
    (-0.5) # take out the expectation
}

f5 = function(x) x[, 1] + x[, 3] + x[, 5] + x[, 7] + x[, 8] + x[, 9] -
    (0.5) # take out the expectation

f6 = function(x) {
    4 * (x[,1]>1) * (x[,3]>0) + 4 * (x[,5]>1) * (x[,7]>0) + 2 * x[,8] * x[,9] - 1 -
    (4*pnorm(-1)-1) # take out the expectation
}

f7 = function(x) {
  ( x[, 1]^2 + 
    x[, 2] + 
    x[, 3]^2 + 
    x[, 4] + 
    x[, 5]^2 + 
    x[, 6] + 
    x[, 7]^2 +
    x[, 8] + 
    x[, 9]^2 ) / 
    sqrt(2) - 5 -
    (7/sqrt(2) - 5) # take out the expectation
}

f8 = function(x) (f4(x) + f5(x)) / sqrt(2) 
    # ((-0.5 - )sqrt(2)) # take out the expectation

f9 = function(x) (x[,1])^2 - 1

create_unbiased_propensity = function(p=0.5) {
    function(x) {
        rep(p, nrow(x))
    }
}

create_biased_propensity = function(mean_fun, effect_fun) {
    function(x) {
        log_p = (mean_fun(x) - effect_fun(x))/2
        exp(log_p) / (1 + exp(log_p))
    }
}

# Scenarios
unbiased_DGPs = list(
    scenario_1 = list(mean_fun=f8, effect_fun=f1, propensity_fun=create_unbiased_propensity(), sigma=1),
    scenario_2 = list(mean_fun=f5, effect_fun=f2, propensity_fun=create_unbiased_propensity(), sigma=1),
    scenario_3 = list(mean_fun=f4, effect_fun=f3, propensity_fun=create_unbiased_propensity(), sigma=1),
    scenario_4 = list(mean_fun=f7, effect_fun=f4, propensity_fun=create_unbiased_propensity(), sigma=1),
    scenario_5 = list(mean_fun=f3, effect_fun=f5, propensity_fun=create_unbiased_propensity(), sigma=1),
    scenario_6 = list(mean_fun=f1, effect_fun=f6, propensity_fun=create_unbiased_propensity(), sigma=1),
    scenario_7 = list(mean_fun=f2, effect_fun=f7, propensity_fun=create_unbiased_propensity(), sigma=1),
    scenario_8 = list(mean_fun=f6, effect_fun=f8, propensity_fun=create_unbiased_propensity(), sigma=1)
)
biased_DGPs = unbiased_DGPs %>%
    map(~list(mean_fun=.$mean_fun, 
              effect_fun=.$effect_fun, 
              propensity_fun=create_biased_propensity(.$mean_fun, .$effect_fun), 
              sigma=.$sigma))
names(biased_DGPs) = names(unbiased_DGPs) %>% map(~str_c("biased", ., sep="_"))
DGPs = c(unbiased_DGPs, biased_DGPs)

data_list_to_df = function(data_list) {
    data_list$covariates %<>% set_names(paste("covariate", 1:ncol(.), sep="_"))
    data_list %$% cbind(subject, outcome, treatment, covariates) %>% data.frame
}

data_df_to_list = function(data_df) {
    data_list = list()
    data_list$subject = data_df %>% pull(subject)
    data_list$outcome = data_df %>% pull(outcome)
    data_list$treatment = data_df %>% pull(treatment)
    data_list$covariates = as.matrix(data_df %>% select(starts_with()))
    return(data_list)
}

generate_covariates = function(n, p=10) {
    X = matrix(rnorm(n*p), nrow = n, ncol = p)
    X[, seq(2, p, by = 2)] = (X[, seq(2, p, by = 2)] > 0)
    return(X)
}

create_data = function(mean_fun, propensity_fun, effect_fun, sigma, n, p=10) {
    X = generate_covariates(n,p)

    mu = mean_fun(X)
    ef = effect_fun(X)
    p = propensity_fun(X)

    W = rbinom(n, 1, p)
    Y = mu + ef * (2*W - 1) / 2 + rnorm(n, 0, sigma)

    data = list()
    data$subject = 1:length(Y)
    data$outcome = Y
    data$treatment = as.logical(W)
    data$covariates = X %>% data.frame

    return(data %>% data_list_to_df)
}