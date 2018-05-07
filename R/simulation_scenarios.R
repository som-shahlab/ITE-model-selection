

#' Generate simulated observational data
#'
#' Generates N tuples from the joint distribution of P(X,W,Y).
#' Returns a list with two elements: a dataframe called data and a dataframe 
#' called aux_data containing true expecations conditional on X.
#' @param covariate_fun a function that randomly takes n samples from P(X)
#' @param mean_fun a vectorized function of a vector-valued arugment x that defines E[Y|X=x]
#' @param effect_fun a vectorized function of a single vector-valued arugment x that defines E[Y|W=1,X=x] - E[Y|W=0,X=x]
#' @param propensity_fun a vectorized function of a single vector-valued arugment x that defines E[W|X=x]
#' @param sigma standard deviation of the normally distributed zero-mean additive noise to Y
#' @param randomized do we know the true propensities?
#' @export
dgp = function(covariate_fun, propensity_fun, mean_fun, effect_fun, sigma, randomized=F) {
    list(covariate_fun=covariate_fun, mean_fun=mean_fun, propensity_fun=propensity_fun, 
         effect_fun=effect_fun, sigma=sigma, randomized=randomized)
}

#' Generate simulated observational data
#'
#' Generates N tuples from the joint distribution of P(X,W,Y).
#' Returns a list with two elements: a dataframe called data and a dataframe 
#' called aux_data containing true expecations conditional on X.
#' @param DGP a data generating process, which is the list output of the dgp() function
#' @param n the number of desired samples
#' @export
create_data = function(DGP, n=1) {
    x = DGP$covariate_fun(n)

    mu = DGP$mean_fun(x)
    tau = DGP$effect_fun(x)
    p = DGP$propensity_fun(x)

    w = rbinom(n, 1, p) %>% as.logical
    y = mu + tau * (2*w - 1) / 2 + rnorm(n, 0, DGP$sigma)

    x %<>% data.frame %>% 
        set_names(paste("covariate", 1:ncol(.), sep="_")) %>%
        make_matrix()

    mu1=mu+0.5*tau 
    mu0=mu-0.5*tau

    return(list(x=x, w=w, y=y, p=p, mu0=mu0, mu1=mu1, tau=tau))
}

x1 = function(n, p=10) {
    X = matrix(rnorm(n*p), nrow = n, ncol = p)
    X[, seq(2, p, by = 2)] = (X[, seq(2, p, by = 2)] > 0)
    return(X)
}

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

create_p_rand = function(p=0.5) {
    function(x) {
        rep(p, nrow(x))
    }
}

create_p_bias = function(mean_fun, effect_fun) {
    function(x) {
        log_p = (mean_fun(x) - effect_fun(x))/2
        exp(log_p) / (1 + exp(log_p))
    }
}

#' Simulations from Powers et al. 2018
#'
#' Creates a set of data generating processes from Powers et al 2018 that can be passed to the
#' create_data function to simulate observational data
#' @export
powers_DGPs = function() {
    unbiased_DGPs = list(
        scenario_1 = dgp(x1, create_p_rand(), f8, f1, 1, randomized=T),
        scenario_2 = dgp(x1, create_p_rand(), f5, f2, 0.25, randomized=T),
        scenario_3 = dgp(x1, create_p_rand(), f4, f3, 1, randomized=T),
        scenario_4 = dgp(x1, create_p_rand(), f7, f4, 0.25, randomized=T),
        scenario_5 = dgp(x1, create_p_rand(), f3, f5, 1, randomized=T),
        scenario_6 = dgp(x1, create_p_rand(), f1, f6, 1, randomized=T),
        scenario_7 = dgp(x1, create_p_rand(), f2, f7, 4, randomized=T),
        scenario_8 = dgp(x1, create_p_rand(), f6, f8, 4, randomized=T)
    )
    biased_DGPs = unbiased_DGPs %>%
        map(~dgp(.$covariate_fun, 
             create_p_bias(.$mean_fun, .$effect_fun), 
             .$mean_fun, .$effect_fun, .$sigma,
             randomized=F))
    names(biased_DGPs) = names(unbiased_DGPs) %>% map(~str_c("biased", ., sep="_"))
    c(unbiased_DGPs, biased_DGPs)
}