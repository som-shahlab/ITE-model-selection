#' @import dplyr
#' @import purrr
#' @import tidyr
#' @import magrittr
#' @import stringr
#' @import caret
#' @import Matching
#' @import distr
#' @import distrEx

## Each patient will now get their own control Y0 and treatment Y1 distributions from 
# package distr, with parameters f(X) for more or less arbitrary f()
# Then we'll calculate tau as E[Y1] - E[Y0] using the E() function
# from distrEx. 
# yi = w_i*r(Y1) + (1-w_i)*r(Y0)
# 

data_list_to_df = function(data_list) {
    data_list$covariates %<>% data.frame %>% set_names(paste("covariate", 1:ncol(.), sep="_"))
    data_list %$% cbind(subject, event, time, treatment, covariates) %>% data.frame
}

data_df_to_list = function(data_df) {
    data_list = list()
    data_list$subject = data_df %>% pull(subject)
    data_list$outcome = data_df %>% pull(outcome) # an indicator: 1 for death, 0 for censoring
    data_list$time = data_df %>% pull(time) # the time of either death or censoring
    data_list$treatment = data_df %>% pull(treatment)
    data_list$covariates = as.matrix(data_df %>% select(starts_with()))
    return(data_list)
}

#' Generate simulated observational data
#'
#' Generates N tuples from the joint distribution of P(X,W,Y).
#' Returns a list with two elements: a dataframe called data and a dataframe 
#' called aux_data containing true expecations conditional on X.
#' @param X a list of distribution objects that generate a covariate vector
#' @param f_W_x a function of x that returns the conditional distribution W|X=x
#' @param f_Y_xw a function of x, w that returns the conditional distribution Y|X=x,W=w
#' @param f_C_xw a function of x, w that returns the conditional distribution C|X=x,W=w
#' @keywords
#' @export
#' @examples
dgp = function() {
    list(X=X, f_W_x=f_W_x, f_Y_xw=f_Y_xw, f_C_xw=f_C_xw)
}
# These functions e.g. f_W_X will likely be generated on-the-fly:
# make_f_W_X_dist = function(f1, f2, ...) {
#       function(x) SomeDist(param1 = f1(x), param2 = f2(x)... )
# }

sample1 = function(dist) {
    r(dist)(1)
}

#' Generate simulated observational data
#'
#' Generates N tuples from the joint distribution of P(X,W,Y).
#' Returns a list with two elements: a dataframe called data and a dataframe 
#' called aux_data containing true expecations conditional on X.
#' @param DGP a data generating process, which is the list output of the dgp() function
#' @param n the number of desired samples
#' @keywords
#' @export
#' @examples
create_data = function(DGP, n=1) {
    x = n %>% rerun(DGP$X %>% map_dbl(sample1)) 

    W = x %>% map(DGP$f_W_x) 
    w = W %>% map_dbl(sample1)
    pw = W %>% map_dbl(E) # propensity scores

    Y1 = x %>% map(~DGP$f_Y_xw(.,1))
    Y0 = x %>% map(~DGP$f_Y_xw(.,0)) 
    mu1 = Y1 %>% map_dbl(E)
    mu0 = Y0 %>% map_dbl(E)
    tau = mu1 - mu0
    y = list(w,Y1,Y0) %>% 
        pmap_dbl(function(w, Y1, Y0) ifelse(w, sample1(Y1), sample1(Y0)))

    C1 = x %>% map(~DGP$f_C_xw(.,0))
    C0 = x %>% map(~DGP$f_C_xw(.,1))
    ce = list(w,C1,C0) %>% 
        pmap_dbl(function(w, C1, C0) ifelse(w, sample1(C1), sample1(C0)))

    t = pmin(y,ce)
    d = y < ce

    pc = list(w, t, C1, C0) %>% # censoring probability at t
         pmap_dbl(function(w,t,C1,C0) ifelse(w, 1-p(C1)(t), 1-p(C0)(t)))

    data = list()
    data$subject = 1:length(w)
    data$time = t
    data$event = d
    data$treatment = as.logical(w)
    data$covariates = x %>% reduce(rbind)
    rownames(data$covariates) = NULL

    aux_data = data.frame(
        subject=data$subject, 
        treated_mean=mu1, 
        control_mean=mu0, 
        effect=tau,
        opt_treatment = tau>0,
        iptw=1/(1-w + 2*w*pw - pw),
        ipcw=1/pc)

    return(list(data=(data %>% data_list_to_df), aux_data=aux_data))
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

#' Simulations 
#'
#' @export
schuler_DGPs = function() {
    X = list(Norm(-1), Norm(1))

    f_W_x = function(x, w) {
        logit_p = x[1] + x[2]
        p = exp(logit_p) / (1 + exp(logit_p))
        Binom(prob=p)
    }

    f_Y_xw = function(x, w) {
        if(w) {
            Weibull(scale=abs(x[1] + x[2]) + 0.7, 
                    shape=1.2)
        } else {
            Weibull(scale= abs(x[1]) + abs(x[2]), 
                    shape=1.5)
        }
    }

    f_C_xw = function(x, w) {
        Weibull(scale=4, 
                shape=1.4)
    }

    list("biased" = dgp(X, function(x,w) 0.5, f_Y_xw, f_C_xw), 
         "unbiased" = dgp(X, f_W_x, f_Y_xw, f_C_xw))
}