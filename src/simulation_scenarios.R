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

setup_simulation = function(scenario) {
  
  if (scenario == 0) {
    mean.function = function(x) f2(x)
    propensity.function = function(x) rep(.5, nrow(x))
    effect.function = function(x) f9(x)
    m = 1000
    n = 200
    p = 400
    sigma = 1
    generate.data = function(m,n,p) {
      x = matrix(rnorm((m + n) * p), nrow = m + n, ncol = p)
      x[, seq(2, p, by = 2)] = x[, seq(2, p, by = 2)] > 0
      x
    }
  }

  if (scenario == 1) {
    mean.function = function(x) f8(x)
    propensity.function = function(x) rep(.5, nrow(x))
    effect.function = function(x) f1(x)
    m = 1000
    n = 200
    p = 400
    sigma = 1
    generate.data = function(m,n,p) {
      x = matrix(rnorm((m + n) * p), nrow = m + n, ncol = p)
      x[, seq(2, p, by = 2)] = x[, seq(2, p, by = 2)] > 0
      x
    }
  }

  if (scenario == 2) {
    mean.function = function(x) f5(x)
    propensity.function = function(x) rep(.5, nrow(x))
    effect.function = function(x) f2(x)
    m = 5000
    n = 200
    p = 400
    sigma = 1/2
    generate.data = function(m,n,p) {
      x = matrix(rnorm((m + n) * p), nrow = m + n, ncol = p)
      x[, seq(2, p, by = 2)] = x[, seq(2, p, by = 2)] > 0
      x
    }
  }

  if (scenario == 3) {
    mean.function = function(x) f4(x)
    propensity.function = function(x) rep(.5, nrow(x))
    effect.function = function(x) f3(x)
    m = 1000
    n = 300
    p = 300
    sigma = 4
    generate.data = function(m,n,p) {
      x = matrix(rnorm((m + n) * p), nrow = m + n, ncol = p)
      x[, seq(2, p, by = 2)] = x[, seq(2, p, by = 2)] > 0
      x
    }
  }

  if (scenario == 4) {
    mean.function = function(x) f7(x)
    propensity.function = function(x) rep(.5, nrow(x))
    effect.function = function(x) f4(x)
    m = 1000
    n = 300
    p = 300
    sigma = 1/2
    generate.data = function(m,n,p) {
      x = matrix(rnorm((m + n) * p), nrow = m + n, ncol = p)
      x[, seq(2, p, by = 2)] = x[, seq(2, p, by = 2)] > 0
      x
    }
  }

  if (scenario == 5) {
    mean.function = function(x) f3(x)
    propensity.function = function(x) rep(.5, nrow(x))
    effect.function = function(x) f5(x)
    m = 1000
    n = 400
    p = 200
    sigma = 1
    generate.data = function(m,n,p) {
      x = matrix(rnorm((m + n) * p), nrow = m + n, ncol = p)
      x[, seq(2, p, by = 2)] = x[, seq(2, p, by = 2)] > 0
      x
    }
  }

  if (scenario == 6) {
    mean.function = function(x) f1(x)
    propensity.function = function(x) rep(.5, nrow(x))
    effect.function = function(x) f6(x)
    m = 1000
    n = 400
    p = 200
    sigma = 1/4
    generate.data = function(m,n,p) {
      x = matrix(rnorm((m + n) * p), nrow = m + n, ncol = p)
      x[, seq(2, p, by = 2)] = x[, seq(2, p, by = 2)] > 0
      x
    }
  }

  if (scenario == 7) {
    mean.function = function(x) f2(x)
    propensity.function = function(x) rep(.5, nrow(x))
    effect.function = function(x) f7(x)
    m = 1000
    n = 1000
    p = 100
    sigma = 2
    generate.data = function(m,n,p) {
      x = matrix(rnorm((m + n) * p), nrow = m + n, ncol = p)
      x[, seq(2, p, by = 2)] = x[, seq(2, p, by = 2)] > 0
      x
    }
  }

  if (scenario == 8) {
    mean.function = function(x) f6(x)
    propensity.function = function(x) rep(.5, nrow(x))
    effect.function = function(x) f8(x)
    m = 1000
    n = 1000
    p = 100
    sigma = 2
    generate.data = function(m,n,p) {
      x = matrix(rnorm((m + n) * p), nrow = m + n, ncol = p)
      x[, seq(2, p, by = 2)] = x[, seq(2, p, by = 2)] > 0
      x
    }
  }

  # like 1 but with biased treatment assignment
  if (scenario == 9) {
    mean.function = function(x) f8(x)
    propensity.function = function(x) {
      f0 = f8(x) - f1(x) / 2
      exp(f0) / (1 + exp(f0))
    }
    effect.function = function(x) f1(x)
    m = 1000
    n = 200
    p = 400
    sigma = 1
    generate.data = function(m,n,p) {
      x = matrix(rnorm((m + n) * p), nrow = m + n, ncol = p)
      x[, seq(2, p, by = 2)] = x[, seq(2, p, by = 2)] > 0
      x
    }
  }

  # like 2 but with biased treatment assignment
  if (scenario == 10) {
    mean.function = function(x) f5(x)
    propensity.function = function(x) {
      f0 = f5(x) - f2(x) / 2
      exp(f0) / (1 + exp(f0))
    }
    effect.function = function(x) f2(x)
    m = 1000
    n = 200
    p = 400
    sigma = 1/2
    generate.data = function(m,n,p) {
      x = matrix(rnorm((m + n) * p), nrow = m + n, ncol = p)
      x[, seq(2, p, by = 2)] = x[, seq(2, p, by = 2)] > 0
      x
    }
  }

  # like 3 but with biased treatment assignment
  if (scenario == 11) {
    mean.function = function(x) f4(x)
    propensity.function = function(x) {
      f0 = f4(x) - f3(x) / 2
      exp(f0) / (1 + exp(f0))
    }
    effect.function = function(x) f3(x)
    m = 1000
    n = 300
    p = 300
    sigma = 4
    generate.data = function(m,n,p) {
      x = matrix(rnorm((m + n) * p), nrow = m + n, ncol = p)
      x[, seq(2, p, by = 2)] = x[, seq(2, p, by = 2)] > 0
      x
    }
  }

  # like 4 but with biased treatment assignment
  if (scenario == 12) {
    mean.function = function(x) f7(x)
    propensity.function = function(x) {
      f0 = f7(x) - f4(x) / 2
      exp(f0) / (1 + exp(f0))
    }
    effect.function = function(x) f4(x)
    m = 1000
    n = 300
    p = 300
    sigma = 1/2
    generate.data = function(m,n,p) {
      x = matrix(rnorm((m + n) * p), nrow = m + n, ncol = p)
      x[, seq(2, p, by = 2)] = x[, seq(2, p, by = 2)] > 0
      x
    }
  }

  # like 5 but with biased treatment assignment
  if (scenario == 13) {
    mean.function = function(x) f3(x)
    propensity.function = function(x) {
      f0 = f3(x) - f5(x) / 2
      exp(f0) / (1 + exp(f0))
    }
    effect.function = function(x) f5(x)
    m = 1000
    n = 400
    p = 200
    sigma = 1
    generate.data = function(m,n,p) {
      x = matrix(rnorm((m + n) * p), nrow = m + n, ncol = p)
      x[, seq(2, p, by = 2)] = x[, seq(2, p, by = 2)] > 0
      x
    }
  }

  # like 6 but with biased treatment assignment
  if (scenario == 14) {
    mean.function = function(x) f1(x)
    propensity.function = function(x) {
      f0 = f1(x) - f6(x) / 2
      exp(f0) / (1 + exp(f0))
    }
    effect.function = function(x) f6(x)
    m = 1000
    n = 400
    p = 200
    sigma = 1/4
    generate.data = function(m,n,p) {
      x = matrix(rnorm((m + n) * p), nrow = m + n, ncol = p)
      x[, seq(2, p, by = 2)] = x[, seq(2, p, by = 2)] > 0
      x
    }
  }

  # like 7 but with biased treatment assignment
  if (scenario == 15) {
    mean.function = function(x) f2(x)
    propensity.function = function(x) {
      f0 = f2(x) - f7(x) / 2
      exp(f0) / (1 + exp(f0))
    }
    effect.function = function(x) f7(x)
    m = 1000
    n = 1000
    p = 100
    sigma = 2
    generate.data = function(m,n,p) {
      x = matrix(rnorm((m + n) * p), nrow = m + n, ncol = p)
      x[, seq(2, p, by = 2)] = x[, seq(2, p, by = 2)] > 0
      x
    }
  }

  # like 8 but with biased treatment assignment
  if (scenario == 16) {
    mean.function = function(x) f6(x)
    propensity.function = function(x) {
      f0 = f6(x) - f8(x) / 2
      exp(f0) / (1 + exp(f0))
    }
    effect.function = function(x) f8(x)
    m = 1000
    n = 1000
    p = 100
    sigma = 2
    generate.data = function(m,n,p) {
      x = matrix(rnorm((m + n) * p), nrow = m + n, ncol = p)
      x[, seq(2, p, by = 2)] = x[, seq(2, p, by = 2)] > 0
      x
    }
  }

    sim_functions = list()
    sim_functions$data = generate.data
    sim_functions$mean = mean.function
    sim_functions$effect = effect.function
    sim_functions$propensity = propensity.function
    sim_functions$sigma = sigma

    sim_params = list()
    sim_params$m = m
    sim_params$n = n
    sim_params$p = p
    sim_params$sigma = sigma

    setup = list()
    setup$functions = sim_functions
    setup$params = sim_params

    return(setup)
}

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

create_data = function(setup) {
    params = setup$params
    m = params$m
    n = params$n
    p = params$p
    sigma = params$sigma
    
    functions = setup$functions    
    X = functions$data(m,n,p)
    mu = functions$mean(X)
    ef = functions$effect(X)
    p = functions$propensity(X)

    W = rbinom(m + n, 1, p)
    Y = mu + ef * (2*W - 1) / 2 + rnorm(m + n, 0, sigma)
    
    data = list()
    data$subject = 1:length(Y)
    data$outcome = Y
    data$treatment = as.logical(W)
    data$covariates = X %>% data.frame

    return(data %>% data_list_to_df)
}

# split_cv = function(data, n_folds) {
#     fold = cut(1:nrow(data),breaks=n_folds,labels=FALSE)
#     data = data %>%
#         arrange(sample(1:nrow(.))) %>%
#         mutate(fold = fold) %>%
#         arrange(subject)
#     return(data)
# }


## Use like:
## data = setup_simulation(sim_num) %>% create_data() %$% cbind(X,W,Y) %>% data.frame