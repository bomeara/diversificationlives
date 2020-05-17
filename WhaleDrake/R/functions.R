
# Fits a model where lambda and mu each have 7 parameters, like the seed plant model in Louca & Pennell 2020.
# Note that drake will not work with multiple cores called inside functions https://github.com/ropensci/drake/issues/675#issuecomment-458222414
param_lambda7p_mu7p <- function(desired_interval = 0.1, tree, condition="crown", ncores=1) {
    root_age = castor::get_tree_span(tree)$max_distance
    rho = 1
    age_grid_param = seq(from=0,to=root_age+desired_interval,by=desired_interval)

    lambda_function = function(ages,params){
        results <- params['p1'] * exp(-params['p2']*ages/max(ages)) +
                   params['p3'] +
                   params['p4'] * (ages/max(ages)) +
                   params['p5'] * ((ages/max(ages))^2) +
                   params['p6'] * ((ages/max(ages))^3) +
                   params['p7'] * ((ages/max(ages))^4)
       results[which(results<0)] <- 0
       return(results)
    }
    mu_function = function(ages,params){
        results <-  params['q1'] * exp(-params['q2']*ages/max(ages)) +
                   params['q3'] +
                   params['q4'] * (ages/max(ages)) +
                   params['q5'] * ((ages/max(ages))^2) +
                   params['q6'] * ((ages/max(ages))^3) +
                   params['q7'] * ((ages/max(ages))^4)
        results[which(results<0)] <- 0
        return(results)
    }

    rho_function = function(params){
        return(rho) # rho does not depend on any of the parameters
    }

    param_values <- rep(NA, 14)
    names(param_values) <- c(paste0("p", sequence(7)), paste0("q", sequence(7)))
    param_guess <- rep(0, 14)
    param_guess[c(3, 3+7)] <- 0.1 #the age independent rate
    fit_param = NA
    try( {
        fit_param <- castor::fit_hbd_model_parametric(	tree,
                                               param_values  = param_values,
                                               param_guess   = param_guess,
                                               param_min     = rep(0, 14),
                                               param_max     = rep(2, 14),
                                               param_scale   = 1, # all params are in the order of 1
                                               lambda        = lambda_function,
                                               mu            = mu_function,
                                               rho0          = rho_function,
                                               age_grid      = age_grid_param,
                                               condition     = condition,
                                               Ntrials       = 10,    # perform 10 fitting trials
                                               Nthreads      = ncores,
                                               fit_control       = list(rel.tol=1e-6)
                                          )
        })
    return(list(fit_param=fit_param, lambda_function=lambda_function, mu_function=mu_function, age_grid_param=age_grid_param))
}


# Fits a model where mu is a constant multiple of lambda's value
param_lambda7p_mu_multiplier_lambda <- function(desired_interval = 0.1, tree, condition="crown", ncores=1) {
    root_age = castor::get_tree_span(tree)$max_distance
    rho = 1
    age_grid_param = seq(from=0,to=root_age+desired_interval,by=desired_interval)

    lambda_function = function(ages,params){
        results <- params['p1'] * exp(-params['p2']*ages/max(ages)) +
                   params['p3'] +
                   params['p4'] * (ages/max(ages)) +
                   params['p5'] * ((ages/max(ages))^2) +
                   params['p6'] * ((ages/max(ages))^3) +
                   params['p7'] * ((ages/max(ages))^4)
       results[which(results<0)] <- 0
       return(results)
    }
    mu_function = function(ages,params){
        results <-  lambda_function(ages, params)
        results <- results*params['qmultiplier']
        results[which(results<0)] <- 0
        return(results)
    }

    rho_function = function(params){
        return(rho) # rho does not depend on any of the parameters
    }

    param_values <- rep(NA, 8)
    names(param_values) <- c(paste0("p", sequence(7)), "qmultiplier")
    param_guess <- rep(0, 8)
    param_guess[c(3)] <- 0.1 #the age independent rate
    fit_param = NA
    try({
        fit_param <- fit_hbd_model_parametric(	tree,
                                               param_values  = param_values,
                                               param_guess   = param_guess,
                                               param_min     = rep(0, length(param_values)),
                                               param_max     = rep(2, length(param_values)),
                                               param_scale   = 1, # all params are in the order of 1
                                               lambda        = lambda_function,
                                               mu            = mu_function,
                                               rho0          = rho_function,
                                               age_grid      = age_grid_param,
                                               condition     = condition,
                                               Ntrials       = 10,    # perform 10 fitting trials
                                               Nthreads      = ncores,
                                               fit_control       = list(rel.tol=1e-6)
                                           )
    })
    return(list(fit_param=fit_param, lambda_function=lambda_function, mu_function=mu_function, age_grid_param=age_grid_param))
}
