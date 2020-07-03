
# Fits a model where lambda and mu each have 7 parameters, like the seed plant model in Louca & Pennell 2020.
# Doing https://books.ropensci.org/drake/hpc.html#advanced-options to allow parallel within
param_lambda7p_mu7p <- function(desired_interval = 0.1, tree, condition="crown", ncores=parallel::detectCores()) {
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
# Doing https://books.ropensci.org/drake/hpc.html#advanced-options to allow parallel within

param_lambda7p_mu_multiplier_lambda <- function(desired_interval = 0.1, tree, condition="crown", ncores=parallel::detectCores()) {
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

# desired_interval: how finely to divide the slices for calculating smoothly changing (or not) lambda and mu
# slice_ages: when to divide regimes
# interpolation_method: same as method in ?approx. Constant or linear
param_lambda_discreteshift_mu_discreteshift <- function(desired_interval = 0.1, tree, condition="crown", ncores=parallel::detectCores(), slice_ages = seq(from=0, to=ceiling(castor::get_tree_span(tree)$max_distance), by=1), interpolation_method="constant") {
    root_age = castor::get_tree_span(tree)$max_distance
    rho = 1
    age_grid_param = seq(from=0,to=root_age+desired_interval,by=desired_interval)

    lambda_params <- rep(NA, length(slice_ages))
    names(lambda_params) <- paste0("lambda", sequence(length(slice_ages))-1)

    mu_params <- rep(NA, length(slice_ages))
    names(mu_params) <- paste0("mu", sequence(length(slice_ages))-1)

    lambda_function = function(ages,params){
        results <- stats::approx(x=slice_ages, y=params[grepl("lambda", names(params))], xout=ages, method=interpolation_method, rule=2)$y

       return(results)
    }

    mu_function = function(ages,params){
        results <- stats::approx(x=slice_ages, y=params[grepl("mu", names(params))], xout=ages, method=interpolation_method, rule=2)$y

       return(results)
    }


    rho_function = function(params){
        return(rho) # rho does not depend on any of the parameters
    }
    param_values <- c(lambda_params, mu_params)
    param_guess <- c(rep(0.1, length(lambda_params)), rep(0.01, length(mu_params)))
    names(param_guess) <- names(param_values)

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

param_lambda_discreteshift_ef_fixed <- function(desired_interval = 0.1, tree, condition="crown", ncores=parallel::detectCores(), slice_ages = seq(from=0, to=ceiling(castor::get_tree_span(tree)$max_distance), by=1), interpolation_method="constant", ef=0.0) {
    root_age = castor::get_tree_span(tree)$max_distance
    rho = 1
    age_grid_param = seq(from=0,to=root_age+desired_interval,by=desired_interval)

    lambda_params <- rep(NA, length(slice_ages))
    names(lambda_params) <- paste0("lambda", sequence(length(slice_ages))-1)


    lambda_function = function(ages,params){
        results <- stats::approx(x=slice_ages, y=params[grepl("lambda", names(params))], xout=ages, method=interpolation_method, rule=2)$y
       return(results)
    }

    mu_function = function(ages,params){
        results <- lambda_function(ages, params)
        results <- ef*results
        return(results)
    }


    rho_function = function(params){
        return(rho) # rho does not depend on any of the parameters
    }
    param_values <- c(lambda_params)
    param_guess <- c(rep(0.1, length(param_values)))
    names(param_guess) <- names(param_values)
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


#' Split a tree into even sized chunks
#'
#' By default, much of Louca and Pennell's work splits a tree into regimes based on equal spacing of ages (i.e., the seed plant work). However, this means regimes close to the root have few data points, while those near the tips can have many. An alternative would be to split so that every regime has the same number of points
#' @param tree The phylo object. Assumed to be ultrametric
#' @param nregimes How many regimes to split the tree into
#' @param minsize Min number of taxa in a regime. If >1, it will decrease nregimes until this is met
#' @param type Either "data" or "time" the way to slice the tree: equal number of data points per regime or equal span of time for each
#' @param minage The end (tipmost) time of the most recent regime (default 0)
#' @param maxage The beginning (rootmost) time of the oldest regime (default root of tree)
#' @examples
#' phy <- ape::rcoal(100)
#' splits <- EvenSplit(phy, nregimes=4)
#' plot(phy)
#' abline(v=splits$time+max(ape::branching.times(phy)), col="red")
#' ltt.plot(phy, log="y")
#' abline(v=splits$time, col="red")
EvenSplit <- function(tree, nregimes, minsize=1, type="data", minage=0, maxage=castor::get_tree_span(tree)$max_distance) {
  if(!is.ultrametric(tree)) {
    stop("This assumes the tree is ultrametric")
  }
  ltt_points <- ape::ltt.plot.coords(tree)
  qualified <- FALSE
  while(!qualified) {
    if(type=="data") {
        splitrows <- round(sequence(nregimes-1)*(nrow(ltt_points)-1)/nregimes) #split after this row
        splittimes <- (ltt_points[splitrows,1]+ltt_points[splitrows+1,1])/2
    } else {
      splittimes <- sequence(nregimes-1)*min(ltt_points[,1])/(nregimes)
      splitrows <- c()
      for (i in seq_along(splittimes)) {
        splitrows[i] <- which(ltt_points[,1]>splittimes[i])[1]-1
      }
    }
    result <- data.frame(ntax.before=ltt_points[splitrows,2], ntax.after=ltt_points[splitrows+1,2], time=splittimes)
    if(min(result$ntax.before)<minsize) {
      nregimes <- nregimes-1
    } else {
      qualified <- TRUE
    }
  }
  return(result)
}



SplitAndLikelihood <- function(tree, nregimes, minsize=1, type="data", interpolation_method="linear") {
    splits <- EvenSplit(tree=tree, nregimes=nregimes, minsize=minsize, type=type)
    desired_interval = min(0.05, 0.2*min(diff(splits$time)))
    results <- param_lambda_discreteshift_mu_discreteshift(desired_interval = desired_interval, tree=tree, condition="crown", ncores=parallel::detectCores(), slice_ages = unique(sort(c(0, abs(splits$time), castor::get_tree_span(tree)$max_distance))), interpolation_method=interpolation_method)
    return(list(splits=splits, results=results, desired_interval=desired_interval))
}

TryManyRegimes <- function(tree, maxregimes=5) {
    conditions <- expand.grid(nregimes=sequence(maxregimes), interpolation_method=c("linear", "constant"))
    full_results <- list()
    summarized_results <- data.frame()
    for(i in sequence(nrow(conditions))) {
        #print(conditions[i,])
        local_result <- SplitAndLikelihood(tree, nregimes=conditions$nregimes[i], interpolation_method=conditions$interpolation_method[i])
        local_result$nregimes=conditions$nregimes[i]
        local_result$interpolation_method=conditions$interpolation_method[i]
        local.df <- data.frame(nregimes=conditions$nregimes[i], interpolation_method=conditions$interpolation_method[i], AIC=local_result$results$fit_param$AIC, loglikelihood=local_result$results$fit_param$loglikelihood
        ,
            lambda0=local_result$results$fit_param$param_fitted['lambda0'],
            lambda1=local_result$results$fit_param$param_fitted['lambda1'],
            lambda2=local_result$results$fit_param$param_fitted['lambda2'],
            lambda3=local_result$results$fit_param$param_fitted['lambda3'],
            lambda4=local_result$results$fit_param$param_fitted['lambda4'],
            lambda5=local_result$results$fit_param$param_fitted['lambda5'],
            lambda6=local_result$results$fit_param$param_fitted['lambda6'],
            lambda7=local_result$results$fit_param$param_fitted['lambda7'],
            mu0=local_result$results$fit_param$param_fitted['mu0'],
            mu1=local_result$results$fit_param$param_fitted['mu1'],
            mu2=local_result$results$fit_param$param_fitted['mu2'],
            mu3=local_result$results$fit_param$param_fitted['mu3'],
            mu4=local_result$results$fit_param$param_fitted['mu4'],
            mu5=local_result$results$fit_param$param_fitted['mu5'],
            mu6=local_result$results$fit_param$param_fitted['mu6'],
            mu7=local_result$results$fit_param$param_fitted['mu7']
        )
        if(i==1) {
            summarized_results <- local.df
        } else {
            summarized_results <- rbind(summarized_results, local.df)
        }
        print(tail(summarized_results,1))
        full_results[[i]] <-local_result
    }
    return(list(summarized_results=summarized_results, full_results=full_results))
}
