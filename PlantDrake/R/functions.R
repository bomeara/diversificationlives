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

SplitAndLikelihood <- function(tree, nregimes, minsize=1, type="data", interpolation_method="linear") {
    splits <- EvenSplit(tree=tree, nregimes=nregimes, minsize=minsize, type=type)

    results <- param_lambda_discreteshift_mu_discreteshift(desired_interval = 0.1, tree=tree, condition="crown", ncores=parallel::detectCores(), slice_ages = sort(c(0, abs(splits$time), castor::get_tree_span(tree)$max_distance)), interpolation_method=interpolation_method)
    return(list(splits=splits, results=results))
}
