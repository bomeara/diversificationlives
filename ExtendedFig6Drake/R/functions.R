
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
                                               fit_control       = list(rel.tol=1e-8)
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
                                               fit_control       = list(rel.tol=1e-8)
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
                                               fit_control       = list(rel.tol=1e-8)
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
                                               fit_control       = list(rel.tol=1e-8)
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



SplitAndLikelihood <- function(tree, nregimes, minsize=1, type="data", interpolation_method="linear", verbose=TRUE) {
    splits <- EvenSplit(tree=tree, nregimes=nregimes, minsize=minsize, type=type)
    desired_interval = min(0.05, 0.2*min(abs(diff(splits$time))))
    results <- param_lambda_discreteshift_mu_discreteshift(desired_interval = desired_interval, tree=tree, condition="crown", ncores=parallel::detectCores(), slice_ages = unique(sort(c(0, abs(splits$time), castor::get_tree_span(tree)$max_distance))), interpolation_method=interpolation_method)
    if(verbose) {
        print(c(nregimes=nregimes, interpolation_method=interpolation_method, AIC=results$fit_param$AIC))
    }
    return(list(splits=splits, results=results, desired_interval=desired_interval, nregimes=nregimes, interpolation_method=interpolation_method, type=type, AIC=results$fit_param$AIC, loglikelihood=results$fit_param$loglikelihood))
}

SummarizeSplitsAndLikelihoods <- function(x) {
    summary.df <- data.frame(nregimes=sapply(x, "[[", "nregimes"), interpolation_method=sapply(x, "[[", "interpolation_method"), AIC=sapply(x, "[[", "AIC"), loglikelihood=sapply(x, "[[", "loglikelihood"), stringsAsFactors=FALSE)
    summary.df$deltaAIC <- summary.df$AIC-min(summary.df$AIC)
    return(summary.df)
}

AdaptiveSampleBestModels <- function(everything, result_summary, tree, deltaAIC_cutoff=20) {
    adaptive_results <- vector(mode="list", length=nrow(result_summary))
    for (i in seq_along(everything)) {
        if(result_summary$deltaAIC[i]<deltaAIC_cutoff) {
            adaptive_results[[i]] <- AdaptiveSupport(fitted.model=everything[[i]], tree=tree)
        }
    }
    return(adaptive_results)
}
# TryManyRegimes <- function(tree, maxregimes=5) {
#     conditions <- expand.grid(nregimes=sequence(maxregimes), interpolation_method=c("linear", "constant"))
#     full_results <- list()
#     summarized_results <- data.frame()
#     for(i in sequence(nrow(conditions))) {
#         #print(conditions[i,])
#         local_result <- SplitAndLikelihood(tree, nregimes=conditions$nregimes[i], interpolation_method=conditions$interpolation_method[i])
#         local_result$nregimes=conditions$nregimes[i]
#         local_result$interpolation_method=conditions$interpolation_method[i]
#         local.df <- data.frame(nregimes=conditions$nregimes[i], interpolation_method=conditions$interpolation_method[i], AIC=local_result$results$fit_param$AIC, loglikelihood=local_result$results$fit_param$loglikelihood
#         ,
#             lambda0=local_result$results$fit_param$param_fitted['lambda0'],
#             lambda1=local_result$results$fit_param$param_fitted['lambda1'],
#             lambda2=local_result$results$fit_param$param_fitted['lambda2'],
#             lambda3=local_result$results$fit_param$param_fitted['lambda3'],
#             lambda4=local_result$results$fit_param$param_fitted['lambda4'],
#             lambda5=local_result$results$fit_param$param_fitted['lambda5'],
#             lambda6=local_result$results$fit_param$param_fitted['lambda6'],
#             lambda7=local_result$results$fit_param$param_fitted['lambda7'],
#             mu0=local_result$results$fit_param$param_fitted['mu0'],
#             mu1=local_result$results$fit_param$param_fitted['mu1'],
#             mu2=local_result$results$fit_param$param_fitted['mu2'],
#             mu3=local_result$results$fit_param$param_fitted['mu3'],
#             mu4=local_result$results$fit_param$param_fitted['mu4'],
#             mu5=local_result$results$fit_param$param_fitted['mu5'],
#             mu6=local_result$results$fit_param$param_fitted['mu6'],
#             mu7=local_result$results$fit_param$param_fitted['mu7']
#         )
#         if(i==1) {
#             summarized_results <- local.df
#         } else {
#             summarized_results <- rbind(summarized_results, local.df)
#         }
#         print(tail(summarized_results,1))
#         full_results[[i]] <-local_result
#     }
#     return(list(summarized_results=summarized_results, full_results=full_results))
# }

ComputeRates <- function(fitted.model, params=NULL) {
    ages <- fitted.model$results$age_grid_param
    if(is.null(params)) {
        params <- fitted.model$results$fit_param$param_fitted
    }
    if(is.null(names(params))) {
        names(params) <- names(fitted.model$results$fit_param$param_fitted)
    }
    lambdas <- fitted.model$results$lambda_function(ages=ages, params=params)
    mus <- fitted.model$results$mu_function(ages=ages, params=params)
    return(data.frame(age=-ages, lambda=lambdas, mu=mus, netdiv = lambdas-mus, turnover=lambdas+mus, ef=mus/lambdas))
}

PlotRates <- function(fitted.model, tree,...) {
    rates <- ComputeRates(fitted.model)
    ape::ltt.plot(tree,...)
    relrates <- rates
    relrates$lambda <- ape::Ntip(tree) * relrates$lambda/diff(range(c(rates$lambda, rates$mu, 0)))
    relrates$mu <- ape::Ntip(tree) * relrates$mu/diff(range(c(rates$lambda, rates$mu, 0)))
    lines(relrates$age, relrates$lambda, col="blue", lwd=2)
    lines(relrates$age, relrates$mu, col="red", lwd=2)
    axis(side=4, at=seq(from=0, to=ape::Ntip(tree), length.out=5), labels=signif(seq(from=min(0, min(c(rates$lambda, rates$mu))), to=max(c(rates$lambda, rates$mu)), length.out=5),2))
}

PlotRateUncertainty <- function(fitted.model, tree, good_adaptive_samples, ...) {
    lambda <- data.frame(matrix(NA, nrow=length(fitted.model$results$age_grid_param), ncol=nrow(good_adaptive_samples)))
    mu <- lambda
    netdiv <- lambda
    ages <- (-1)*fitted.model$results$age_grid_param

    for (i in sequence(nrow(good_adaptive_samples))) {
        rates_local <- ComputeRates(fitted.model, good_adaptive_samples[i,-1])
        lambda[,i] <- rates_local$lambda
        mu[,i] <- rates_local$mu
        netdiv[,i] <- rates_local$netdiv
    }
    # ltt_data <- ape::ltt.plot.coords(tree)
    # ltt_data$logN <- log(ltt_data$N)
    par(mfcol=c(1,3))
    ylimits <- range(c(min(max(mu), max(lambda),0.5), max(min(mu), min(lambda), -0.5)))
    for (i in sequence(3)) {
        rates <- lambda
        title <- "Speciation rate"
        if(i==2) {
            rates <- mu
            title <- "Extinction rate"
        }
        if(i==3) {
            rates <- netdiv
            title <- "Net diversification rate"
            ylimits <- range(c(min(max(netdiv),0.5), max(min(netdiv), -0.5)))
        }
        plot(x=ages, y=rates[,1], type="n", bty="n", ylim=ylimits, ylab=title, xlab="Time")
        polygon(x=c(ages, rev(ages)), y=c(apply(rates,1,min), rev(apply(rates,1,max))), col="gray", border=NA)
        # for(j in sequence(ncol(rates))) {
        #     lines(ages, rates[,j], col=rgb(0,0,0,0.1))
        # }
        lines(ages, rates[,1], col=c("blue", "red", "purple")[i], lwd=2)
        axis(side=3, at=c(0, fitted.model$splits$time, min(ages)), labels=c(ape::Ntip(tree), fitted.model$splits$ntax.after, 2))
        mtext("Number of taxa", side=3, line=3, cex=0.6)
    }
}

PlotAllUncertainty <- function(x, tree, adaptive_list, file="uncertainty.pdf", desired_delta=2) {
    pdf(file=file, width=10, height=5)
    for(i in seq_along(adaptive_list)) {
        if(class(adaptive_list[[i]])=="data.frame") {
            good_enough <- subset(adaptive_list[[i]], loglikelihood>max(loglikelihood)-desired_delta)
            fitted.model <- x[[i]]
            PlotRateUncertainty(fitted.model, tree, good_enough)
        }
    }
    dev.off()
}

PlotAll <- function(x, tree, file="plot.pdf") {
    summaries <- SummarizeSplitsAndLikelihoods(x)
    pdf(file=file)
    for (i in seq_along(x)) {
        PlotRates(x[[i]], tree, main=paste0("Regimes: ", summaries$nregimes[i], " deltaAICc: ", round(summaries$deltaAIC[i],2)))
    }
    dev.off()
}

# delta is desired ∆lnL to sample along
AdaptiveSupport <- function(fitted.model, tree, delta=2, n_per_rep=12, n_per_good=36) {
    original_params <- fitted.model$results$fit_param$param_fitted
    best_loglikelihood <- fitted.model$loglikelihood
    results <- likelihood_lambda_discreteshift_mu_discreteshift(fitted.model=fitted.model, tree=tree, randomize=FALSE)

    # univariate
    for(focal_param in seq_along(original_params)) {
        multiplier <- -2
        #local_result_aggregate <- data.frame()
        good_sample <- FALSE
        run_num <- 0
        run_max <- 10
        while(!good_sample & run_num < run_max) {
            print(paste(names(original_params)[focal_param], multiplier))
            run_num <- run_num+1
            min_value <- (1-10^multiplier)*original_params[focal_param]
            max_value <- (1+(10^(.5*multiplier)))*original_params[focal_param]
            if(max_value==0) {
                max_value <- (10^(.5*multiplier))*1e-2
            }
            if(multiplier>=0) {
                min_value <- 0
            }
            param_min <- original_params
            param_min[focal_param] <- min_value
            param_max <- original_params
            param_max[focal_param] <- max_value

            # repeated below
            local_results <- do.call(rbind, parallel::mclapply(rep(list(fitted.model),n_per_rep), likelihood_lambda_discreteshift_mu_discreteshift, param_min=param_min, param_max=param_max, tree=tree, randomize=TRUE, mc.cores=parallel::detectCores()))
            results <- rbind(results, local_results)

            if(best_loglikelihood > max(unlist(local_results[,'loglikelihood']))+delta) {
                print("too wide")
                multiplier <- multiplier-.5
            } else if(best_loglikelihood < min(unlist(local_results[,'loglikelihood']))+delta) {
                print("too narrow")
                multiplier <- multiplier+0.5

            } else {
                good_sample <- TRUE
                local_results_good <- do.call(rbind, parallel::mclapply(rep(list(fitted.model),n_per_good), likelihood_lambda_discreteshift_mu_discreteshift, param_min=param_min, param_max=param_max, tree=tree, randomize=TRUE, mc.cores=parallel::detectCores()))
                results <- rbind(results, local_results_good)
            }

            #print(local_results[,1]-best_loglikelihood)
            #plot(local_results[,1+focal_param], local_results[,1]-best_loglikelihood, pch=21, main=colnames(local_results)[1+focal_param])
            #print(cbind(local_results[,1]-best_loglikelihood, local_results[,1+focal_param]))
            #print(range(local_results[,1+focal_param]))

        }
    }

    # bivariate
    indices <- sort(as.numeric(unique(gsub("mu", "", gsub("lambda", "", names(original_params))))))
    good_enough_univariate <- subset(results, loglikelihood+delta>=best_loglikelihood)
    for(focal_pair in seq_along(indices)) {
        print(paste0("focal_pair ", focal_pair))
        param_min <- original_params
        param_max <- original_params
        multiplier <- -2
        good_sample <- FALSE
        run_num <- 0
        run_max <- 10
        while(!good_sample & run_num < run_max) {
            print(multiplier)
            run_num <- run_num+1
            param_min[paste0("mu",indices[focal_pair])] <- max(0,(1-10^multiplier))*min(good_enough_univariate[paste0("mu",indices[focal_pair])])
            param_min[paste0("lambda",indices[focal_pair])] <- max(0,(1-10^multiplier))*min(good_enough_univariate[paste0("lambda",indices[focal_pair])])


            param_max[paste0("mu",indices[focal_pair])] <- (1+(10^(.5*multiplier)))*max(good_enough_univariate[paste0("mu",indices[focal_pair])])
            param_max[paste0("lambda",indices[focal_pair])] <- (1+(10^(.5*multiplier)))*max(good_enough_univariate[paste0("lambda",indices[focal_pair])])

            local_results <- do.call(rbind, parallel::mclapply(rep(list(fitted.model),2*n_per_rep), likelihood_lambda_discreteshift_mu_discreteshift, param_min=param_min, param_max=param_max, tree=tree, randomize=TRUE, mc.cores=parallel::detectCores()))
            results <- rbind(results, local_results)
            if(best_loglikelihood > max(unlist(local_results[,'loglikelihood']))+delta) {
                print("too wide")
                multiplier <- multiplier-.5
            } else if(best_loglikelihood < min(unlist(local_results[,'loglikelihood']))+delta) {
                print("too narrow")
                multiplier <- multiplier+0.5

            } else {
                good_sample <- TRUE
                local_results_good <- do.call(rbind, parallel::mclapply(rep(list(fitted.model),n_per_good), likelihood_lambda_discreteshift_mu_discreteshift, param_min=param_min, param_max=param_max, tree=tree, randomize=TRUE, mc.cores=parallel::detectCores()))
                results <- rbind(results, local_results_good)
            }

        }

    }

    # multivariate

    # multivariate using mcmc
    print("starting mcmc")

    good_enough_already <- subset(results, loglikelihood+delta>=best_loglikelihood)

    print(paste("already have", nrow(good_enough_already), "results in close enough range in ", nrow(results), "trials of parameters"))
    good_enough_ranges <- abs(apply(good_enough_already[,-1], 2, max) - apply(good_enough_already[,-1], 2, min))
    if(any(good_enough_ranges==0)) {
        good_enough_ranges[which(good_enough_ranges==0)] <- 1e-8
    }
    w <- 0.1*abs(log(good_enough_ranges))
    print("w")
    print(w)

    # lik <- function(x, fitted.model, tree) {
    #
    #         return(likelihood_lambda_discreteshift_mu_discreteshift_for_mcmc(x,fitted.model=fitted.model, tree=tree))
    #
    # }
    original_params[original_params==0] <- 1e-8
    #save(list=ls(), file="before_mcmc.rda")
    mcmc_results <- sample_ridge(obj=likelihood_lambda_discreteshift_mu_discreteshift_for_mcmc, initial=log(original_params), nsteps=100, scale=w, fitted.model=fitted.model, tree=tree)
    #mcmc_results <- mcmc::metrop(obj=likelihood_lambda_discreteshift_mu_discreteshift_for_mcmc, initial=log(original_params), nbatch=10000, blen=1, scale=w, fitted.model=fitted.model, tree=tree, debug=TRUE)
    #likelihoods <- apply(mcmc_results$batch, 1, likelihood_lambda_discreteshift_mu_discreteshift_for_mcmc, fitted.model=fitted.model, tree=tree)
    mcmc_params <- exp(mcmc_results$parameters)
    colnames(mcmc_params) <- names(original_params)
    mcmc_results_fixed = data.frame(loglikelihood=mcmc_results$loglikelihoods, mcmc_params)

    results <- rbind(results, mcmc_results_fixed)

    # mcmc_results <- MCMCpack::MCMCmetrop1R(fun=likelihood_lambda_discreteshift_mu_discreteshift_for_mcmc, theta.init=log(original_params), mcmc=10, burnin=2, verbose=TRUE, fitted.model=fitted.model, tree=tree, force.samp=TRUE)

    # mcmc_results <- diversitree::mcmc(lik=likelihood_lambda_discreteshift_mu_discreteshift_for_mcmc, x.init=log(original_params), w=w, nsteps=1, fitted.model=fitted.model, tree=tree)


    # print("beginning multivariate")
    # good_enough_already <- subset(results, loglikelihood+delta>=best_loglikelihood)

    #
    # param_min_original <- apply(good_enough_already, 2, min)[-1]
    # param_max_original <- apply(good_enough_already, 2, max)[-1]
    # multiplier <- -6
    # good_sample <- FALSE
    # run_num <- 0
    # run_max <- 10
    # while(!good_sample & run_num < run_max) {
    #     print(multiplier)
    #     run_num <- run_num+1
    #     param_min <- max(0,(1-10^multiplier))*param_min_original
    #     param_max <- (1+(10^(.5*multiplier)))*param_max_original
    #
    #     local_results <- do.call(rbind, parallel::mclapply(rep(list(fitted.model),8*n_per_rep), likelihood_lambda_discreteshift_mu_discreteshift, param_min=param_min, param_max=param_max, tree=tree, randomize=TRUE, mc.cores=parallel::detectCores()))
    #     results <- rbind(results, local_results)
    #     if(best_loglikelihood > max(unlist(local_results[,'loglikelihood']))+delta) {
    #         print("too wide")
    #         multiplier <- multiplier-.25
    #     } else if(best_loglikelihood < min(unlist(local_results[,'loglikelihood']))+delta) {
    #         print("too narrow")
    #         multiplier <- multiplier+0.5
    #
    #     } else {
    #         good_sample <- TRUE
    #         local_results_good <- do.call(rbind, parallel::mclapply(rep(list(fitted.model),8*n_per_good), likelihood_lambda_discreteshift_mu_discreteshift, param_min=param_min, param_max=param_max, tree=tree, randomize=TRUE, mc.cores=parallel::detectCores()))
    #         results <- rbind(results, local_results_good)
    #     }
    #
    # }


    return(results)
}


# desired_interval: how finely to divide the slices for calculating smoothly changing (or not) lambda and mu
# slice_ages: when to divide regimes
# interpolation_method: same as method in ?approx. Constant or linear
likelihood_lambda_discreteshift_mu_discreteshift <- function(fitted.model, tree, param_min=0*fitted.model$results$fit_param$param_fitted, param_max=4*(0.05+fitted.model$results$fit_param$param_fitted), randomize=TRUE) {
    original_params <- fitted.model$results$fit_param$param_fitted
    seed_params <- original_params
    # if any are exactly zero, won't get variation. So allow change
    new_params <- original_params
    if(randomize) {
        new_params <- runif(n=length(original_params), min=param_min, max=param_max)
        names(new_params) <- names(original_params)
    }
    mu_values <- fitted.model$results$mu_function(fitted.model$results$age_grid_param, new_params)

    lambda_values <- fitted.model$results$lambda_function(fitted.model$results$age_grid_param, new_params)

    loglikelihood_result <- castor::loglikelihood_hbd(
        tree=tree,
        age_grid = fitted.model$results$age_grid_param,
        lambda=lambda_values,
        mu=mu_values,
        rho0=1,
        splines_degree=1
    )
    result <- c(loglikelihood=ifelse(is.finite(loglikelihood_result$loglikelihood), loglikelihood_result$loglikelihood, -1e12), new_params)
    return(data.frame(t(result)))
}


likelihood_lambda_discreteshift_mu_discreteshift_for_mcmc <- function(log_params, fitted.model, tree) {
    params <- exp(log_params)
    names(params) <- names(fitted.model$results$fit_param$param_fitted)
    #print(params)
    mu_values <- fitted.model$results$mu_function(fitted.model$results$age_grid_param, params)

    lambda_values <- fitted.model$results$lambda_function(fitted.model$results$age_grid_param, params)

    loglikelihood_result <- castor::loglikelihood_hbd(
        tree=tree,
        age_grid = fitted.model$results$age_grid_param,
        lambda=lambda_values,
        mu=mu_values,
        rho0=1,
        splines_degree=1
    )
    #print(loglikelihood_result$loglikelihood)
    return(ifelse(is.finite(loglikelihood_result$loglikelihood), loglikelihood_result$loglikelihood, -Inf))
}

#(obj=likelihood_lambda_discreteshift_mu_discreteshift_for_mcmc, initial=log(original_params), nbatch=10000, blen=1, scale=w, fitted.model=fitted.model, tree=tree, debug=TRUE)
#this assumes obj gives back the loglikelihood, NOT the negative log likelihood. Bigger is better, the MLE should be the max value for likelihood.
sample_ridge <- function(obj, initial,  scale, nsteps=1000, desired_delta=2, ...) {
    #save(list=ls(), file="within_sample_ridge.rda")
    loglikelihoods <- rep(NA, nsteps+1)
    parameters <- data.frame(matrix(NA, nrow=nsteps+1, ncol=length(initial)))
    parameters[1,] <- initial
    colnames(parameters) <- names(initial)
    loglikelihoods[1] <- obj(parameters[1,], ...)
    generating_params <- initial
    target_loglikelihood <- loglikelihoods[1]-desired_delta
    print(paste("target_loglikelihood", target_loglikelihood))
    for (i in sequence(nsteps)) {
        new_params <- generating_params
        param_to_tweak <- sample.int(length(initial), 1)
        new_params[param_to_tweak] <- stats::rnorm(1, mean=new_params[param_to_tweak], sd=scale[param_to_tweak])
        loglikelihoods[i+1] <- obj(new_params, ...)
        parameters[i+1,] <- new_params
        loglikdiff <- target_loglikelihood-loglikelihoods[i+1]
        # if 0, we're hitting exactly; if <0, the new params are "too good", if >0, the new params are too far away from a good region
        if(abs(loglikdiff)<runif(1,0,1)) {
            generating_params <- new_params # we're within ∆1 lnL of the desired value, so pretty good chance we want to accept this change, esp as we get closer
        }
        if(loglikdiff<0) {
            scale[param_to_tweak] <- 1.1 * scale[param_to_tweak] # get a bit more variance in proposals
        }
        if(loglikdiff>0) {
            scale[param_to_tweak] <- 0.9 * scale[param_to_tweak] # let's not be so bold
        }
        if(loglikediff>10) {
            scale[param_to_tweak] <- 0.1 * scale[param_to_tweak] # we're very far away from where we want to be
        }
        print(paste("mcmc step", i, "difference from targeted likelihood is", loglikdiff))
    }
    return(list(loglikelihood=loglikelihood, parameters=parameters))
}
