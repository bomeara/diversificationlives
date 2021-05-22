source("R/packages.R")
source("R/functions.R")

tree <- ape::read.tree("data/tree_Extended_Data_Fig_6.tre")

#best <- SplitAndLikelihood(tree, nregimes=12, interpolation_method="linear", type="time", Ntrials=1, ncores=1, seed_to_use=42)
#save(best, tree, file="best.rda")
#load("best.rda")
#tree <- rcoal(500)
#tree  <- geiger::drop.extinct(geiger::sim.bdtree(stop="taxa", n=100, seed=42, extinct=FALSE))
best <- SplitAndLikelihoodPDRdiscreteshift(tree=tree,nregimes=10, seed_to_use=42, ncores=2)
best_neglnL <- likelihood_pdr_discreteshift_for_mcmc(best$results$fit_param$param_fitted, best, tree=tree, return_neg=TRUE, params_are_log_transformed=FALSE)
# tree <- rcoal(500)
# nregimes=10
# minsize=1
# type="time"
# interpolation_method="linear"
# verbose=TRUE
# Ntrials=3
# ncores=1
# instance=1
# seed_to_use=42
# set.seed(seed_to_use)
# iterate_on_seed <- runif(seed_to_use)
# rm(iterate_on_seed)
# splits <- EvenSplit(tree=tree, nregimes=nregimes, minsize=minsize, type=type)
# desired_interval = min(0.05, 0.2*min(abs(diff(splits$time))))
# results <- NA
# condition="crown"
# slice_ages = unique(sort(c(0, abs(splits$time), castor::get_tree_span(tree)$max_distance)))
# root_age = castor::get_tree_span(tree)$max_distance
# rho = 1
# age_grid_param = seq(from=0,to=root_age+desired_interval,by=desired_interval)

# pdr_params <- rep(NA, length(slice_ages))
# names(pdr_params) <- paste0("pdr", sequence(length(slice_ages))-1)


# rholambda0_function = function(params){
# 		return(params['rholambda0'])
# }
# PDR_function = function(ages,params){
# results <- stats::approx(x=slice_ages, y=params[grepl("pdr", names(params))], xout=ages, method=interpolation_method, rule=2)$y
# return(results)
# }
# param_values  = c(pdr_params, rholambda0=NA)

# ape_estimate <- ape::birthdeath(ape::multi2di(tree))

# netdiv_range <- unname(ape_estimate$CI['b-d',])
# #birth_range <- abs(range(c(netdiv_range,rev(netdiv_range)) / (1-ef_range)))
# # death_range <- abs(range(ef_range * c(netdiv_range, rev(netdiv_range)) / (1-ef_range)))


# param_guess <- c(runif(n=length(pdr_params), min=min(netdiv_range), max=max(netdiv_range)),.9)

# names(param_guess) <- names(param_values)

# fit_param = NA


# fit_param <- castor::fit_hbd_pdr_parametric(	tree,
# 	param_values  = param_values,
# 	param_guess   = param_guess,
# 	param_min     = c(rep(-10,length(pdr_params)),0),
# 	param_max     = c(rep(10,length(pdr_params)),1),
# 	param_scale   = 1, # all params are in the order of 1
# 	PDR = PDR_function, 
# 	rholambda0 = rholambda0_function,
# 	age_grid      = age_grid_param,
# 	condition     = condition,
# 	Ntrials       = Ntrials,    # perform 10 fitting trials
# 	Nthreads      = ncores,
# 	fit_control       = list(rel.tol=1e-8, trace=1)
# )

params <- best$results$fit_param$param_fitted
names(params) <- names(best$results$fit_param$param_fitted)


#best_neglnL <- likelihood_lambda_discreteshift_mu_discreteshift_for_mcmc(params, tree=tree, fitted.model=best, return_neg=TRUE, params_are_log_transformed=FALSE)

#param_numbers <- as.numeric(gsub("pdr", "", gsub("mu", "", names(params))))
#sd_vector <- 0.02*params*sqrt(1+param_numbers)
sd_vector <- 0.02*abs(params)

likelihood_pdr_discreteshift_for_mcmc(dentist:::dent_pr, best, tree=tree, return_neg=TRUE, params_are_log_transformed=FALSE)
results <- dent_walk(par=params, fn=likelihood_pdr_discreteshift_for_mcmc, best_neglnL=best_neglnL, nsteps=100, tree=tree, fitted.model=best, return_neg=TRUE, print_freq=5, debug=TRUE, delta=2, params_are_log_transformed=FALSE, sd_vector=sd_vector, lower_bound=c(rep(-100, length(params)-1), 0), upper_bound=c(rep(100,length(params)-1), 100))

# all_results <- list()
# for (i in sequence(1000)) {
# 	param_numbers <- as.numeric(gsub("lambda", "", gsub("mu", "", names(params))))
# 	sd_vector <- runif(1, min=0.001, max=0.05)*params*sqrt(1+param_numbers)
# 	results <- dent_walk(par=params, fn=likelihood_lambda_discreteshift_mu_discreteshift_for_mcmc, best_neglnL=best_neglnL, nsteps=1000, tree=tree, fitted.model=best, return_neg=TRUE, print_freq=5, debug=TRUE, delta=2, params_are_log_transformed=FALSE, sd_vector=sd_vector, restart_after=15)
# 	all_results[[i]] <- results
# 	save(list=ls(), file="best_and_results.rda")
# }

# load("best_and_results.rda")
# source("R/functions.R")
# focal_results <- do.call("rbind", lapply(all_results, "[[", "results"))

# good_adaptive_samples=focal_results[which(focal_results$neglnL<2+min(focal_results$neglnL)),]

# PlotRates(fitted.model=best, tree=tree)
# PlotRateUncertainty(fitted.model=best, tree=tree, good_adaptive_samples, oldest_age=-100)

