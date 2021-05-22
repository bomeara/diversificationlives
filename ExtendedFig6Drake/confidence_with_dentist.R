source("R/packages.R")
source("R/functions.R")

tree <- ape::read.tree("data/tree_Extended_Data_Fig_6.tre")

#best <- SplitAndLikelihood(tree, nregimes=12, interpolation_method="linear", type="time", Ntrials=1, ncores=1, seed_to_use=42)
#save(best, tree, file="best.rda")
load("best.rda")

params <- best$results$fit_param$param_fitted
names(params) <- names(best$results$fit_param$param_fitted)


best_neglnL <- likelihood_lambda_discreteshift_mu_discreteshift_for_mcmc(params, tree=tree, fitted.model=best, return_neg=TRUE, params_are_log_transformed=FALSE)

param_numbers <- as.numeric(gsub("lambda", "", gsub("mu", "", names(params))))
sd_vector <- 0.02*params*sqrt(1+param_numbers)

all_results <- list()
for (i in sequence(1000)) {
	param_numbers <- as.numeric(gsub("lambda", "", gsub("mu", "", names(params))))
	sd_vector <- runif(1, min=0.001, max=0.05)*params*sqrt(1+param_numbers)
	results <- dent_walk(par=params, fn=likelihood_lambda_discreteshift_mu_discreteshift_for_mcmc, best_neglnL=best_neglnL, nsteps=1000, tree=tree, fitted.model=best, return_neg=TRUE, print_freq=5, debug=TRUE, delta=2, params_are_log_transformed=FALSE, sd_vector=sd_vector, restart_after=15)
	all_results[[i]] <- results
	save(list=ls(), file="best_and_results.rda")
}

load("best_and_results.rda")
source("R/functions.R")
focal_results <- do.call("rbind", lapply(all_results, "[[", "results"))

good_adaptive_samples=focal_results[which(focal_results$neglnL<2+min(focal_results$neglnL)),]

PlotRates(fitted.model=best, tree=tree)
PlotRateUncertainty(fitted.model=best, tree=tree, good_adaptive_samples, oldest_age=-100)

