source("R/packages.R")
source("R/functions.R")

tree <- ape::read.tree("data/tree_Extended_Data_Fig_6.tre")

#best <- SplitAndLikelihood(tree, nregimes=12, interpolation_method="linear", type="time", Ntrials=1, ncores=1, seed_to_use=42)
#save(best, tree, file="best.rda")
#load("best.rda")
#tree <- rcoal(500)
#tree  <- geiger::drop.extinct(geiger::sim.bdtree(stop="taxa", n=100, seed=42, extinct=FALSE))
best <- SplitAndLikelihoodPDRdiscreteshift(tree=tree,nregimes=12, seed_to_use=42, ncores=2, type="time")
best_neglnL <- likelihood_pdr_discreteshift_for_mcmc(best$results$fit_param$param_fitted, best, tree=tree, return_neg=TRUE, params_are_log_transformed=FALSE)
save(list=ls(), file="best_pdr.rda")


params <- best$results$fit_param$param_fitted
names(params) <- names(best$results$fit_param$param_fitted)

all_results <- list()
for (i in sequence(1000)) {
	param_numbers <- as.numeric(gsub("lambda", "", gsub("mu", "", names(params))))
	sd_vector <- abs(runif(1, min=0.001, max=0.05)*params)
	results <- dent_walk(par=params, fn=likelihood_pdr_discreteshift_for_mcmc, best_neglnL=best_neglnL, nsteps=runif(1, min=200, max=10000), tree=tree, fitted.model=best, return_neg=TRUE, print_freq=50, debug=TRUE, delta=2, params_are_log_transformed=FALSE, sd_vector=sd_vector, lower_bound=c(rep(-100, length(params)-1), 0), upper_bound=c(rep(100,length(params)-1), 100))
	all_results[[i]] <- results
	save(list=ls(), file="best_pdr_and_results.rda")
	focal_results <- do.call("rbind", lapply(all_results, "[[", "results"))
	good_adaptive_samples=focal_results[which(focal_results$neglnL<2+min(focal_results$neglnL)),]
	pdf(file="~/Dropbox/pdr_dentist.pdf", width=7, height=7)
	PlotRateUncertaintyPDR(fitted.model=best, tree=tree, good_adaptive_samples, oldest_age=-100)
	dev.off()
}

# load("best_and_results.rda")
# source("R/functions.R")
#focal_results <- do.call("rbind", lapply(all_results, "[[", "results"))

#good_adaptive_samples=focal_results[which(focal_results$neglnL<2+min(focal_results$neglnL)),]

# PlotRates(fitted.model=best, tree=tree)
# PlotRateUncertainty(fitted.model=best, tree=tree, good_adaptive_samples, oldest_age=-100)

