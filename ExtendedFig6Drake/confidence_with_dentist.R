source("R/packages.R")
source("R/functions.R")

tree <- ape::read.tree("data/tree_Extended_Data_Fig_6.tre")

#best <- SplitAndLikelihood(tree, nregimes=12, interpolation_method="linear", type="time", Ntrials=1, ncores=1, seed_to_use=42)
#save(best, tree, file="best.rda")
load("best.rda")

params <- best$results$fit_param$param_fitted
names(params) <- names(best$results$fit_param$param_fitted)
results <- dent_walk(par=params, fn=likelihood_lambda_discreteshift_mu_discreteshift_for_mcmc, best_neglnL=-1*best$loglikelihood, nsteps=1000, tree=tree, fitted.model=best, return_neg=TRUE, print_freq=5, adjust_width_interval=20)
