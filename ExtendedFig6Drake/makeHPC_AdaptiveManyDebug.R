# Doing https://books.ropensci.org/drake/hpc.html#advanced-options to allow parallel within

Sys.setenv('R_MAX_VSIZE'=32000000000)

source("R/packages.R")
source("R/functions.R")
source("R/planHPC_AdaptiveMany.R")


#future::plan(cluster, workers=cl)


#make(plan_hpc_adaptivemany)





# now stop with drake given its issues with parallel

workers <- sample(c(rep(c("omearaclustera.local", "omearaclusterb.local", "omearaclusterh.local", "omearaclusterl.local", "omearaclusterk.local"), 2), rep(c("omearatc1.local", "omearatc2.local"),1)))

print("making cluster")
cl <- makeClusterPSOCK(workers)

print("setting stream")
parallel::clusterSetRNGStream(cl=cl)

loadd(subset_list)
loadd(tree)

print("trees and subset loaded")

GetAIC <- function(x) {
	return(x$AIC)
}

gethostname <- function(...) {
	return(system("hostname", intern=TRUE))
}

sl <-subset_list[1:2]

print("subsetted")

clusterExport(cl, list(ls()))


#results <- parallel::parLapplyLB(cl=cl, sl, fun=GetAIC) # AdaptiveSupport <- function(fitted.model, tree, delta=2
#results <- parallel::parLapplyLB(cl=cl, sl, fun=AdaptiveSupport, tree=tree, delta=10, n_per_rep=1, n_per_good=1) # AdaptiveSupport <- function(fitted.model, tree, delta=2

results <- parallel::parLapplyLB(cl=cl, AdaptiveSupport, fun=gethostname) # AdaptiveSupport <- function(fitted.model, tree, delta=2


save(results, file="length.rda")

stopCluster(cl)