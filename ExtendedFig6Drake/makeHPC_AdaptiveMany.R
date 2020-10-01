# Doing https://books.ropensci.org/drake/hpc.html#advanced-options to allow parallel within

Sys.setenv('R_MAX_VSIZE'=32000000000)

source("R/packages.R")
source("R/functions.R")
source("R/planHPC_AdaptiveMany.R")


#future::plan(cluster, workers=cl)


make(plan_hpc_adaptivemany)





# now stop with drake given its issues with parallel

workers <- sample(c(rep(c("omearaclustera.local", "omearaclusterb.local", "omearaclusterh.local", "omearaclusterl.local", "omearaclusterk.local"), 2), rep(c("omearatc1.local", "omearatc2.local"),1)))
cl <- makeClusterPSOCK(workers)

parallel::clusterSetRNGStream(cl=cl)

clusterExport(cl, list(ls()))


loadd(subset_list)
loadd(tree)

results <- parallel::parLapplyLB(cl=cl, subset_list, fun=AdaptiveSupport, tree=tree, delta=10) # AdaptiveSupport <- function(fitted.model, tree, delta=2
save(results, file="adaptive.rda")

stopCluster(cl)