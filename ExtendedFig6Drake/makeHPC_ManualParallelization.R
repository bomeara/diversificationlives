# Doing https://books.ropensci.org/drake/hpc.html#advanced-options to allow parallel within

Sys.setenv('R_MAX_VSIZE'=32000000000)

set.seed(as.integer(sqrt(as.numeric(gsub("\\.", "", as.character(ipify::get_ip()))))))

source("R/packages.R")
source("R/functions.R")


#future::plan(cluster, workers=cl)







# now stop with drake given its issues with parallel

#workers <- sample(c(rep(c("omearaclustera.local", "omearaclusterb.local", "omearaclusterh.local", "omearaclusterl.local", "omearaclusterk.local"), 2), rep(c("omearatc1.local", "omearatc2.local"),1)))
#cl <- makeClusterPSOCK(workers)

#parallel::clusterSetRNGStream(cl=cl)

#clusterExport(cl, list(ls()))


load("InputForAdaptive") #contains subset_list and tree

results <- parallel::parLapplyLB(cl=cl, AdaptiveSupport, fun=gethostname) # AdaptiveSupport <- function(fitted.model, tree, delta=2

subset_list <- sample(subset_list, replace=FALSE) #randomize order

results <- list()

for (i in seq_along(subset_list)) {
	results[[i]] <- AdaptiveSupport(subset_list[[i]], tree=tree, delta=10, n_per_rep=parallel::detectCores(), n_per_good=3*parallel::detectCores(), cache_file=paste0("temp_", Sys.info()['nodename'], ".rda"))
	save(results, file=paste0("adaptive_", Sys.info()['nodename'], ".rda"))
}
