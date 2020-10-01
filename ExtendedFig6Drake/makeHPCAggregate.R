# Doing https://books.ropensci.org/drake/hpc.html#advanced-options to allow parallel within

Sys.setenv('R_MAX_VSIZE'=32000000000)

source("R/packages.R")
source("R/functions.R")
source("R/planHPCAggregate.R")

#workers <- c(rep(c("omearaclustera.local", "omearaclusterb.local", "omearaclusterc.local"), 24), rep(c("omearatc1.local", "omearatc2.local"),12))
#cl <- makeClusterPSOCK(workers)

#future::plan(cluster, workers=cl)


make(plan_hpc_aggregate)
#
# source("R/packages.R")  # loads packages
# source("R/functions.R")
# source("R/plan.R")      # creates the drake plan
#
# make(plan, parallelism = "future", jobs = parallel::detectCores())
