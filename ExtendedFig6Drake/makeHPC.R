# Doing https://books.ropensci.org/drake/hpc.html#advanced-options to allow parallel within


source("R/packages.R")
source("R/functions.R")
source("R/planHPC.R")

workers <- rep(c("omearatc1.local", "omearatc2.local"),24)
cl <- makeClusterPSOCK(workers)

future::plan(cluster, workers=cl)


make(planHPC, parallelism = "future", jobs = 24)
#
# source("R/packages.R")  # loads packages
# source("R/functions.R")
# source("R/plan.R")      # creates the drake plan
#
# make(plan, parallelism = "future", jobs = parallel::detectCores())
