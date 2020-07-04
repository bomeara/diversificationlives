# Doing https://books.ropensci.org/drake/hpc.html#advanced-options to allow parallel within


envir <- new.env(parent = globalenv())
source("R/packages.R", local = envir)
source("R/functions.R", local = envir)
source("R/plan.R", local = envir)
make(envir$planmany, envir = envir)
#
# source("R/packages.R")  # loads packages
# source("R/functions.R")
# source("R/plan.R")      # creates the drake plan
#
# make(plan, parallelism = "future", jobs = parallel::detectCores())
