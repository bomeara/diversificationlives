# Doing https://books.ropensci.org/drake/hpc.html#advanced-options to allow parallel within

#
# envir <- new.env(parent = globalenv())
# source("R/packages.R", local = envir)
# source("R/functions.R", local = envir)
# source("R/plan.R", local = envir)
# make(envir$plan_original, envir = envir, parallelism = "future", jobs = parallel::detectCores())



source("R/packages.R")  # loads packages
source("R/functions.R")
source("R/plan.R")      # creates the drake plan
drake::drake_cache("/home/bomeara/Documents/diversificationlives/ExtendedFig6Drake/.drake")$unlock()
drake::clean()
options(clustermq.scheduler = "multicore")

make(plan_manystart, parallelism = "clustermq", jobs = parallel::detectCores())
