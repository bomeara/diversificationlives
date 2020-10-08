# Doing https://books.ropensci.org/drake/hpc.html#advanced-options to allow parallel within


future::plan(future::multiprocess)

setwd("/share/diversificationlives/PlantDrake")

# envir <- new.env(parent = globalenv())
# source("R/packages.R", local = envir)
# source("R/functions.R", local = envir)
# source("R/plan.R", local = envir)
# make(envir$plan, envir = envir, jobs = 12)

source("R/packages.R")  # loads packages
source("R/functions.R")
source("R/plan.R")      # creates the drake plan

try(drake::drake_cache("/share/diversificationlives/PlantDrake/.drake")$unlock())


make(plan, parallelism = "future", jobs = parallel::detectCores(), cache=drake::new_cache("~/Documents/localcache"), keep_going=TRUE) # cache so that they don't all try writing to same cache


# ansible linux -a 'nohup Rscript /share/diversificationlives/PlantDrake/make.R &'