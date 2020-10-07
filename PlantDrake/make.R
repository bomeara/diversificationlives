# Doing https://books.ropensci.org/drake/hpc.html#advanced-options to allow parallel within


future::plan(future::multiprocess)

setwd("/share/diversificationlives/PlantDrake")
try(drake::drake_cache("/share/diversificationlives/PlantDrake/.drake")$unlock())

# envir <- new.env(parent = globalenv())
# source("R/packages.R", local = envir)
# source("R/functions.R", local = envir)
# source("R/plan.R", local = envir)
# make(envir$plan, envir = envir, jobs = 12)

source("R/packages.R")  # loads packages
source("R/functions.R")
source("R/plan.R")      # creates the drake plan

make(plan, parallelism = "future", jobs = parallel::detectCores())


# ansible linux -a 'nohup Rscript /share/diversificationlives/PlantDrake/make.R &'