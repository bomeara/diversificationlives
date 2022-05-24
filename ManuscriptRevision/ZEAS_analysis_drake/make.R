# Doing https://books.ropensci.org/drake/hpc.html#advanced-options to allow parallel within


#future::plan(future::multiprocess)

try(drake::drake_cache("/home/bomeara/Documents/localcache")$unlock())


# envir <- new.env(parent = globalenv())
# source("R/packages.R", local = envir)
# source("R/functions.R", local = envir)
# source("R/plan.R", local = envir)
# make(envir$plan, envir = envir, jobs = 12)

source("R/packages.R")  # loads packages
source("R/functions.R")
source("R/plan.R")      # creates the drake plan



#make(plan, cache=drake::new_cache("~/Documents/localcache")) # cache so that they don't all try writing to same cache
results = RunMany(nrep=100, ncore=parallel::detectCores())
write(results, file=paste0(system("hostname", intern=TRUE), "_final_results.rda"))

# ansible linux -a 'nohup Rscript /share/diversificationlives/PlantDrake/make.R &'