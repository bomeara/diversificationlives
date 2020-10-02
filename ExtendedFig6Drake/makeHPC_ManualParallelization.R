# Doing https://books.ropensci.org/drake/hpc.html#advanced-options to allow parallel within

Sys.setenv('R_MAX_VSIZE'=32000000000)

set.seed(as.integer(round(runif(1, min=1, max=1e5)) + sqrt(as.numeric(gsub("\\.", "", as.character(ipify::get_ip()))))))

setwd("/share/diversificationlives/ExtendedFig6Drake")

# ansible linux -a 'nohup Rscript /share/diversificationlives/ExtendedFig6Drake/makeHPC_ManualParallelization.R &'


source("R/packages.R")
source("R/functions.R")





load("InputForAdaptive.rda") #contains subset_list and tree

subset_list <- sample(subset_list, replace=FALSE) #randomize order

results <- list()

for (i in seq_along(subset_list)) {
	results[[i]] <- AdaptiveSupport(subset_list[[i]], tree=tree, delta=10, n_per_rep=parallel::detectCores(), n_per_good=3*parallel::detectCores(), cache_file=paste0("temp_", Sys.info()['nodename'], ".rda"))
	save(results, file=paste0("adaptive_", Sys.info()['nodename'], ".rda"))
}

#rando <- runif(1)
#save(rando, file=paste0("rando_", Sys.info()['nodename'], ".rda"))
#Sys.sleep(60)
