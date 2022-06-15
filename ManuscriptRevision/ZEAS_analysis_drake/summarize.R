library(hisse)


all_results_model_average <- data.frame()
all_results_best <- data.frame()
all_results_overview <- data.frame()

inputs <- list.files(pattern="omeara.*results.rda", full.names=TRUE)

for (i in sequence(length(inputs))) {
  load(inputs[i])
  for(j in sequence(length(results))) {
	local <- hisse::SummarizeMiSSEGreedy(results[[j]]$greedy)
	local_model_average <- as.data.frame(local$rates[,,"model_average"])
	local_model_average$input_file <- inputs[i]
	local_model_average$rep_number <- j
	local_model_average$Ntip <- ape::Ntip(results[[j]]$phy)
	all_results_model_average <- rbind(all_results_model_average, local_model_average)
	
	local_model_best <- as.data.frame(local$rates[,,"best"])
	local_model_best$input_file <- inputs[i]
	local_model_best$rep_number <- j
	local_model_best$Ntip <- ape::Ntip(results[[j]]$phy)

	all_results_best <- rbind(all_results_best, local_model_best)
	
	local_overview <- local$overview
	local_overview$input_file <- inputs[i]
	local_overview$rep_number <- j
	local_overview$Ntip <- ape::Ntip(results[[j]]$phy)

	all_results_overview <- rbind(all_results_overview, local_overview)
	save(all_results_best, file="all_results_best.rda")
  	save(all_results_model_average, file="all_results_model_average.rda")
  	save(all_results_overview, file="all_results_overview.rda")
  	print(paste0("Finished ", inputs[i], " rep ", j))
  }

}