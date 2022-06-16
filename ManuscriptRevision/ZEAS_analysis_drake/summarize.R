library(hisse)
library(Metrics)


all_results_best <- data.frame()
aggregate_results_best <- data.frame()

inputs <- list.files(pattern="omeara.*results.rda", full.names=TRUE)

for (i in sequence(length(inputs))) {
  load(inputs[i])
  for(j in sequence(length(results))) {
	local_best <- as.data.frame(results[[j]]$summaries$rates[,,"best"])
	local_best$speciation_true <- results[[j]]$lambda[length(results[[j]]$lambda)]
	local_best$extinction_true <- 0
	local_best$final_time <- results[[j]]$final_time
	local_best$crown_time <- results[[j]]$root_time
	local_best$turnover_true <- local_best$speciation_true + local_best$extinction_true
	local_best$extinction.fraction_true <- local_best$extinction_true/local_best$speciation_true 
	local_best$net.div_true <- local_best$speciation_true - local_best$extinction_true
	local_best$input_file <- inputs[i]
	local_best$rep_number <- j
	local_best <- local_best[sequence(ceiling(nrow(local_best)/2)),] # to get only the tips
	local_best$taxon_number <- sequence(nrow(local_best))
	all_results_best <- rbind(all_results_best, local_best)
	aggregate_results_best <- rbind(aggregate_results_best, data.frame(net.div_range=diff(range(local_best$net.div)), net.div_rmse=Metrics::rmse(local_best$net.div, local_best$net.div_true), net.div_median = median(local_best$net.div), net.div_true=local_best$net.div_true[1],input_file = inputs[i], rep_number = j))
	
  }
  save(all_results_best, file="all_results_best.rda")
  save(aggregate_results_best, file="aggregate_results_best.rda")
  print(paste0("Finished ", inputs[i], " rep ", j))
}

