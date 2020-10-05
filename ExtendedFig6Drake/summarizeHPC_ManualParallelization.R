library(dplyr)
all_results <- data.frame()
finished <- list.files(pattern="adaptive.*rda")
for (i in seq_along(finished)) {
	print(finished[i])
	load(finished[i])
	print(length(results))
	for (j in seq_along(results)) {
		all_results <- dplyr::bind_rows(all_results, results[[j]])
	}
	rm(results)
}
all_results <- subset(all_results, is.finite(all_results$loglikelihood))
all_results$negloglikelihood <- -1*all_results$loglikelihood
all_results$deltalnl <- all_results$negloglikelihood - min(all_results$negloglikelihood)

param_range <- as.numeric(gsub("lambda", "", colnames(all_results)[grepl("lambda", colnames(all_results))]))
for (param_index in seq_along(param_range)) {
	all_results$divNEW <- all_results[,paste0("lambda", param_range[param_index])] - all_results[,paste0("mu", param_range[param_index])]
	colnames(all_results)[ncol(all_results)] <- paste0("netdiv", param_range[param_index])
}

for (param_index in seq_along(param_range)) {
	all_results$turnNEW <- all_results[,paste0("lambda", param_range[param_index])] + all_results[,paste0("mu", param_range[param_index])]
	colnames(all_results)[ncol(all_results)] <- paste0("turnover", param_range[param_index])
}

for (param_index in seq_along(param_range)) {
	all_results$efNEW <- all_results[,paste0("mu", param_range[param_index])] / all_results[,paste0("lambda", param_range[param_index])]
	colnames(all_results)[ncol(all_results)] <- paste0("ef", param_range[param_index])
}


all_results_delta2 <- subset(all_results, deltalnl<=2)
write.csv(apply(all_results_delta2, 2, range), file="summarydelta2.csv")
save(all_results, file="allruns.rda")