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
	aggregate_results_best <- rbind(aggregate_results_best, data.frame(net.div_range=diff(range(local_best$net.div)), net.div_rmse=Metrics::rmse(local_best$net.div, local_best$net.div_true), net.div_median = median(local_best$net.div), net.div_true=local_best$net.div_true[1], net.div_sd=sd(local_best$net.div), speciation_median = median(local_best$speciation), speciation_sd=sd(local_best$speciation), extinction_median=median(local_best$extinction), extinction_sd=sd(local_best$extinction), input_file = inputs[i], rep_number = j, final_time=results[[j]]$final_time))
	
  }
  save(all_results_best, file="all_results_best.rda")
  save(aggregate_results_best, file="aggregate_results_best.rda")
  print(paste0("Finished ", inputs[i], " rep ", j))
}

time_grid = seq(0,500,length.out=10000)
time_grid <- time_grid[which(time_grid<=max(aggregate_results_best$final_time))]

     # specify the time-dependent spexciation rate lambda on the time-grid
     #lambda_grid = 0.01 + .02*exp(0.01*time_grid)
offset <- 0.05
lambda_grid <- -offset + offset*exp(0.01*time_grid)

pdf(file="ZEAS.pdf", width=10, height=10)
par(mfrow=c(2,2))
plot(time_grid, lambda_grid, type="l", ylim=c(0, max(quantile(aggregate_results_best$extinction_median, 0.999), quantile(aggregate_results_best$speciation_median, 0.999))), xlab="Time", ylab="Speciation Rate", col="red", main="True speciation rate at tips versus estimates")
points(aggregate_results_best$final_time, aggregate_results_best$speciation_median, pch='.')
for (i in sequence(nrow(aggregate_results_best))) {
	lines(x=rep(aggregate_results_best$final_time[i],2), c(aggregate_results_best$speciation_median[i]+aggregate_results_best$speciation_sd[i], aggregate_results_best$speciation_median[i]-aggregate_results_best$speciation_sd[i]), col="gray")
	
}

plot(time_grid, 0*lambda_grid, type="l", ylim=c(0, max(quantile(aggregate_results_best$extinction_median, 0.999), quantile(aggregate_results_best$speciation_median, 0.999))), xlab="Time", ylab="Extinction Rate", col="red", main="True extinction rate at tips versus estimates")
points(aggregate_results_best$final_time, aggregate_results_best$extinction_median, pch='.')
for (i in sequence(nrow(aggregate_results_best))) {
	lines(x=rep(aggregate_results_best$final_time[i],2), c(aggregate_results_best$extinction_median[i]+aggregate_results_best$extinction_sd[i], aggregate_results_best$extinction_median[i]-aggregate_results_best$extinction_sd[i]), col="gray")
	
}


plot(time_grid, lambda_grid, type="l", ylim=c(0, max(quantile(aggregate_results_best$extinction_median, 0.999), quantile(aggregate_results_best$speciation_median, 0.999))), xlab="Time", ylab="Net Diversification Rate", col="red", main="True net diversification rate at tips versus estimates")
points(aggregate_results_best$final_time, aggregate_results_best$net.div_median, pch='.')
for (i in sequence(nrow(aggregate_results_best))) {
	lines(x=rep(aggregate_results_best$final_time[i],2), c(aggregate_results_best$net.div_median[i]+aggregate_results_best$net.div_sd[i], aggregate_results_best$net.div_median[i]-aggregate_results_best$net.div_sd[i]), col="gray")
	
}

hist(aggregate_results_best$net.div_range[which(aggregate_results_best$net.div_range<max(aggregate_results_best$net.div_range))], main="Range of net diversification rate tip estimates\n(not shown: a single simulation with a difference of 6.7)", breaks=50, freq=TRUE, xlab="Maximum difference between tip estimates of rate")
print(paste0(sum(aggregate_results_best$net.div_range==0), " of ", nrow(aggregate_results_best), " differences are exactly zero"))
dev.off()