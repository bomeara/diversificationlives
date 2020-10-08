setwd("/share/diversificationlives/PlantDrake")
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
summaries <- list.files(pattern="ef_*inprocess", full.names=TRUE)
summary_df <- data.frame()
for (i in seq_along(summaries)) {
	print(summaries[i])
	load(summaries[i])
	trials <-  ls(pattern="ef_")
	for (j in seq_along(trials)) {
		summary_df <- dplyr::bind_rows(summary_df, data.frame(name=trials[j], AIC=get(trials[j])$AIC, type=get(trials[j])$type, nregimes=get(trials[j])$nregimes, source=summaries[i]))	
		print(tail(summary_df, 1))
	}
}
summary_df$deltaAIC <- summary_df$AIC - min(summary_df$AIC)
save(summary_df, file="allsummaries.rda")
g <- ggplot(data=summary_df, aes(x=nregimes, y=deltaAIC, group=type)) + geom_line(aes(color=type)) + geom_beeswarm(aes(color=type))
pdf(file="summary.pdf")
print(g)
dev.off()