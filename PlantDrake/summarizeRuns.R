library(drake)
setwd("/share/diversificationlives/PlantDrake")
loadd(path="~/Documents/localcache")
print(ls())
save(list=ls(), file=paste0("ef_", system("hostname", intern=TRUE), "_inprocess_everything.rda"))

# ansible linux -a 'nohup Rscript /share/diversificationlives/PlantDrake/summarizeRuns.R &'

# then on one computer, summarizeRunsStep2.R