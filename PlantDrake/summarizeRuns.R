library(drake)
setwd("/share/diversificationlives/PlantDrake")
loadd(path="~/Documents/localcache")
save(list=ls(), file=paste0("ef_", system("hostname", intern=TRUE), "_inprocess_everything.rda"))