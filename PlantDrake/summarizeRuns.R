library(drake)
setwd("/share/diversificationlives/PlantDrake")
loadd(path="~/Documents/localcache")
save(list=ls(), paste0(system("hostname", intern=TRUE), "_inprocess_everything.rda"))