source("packages.R")  # loads packages
source("functions.R")



results = RunMany(nrep=100, ncore=parallel::detectCores())
write(results, file=paste0(system("hostname", intern=TRUE), "_final_results.rda"))

# ansible linux -a 'nohup Rscript /share/diversificationlives/ManuscriptRevision/ZEAS_analysis_drake/make.R &'