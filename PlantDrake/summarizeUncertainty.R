library(drake)
setwd("/share/diversificationlives/PlantDrake")
files <- list.files(pattern="ef_adaptive_new_.*.rda")
for(i in seq_along(files)) {
	load(files[i])
	everything_adaptive_compiled$source <- files[i]
	if(i==1) {
		global_everything_adaptive_compiled <- everything_adaptive_compiled
	} else {
		for(run_index in seq_along(global_everything_adaptive_compiled)) {
			global_everything_adaptive_compiled[[run_index]] <- rbind(global_everything_adaptive_compiled[[run_index]], everything_adaptive_compiled[[run_index]])
		}
	}
	rm(everything_adaptive_compiled)
}
save(list=ls(), file="global_adaptive_new.rda")
PlotAllUncertainty(everything_adaptive, tree, global_everything_adaptive_compiled,  file=paste0("global_adaptive_new_UNCERTAINTY_2_plot.pdf"), desired_delta=2)
PlotAllUncertainty(everything_adaptive, tree, global_everything_adaptive_compiled,  file=paste0("global_adaptive_new_UNCERTAINTY_5_plot.pdf"), desired_delta=5)
PlotAllUncertainty(everything_adaptive, tree, global_everything_adaptive_compiled,  file=paste0("global_adaptive_new_UNCERTAINTY_10_plot.pdf"), desired_delta=10)

