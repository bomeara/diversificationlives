source("../WhaleDrake/R/functions.R")

CompileEverythingAdaptiveNew <- function(everything_adaptive_new) {
	everything_adaptive_compiled <- everything_adaptive_new[[1]]
	for (i in seq(from=2, to=length(everything_adaptive_new), by=1)) {
		for (j in seq_along(everything_adaptive_new[[i]])) {
			everything_adaptive_compiled[[j]] <- rbind(everything_adaptive_compiled[[j]], everything_adaptive_new[[i]][[j]])
		}
	}
	return(everything_adaptive_compiled)
}