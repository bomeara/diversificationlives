

SimulateZEAS <- function(ntax=100) {
	# define time grid on which lambda, mu and psi will be specified
     time_grid = seq(0,500,length.out=10000)

     # specify the time-dependent spexciation rate lambda on the time-grid
     #lambda_grid = 0.01 + .02*exp(0.01*time_grid)
	 offset <- 0.05
	 lambda_grid=-offset + offset*exp(0.01*time_grid)

    simul <- list(success=FALSE)
	while(!simul$success) {
     # generate tree with a constant speciation & sampling rate,
     # time-variable extinction rate and additional discrete sampling points
     # assuming that all continuously sampled lineages are removed from the pool
     simul = generate_tree_hbds( max_extant_tips = ntax+1,
                                 include_extant  = TRUE,
                                 include_extinct = FALSE,
                                 time_grid       = time_grid,
                                 lambda          = lambda_grid,
                                 mu              = 0,
                                 psi             = 0,
                                 kappa           = 0)
	}
	tip_lengths <- simul$tree$edge.length[which(simul$tree$edge[,2]<=ntax+1)]
	nonzero_tip_lengths <- tip_lengths[tip_lengths>0]
	phy <- TreeSim::cuttree(simul$tree, 0.5*min(nonzero_tip_lengths))
	return(list(phy=phy, final_time=simul$final_time, root_time=simul$root_time, lambda = lambda_grid[which(time_grid<=simul$final_time)]))	
}

RunMiSSe <- function(zeas_sim, ncore=3) {
	greedy <- hisse::MiSSEGreedy(zeas_sim$phy, n.cores=ncore)
	greedy_summaries <- hisse::SummarizeMiSSEGreedy(greedy,n.cores=ncore) 
	return(list(greedy=greedy, summaries=greedy_summaries, phy=zeas_sim$phy, final_time=zeas_sim$final_time, root_time=zeas_sim$root_time, lambda=zeas_sim$lambda))	
}

RunMany <- function(nrep=100, ncore=parallel::detectCores()) {
	results <- list()
	for (i in sequence(nrep)) {
		print(i)
		zeas_sim <- SimulateZEAS()
		try(results[[i]] <- RunMiSSe(zeas_sim, ncore=ncore))
		save(results, file=paste0(system("hostname", intern=TRUE), "_results.rda"))
	}
	return(results)
}