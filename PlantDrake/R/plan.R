# Note that drake will not work with multiple cores called inside functions https://github.com/ropensci/drake/issues/675#issuecomment-458222414

# plan <- drake_plan(
#     tree = ape::read.tree("data/GBMB.ultra.tre"),
#     try_many = target(
#         SplitAndLikelihood(tree, nregimes=nregimes, interpolation_method=interpolation_method, type=type, ncores=1),
#         transform = cross(
#             nregimes=c(1,2,3,5,7,10,15,20,25,30,35,40),
#             interpolation_method=c("linear"),
# 			type=c("data", "time")
#         )
#     ),
#     everything = target(
#         list(try_many),
#         transform=combine(try_many)
#     ),
#     result_summary = SummarizeSplitsAndLikelihoods(everything),
#     plot_all = PlotAll(everything, tree, file=paste0(system("hostname", intern=TRUE), "_plot.pdf")),
#     save_all = save(everything, result_summary, file=paste0(system("hostname", intern=TRUE), "_everything.rda"))
# )


plan <- drake_plan(
    tree = ape::read.tree("data/GBMB.ultra.tre"),
    ef_fixed = target(
        SplitAndLikelihoodEFFixed(tree, nregimes=nregimes, interpolation_method=interpolation_method, type=type, ncores=1),
        transform = cross(
            nregimes=c(1,2,3,5,7,10,15,20,25,30,35,40),
            interpolation_method=c("linear"),
			type=c("data", "time")
        )
    ),
    everything = target(
        list(ef_fixed),
        transform=combine(ef_fixed)
    ),
    result_summary = SummarizeSplitsAndLikelihoods(everything),
    plot_all = PlotAll(everything, tree, file=paste0("ef_", system("hostname", intern=TRUE), "_plot.pdf")),
    save_all = save(everything, result_summary, file=paste0("ef_", system("hostname", intern=TRUE), "_everything.rda"))
)

plan_adaptive <- drake_plan(
    tree = ape::read.tree("data/GBMB.ultra.tre"),
    ef_fixed = target(
        SplitAndLikelihoodEFFixed(tree, nregimes=nregimes, interpolation_method=interpolation_method, type=type, ncores=1),
        transform = cross(
            nregimes=c(31,32,33,34,35,36,37,38,39),
            interpolation_method=c("linear"),
			type=c("data")
        )
    ),
    everything = target(
        list(ef_fixed),
        transform=combine(ef_fixed)
    ),
    result_summary = SummarizeSplitsAndLikelihoods(everything),
    print_result_summary = print(result_summary),
	adaptive_list = target(
		AdaptiveSampleBestModels(everything, result_summary, tree, deltaAIC_cutoff=deltaAIC_cutoff),
		transform = cross(
			deltaAIC_cutoff=c(2,3,5,7,10)
		)
	),
	everything2 = target(
        list(adaptive_list),
        transform=combine(adaptive_list)
    ),
    save_adaptive = save(tree, session, result_summary, everything, everything2, file=paste0("ef_adaptive_", system("hostname", intern=TRUE), ".rda"))
)

