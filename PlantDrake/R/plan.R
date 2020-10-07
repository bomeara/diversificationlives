# Note that drake will not work with multiple cores called inside functions https://github.com/ropensci/drake/issues/675#issuecomment-458222414

plan <- drake_plan(
    tree = ape::read.tree("data/GBMB.ultra.tre"),
    #many_regimes = TryManyRegimes(tree, maxregimes=13),
    #save(many_regimes, file=file_out("results.rda"))
    try_many = target(
        SplitAndLikelihood(tree, nregimes=nregimes, interpolation_method=interpolation_method, type=type, ncores=1),
        transform = cross(
            nregimes=c(1,2,3,5,7,10,15,20,25,30,35,40),
            interpolation_method=c("linear"),
			type=c("data", "time")
        )
    ),
    everything = target(
        list(try_many),
        transform=combine(try_many)
    ),
    result_summary = SummarizeSplitsAndLikelihoods(everything),
    plot_all = PlotAll(everything, tree, file=paste0(system("hostname", intern=TRUE), "_plot.pdf")),
    save_all = save(everything, result_summary, file=paste0(system("hostname", intern=TRUE), "_everything.rda"))
)
