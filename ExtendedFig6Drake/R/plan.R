# Note that drake will not work with multiple cores called inside functions https://github.com/ropensci/drake/issues/675#issuecomment-458222414

future::plan(future::multiprocess)
plan <- drake_plan(
    session = utils::sessionInfo(),
    tree = ape::read.tree("data/tree_Extended_Data_Fig_6.tre"),
    #many_regimes = TryManyRegimes(tree, maxregimes=13),
    #save(many_regimes, file=file_out("results.rda"))
    try_many = target(
        SplitAndLikelihood(tree, nregimes=nregimes, interpolation_method=interpolation_method, type=type),
        transform = cross(
            nregimes=c(7,10,11,12),
            interpolation_method=c("linear"),
            #type=c("data", "time")
            type=c("time")
        )
    ),
    everything = target(
        list(try_many),
        transform=combine(try_many)
    ),
    result_summary = SummarizeSplitsAndLikelihoods(everything),
    adaptive_list = AdaptiveSampleBestModels(everything, result_summary, tree, deltaAIC_cutoff=10),
    plot_all = PlotAll(everything, tree, file=file_out("plot.pdf")),
    plot_uncertainty = PlotAllUncertainty(everything, tree, adaptive_list, file=file_out("uncertainty.pdf")),
    save_all = save(everything, result_summary, adaptive_list,tree, session, file=file_out("everything.rda"))
)
