# Note that drake will not work with multiple cores called inside functions https://github.com/ropensci/drake/issues/675#issuecomment-458222414

future::plan(future::multiprocess)
plan <- drake_plan(
    tree = ape::read.tree("data/tree_Extended_Data_Fig_6.tre"),
    #many_regimes = TryManyRegimes(tree, maxregimes=13),
    #save(many_regimes, file=file_out("results.rda"))
    try_many = target(
        SplitAndLikelihood(tree, nregimes=nregimes, interpolation_method=interpolation_method, type=type),
        transform = cross(
            nregimes=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14),
            interpolation_method=c("linear"),
            type=c("data", "time")
        )
    ),
    everything = target(
        list(try_many),
        transform=combine(try_many)
    ),
    result_summary = SummarizeSplitsAndLikelihoods(everything),
    plot_all = PlotAll(everything, tree, file=file_out("plot.pdf")),
    save_all = save(everything, result_summary, file=file_out("everything.rda"))
)
