# Note that drake will not work with multiple cores called inside functions https://github.com/ropensci/drake/issues/675#issuecomment-458222414



plan_hpc <- drake_plan(
    session = utils::sessionInfo(),
    tree = ape::read.tree("~/Documents/diversificationlives/ExtendedFig6Drake/data/tree_Extended_Data_Fig_6.tre"),
    #many_regimes = TryManyRegimes(tree, maxregimes=13),
    #save(many_regimes, file=file_out("results.rda"))
    try_many = target(
        SplitAndLikelihood(tree, nregimes=nregimes, interpolation_method=interpolation_method, type=type, Ntrials=4, ncores=1),
        transform = cross(
            nregimes=c(10,11,12),
            interpolation_method=c("linear"),
            #type=c("data", "time")
            type=c("time"),
            instance=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32)
        )
    ),
    optimize_many = target(
        OptimizeLogSpace(try_many, tree=tree),
        transform = map(try_many)
    ),
    everything = target(
        list(optimize_many),
        transform=combine(optimize_many)
    ),
    result_summary = SummarizeSplitsAndLikelihoods(everything),
    print_result_summary = print(result_summary),
    adaptive_list = AdaptiveSampleBestModels(everything, result_summary, tree, deltaAIC_cutoff=10),
    plot_all = PlotAll(everything, tree, file=file_out("plot.pdf")),
    plot_uncertainty = PlotAllUncertainty(everything, tree, adaptive_list, file=file_out("uncertainty.pdf")),
    save_all = save(everything, result_summary, adaptive_list,tree, session, file=file_out("everything.rda"))
)
