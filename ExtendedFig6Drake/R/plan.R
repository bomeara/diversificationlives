# Note that drake will not work with multiple cores called inside functions https://github.com/ropensci/drake/issues/675#issuecomment-458222414

#future::plan(future::multiprocess)

plan_original <- drake_plan(
    session = utils::sessionInfo(),
    tree = ape::read.tree("data/tree_Extended_Data_Fig_6.tre"),
    #many_regimes = TryManyRegimes(tree, maxregimes=13),
    #save(many_regimes, file=file_out("results.rda"))
    try_many = target(
        SplitAndLikelihood(tree, nregimes=nregimes, interpolation_method=interpolation_method, type=type, Ntrials=1, ncores=1),
        transform = cross(
            nregimes=c(10,11,12),
            interpolation_method=c("linear"),
            #type=c("data", "time")
            type=c("time"),
            instance=c(1,2,3,4,5,6,7,8,9,10,11,12)
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

instances <- sequence(parallel::detectCores())

plan_manystart <- drake_plan(
    session = utils::sessionInfo(),
    tree = ape::read.tree("data/tree_Extended_Data_Fig_6.tre"),
    #many_regimes = TryManyRegimes(tree, maxregimes=13),
    #save(many_regimes, file=file_out("results.rda"))
    hostname = system("hostname -I", intern=TRUE),
    try_many = target(
        SplitAndLikelihood(tree, nregimes=nregimes, interpolation_method=interpolation_method, type=type, Ntrials=1, ncores=1),
        transform = cross(
            nregimes=c(11),
            interpolation_method=c("linear"),
            #type=c("data", "time")
            type=c("time"),
            instance=instances
        )
    ),
    save_try = save(tree, session, try_many, file=file_out("trymany.rda")),
    optimize_many = target(
        OptimizeLogSpace(try_many, tree=tree),
        transform = map(try_many)
    ),
    save_optim = save(tree, session, optimize_many, file=file_out("optimmany.rda")),
    everything = target(
        list(optimize_many),
        transform=combine(optimize_many)
    ),
    result_summary = SummarizeSplitsAndLikelihoods(everything),
    save_summary = save(tree, session, result_summary, everything, file=file_out("summary.rda")),
    print_result_summary = print(result_summary),
    adaptive_list = AdaptiveSampleBestModels(everything, result_summary, tree, deltaAIC_cutoff=10),
    save_adaptive = save(tree, session, result_summary, everything, adaptive_list, file=file_out("adaptive.rda")),
    plot_all = PlotAll(everything, tree, file=file_out(paste0(hostname,"_plot.pdf"))),
    plot_uncertainty = PlotAllUncertainty(everything, tree, adaptive_list, file=file_out("uncertainty.pdf")),
    save_all = save(everything, result_summary, adaptive_list,tree, session, file=file_out("everything.rda"))
)

# Note that drake will not work with multiple cores called inside functions https://github.com/ropensci/drake/issues/675#issuecomment-458222414

#future::plan(future::multiprocess)
my_plan <- drake_plan(
    session = utils::sessionInfo(),
    tree = ape::read.tree("data/tree_Extended_Data_Fig_6.tre"),
    #many_regimes = TryManyRegimes(tree, maxregimes=13),
    #save(many_regimes, file=file_out("results.rda"))
    try_many = target(
        SplitAndLikelihood(tree, nregimes=nregimes, interpolation_method=interpolation_method, type=type, Ntrials=3, ncores=1),
        transform = cross(
            nregimes=c(10,11,12),
            interpolation_method=c("linear"),
            #type=c("data", "time")
            type=c("time"),
            instance=c(1,2,3,4)
        )
    ),
    optimize_many = target(
        OptimizeLogSpace(try_many, tree=tree),
        transform = map(try_many)
    ),
    everything_preoptimize = target(
        list(try_many),
        transform=combine(try_many)
    ),
    everything = target(
        list(optimize_many),
        transform=combine(optimize_many)
    ),
    result_summary = SummarizeSplitsAndLikelihoods(everything),
    result_summary_preoptimize = SummarizeSplitsAndLikelihoods(everything_preoptimize),
    adaptive_list = AdaptiveSampleBestModels(everything, result_summary, tree, deltaAIC_cutoff=10),
    plot_all = PlotAll(everything, tree, file=file_out("plot.pdf")),
    plot_uncertainty = PlotAllUncertainty(everything, tree, adaptive_list, file=file_out("uncertainty.pdf")),
    save_all = save(everything, everything_preoptimize, result_summary, adaptive_list,tree, session, file=file_out("everything.rda"))
)
