# Note that drake will not work with multiple cores called inside functions https://github.com/ropensci/drake/issues/675#issuecomment-458222414



plan_hpc_adaptivemany <- drake_plan(
    session = utils::sessionInfo(),
    tree = ape::read.tree("/Users/bomeara/Documents/MyDocuments/GitClones/diversificationlives/ExtendedFig6Drake/data/tree_Extended_Data_Fig_6.tre"),
    #many_regimes = TryManyRegimes(tree, maxregimes=13),
    #save(many_regimes, file=file_out("results.rda"))
    everything = AggregateEverythingFromSavedFies(),
    result_summary = SummarizeSplitsAndLikelihoods(everything),
    print_result_summary = print(result_summary),
    subset_list = SubsetBestModels(everything, result_summary, tree, deltaAIC_cutoff=10)
)
