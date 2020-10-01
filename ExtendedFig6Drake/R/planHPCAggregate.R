# Note that drake will not work with multiple cores called inside functions https://github.com/ropensci/drake/issues/675#issuecomment-458222414



plan_hpc_aggregate <- drake_plan(
    session = utils::sessionInfo(),
    tree = ape::read.tree("/Users/bomeara/Documents/MyDocuments/GitClones/diversificationlives/ExtendedFig6Drake/data/tree_Extended_Data_Fig_6.tre"),
    #many_regimes = TryManyRegimes(tree, maxregimes=13),
    #save(many_regimes, file=file_out("results.rda"))
    everything = AggregateEverythingFromSavedFies(),
    result_summary = SummarizeSplitsAndLikelihoods(everything),
    print_result_summary = print(result_summary),
    subset_list = SubsetBestModels(everything, result_summary, tree, deltaAIC_cutoff=10),
    #plot_all = PlotAll(everything, tree, file=file_out("plot.pdf")),
    #plot_uncertainty = PlotAllUncertainty(everything, tree, subset_list, file=file_out("uncertainty.pdf")),
    #save_all = save(everything, result_summary, subset_list,tree, session, file=file_out("everything.rda"))
	params_all = GetParamsFromAggregation(subset_list),
	params_quantile = apply(params_all, 2, quantile, probs=seq(from=0, to=1, length.out=11)),
	params_all_csv = write.csv(params_all, file=file_out("params_all.csv")),
	params_quantile_csv = write.csv(params_quantile, file=file_out("params_quantile.csv"))

)
