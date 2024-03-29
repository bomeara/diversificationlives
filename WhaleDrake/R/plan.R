# Note that drake will not work with multiple cores called inside functions https://github.com/ropensci/drake/issues/675#issuecomment-458222414

future::plan(future::multiprocess)
plan <- drake_plan(
    tree = ape::read.tree("data/whales_Steemanetal2009.tre"),
    fit_param_lambda7p_mu7p_crown = param_lambda7p_mu7p(desired_interval = 0.1, tree=tree, condition="crown"),

    fit_param_lambda7p_mu7p_stem = param_lambda7p_mu7p(desired_interval = 0.1, tree=tree, condition="stem"),

    fit_param_lambda7p_mu_multiplier_lambda_crown = param_lambda7p_mu_multiplier_lambda(desired_interval = 0.1, tree=tree, condition="crown"),

    fit_param_lambda7p_mu_multiplier_lambda_stem = param_lambda7p_mu_multiplier_lambda(desired_interval = 0.1, tree=tree, condition="stem"),

    fit_param_lambda_discreteshift_mu_discreteshift_crown_1my_constant = param_lambda_discreteshift_mu_discreteshift(desired_interval = 0.1, tree=tree, condition="crown", interpolation_method="constant", slice_ages = seq(from=0, to=floor(castor::get_tree_span(tree)$max_distance), by=1)),

    fit_param_lambda_discreteshift_mu_discreteshift_crown_1my_linear = param_lambda_discreteshift_mu_discreteshift(desired_interval = 0.1, tree=tree, condition="crown", interpolation_method="linear", slice_ages = seq(from=0, to=floor(castor::get_tree_span(tree)$max_distance), by=1)),

    fit_param_lambda_discreteshift_ef_0_crown_1my_linear = param_lambda_discreteshift_ef_fixed(desired_interval = 0.1, tree=tree, condition="crown", interpolation_method="linear", slice_ages = seq(from=0, to=floor(castor::get_tree_span(tree)$max_distance), by=1), ef=0.0),

    fit_param_lambda_discreteshift_ef_0_9_crown_1my_linear = param_lambda_discreteshift_ef_fixed(desired_interval = 0.1, tree=tree, condition="crown", interpolation_method="linear", slice_ages = seq(from=0, to=floor(castor::get_tree_span(tree)$max_distance), by=1), ef=0.9),


    all_results = list(
        fit_param_lambda7p_mu7p_crown=fit_param_lambda7p_mu7p_crown,
        fit_param_lambda7p_mu7p_stem=fit_param_lambda7p_mu7p_stem,
        fit_param_lambda7p_mu_multiplier_lambda_crown=fit_param_lambda7p_mu_multiplier_lambda_crown,
        fit_param_lambda7p_mu_multiplier_lambda_stem=fit_param_lambda7p_mu_multiplier_lambda_stem,
        fit_param_lambda_discreteshift_mu_discreteshift_crown_1my_constant=fit_param_lambda_discreteshift_mu_discreteshift_crown_1my_constant,
        fit_param_lambda_discreteshift_mu_discreteshift_crown_1my_linear=fit_param_lambda_discreteshift_mu_discreteshift_crown_1my_linear,
        fit_param_lambda_discreteshift_ef_0_crown_1my_linear=fit_param_lambda_discreteshift_ef_0_crown_1my_linear,
        fit_param_lambda_discreteshift_ef_0_9_crown_1my_linear=fit_param_lambda_discreteshift_ef_0_9_crown_1my_linear
    ),
    treeultra = phytools::force.ultrametric(tree, method="extend"),

    #splits7datalinear = SplitAndLikelihood(tree=treeultra, nregimes=7, type="data", interpolation_method="linear"),

    #many_regimes = TryManyRegimes(treeultra, maxregimes=5),
    all_results_saved = save(all_results, file=file_out("results.rda"))
)

planmany <- drake_plan(
    tree = phytools::force.ultrametric(ape::read.tree("data/whales_Steemanetal2009.tre"),method="extend"),

    # t3 = SplitAndLikelihood(tree, nregimes=3, interpolation_method="constant"),
    try_many = target(
        SplitAndLikelihood(tree, nregimes=nregimes, interpolation_method=interpolation_method),
        transform = cross(
            nregimes=c(1,2,3,4,5,6,7),
            interpolation_method=c("linear")
        )
    ),
    everything = target(
        list(try_many),
        transform=combine(try_many)
    ),
    result_summary = SummarizeSplitsAndLikelihoods(everything),
    plot_all = PlotAll(everything, tree, file=file_out("plot.pdf"))

    #many_results_saved = save(try_many, file=file_out("results_many.rda"))
)
