# Note that drake will not work with multiple cores called inside functions https://github.com/ropensci/drake/issues/675#issuecomment-458222414

plan <- drake_plan(
    tree = ape::read.tree("data/whales_Steemanetal2009.tre"),
    fit_param_lambda7p_mu7p_crown = param_lambda7p_mu7p(desired_interval = 0.1, tree=tree, condition="crown"),

    fit_param_lambda7p_mu7p_stem = param_lambda7p_mu7p(desired_interval = 0.1, tree=tree, condition="stem"),

    fit_param_lambda7p_mu_multiplier_lambda_crown = param_lambda7p_mu_multiplier_lambda(desired_interval = 0.1, tree=tree, condition="crown"),

    fit_param_lambda7p_mu_multiplier_lambda_stem = param_lambda7p_mu_multiplier_lambda(desired_interval = 0.1, tree=tree, condition="stem"),

    all_results = list(
        fit_param_lambda7p_mu7p_crown=fit_param_lambda7p_mu7p_crown,
        fit_param_lambda7p_mu7p_stem=fit_param_lambda7p_mu7p_stem,
        fit_param_lambda7p_mu_multiplier_lambda_crown=fit_param_lambda7p_mu_multiplier_lambda_crown,
        fit_param_lambda7p_mu_multiplier_lambda_stem=fit_param_lambda7p_mu_multiplier_lambda_stem
    ),

    all_results_saved = save(all_results, file=file.out("results.rda"))
)
