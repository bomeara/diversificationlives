# Note that drake will not work with multiple cores called inside functions https://github.com/ropensci/drake/issues/675#issuecomment-458222414

future::plan(future::multiprocess)
plan <- drake_plan(
    tree = ape::read.tree("data/GBMB.ultra.tre"),
    splits7datalinear = SplitAndLikelihood(tree=tree, nregimes=7, type="data", interpolation_method="linear")
)
