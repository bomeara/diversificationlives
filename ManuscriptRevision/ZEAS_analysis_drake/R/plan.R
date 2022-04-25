
plan <- drake_plan(
  simulated_trees = SimulateZEAS(),
  misse_fit = MisseFit(simulated_trees)
)
