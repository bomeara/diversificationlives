library(ape)
library(castor)

rm(list=ls())
tree <- rcoal(200)
tree$edge.length <- tree$edge.length/max(branching.times(tree))

#models are congruent if the rp_a == rp_b and lambda_p_a == lambda_p_b
#rp(t) = lambda(t) - mu(t) + (1/lambda(t)) * dlambda_dt
#lambda_p = lambda_0 * rho_0

lambda_0_a = 0.5
lambda_0_b = lambda_0_a

mu_0_a = 0.2

dlambda_dt_a = 0
dlambda_dt_b = 0
dmu_dt_a = 0

rho_0_a = 1
rho_0_b = rho_0_a

#  so lambda_a(t) = lambda_0_a + t*dlambda_dt_a
# and lambda_b(t) = lambda_0_b + t*dlambda_dt_b

lambda_a_fn <- function(t) {
  results <- lambda_0_a + t*dlambda_dt_a
  if(any(results<0)) {
    stop("lambda_a < 0")
  }
  return(results)
}

lambda_b_fn <- function(t) {
  results <- lambda_0_b + t*dlambda_dt_b
  if(any(results<0)) {
    stop("lambda_b < 0")
  }
  return(results)
}



#thus, lambda_p_a == lambda_p_b

#now derive mu_b if we know mu_a, dlambda_dt_a, and dlambda_dt_b



mu_a_fn <- function(t) {
 results <- mu_0_a+dmu_dt_a*t
  if(any(results<0)) {
    stop("mu_a < 0")
  }
  return(results)
}

# lambda_b(t) - mu_b(t) + (1/lambda_b(t)) * dlambda_dt_b == lambda_a(t) - mu_a(t) + (1/lambda_a(t)) * dlambda_dt_a
# substituting
# lambda_0_b + t*dlambda_dt_b -  mu_b(t) + (1/(lambda_0_b + t*dlambda_dt_b)) * dlambda_dt_b == lambda_0_a + t*dlambda_dt_a - mu_0_a + (1/(lambda_0_a + t*dlambda_dt_a)) * dlambda_dt_a
# since dlambda_dt_a == 0
# lambda_0_b + t*dlambda_dt_b -  mu_b(t) + (1/(lambda_0_b + t*dlambda_dt_b)) * dlambda_dt_b == lambda_0_a - mu_0_a
#  mu_b(t) = lambda_0_b + t*dlambda_dt_b + (1/(lambda_0_b + t*dlambda_dt_b)) * dlambda_dt_b - lambda_0_a + mu_0_a

mu_b_fn <- function(t) {
   results <- lambda_0_b + t*dlambda_dt_b + (1/(lambda_0_b + t*dlambda_dt_b)) * dlambda_dt_b - lambda_0_a + mu_0_a
   if(any(results<0)) {
     stop("mu_b < 0")
   }
   return(results)
}

age_grid = seq(from=0, to=get_tree_span(tree)$max_distance, length.out=1001)

lik_a = loglikelihood_hbd(
  tree=tree,
  rho0=rho_0_a,
  age_grid=age_grid,
  lambda=lambda_a_fn(age_grid),
  mu=mu_a_fn(age_grid)
)

lik_b = loglikelihood_hbd(
  tree=tree,
  rho0=rho_0_b,
  age_grid=age_grid,
  lambda=lambda_b_fn(age_grid),
  mu=mu_b_fn(age_grid)
)



sim_a <- simulate_deterministic_hbd(
  LTT0 = ape::Ntip(tree),
  oldest_age = get_tree_span(tree)$max_distance,
  age0 = 0,
  rho0 = rho_0_a,
  age_grid = age_grid,
  lambda=lambda_a_fn(age_grid),
  mu=mu_a_fn(age_grid)
)

sim_b <- simulate_deterministic_hbd(
  LTT0 = ape::Ntip(tree),
  oldest_age = get_tree_span(tree)$max_distance,
  age0 = 0,
  rho0 = rho_0_b,
  age_grid = age_grid,
  lambda=lambda_b_fn(age_grid),
  mu=mu_b_fn(age_grid)
)

par(mfcol=c(1,3))
plot(age_grid, lambda_a_fn(age_grid), type="l", lwd=2, ylim=range(c(lambda_a_fn(age_grid),lambda_b_fn(age_grid))), ylab="lambda")
lines(age_grid, lambda_b_fn(age_grid), col="red")


plot(age_grid, mu_a_fn(age_grid), type="l", lwd=2, ylim=range(c(mu_a_fn(age_grid),mu_b_fn(age_grid))), ylab="mu")
lines(age_grid, mu_b_fn(age_grid), col="red")

plot(sim_a$ages, sim_a$LTT, type="l", lwd=2, ylab="dLTT")
lines(sim_b$ages, sim_b$LTT, col="red")

print(lik_a)
print(lik_b)