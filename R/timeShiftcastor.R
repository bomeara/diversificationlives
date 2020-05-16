
library(castor)
library(ape)




######################################################################################################################################
######################################################################################################################################
### Running analyses and Plotting code
######################################################################################################################################
######################################################################################################################################

tree <- read.tree("../data/whales_Steemanetal2009.tre")
root_age <- castor::get_tree_span(tree)$max_distance


### set up done, now do grid

desired_interval=1
age_grid  = seq(from=0,to=root_age+desired_interval,by=desired_interval)
fit_yule_splines0 = fit_hbd_model_on_grid(tree, 
                            age_grid    = age_grid,
                            max_mu      = 1,
                            fixed_rho0   = 1,
                            fixed_mu = 0,
                            condition="crown",
                            Ntrials=10,
                            Nthreads=parallel::detectCores(),
                            splines_degree = 0)
save(list=ls(), file="yule0.rda")

fit_yule_splines1 = fit_hbd_model_on_grid(tree, 
                                          age_grid    = age_grid,
                                          max_mu      = 1,
                                          fixed_rho0   = 1,
                                          fixed_mu = 0,
                                          condition="crown",
                                          Ntrials=10,
                                          Nthreads=parallel::detectCores(),
                                          splines_degree = 1)
save(list=ls(), file="yule1.rda")

fit_bd_splines0 = fit_hbd_model_on_grid(tree, 
                                          age_grid    = age_grid,
                                          max_mu      = 1,
                                          fixed_rho0   = 1,
                                          condition="crown",
                                          Ntrials=10,
                                          Nthreads=parallel::detectCores(),
                                          splines_degree = 0)
save(list=ls(), file="bd0.rda")

fit_bd_splines1 = fit_hbd_model_on_grid(tree, 
                                        age_grid    = age_grid,
                                        max_mu      = 1,
                                        fixed_rho0   = 1,
                                        condition="crown",
                                        Ntrials=10,
                                        Nthreads=parallel::detectCores(),
                                        splines_degree = 1)
save(list=ls(), file="bd1.rda")


#plot(x=-1*fit$age_grid, y=fit$fitted_lambda, type="l")

fits_grid <- list(fit_yule_splines0=fit_yule_splines0, fit_yule_splines1=fit_yule_splines1, fit_bd_splines0=fit_bd_splines0, fit_bd_splines1=fit_bd_splines1)

### Now for parametric

rho = 1
desired_interval = 0.1
age_grid_param = seq(from=0,to=root_age+desired_interval,by=desired_interval)

lambda_function = function(ages,params){
    return(params['lambda0']*exp(-params['alpha']*ages));
}
mu_function = function(ages,params){
    return(params['mu0']*exp(-params['beta']*ages));
}
rho_function = function(params){
    return(rho) # rho does not depend on any of the parameters
}

fit_param_yule = fit_hbd_model_parametric(	tree, 
                                param_values  = c(lambda0=NA, mu0=0, alpha=NA, beta=NA),
                                param_guess   = c(0.1,0,0,0),
                                param_min     = c(0,0,-1,-1),
                                param_max     = c(2,2,1,1),
                                param_scale   = 1, # all params are in the order of 1
                                lambda        = lambda_function,
                                mu            = mu_function,
                                rho0          = rho_function,
                                age_grid      = age_grid_param,
                                condition     = "crown",
                                Ntrials       = 10,    # perform 10 fitting trials
                                Nthreads      = parallel::detectCores(),     # use 2 CPUs
                                max_model_runtime = 1, # limit model evaluation to 1 second
                                fit_control       = list(rel.tol=1e-6)
                            )


save(list=ls(), file="runparam.rda")


# plot(range(age_grid_param), range(c(lambda_function(age_grid_param, fit_param_yule$param_fitted), mu_function(age_grid_param, fit_param_yule$param_fitted))), type="n", bty="n", ylab="rate")
# lines(age_grid_param, lambda_function(age_grid_param, fit_param_yule$param_fitted))
# lines(age_grid_param, mu_function(age_grid_param, fit_param_yule$param_fitted), lty="dotted")

# parametric with seven params each rate

rho = 1
desired_interval = 0.1
age_grid_param = seq(from=0,to=root_age+desired_interval,by=desired_interval)

lambda_function = function(ages,params){
    results <- params['p1'] * exp(-params['p2']*ages/max(ages)) + 
               params['p3'] +
               params['p4'] * (ages/max(ages)) +
               params['p5'] * ((ages/max(ages))^2) +
               params['p6'] * ((ages/max(ages))^3) +
               params['p7'] * ((ages/max(ages))^4)
   results[which(results<0)] <- 0
   return(results)
}
mu_function = function(ages,params){
    results <-  params['q1'] * exp(-params['q2']*ages/max(ages)) + 
               params['q3'] +
               params['q4'] * (ages/max(ages)) +
               params['q5'] * ((ages/max(ages))^2) +
               params['q6'] * ((ages/max(ages))^3) +
               params['q7'] * ((ages/max(ages))^4)
    results[which(results<0)] <- 0
    return(results)
}

rho_function = function(params){
    return(rho) # rho does not depend on any of the parameters
}

param_values <- rep(NA, 14)
names(param_values) <- c(paste0("p", sequence(7)), paste0("q", sequence(7)))
param_guess <- rep(0, 14)
param_guess[c(3, 3+7)] <- 0.1 #the age independent rate

fit_param_yule_7 = fit_hbd_model_parametric(	tree, 
                                           param_values  = param_values,
                                           param_guess   = param_guess,
                                           param_min     = rep(0, 14),
                                           param_max     = rep(2, 14),
                                           param_scale   = 1, # all params are in the order of 1
                                           lambda        = lambda_function,
                                           mu            = mu_function,
                                           rho0          = rho_function,
                                           age_grid      = age_grid_param,
                                           condition     = "crown",
                                           Ntrials       = 10,    # perform 10 fitting trials
                                           Nthreads      = parallel::detectCores(),     # use 2 CPUs
                                           max_model_runtime = 1, # limit model evaluation to 1 second
                                           fit_control       = list(rel.tol=1e-6)
)
save(list=ls(), file="runall.rda")
# plot(range(age_grid_param), range(c(lambda_function(age_grid_param, fit_param_yule_7$param_fitted), mu_function(age_grid_param, fit_param_yule_7$param_fitted))), type="n", bty="n", ylab="rate")
# lines(age_grid_param, lambda_function(age_grid_param, fit_param_yule_7$param_fitted))
# lines(age_grid_param, mu_function(age_grid_param, fit_param_yule_7$param_fitted), lty="dotted")
# 
