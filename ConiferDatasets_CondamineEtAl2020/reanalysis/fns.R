
fit_env2 <- function (phylo, env_data, tot_time, f.lamb, f.mu, lamb_par,
    mu_par, df = NULL, f = 1, meth = "Nelder-Mead", cst.lamb = FALSE,
    cst.mu = FALSE, expo.lamb = FALSE, expo.mu = FALSE, fix.mu = FALSE,
    dt = 0, cond = "crown")
{
    if (tot_time > max(env_data[, 1])) {
        stop("The environmental data does not cover the time span of the phylogeny: either enter data that covers the full time span or run analyses on younger clades")
    }
    if (is.null(df)) {
        df <- smooth.spline(x = env_data[, 1], env_data[, 2])$df
    }
    spline_result <- pspline::sm.spline(env_data[, 1], env_data[, 2],
        df = df)
    env_func <- function(t) {
        predict(spline_result, t)
    }
    lower_bound_control <- 0.1
    upper_bound_control <- 0.1
    lower_bound <- min(env_data[, 1])
    upper_bound <- max(env_data[, 1])
    time_tabulated <- seq(from = lower_bound * (1 - lower_bound_control),
        to = upper_bound * (1 + upper_bound_control), length.out = 1 +
            1e+06)
    env_tabulated <- env_func(time_tabulated)
    env_func_tab <- function(t) {
        b <- upper_bound * (1 + upper_bound_control)
        a <- lower_bound * (1 - lower_bound_control)
        n <- length(env_tabulated) - 1
        index <- 1 + as.integer((t - a) * n/(b - a))
        return(env_tabulated[index])
    }
    f.lamb.env <- function(t, y) {
        f.lamb(t, env_func_tab(t), y)
    }
    f.mu.env <- function(t, y) {
        f.mu(t, env_func_tab(t), y)
    }
    res <- fit_bd2(phylo, tot_time, f.lamb.env, f.mu.env, lamb_par,
        mu_par, f, meth, cst.lamb, cst.mu, expo.lamb, expo.mu,
        fix.mu, dt, cond)
    # res$model <- "environmental birth death"
    # res$f.lamb <- function(t) {
    #     abs(f.lamb(t, env_func_tab(t), res$lamb_par))
    # }
    # if (fix.mu == FALSE) {
    #     res$f.mu <- function(t) {
    #         abs(f.mu(t, env_func_tab(t), res$mu_par))
    #     }
    # }
    # class(res) <- "fit.env"
    return(res)
}

fit_bd2 <- function (phylo, tot_time, f.lamb, f.mu, lamb_par, mu_par, f = 1,
    meth = "Nelder-Mead", cst.lamb = FALSE, cst.mu = FALSE, expo.lamb = FALSE,
    expo.mu = FALSE, fix.mu = FALSE, dt = 0, cond = "crown")
{
    if (!inherits(phylo, "phylo"))
        stop("object \"phylo\" is not of class \"phylo\"")
    nobs <- Ntip(phylo)
    if (fix.mu == FALSE) {
        init <- c(lamb_par, mu_par)
        p <- length(init)
        optimLH <- function(init) {
            lamb_par <- init[1:length(lamb_par)]
            mu_par <- init[(1 + length(lamb_par)):length(init)]
            f.lamb.par <- function(t) {
                abs(f.lamb(t, lamb_par))
            }
            f.mu.par <- function(t) {
                abs(f.mu(t, mu_par))
            }
            LH <- likelihood_bd(phylo, tot_time, f.lamb.par,
                f.mu.par, f, cst.lamb = cst.lamb, cst.mu = cst.mu,
                expo.lamb = expo.lamb, expo.mu = expo.mu, dt = dt,
                cond = cond)
            return(-LH)
        }
        temp <- suppressWarnings(optimLH(init))
        # lamb.par <- temp$par[1:length(lamb_par)]
        # mu.par <- temp$par[(1 + length(lamb_par)):length(init)]
        # f.lamb.par <- function(t) {
        #     abs(f.lamb(t, lamb.par))
        # }
        # f.mu.par <- function(t) {
        #     abs(f.mu(t, mu.par))
        # }
		res <- temp
        res <- list(model = "birth death", LH = -temp,
            aicc = 2 * (-temp) + 2 * p + (2 * p * (p + 1))/(nobs -
                p - 1), lamb_par = NA, mu_par = NA,
            f.lamb = NA, f.mu = NA)
    }
    else {
        init <- c(lamb_par)
        p <- length(init)
        optimLH <- function(init) {
            lamb_par <- init[1:length(lamb_par)]
            f.lamb.par <- function(t) {
                abs(f.lamb(t, lamb_par))
            }
            f.mu.par <- function(t) {
                abs(f.mu(t, mu_par))
            }
            LH <- likelihood_bd(phylo, tot_time, f.lamb.par,
                f.mu.par, f, cst.lamb = cst.lamb, cst.mu = TRUE,
                expo.lamb = expo.lamb, dt = dt, cond = cond)
            return(-LH)
        }
        temp <- suppressWarnings(optimLH(init))
        #lamb.par <- temp$par[1:length(lamb_par)]
        # f.lamb.par <- function(t) {
        #     abs(f.lamb(t, lamb.par))
        # }
        # f.mu.par <- function(t) {
        #     abs(f.mu(t, mu_par))
        # }
		res <- temp
        res <- list(model = "birth death", LH = -temp,
            aicc = 2 * (-temp) + 2 * p + (2 * p * (p + 1))/(nobs -
                p - 1), lamb_par = NA, mu_par = NA,
            f.lamb = NA, f.mu = NA)
    }
    class(res) <- "fit.bd"
    return(res)
}
