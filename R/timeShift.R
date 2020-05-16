
library(TreeSim)
library(TreePar)
library(nloptr)


######################################################################################################################################
######################################################################################################################################
### Likelihood function and the optimization function
######################################################################################################################################
######################################################################################################################################

GetLikelihoodShifts <- function(p, phy, times, ef.fixed=NULL){
    new.p <- exp(p)
    lambda <- new.p[1:length(times)]
    if(is.null(ef.fixed)){
        mu <- lambda * new.p[length(p)]
    }else{
        mu <- lambda * ef.fixed
    }
    frac <- rep(1, length(lambda))
    loglik <- LikShifts(getx(phy), t=times, lambda=lambda, mu=mu, sampling=frac)
    return(loglik)
}


TimeSliceModel <- function(phy, times, ef.fixed=NULL){
    opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"=100000, "ftol_rel"=.Machine$double.eps^.5)
    out <- nloptr(x0=log(c(rep(.1, length(times)),0.5)), eval_f=GetLikelihoodShifts, ub=log(c(rep(10, length(times)),0.99)), lb=rep(-21, length(times)+1), opts=opts, phy=phy, times=times, ef.fixed=ef.fixed)
    if(is.null(ef.fixed)){
        out <- nloptr(x0=log(c(rep(.1, length(times)),0.5)), eval_f=GetLikelihoodShifts, ub=log(c(rep(10, length(times)),0.99)), lb=rep(-21, length(times)+1), opts=opts, phy=phy, times=times, ef.fixed=ef.fixed)
        fit <- list(logLik=-out$objective, aic=(2*out$objective)+(2*length(out$solution)), times=times, lambda=exp(out$solution[1:length(times)]), mu=exp(out$solution[1:length(times)])*exp(out$solution[length(out$solution)]), ef=exp(out$solution[length(out$solution)]))
    }else{
        out <- nloptr(x0=log(c(rep(.1, length(times)))), eval_f=GetLikelihoodShifts, ub=log(c(rep(10, length(times)))), lb=rep(-21, length(times)), opts=opts, phy=phy, times=times, ef.fixed=ef.fixed)
        fit <- list(logLik=-out$objective, aic=(2*out$objective)+(2*length(out$solution)), times=times, lambda=exp(out$solution[1:length(times)]), mu=exp(out$solution[1:length(times)])*ef.fixed, ef=ef.fixed)
    }
    return(fit)
}



######################################################################################################################################
######################################################################################################################################
### Running analyses and Plotting code
######################################################################################################################################
######################################################################################################################################

tree <- read.tree("../data/whales_Steemanetal2009.tre")
times <- seq(0,30, by=10)

pdf("time.shift.10myrbin.pdf", height=3, width=7)
par(mfcol=c(1,3))

#### PANEL A ####
pp <- TimeSliceModel(phy=tree, times=times, ef.fixed=NULL)
times.plot <- seq(0.1,max(branching.times(tree)),0.1)
lambdas.plot <- c()
mus.plot <- c()
for(index in 1:length(times.plot)){
    lambdas.plot <- c(lambdas.plot, pp$lambda[max(which(times.plot[index] > times))])
    mus.plot <- c(mus.plot, pp$ef*pp$lambda[max(which(times.plot[index] > times))])
}
plot(times.plot, lambdas.plot, col=0, ylim=c(0,1.2), xlim=c(35,0), xlab="age (Myr)", ylab="rate (1/Myr)", main=expression(epsilon~est.))
lines(times.plot, lambdas.plot, lty=2, col="blue")
lines(times.plot, mus.plot, lty=3, col="darkorange")
mtext(paste("logLik=", round(pp$logLik,2)), cex=.75)


#### PANEL B ####
pp <- TimeSliceModel(phy=tree, times=times, ef.fixed=0)
times.plot <- seq(0.1,max(branching.times(tree)),0.1)
lambdas.plot <- c()
mus.plot <- c()
for(index in 1:length(times.plot)){
    lambdas.plot <- c(lambdas.plot, pp$lambda[max(which(times.plot[index] > times))])
    mus.plot <- c(mus.plot, pp$ef*pp$lambda[max(which(times.plot[index] > times))])
}
plot(times.plot, lambdas.plot, col=0, ylim=c(0,1.2), xlim=c(35,0), xlab="age (Myr)", ylab="rate (1/Myr)", main=expression(epsilon==0))
lines(times.plot, lambdas.plot, lty=2, col="blue")
lines(times.plot, mus.plot, lty=3, col="darkorange")
mtext(paste("logLik=", round(pp$logLik,2)), cex=.75)


#### PANEL C ####
pp <- TimeSliceModel(phy=tree, times=times, ef.fixed=0.9)
times.plot <- seq(0.1,max(branching.times(tree)),0.1)
lambdas.plot <- c()
mus.plot <- c()
for(index in 1:length(times.plot)){
    lambdas.plot <- c(lambdas.plot, pp$lambda[max(which(times.plot[index] > times))])
    mus.plot <- c(mus.plot, pp$ef*pp$lambda[max(which(times.plot[index] > times))])
}
plot(times.plot, lambdas.plot, col=0, ylim=c(0,1.2), xlim=c(35,0), xlab="age (Myr)", ylab="rate (1/Myr)", main=expression(epsilon==0.9))
lines(times.plot, lambdas.plot, lty=2, col="blue")
lines(times.plot, mus.plot, lty=3, col="darkorange")
mtext(paste("logLik=", round(pp$logLik,2)), cex=.75)

dev.off()


######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

