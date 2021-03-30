GetX <- function(phy) {
	return(TreeSim::getx(phy))
}

#magnitude is a fraction: 0.5 mean shrink the focal interval by half, 1.4 increase by 40%
ChangeX <- function(magnitude, interval, X) {
	offsets <- c(X, X[length(X)])-c(0,X)
	new_length <- offsets[interval]*magnitude
	X[seq(from=interval, to=length(X), by=1)] <- X[seq(from=interval, to=length(X), by=1)] - offsets[interval] + new_length
	return(X)
}

Optimize <- function(X, phy) {
	opts <- list("algorithm" = "NLOPT_LN_SBPLX", "maxeval" = 100000, "ftol_rel" = .Machine$double.eps^.5)
	condition.type="survival"
	result <- nloptr(x0=log(c(.1,.05)), eval_f=GetLikelihood, ub=log(c(10,.99)), lb=c(-6,-6), opts=opts, tree=phy, X=X, condition.type=condition.type, verbose=FALSE)
	
	
	#result <- bd.shifts.optim(X, sampling=1, 0.2, 0.1)[[1]][[1]][[1]]
	return(c(lnL=result$objective, lambda=exp(result$solution[1]), mu=exp(result$solution[2])))
}

TryAllIntervals <- function(magnitude, X, phy) {
	results <- Optimize(X, phy)
	results$interval <- NA
	results$magnitude <- NA
	for(interval in sequence(length(X)-1)) {
		result <- Optimize(ChangeX(magnitude, interval, X), phy)
		result$interval <- interval
		result$magnitude <- magnitude
		results <- rbind(results, result)
	}
	rownames(results) <- NULL
	results <- as.data.frame(results)
	results <- rbind(results, rep(NA, ncol(results)))
	#results[nrow(results)+1,sequence(ncol(results))] <- NA
	results$X <- c(NA, unname(X))
	results$Ntip <- c(ape::Ntip(phy), seq(from=ape::Ntip(phy), to=2, by=-1))
	return(results)
}



#### Code from earlier work

p0 <- function(t,l,m,rho){
    1- rho*(l-m)/(rho*l+(l*(1-rho)-m)*exp(-(l-m)*t))
}
​
​
p1 <- function(t,l,m,rho){
    rho*(l-m)^2 * exp(-(l-m)*t)/(rho*l+(l*(1-rho)-m)*exp(-(l-m)*t))^2
}
​
​
qhelp <- function(t,l,m,rho){
    rho*l*(1-exp(-(l-m)*t))/(rho*l+(l*(1-rho)-m)*exp(-(l-m)*t))
}
​
​
stadler.pn <- function(birth.rate, death.rate, time, n){
    rho=1
    l <- birth.rate
    m <- death.rate
    t <- time
    part1 <- p1(t,l,m,rho)
    part2 <- qhelp(t,l,m,rho)
    pn.t <- part1 * part2^(n-1)
    return(pn.t)
}
​
​
AscertainBiasNew <- function(k,l,m,t,rho, end=10000000, get.partials=FALSE) {
    part1 <- p1(t,l,m,rho)
    part2 <- qhelp(t,l,m,rho)
    res.res <- numeric(10000000)
    for(n.index in k:end){
        res.res[n.index] <- (n.index-1)*(part1^2)*(part2^(n.index-2))
    }
    if(get.partials == TRUE){
        return(res.res)
    }else{
        prob.nplus <- sum(res.res)
        return(prob.nplus)
    }
}
​
​
AscertainBiasOld <- function(k, l, m, t, rho) {
    net.diver.rate <- l - m
    #Magallon and Sanderson 2001 -- Eq. 2a:
    exprt <- exp(net.diver.rate * t)
    beta <- (exprt - 1) / (exprt - (m/l))
    #Magallon and Sanderson 2001 -- Eq. 2b:
    alpha <- (m/l) * beta
    #Magallon and Sanderson 2001 -- Eq. 10a:
    probNgeNtax <- beta^(k-1)
    survival.prob <- (1 - alpha)^2
    probNgeNtax <- (beta^(k-2))*(k*(1 - alpha - beta + alpha*beta) + alpha + 2*beta-1)/(1 - alpha + 2*alpha*beta)
    combined.prob <- probNgeNtax * survival.prob
    return(combined.prob)
}
​
​
#conditional probability based on survival only
survival.conditional.p <- function(time, turn, eps){
    l <- turn/(1+eps)
    m <- (turn * eps) / (1 + eps)
    conditional.p <- (1-p0(time,l,m,rho=1))^2
    return(conditional.p)
}
​
#conditional probability based on survival and n taxa
survival.n.conditional.p <- function(time, turn, eps, n){
    l <- turn/(1+eps)
    m <- (turn * eps) / (1 + eps)
    part1 <- p1(time,l,m,rho=1)
    part2 <- qhelp(time,l,m,rho=1)
    conditional.p <- (n-1)*(part1^2)*(part2)^(n-2)
    return(conditional.p)
}
​
​
GetLikelihood <- function(p, X, tree, condition.type="survival", verbose=FALSE){
    p <- exp(p)
    n <- Ntip(tree)
    rho <- 1
    l <- p[1] / (1 + p[2])
    m <- (p[1] * p[2]) / (1 + p[2])
    #x <- getx(tree)
    x <- sort(X,decreasing=TRUE)
    t <- x
​
    lik <- 0
    for (i in 2:length(t)){
        lik <- lik+log(l*p1(t[i],l,m,rho))
    }
    if(condition.type == "survival"){
        lik.unconditioned <- lik + 2 * log(p1(t[1],l,m,rho))
        condition.p <- 2 * log(1-p0(t[1],l,m,rho))
    }
    if(condition.type == "exactlyn"){
        lik.unconditioned <- lik + 2*log(p1(t[1],l,m,rho))
        part1 <- p1(t[1],l,m,rho)
        part2 <- qhelp(t[1],l,m,rho)
        condition.p <- log((n-1)*(part1^2)*(part2)^(n-2))
    }
    if(condition.type == "morethann"){
        lik.unconditioned <- lik + 2*log(p1(t[1],l,m,rho))
        condition.p <- log(AscertainBiasNew(n,l,m,t[1],rho=1))
    }
    if(condition.type == "magsan"){
        lik.unconditioned <- lik + 2*log(p1(t[1],l,m,rho))
        condition.p <- log(AscertainBiasOld(n,l,m,t[1],rho=1))
    }
    if(condition.type == "exactlyn.stem0"){
        lik.unconditioned <- lik + log(l*p1(t[1],l,m,rho))
        condition.p <- log(n * (p1(t[1],l,m,rho) / (1-p0(t[1],l,m,rho))))
    }
    if(condition.type == "exactlyn.stem1"){
        lik.unconditioned <- lik + log(l*p1(t[i],l,m,rho))
        condition.p <- log((qhelp(t=t[1], l=l, m=m, rho=rho))^(1-n))
    }
    if(condition.type == "none") {
      lik.unconditioned <- lik + log(l*p1(t[i],l,m,rho))
      condition.p <- log(1)
    }
    if(verbose==TRUE){
        print(paste("unconditional prob", lik.unconditioned))
        print(paste("conditional prob", condition.p))
    }
​
    if(condition.type == "exactlyn.stem0" | condition.type == "exactlyn.stem1" | condition.type == "exactlyn" | condition.type == "magsan" | condition.type == "morethann"){
        #tmp <- c(l+m, m/l, lik.unconditioned, log(condition.p))
        lik.conditioned.both.check <- lik.unconditioned - condition.p
        return(-lik.conditioned.both.check)
    }else{
        lik.conditioned.both.check <- lik.unconditioned - condition.p
        return(-lik.conditioned.both.check)
    }
}
​
​