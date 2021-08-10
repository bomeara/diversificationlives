library(dplyr)
library(ape)
library(RPANDA)
library(ggplot2)
library(ggpubr)
#devtools::install_github("willgearty/deeptime")
library(deeptime)
library(reshape2)
library(truncnorm)
library("future.apply")
library(zoo)
fossils <- read.delim("../Fossil record/Fossil record of conifers.txt")
fossil_by_clade <- fossils %>% group_by(Species) %>% mutate(min=min(min_age), max=max(max_age)) %>% select(Species, min, max) %>% distinct()
phy <- ape::read.nexus("../Phylogenies/Phylogeny of conifers 2018.tre")
angios <- read.delim("../Paleo-environmental variables/Angiosperms.txt")
f.constant <- function(t,x,y){y[1]}
f.linear <- function(t,x,y){y[1] + (y[2] * x)}
f.exp <-function(t,x,y){y[1] * exp(y[2] * x)}
dt <- 1e-3

tmp <- tempfile()
download.file("https://datasets.imdbws.com/title.episode.tsv.gz",tmp) #downloaded Nov 15, 2020
imdb_episodes <- read.table(tmp, header=TRUE)
download.file("https://datasets.imdbws.com/title.ratings.tsv.gz",tmp)  #downloaded Nov 15, 2020
imdb_ratings <- read.table(tmp, header=TRUE)
rownames(imdb_ratings) <- imdb_ratings$tconst

series <- c(Simpsons = "tt0096697", My_Little_Pony = "tt1751105", Star_Trek_DS9 = "tt0106145", Daily_Show = "tt0115147")

GetTVRankings <- function(parent_tconst="tt0096697") {
	tvseries_episodes <- subset(imdb_episodes, parentTconst==parent_tconst)
	tvseries_episodes$seasonNumber <- as.numeric(tvseries_episodes$seasonNumber)
	tvseries_episodes$episodeNumber <- as.numeric(tvseries_episodes$episodeNumber)
	tvseries_episodes <- tvseries_episodes[!is.na(tvseries_episodes$seasonNumber),]
	tvseries_episodes <- tvseries_episodes[!is.na(tvseries_episodes$episodeNumber),]
	tvseries_episodes <- tvseries_episodes[order(tvseries_episodes$seasonNumber, tvseries_episodes$episodeNumber), ]
	tvseries_episodes$averageRating <- NA
	for (i in sequence(nrow(tvseries_episodes))) {
		try(tvseries_episodes$averageRating[i] <- imdb_ratings$averageRating[which(imdb_ratings$tconst==tvseries_episodes$tconst[i])])
	}
	tvseries_episodes <- tvseries_episodes[!is.na(tvseries_episodes$averageRating),]
	tvseries_episodes$order <- sequence(nrow(tvseries_episodes))
	tvseries_episodes$averageProportion <- tvseries_episodes$averageRating / 10
	tvseries_episodes$age_from_recent <- seq(from=nrow(tvseries_episodes)-1, by=-1, length=nrow(tvseries_episodes))
	tvseries_as_predictor <- data.frame(Age=angios$Age, TV_seriesRating=approx(x=seq(from=0, to=1, length.out = nrow(tvseries_episodes)), y=tvseries_episodes$averageProportion, xout=angios$Age/max(angios$Age))$y) #so we go from the present to the past, interpolating ratings at the same point in the series as in our angiosperm fossils
	return(tvseries_as_predictor)
}

ModelPerSeries <- function(chosen_series) {
	rankings <- GetTVRankings(chosen_series)
	series_fit <- fit_env(phy, env_data=rankings, tot_time=max(node.age(phy)$ages), f.lamb=f.constant, f.mu=f.exp, lamb_par=c(0.2161), mu_par= c(0.2023, 0.0772), f=0.9, dt=dt)
	return(series_fit)
}


TVSeries_rankings <- parallel::mclapply(series, GetTVRankings, mc.cores=min(4, length(series)))

TVSeries_models <- parallel::mclapply(series, ModelPerSeries, mc.cores=min(4, length(series)))

# TVSeries_model_1 <- ModelPerSeries(series[1])
# TVSeries_model_2 <- ModelPerSeries(series[2])
# TVSeries_model_3 <- ModelPerSeries(series[3])
# TVSeries_model_4 <- ModelPerSeries(series[4])

AllPredictors <- rbind(
	data.frame(Age=TVSeries_rankings[[1]]$Age, Predictor=TVSeries_rankings[[1]]$TV_seriesRating, Set=names(series)[1]),
	data.frame(Age=TVSeries_rankings[[2]]$Age, Predictor=TVSeries_rankings[[2]]$TV_seriesRating, Set=names(series)[2]),
	data.frame(Age=TVSeries_rankings[[3]]$Age, Predictor=TVSeries_rankings[[3]]$TV_seriesRating, Set=names(series)[3]),
	data.frame(Age=TVSeries_rankings[[4]]$Age, Predictor=TVSeries_rankings[[4]]$TV_seriesRating, Set=names(series)[4]),
	data.frame(Age=angios$Age, Predictor=angios$Angiosperms, Set="Angiosperms")

)
pdf(file="PredictorsPlot.pdf")
all_predictors_plot <- ggplot(AllPredictors, mapping=aes(x=Age, y=Predictor, group=Set)) + geom_line(aes(colour=Set)) + scale_x_reverse() + coord_geo() 
print(all_predictors_plot)
dev.off()


# just to check getting similar results
BcstDAngioVar_EXPO <- fit_env(phy, env_data=angios, tot_time=max(node.age(phy)$ages), f.lamb=f.constant, f.mu=f.exp, lamb_par=c(0.2161), mu_par= c(0.2023, 0.0772), f=0.9, dt=dt)
BcstDAngioVar_LIN <- fit_env(phy, env_data=angios, tot_time=max(node.age(phy)$ages), f.lamb=f.constant, f.mu=f.linear, lamb_par=c(0.2161), mu_par= c(0.2022, 0.0162), f=0.9, dt=dt)
BAngioVarDcst_LIN <- fit_env(phy, env_data=angios, tot_time=max(node.age(phy)$ages), f.lamb=f.linear, f.mu=f.constant, lamb_par=c(0.2305, -0.0159), mu_par= c(0.2171), f=0.9, dt=dt)
BAngioVarDcst_EXPO <- fit_env(phy, env_data=angios, tot_time=max(node.age(phy)$ages), f.lamb=f.exp, f.mu=f.constant, lamb_par=c(0.2305, -0.0711), mu_par= c(0.2171), f=0.9, dt=dt)

#BcstDTV_seriesVar_EXPO <- fit_env(phy, env_data=tvseries_as_predictor, tot_time=max(node.age(phy)$ages), f.lamb=f.constant, f.mu=f.exp, lamb_par=c(0.2161), mu_par= c(0.2023, 0.0772), f=0.9, dt=dt)



models <- list(BcstDAngioVar_EXPO=BcstDAngioVar_EXPO, BcstDAngioVar_LIN=BcstDAngioVar_LIN, BAngioVarDcst_LIN=BAngioVarDcst_LIN, BAngioVarDcst_EXPO=BAngioVarDcst_EXPO)

print(unlist(lapply(models, "[[", "aicc")) - min(unlist(lapply(models, "[[", "aicc"))))
print(lapply(models, "[[", "lamb_par"))
print(lapply(models, "[[", "mu_par"))
print(lapply(models, "[[", "LH"))

ltt_data <- as.data.frame(ape::ltt.plot.coords(phy))
ltt_data$age <- abs(ltt_data$time)
midpoint <- abs(ltt_data[round(Ntip(phy)/2),1])

all_times <- data.frame(age=seq(from=0, to=max(abs(ltt_data[,1])), by=0.1), surviving=NA, fossil_origins=NA, fossil_extinctions=NA)
for (i in sequence(nrow(all_times))) {
	all_times$surviving[i] <- ltt_data$N[tail(which(ltt_data$age>=all_times$age[i]),1)]
	all_times$fossil_origins[i] <- sum(fossil_by_clade$max>all_times$age[i])
	all_times$fossil_extinctions[i] <- sum(fossil_by_clade$min>all_times$age[i])
}

all_times_tall  <- melt(all_times, id = "age")

# angios_not_decreasing <- angios
# angios_not_decreasing$Angiosperms[1:44] <- 1
# angios_linear_decrease <- data.frame(Age=angios$Age, Angiosperms=angios$Age/max(angios$Age) )
# angios_linear_increase <- data.frame(Age=angios$Age, Angiosperms=1 - angios$Age/max(angios$Age) )
angios_delayed_rise_50 <- angios
angios_delayed_rise_50$Angiosperms[which(angios_delayed_rise_50$Age>50)] <- 0
angios_delayed_rise_10 <- angios
angios_delayed_rise_10$Angiosperms[which(angios_delayed_rise_10$Age>10)] <- 0
angios_delayed_rise_1 <- angios
angios_delayed_rise_1$Angiosperms[which(angios_delayed_rise_1$Age>1)] <- 0
angios_instant_rise_data_mid <- angios
angios_instant_rise_data_mid$Angiosperms[which(angios_instant_rise_data_mid$Age>midpoint)] <- 0
angios_instant_rise_data_mid$Angiosperms[which(angios_instant_rise_data_mid$Age<=midpoint)] <- 0.86

angios_instant_rise_20 <- angios
angios_instant_rise_20$Angiosperms[which(angios_instant_rise_data_mid$Age>20)] <- 0
angios_instant_rise_20$Angiosperms[which(angios_instant_rise_data_mid$Age<=20)] <- 0.86


angios_delayed_rise_max42 <- angios
angios_delayed_rise_max42$Angiosperms[which(angios_delayed_rise_max42$Age>42)] <- 0

angios_delayed_rise_max66 <- angios
angios_delayed_rise_max66$Angiosperms[which(angios_delayed_rise_max66$Age>66)] <- 0

angios_early_rise_max66 <- angios
angios_early_rise_max66$Angiosperms[which(angios_early_rise_max66$Age>66)] <- 0.86

#angios_delayed_rise_min42 <- angios
#angios_delayed_rise_min42$Angiosperms[which(angios_delayed_rise_min42$Age<=42)] <- 0



# BcstDAngioVar_EXPO_angios_not_decreasing <- fit_env(phy, env_data=angios_not_decreasing, tot_time=max(node.age(phy)$ages), f.lamb=f.constant, f.mu=f.exp, lamb_par=c(0.2161), mu_par= c(0.2023, 0.0772), f=0.9, dt=dt)
# BcstDAngioVar_EXPO_angios_linear_decrease <- fit_env(phy, env_data=angios_linear_decrease, tot_time=max(node.age(phy)$ages), f.lamb=f.constant, f.mu=f.exp, lamb_par=c(0.2161), mu_par= c(0.2023, 0.0772), f=0.9, dt=dt)
# BcstDAngioVar_EXPO_angios_linear_increase <- fit_env(phy, env_data=angios_linear_increase, tot_time=max(node.age(phy)$ages), f.lamb=f.constant, f.mu=f.exp, lamb_par=c(0.2161), mu_par= c(0.2023, 0.0772), f=0.9, dt=dt)
# BcstDAngioVar_EXPO_angios_delayed_rise_50 <- fit_env(phy, env_data=angios_delayed_rise_50, tot_time=max(node.age(phy)$ages), f.lamb=f.constant, f.mu=f.exp, lamb_par=c(0.2161), mu_par= c(0.2023, 0.0772), f=0.9, dt=dt)
# BcstDAngioVar_EXPO_angios_delayed_rise_10 <- fit_env(phy, env_data=angios_delayed_rise_10, tot_time=max(node.age(phy)$ages), f.lamb=f.constant, f.mu=f.exp, lamb_par=c(0.2161), mu_par= c(0.2023, 0.0772), f=0.9, dt=dt)
# BcstDAngioVar_EXPO_angios_delayed_rise_1 <- fit_env(phy, env_data=angios_delayed_rise_1, tot_time=max(node.age(phy)$ages), f.lamb=f.constant, f.mu=f.exp, lamb_par=c(0.2161), mu_par= c(0.2023, 0.0772), f=0.9, dt=dt)
BcstDAngioVar_LIN_angios_delayed_rise_max42 <- fit_env(phy, env_data=angios_delayed_rise_max42, tot_time=max(node.age(phy)$ages), f.lamb=f.constant, f.mu=f.linear, lamb_par=c(0.2161), mu_par= c(0.2023, 0.0772), f=0.9, dt=dt)
BcstDAngioVar_LIN_angios_delayed_rise_max66 <- fit_env(phy, env_data=angios_delayed_rise_max66, tot_time=max(node.age(phy)$ages), f.lamb=f.constant, f.mu=f.linear, lamb_par=c(0.2161), mu_par= c(0.2023, 0.0772), f=0.9, dt=dt)
BcstDAngioVar_LIN_angios_early_rise_max66 <- fit_env(phy, env_data=angios_early_rise_max66, tot_time=max(node.age(phy)$ages), f.lamb=f.constant, f.mu=f.linear, lamb_par=c(0.2161), mu_par= c(0.2023, 0.0772), f=0.9, dt=dt)

threshold_models <- list()
preK_angios <- seq(from=.0, to=1, by=0.1)
for(preK_index in seq_along(preK_angios)){
	angios_66_slice <- angios
	angios_66_slice$Angiosperms[which(angios_66_slice$Age>66)] <- preK_angios[preK_index]
	threshold_models[[preK_index]] <- fit_env(phy, env_data=angios_66_slice, tot_time=max(node.age(phy)$ages), f.lamb=f.constant, f.mu=f.linear, lamb_par=c(0.2161), mu_par= c(0.2023, 0.0772), f=0.9, dt=dt)
	print(c(preK_angios[preK_index], threshold_models[[preK_index]]$aicc))
}



mu_par <- c(0.1, 0.01)
f.mu <- function(t){lamb_par[1] * exp(lamb_par[2] * t)}
f.lamb <- function(t){BcstDAngioVar_LIN$lamb_par[1]}

f.mu<- function(t){BcstDAngioVar_LIN$mu_par[1] + exp(BcstDAngioVar_LIN$mu_par[2] * t)}
f.lamb <- function(t){BcstDAngioVar_LIN$lamb_par[1]}

fitted <- likelihood_bd(phy, tot_time=max(node.age(phy)$ages), f.lamb=f.lamb, f.mu=f.mu, f=0.9, cst.lamb=TRUE, cst.mu=FALSE, expo.lamb=FALSE, expo.mu=FALSE, dt=dt)

source('fns.R')
#BcstDAngioVar_LIN2 <- fit_env2(phy, env_data=angios, tot_time=max(node.age(phy)$ages), f.lamb=f.constant, f.mu=f.linear, lamb_par=BcstDAngioVar_LIN$lamb_par, mu_par= BcstDAngioVar_LIN$mu_par, f=0.9, dt=dt)
current_aicc <- 0
while(abs(current_aicc-BcstDAngioVar_LIN2$aicc)>0.5) {
	angios_wiggle <- angios
	#angios_wiggle$Angiosperms <- abs(angios_wiggle$Angiosperms+runif(nrow(angios_wiggle), min=-0.2, max=0.2))
	random_vals <- c(0.86, runif(runif(1, 2, 15), min=0, max=1))
	angios_wiggle$Angiosperms <- spline(x=c(0, sort(sample(angios$Age[2:(nrow(angios)-1)], size=length(random_vals)-2, replace=FALSE)), max(angios$Age)), y=random_vals, xout=angios_wiggle$Age)$y
	#angios_wiggle$Angiosperms <- runif(nrow(angios_wiggle), min=0, max=1)

	BcstDAngioVar_LIN3 <- fit_env2(phy, env_data=angios_wiggle, tot_time=max(node.age(phy)$ages), f.lamb=f.constant, f.mu=f.linear, lamb_par=BcstDAngioVar_LIN$lamb_par, mu_par= BcstDAngioVar_LIN$mu_par, f=0.9, dt=dt)
	print(BcstDAngioVar_LIN2)
	print(BcstDAngioVar_LIN3)
	current_aicc <- BcstDAngioVar_LIN3$aicc
	plot(angios_wiggle$Age, angios_wiggle$Angiosperms, col="red", type="l", ylim=c(0, max(c(angios$Angiosperms, angios_wiggle$Angiosperms))), main=paste0("∆AICc = ", current_aicc-BcstDAngioVar_LIN2$aicc))
	lines(angios$Age, angios$Angiosperms, col="blue")
}

phylo=phy
cst.lamb = FALSE
 cst.mu = FALSE
 expo.lamb = FALSE
 expo.mu = FALSE
  dt = 1e-3
cond = "crown"

    
    nbtips <- Ntip(phylo)
    log_indLikelihood <- c()
    from_past <- cbind(phylo$edge, node.age(phylo)$ages)
    ages <- rbind(from_past[, 2:3], c(nbtips + 1, 0))
    ages <- ages[order(ages[, 1]), ]
    age <- max(ages[, 2])
    for (j in 1:(nbtips - 1)) {
        node <- (nbtips + j)
        edges <- phylo$edge[phylo$edge[, 1] == node, ]
        tj <- age - ages[edges[1, 1], 2]
        Psi_timevar_errap_tj <- RPANDA:::.Psi(0, tj, f.lamb, f.mu, f,
            cst.lamb = cst.lamb, cst.mu = cst.mu, expo.lamb = expo.lamb,
            expo.mu = expo.mu, dt = dt)
			print(c(tj, Psi_timevar_errap_tj))
        log_lik_tj <- log(f.lamb(tj)) + log(Psi_timevar_errap_tj)
        log_indLikelihood <- c(log_indLikelihood, log_lik_tj)
    }
    log_indLikelihood <- c(log_indLikelihood, log(RPANDA:::.Psi(0, tot_time,
        f.lamb, f.mu, f, cst.lamb = cst.lamb, cst.mu = cst.mu,
        expo.lamb = expo.lamb, expo.mu = expo.mu, dt = dt)))
    log_data_lik <- sum(log_indLikelihood) + nbtips * log(f)
    if (cond == FALSE) {
        log_final_lik <- log_data_lik
    }
    else if (cond == "stem") {
        Phi <- RPANDA:::.Phi(tot_time, f.lamb, f.mu, f, cst.lamb = cst.lamb,
            cst.mu = cst.mu, expo.lamb = expo.lamb, expo.mu = expo.mu,
            dt = dt)
        log_final_lik <- log_data_lik - log(1 - Phi)
    }
    else if (cond == "crown") {
        Phi <- RPANDA:::.Phi(tot_time, f.lamb, f.mu, f, cst.lamb = cst.lamb,
            cst.mu = cst.mu, expo.lamb = expo.lamb, expo.mu = expo.mu,
            dt = dt)
        log_final_lik <- log_data_lik - log(f.lamb(tot_time)) -
            2 * log(1 - Phi)
    }



#BcstDAngioVar_EXPO_angios_instant_rise_data_mid <- fit_env(phy, env_data=angios_instant_rise_data_mid, tot_time=max(node.age(phy)$ages), f.lamb=f.constant, f.mu=f.exp, lamb_par=c(0.2161), mu_par= c(0.2023, 0.0772), f=0.9, dt=dt)
BcstDAngioVar_LIN_angios_instant_rise_20 <- fit_env(phy, env_data=angios_instant_rise_20, tot_time=max(node.age(phy)$ages), f.lamb=f.constant, f.mu=f.linear, lamb_par=c(0.2161), mu_par= c(0.2023, 0.0772), f=0.9, dt=dt)


angios_models <- list(BcstDAngioVar_LIN_angios_true=BcstDAngioVar_LIN, BcstDAngioVar_LIN_angios_delayed_rise_max42=BcstDAngioVar_LIN_angios_delayed_rise_max42, 
BcstDAngioVar_LIN_angios_delayed_rise_max66=BcstDAngioVar_LIN_angios_delayed_rise_max66,
BcstDAngioVar_LIN_angios_early_rise_max66, BcstDAngioVar_LIN_angios_early_rise_max66
BcstDAngioVar_LIN_angios_instant_rise_20=BcstDAngioVar_LIN_angios_instant_rise_20)

print(unlist(lapply(angios_models, "[[", "aicc")) - min(unlist(lapply(angios_models, "[[", "aicc"))))
print(lapply(angios_models, "[[", "lamb_par"))
print(lapply(angios_models, "[[", "mu_par"))
print(lapply(angios_models, "[[", "LH"))


all_models <- c(list(BcstDAngioVar_EXPO=BcstDAngioVar_EXPO, BcstDAngioVar_LIN=BcstDAngioVar_LIN, BAngioVarDcst_LIN=BAngioVarDcst_LIN, BAngioVarDcst_EXPO=BAngioVarDcst_EXPO, BcstDAngioVar_EXPO_angios_true=BcstDAngioVar_EXPO, BcstDAngioVar_LIN_angios_delayed_rise_max42=BcstDAngioVar_LIN_angios_delayed_rise_max42, 
BcstDAngioVar_LIN_angios_delayed_rise_max66=BcstDAngioVar_LIN_angios_delayed_rise_max66,
BcstDAngioVar_LIN_angios_instant_rise_20=BcstDAngioVar_LIN_angios_instant_rise_20, BcstDAngioVar_LIN_angios_early_rise_max66=BcstDAngioVar_LIN_angios_early_rise_max66),TVSeries_models)

deltaAICc <- unlist(lapply(all_models, "[[", "aicc")) - min(unlist(lapply(all_models, "[[", "aicc")))
print(deltaAICc)
print(lapply(all_models, "[[", "lamb_par"))
print(lapply(all_models, "[[", "mu_par"))
print(lapply(all_models, "[[", "LH"))

all_models_df <- data.frame(name=names(deltaAICc), deltaAICc=unname(unlist(deltaAICc)), likelihood=unname(unlist(lapply(all_models, "[[", "LH"))), aicc=unname(unlist(lapply(all_models, "[[", "aicc"))))
lambda_par_df <- t(plyr::rbind.fill(as.data.frame(lapply(all_models, "[[", "lamb_par"))))
mu_par_df <- t(plyr::rbind.fill(as.data.frame(lapply(all_models, "[[", "mu_par"))))
all_models_df$lambda_1 <- NA
all_models_df$lambda_2 <- NA
all_models_df$mu_1 <- NA
all_models_df$mu_2 <- NA
for (i in sequence(nrow(lambda_par_df))) {
  all_models_df$lambda_1[i] <- lambda_par_df[i,1]
  if(lambda_par_df[i,1]!=lambda_par_df[i,2]) {
    all_models_df$lambda_2[i] <- lambda_par_df[i,2]
  }
  all_models_df$mu_1[i] <- mu_par_df[i,1]
  if(mu_par_df[i,1]!=mu_par_df[i,2]) {
    all_models_df$mu_2[i] <- mu_par_df[i,2]
  }
}

all_models_df <- all_models_df[order(all_models_df$deltaAICc),]
write.csv(all_models_df, file="AllModels.csv")

mu_original <- f.exp(angios$Age, angios$Angiosperms, BcstDAngioVar_EXPO$mu_par)
mu_max42 <- f.exp(angios_delayed_rise_max42$Age, angios_delayed_rise_max42$Angiosperms, BcstDAngioVar_EXPO_angios_delayed_rise_max42$mu_par)
mu_instant_rise_20 <- f.exp(angios_instant_rise_20$Age, angios_instant_rise_20$Angiosperms, BcstDAngioVar_EXPO_angios_instant_rise_20$mu_par)

lambda_original <- rep(f.constant(angios$Age, angios$Angiosperms, BcstDAngioVar_EXPO$lamb_par), length(mu_original))
lambda_max42 <- rep(f.constant(angios_delayed_rise_max42$Age, angios_delayed_rise_max42$Angiosperms, BcstDAngioVar_EXPO_angios_delayed_rise_max42$lamb_par), length(mu_max42))
lambda_instant_rise_20 <- rep(f.constant(angios_instant_rise_20$Age, angios_instant_rise_20$Angiosperms, BcstDAngioVar_EXPO_angios_instant_rise_20$lamb_par), length(mu_instant_rise_20))


div_data <- data.frame(mu=c(mu_original, mu_max42, mu_instant_rise_20), lambda=c(lambda_original, lambda_max42, lambda_instant_rise_20), angiosperms=c(angios$Angiosperms, angios_delayed_rise_max42$Angiosperms, angios_instant_rise_20$Angiosperms),  age=rep(angios$Age, 3), category=c(rep("All", nrow(angios)), rep("No flowers before 42 MYA", nrow(angios)), rep("Instant domination 20 MYA", nrow(angios)))) 



mu_plot <- ggplot(div_data, mapping=aes(x=age, y=mu, group=category)) + geom_line(aes(colour=category)) + scale_x_reverse() + coord_geo()
lambda_plot <- ggplot(div_data, mapping=aes(x=age, y=lambda, group=category)) + geom_line(aes(colour=category)) + scale_x_reverse() + coord_geo()
angio_plot <- ggplot(div_data, mapping=aes(x=age, y=angiosperms, group=category)) + geom_line(aes(colour=category)) + scale_x_reverse() + coord_geo()

ltt_plot <- ggplot(all_times_tall, mapping=aes(x=age, y=value, group=variable)) + geom_line(aes(colour=variable)) + scale_x_reverse() + coord_geo() + scale_y_log10()

pdf(file="LambaMuAngio.pdf", width=15, height=15)
print(ggarrange(lambda_plot, mu_plot, angio_plot, ltt_plot, ncol=2, nrow=2))
dev.off
rise_of_angios <- (apply(subset(all_times, age>=80 & age<=150), 2, range))
rise_of_angios_diff <- diff(rise_of_angios)
print("during the rise of angiosperms")
print(rise_of_angios)
print(rise_of_angios_diff)

mid_cretaceous <- (apply(subset(all_times, age>=100 & age<=110), 2, range))
mid_cretaceous_diff <- diff(mid_cretaceous)
print("during the mid cretaceous supposed extinction pulse")
print(mid_cretaceous)
print(mid_cretaceous_diff)

try_random <- function(angios) {
	sd <- runif(1, 0, 1)
	k <- sample(5)[1]
	angios_random <- angios
	angios_random$Angiosperms <- rtruncnorm(n=length(angios_random$Angiosperms), a=0, b=1, mean=angios_random$Angiosperms, sd=sd)
	values <- rollmean(angios_random$Angiosperms, k=k, fill="extend")
	values[c(1,2)] <- angios_random$Angiosperms[c(1,2)]
	values[c(0,-1)+length(values)] <- 0
	angios_random$Angiosperms <- values
	print(angios_random$Angiosperms)
	result <- fit_env(phy, env_data=angios_random, tot_time=max(node.age(phy)$ages), f.lamb=f.constant, f.mu=f.exp, lamb_par=c(0.2161), mu_par= c(0.2023, 0.0772), f=0.9, dt=dt)
	return(list(angios_random=angios_random, result=result))
}

plan(multisession)
random_results <- future_replicate(100, try_random(angios))
save(random_results, file="ConiferRandomValues.rda")
AICC <- unlist(lapply(random_results[2,], "[[", "aicc"))
deltaAICC <- AICC-BcstDAngioVar_EXPO$aicc
good <- which(deltaAICC<2)
random_results[1,good]
random_df <- data.frame()
for (i in sequence(ncol(random_results))) {
	local_df <- data.frame(random_results[1,i])
	local_df$mu <- f.exp(random_results[1,i]$angios_random$Age, random_results[1,i]$angios_random$Angiosperms, random_results[2,i]$result$mu_par)
	local_df$lambda <- rep(f.constant(random_results[1,i]$angios_random$Age, random_results[1,i]$angios_random$Angiosperms, random_results[2,i]$result$mu_par), length(random_results[1,i]$angios_random$Age))
	local_df$deltaAICC = deltaAICC[i]
	local_df$rep = i
	if(i==1) {
		random_df <- local_df
	}  else {
		random_df <- rbind(random_df, local_df)
	}
}
random_df_good <- subset(random_df, deltaAICC<=2)
random_angios <- ggplot(random_df_good, mapping=aes(x=angios_random.Age, y=angios_random.Angiosperms, group=rep)) + geom_line() + scale_x_reverse() + coord_geo()
# random_angios
# mu_pars_good  <- (lapply(random_results[2,good], "[[", "mu_par"))
# mu_pars_good_v <- unlist(mu_pars_good)
# mu_pars_df <- data.frame(mu=mu_pars_good_v[gtools::odd(seq_along(mu_pars_good))], a=mu_pars_good_v[gtools::even(seq_along(mu_pars_good))])
# plot(mu_pars_df$mu, mu_pars_df$a)
save(list=ls(), file="simple.rda")