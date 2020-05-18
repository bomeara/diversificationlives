#!/usr/local/bin/Rscript
#
# This R script is provided as a Supplemental code to the paper:
#   Louca and Pennell (2020). Extant timetrees are consistent with a myriad of diversification histories. Nature. 580:502-505.
#
# Included are:
#	Analyses of the Cetacea timetree (Steeman et al. 2009)
#	Analyses of the seed plant timetree (Smith et al. 2018)
#	Analyses of fossil-based origination/extinction rate estimates for marine invertebrates (Alroy 2008)
#	Simulation of various hypothetical diversification scenarios and analysis of the generated timetrees
#
# To run this script, type the following in your terminal:
#	Rscript workflow.R
# Note that the script will automatically attempt to download any required packages!
# The script will save all output in a newly generated output folder, each time it is ran (i.e., previous output is not replaced).
#
# Tested on R 3.6.0, MacOS 10.13.6, castor v1.5.7.
# On a 2015 MacBook Pro, the full script takes about 20 hours to finish.
#
# LICENSE AGREEMENT
# - - - - - - - - -
# Use and redistributions of this code is permitted, provided that the following
# conditions are always met:
#
#    * Redistributions must retain the above copyright notice, this list of
#      conditions and the following disclaimer in the code itself, as well
#      as in documentation and/or other materials provided with the code.
#    * Neither the name of the original author (Stilianos Louca), nor the names
#      of its contributors may be used to endorse or promote products derived
#      from this code without specific prior written permission.
#
# THIS CODE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
# IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS CODE,
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# - - - - - - - - -
#
# Stilianos Louca
# May 16, 2020

###################################
# OPTIONS

REQUIRED_PACKAGES	= c("castor", "msir")
OUTPUT_DIR			= "output"

# NUMBER_OF_PARALLEL_THREADS	= 4 #BCO
NUMBER_OF_PARALLEL_THREADS	= parallel::detectCores() #BCO
INTEGRATION_RELATIVE_DT		= 1e-3

# options for fitting models to trees
# DEFAULT_FITTING_NTRIALS			= 4 * NUMBER_OF_PARALLEL_THREADS # this can be increased, but fitting will take longer #BCO
DEFAULT_FITTING_NTRIALS			= 4 * 4 # this can be increased, but fitting will take longer #BCO
FITTING_NITERATIONS 			= 500
FITTING_NEVALUATIONS 			= 500
FITTING_REL_TOLERANCE 			= 1e-8
FITTING_TIPS_PER_RUNTIME_SECOND	= 5e4 		# number of tips for which to allocate one second max runtime per likelihood evaluation. Reducing this number will increase computation time.
FITTING_GRID_SPLINES_DEGREE		= 1 		# splines degree to assume when fitting HBD model params (lambda, mu or PDR) on a temporal grid


# plot styles for specific data types
PLOT_COLOR_LTT="#000000"
PLOT_COLOR_FIT_DLTT="#0071b8"
PLOT_LW_LTT=1
PLOT_LW_FIT_DLTT=2
PLOT_COLOR_TOTAL_DIVERSITY="#9b00a1"
PLOT_COLOR_PDR="#1c568e"
PLOT_COLOR_DIVERSIFICATION_RATE="#202020"
PLOT_COLOR_LAMBDA="#808080"
PLOT_COLOR_MU="#202020"
PLOT_LW_LAMBDA=2
PLOT_LW_MU=1
PLOT_COLOR_PSR="#1c568e"

PLOT_LINE_TYPE_DETERMINISTIC_TRUE=1
PLOT_LINE_TYPE_EMPIRICAL=2
PLOT_LINE_TYPE_DETERMINISTIC_FITTED=3

# default plot style palettes
PLOT_COLOR_PALETTE=c("black", "red", "blue", "brown", "darkgreen", "#D3762A", "#618ACC", "#B80DC4", "#565656")
PLOT_LINE_TYPE_PALETTE=rep(c(1,1,2,2,3,3),times=5)
PLOT_LINE_WIDTH_PALETTE=rep(c(2,1),times=5)


DEFAULT_PLOT_WIDTH	= 5 	# inches
DEFAULT_PLOT_HEIGHT	= 5  	# inches
PLOT_DOWNSAMPLING_RESOLUTION = 1000

# only perform a subset of analyses
INCLUDE_ALROY2008_FOSSIL_DATA=FALSE
INCLUDE_SEED_PLANT_TREE=FALSE
INCLUDE_CETACEA_TREE=TRUE
INCLUDE_SIM_MODEL_ANALYSIS=FALSE


# define cladogenic models for simulating trees, as well as various diversification analyses to perform on them
SIM_MODELS = list(list(	name 			= "radiation_followed_by_extinction",
						base_lambda		= 1,
						base_mu 		= 0.5,
						added_lambda 	= function(time){ return(6*exp(-((time-10)^2)/(2*0.5^2))); },
						added_mu 		= function(time){ return(1.5*exp(-((time-13)^2)/(2*0.5^2))); },
						congruent_added_mus	= list(function(time){ return(10*exp(-((time-9)^2)/(2*0.5^2))); }),
						rho0			= 0.5,				# present-day sampling fraction
						max_time		= 15,				# max simulation time for generating the tree
						max_tips		= 10000000,			# maximum number of allowed tips, before halting the simulation of the tree
						Nfitting_ages	= c(5,10,15), 		# number of discrete ages at which to fit PDR/lambda/mu during model or class fitting. Can be a vector, in which case multiple fits are performed, each at a different age_grid size
						random_seed		= 165992,			# optional fixed randomization seed to use for the generation & analysis of this tree. Can be NULL
						fitting_min_lineages 		= 500, 	# only fit a model during ages where the number of lineages in the tree is at or above this treshold
						include_model_comparison	= TRUE,
						include_PDR_grid_fitting	= TRUE,
						include_PSR_grid_fitting	= TRUE,
						include_grid_model_fitting	= TRUE,
						include_grid_lambda_fitting = TRUE),
				list(	name 			= "extinction_event",
						base_lambda		= 1,
						base_mu			= 0.1,
						added_lambda	= function(time){ return(0*time); },
						added_mu		= function(time){ return(1.5*exp(-((time-12)^2)/(2*0.5^2))); },
						congruent_added_mus	= list( function(time){ return(0*time); },
													function(time){ return(1.5*exp(-((time-7)^2)/(2*0.5^2))); },
													function(time){ return(3*exp(-time/6)); } ),
						rho0			= 0.5,
						max_time		= 17,
						max_tips		= 1000000,
						Nfitting_ages	= c(5,10,15),
						random_seed		= 733501,
						fitting_min_lineages 		= 500,
						include_model_comparison	= TRUE,
						include_PDR_grid_fitting	= TRUE,
						include_PSR_grid_fitting	= TRUE,
						include_grid_model_fitting	= TRUE,
						include_grid_lambda_fitting = TRUE));


###################################
# AUXILIARY FUNCTIONS


check_output_file = function(file_path,force_replace,verbose,verbose_prefix){
	if(file.exists(file_path)){
		if(force_replace){
			cat(sprintf("%sNote: Replacing output file '%s'.\n",verbose_prefix,file_path))
			file.remove(file_path);
		}else{
			stop(sprintf("Output file '%s' already exists. Cowardly refusing to continue.",file_path), call.=FALSE)
		}
	}
	dir.create(dirname(file_path), showWarnings = FALSE, recursive=TRUE);
}


# reverse the sign of the linear trend of a time series, by subtracting twice the linear regression line
reverse_linear_trend = function(X,Y){
	fit = lm(Y~X, data=data.frame(X,Y))
	Yp = fit$coefficients[[1]] + fit$coefficients[[2]] * X
	return(Y-2*Yp);
}


loess_nansafe = function(X, Y, span, degree, minX=-Inf, maxX=Inf){
	valids	= which(!(is.na(Y) | is.infinite(Y) | is.nan(Y) | is.na(X) | is.infinite(X) | is.nan(X) | (X<minX) | (X>maxX)))
	S		= rep(NA, length(X))
	nsigma	= stats::qnorm(0.975)

	smoothing	= msir::loess.sd(X[valids],Y[valids], nsigma=nsigma, span=span, degree=degree)
	S[valids]	= smoothing$y
	return(S)
}


get_non_existent_dir = function(parent_path, child_basename, digits=3){
	counter = 1
	child_path = file.path(parent_path,sprintf(sprintf("%%s%%0%dd",digits),child_basename,counter))
	while(file.exists(child_path)){
		counter = counter + 1;
		child_path = file.path(parent_path,sprintf(sprintf("%%s%%0%dd",digits),child_basename,counter))
	}
	return(child_path);
}


get_PDR_on_grid = function(age_grid, lambdas, mus){
	NG				= length(age_grid)
	lambda_slopes 	= diff(lambdas)/diff(age_grid);
	lambda_slopes 	= c(lambda_slopes[1],lambda_slopes,tail(lambda_slopes,1)); # add dummy slope to the left & right, so that lambda_slopes becomes of length NG+1
	lambda_slopes 	= 0.5*(lambda_slopes[1:NG]+lambda_slopes[2:(NG+1)]);
	PDRs			= lambdas - mus + lambda_slopes/lambdas;
	return(PDRs);
}


# save & plot various quantities vs age
plot_age_curves = function(	file_basepath, 		# e.g. 'output/SILVA_curves_last_1000years'
							data_type,			# e.g. 'curves' or 'birth_rates'
							case_tag,			# e.g. 'model 2'
							curves, 			# list of size Ncurves, each entry of which is a sub-list with two elements (ages, values) specifying a separate curve to be plotted
							curve_names,		# 1D character vector of size >=Ncurves
							max_age				= Inf,			# maximum age to plot. If Inf, the full available age range is shown
							miny				= NA,
							maxy				= NA,
							Nages				= NULL,		# either NULL or an integer, specifying the temporal resolution for downsampling curves. May be needed in order to avoid excessively large data & plot files.
							plot_curves			= NULL,		# either a 1D vector of booleans of size Ncurves, indicating whether a curve should be plotted, or NULL (plot all curves)
							curve_colors		= "#000000",		# either NULL or a 1D vector of size >=Ncurves, specifying the color of a curve. If NULL, colors are picked automatically.
							curve_line_types	= 1,	# either NULL or a 1D vector of size >=Ncurves, specifying the line type of a curve. If NULL, line types are picked automatically.
							curve_widths		= 1,		# either NULL or a 1D vector of size >=Ncurves, specifying the line width of a curve. If NULL, line widths are picked automatically.
							age_label			= "age",	# e.g. 'age (years)'
							value_label			= "value",	# e.g. 'number of lineages'
							plot_log_values		= FALSE,
							legend_pos			= "none",	# (string) if "none", no legend will be shown
							plot_title			= "",
							data_file_comments 	= "",
							verbose 			= FALSE,
							verbose_prefix 		= "  "){
	curves  = curves[sapply(curves,FUN = function(l) !is.null(l))] # remove NULL elements
	Ncurves = length(curves);
	if(is.null(curve_colors)){ curve_colors = PLOT_COLOR_PALETTE[1:Ncurves]; }
	else if(length(curve_colors)==1){ curve_colors = rep(curve_colors,times=Ncurves); }
	if(is.null(curve_line_types)){ curve_line_types = PLOT_LINE_TYPE_PALETTE[1:Ncurves] }
	else if(length(curve_line_types)==1){ curve_line_types = rep(curve_line_types,times=Ncurves); }
	if(is.null(curve_widths)){ curve_widths = PLOT_LINE_WIDTH_PALETTE[1:Ncurves] }
	else if(length(curve_widths)==1){ curve_widths = rep(curve_widths,times=Ncurves); }
	if(is.null(plot_curves)) plot_curves = rep(TRUE,Ncurves)

	if(!is.null(Nages)){
		if(verbose) cat(sprintf("%sDownsampling curves to %d time points..\n",Nages))
		for(n in 1:Ncurves){
			if(sum(!is.na(curves[[n]][[2]]))<2) next;
			X = curves[[n]][[1]];
			Xnew = seq(from=X[1], to=tail(X,1), length.out=Nages);
			curves[[n]][[1]] = Xnew;
			curves[[n]][[2]] = approx(x=X, y=curves[[n]][[2]], xout=Xnew, method="linear", yleft=NaN, yright=NaN, rule = 1, f = 0, ties = mean)$y;
		}
	}

	# save data as text file
	if(verbose) cat(sprintf("%sSaving %s to table (%s)..\n",verbose_prefix,data_type,case_tag))
	output_table=sprintf("%s.tsv",file_basepath)
	check_output_file(output_table,TRUE,TRUE,"  ")
	cat(sprintf("# %s over age (%s)\n%s\n#\n",data_type,case_tag,data_file_comments), file=output_table, append=FALSE)
	for(n in 1:Ncurves){
		cat(sprintf("# age\t%s\n",curve_names[[n]]), file=output_table, append=TRUE);
		write.table(x=cbind(curves[[n]][[1]],curves[[n]][[2]]), file=output_table, append=TRUE, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE);
		cat(sprintf("\n\n"), file=output_table, append=TRUE);
	}

	# determine X & Y ranges
	minx = NA
	maxx = NA
	for(n in 1:Ncurves){
		X = curves[[n]][[1]];
		Y = curves[[n]][[2]];
		if(plot_log_values){
			Y = Y[is.finite(log(Y)) & (X<=max_age)];
		}else{
			Y = Y[is.finite(Y) & (X<max_age)];
		}
		X = X[is.finite(X) & (X<max_age)];
		if(plot_log_values) Y = Y[Y>0]
		if(length(Y)>0){
			miny = (if(is.na(miny)) min(Y) else min(miny,min(Y)));
			maxy = (if(is.na(maxy)) max(Y) else max(maxy,max(Y)));
			minx = (if(is.na(minx)) min(X) else min(minx,min(X)));
			maxx = (if(is.na(maxx)) max(X) else max(maxx,max(X)));
		}
	}
	xrange = maxx-minx;
	yrange = maxy-miny;
	if(is.na(miny) || is.na(maxy) || is.na(minx) || is.na(maxx)){
		if(verbose) cat(sprintf("%sWARNING: No valid points for plotting\n",verbose_prefix))
		return();
	}
	if(is.infinite(max_age)) max_age = maxx

	# plot to PDF
	if(verbose) cat(sprintf("%sPlotting %s over age (%s)..\n",verbose_prefix,data_type,case_tag))
	curves 				= curves[plot_curves]
	curve_names			= curve_names[plot_curves]
	curve_colors		= curve_colors[plot_curves]
	curve_line_types	= curve_line_types[plot_curves]
	curve_widths		= curve_widths[plot_curves]
	plot_file 			= sprintf("%s.pdf",file_basepath)
	check_output_file(plot_file,TRUE,TRUE,"  ")
	pdf(file=plot_file, width=DEFAULT_PLOT_WIDTH+(if(legend_pos=="outside") 2.5 else 0), height=DEFAULT_PLOT_HEIGHT);
	if(legend_pos=="outside") par(mar=c(5,5,5,14))
	if(plot_log_values){
		valids = which(is.finite(log(curves[[1]][[2]])));
	}else{
		valids = which(is.finite(curves[[1]][[2]]));
	}
	plot(	x		= curves[[1]][[1]][valids],
			y		= curves[[1]][[2]][valids],
			col		= curve_colors[1],
			lty		= curve_line_types[1],
			lwd		= curve_widths[1],
			las		= 1,
			type	= "l",
			main	= plot_title,
			xlab	= age_label,
			ylab	= NA,
			cex		= 1.1,
			log		= (if(plot_log_values) "y" else ""),
			yaxt	= (if(plot_log_values) "n" else NULL),
			xaxs	= 'i',
			yaxs	= 'i',
			xlim	= c(max_age,0),
			ylim	= c(miny, maxy+yrange*0.1));
	title(ylab=value_label, line=4)
	if(Ncurves>1){
		for(n in 2:Ncurves){
			if(plot_log_values){
				valids = which(is.finite(log(curves[[n]][[2]])));
			}else{
				valids = which(is.finite(curves[[n]][[2]]));
			}
			lines(	x	= curves[[n]][[1]][valids],
					y	= curves[[n]][[2]][valids],
					type= "l",
					col	= curve_colors[n],
					lty	= curve_line_types[n],
					lwd	= curve_widths[n]);

		}
	}
	if(legend_pos!="none"){
		if(legend_pos=="outside"){
			legendx = minx - 0.2*xrange;
			legendy = maxy;
			legend(x=legendx, y=legendy, legend = curve_names, col=curve_colors, lty=curve_line_types, lwd=curve_widths, xpd=NA);
		}else{
			legend(legend_pos, legend = curve_names, col=curve_colors, lty=curve_line_types, lwd=curve_widths);
		}
	}
	if(plot_log_values){
		# improve appearance of y-axis ticks if logarithmic
		aty = axTicks(2)
		axis(2,at=aty,labels=sapply(aty, function(x) sprintf("%g",x)), las=1)
	}
	invisible(dev.off());
}


save_object_to_file = function(object, filepath){
	dir.create(dirname(filepath), showWarnings = FALSE, recursive=TRUE);
	sink(file=filepath);
	print(fit);
	sink();
}


###############################
# PREPARATIONS

# print warnings as they occur
options(warn=1)

cat(sprintf("Loading %d required packages (%s)..\n",length(REQUIRED_PACKAGES),paste(REQUIRED_PACKAGES,collapse=", ")));
for(p in 1:length(REQUIRED_PACKAGES)){
	if(!suppressMessages(suppressPackageStartupMessages(require(REQUIRED_PACKAGES[p], quietly=TRUE, character.only=TRUE)))){
		cat(sprintf("  Note: Installing missing package '%s'..\n",REQUIRED_PACKAGES[p]))
		install.packages(REQUIRED_PACKAGES[p], dependencies=TRUE, repos="http://cran.r-project.org/");
		suppressMessages(suppressPackageStartupMessages(require(REQUIRED_PACKAGES[p],quietly=TRUE,character.only=TRUE,warn.conflicts = FALSE)))
	}
}

# seed random number generator
set.seed(NULL)
global_random_seed = sample.int(n=1000000,size=1)
set.seed(global_random_seed)
cat(sprintf("Note: Seeding global random generator at %d\n",global_random_seed));


# prepare output dir
output_dir = get_non_existent_dir("output", "run_", 3);
dir.create(output_dir, showWarnings = FALSE, recursive=TRUE);
cat(sprintf("All output will be written to '%s'..\n",output_dir))


# prepare log file
logfile = sprintf("%s/log.txt",output_dir);
cat2 = function(message, file=logfile, append=TRUE){
	cat(message);
	cat(message, file=file, append=append);
}




###############################################
# analysis of marine invertebrate fossil-based speciation/extinction rates (Alroy 2008)
# Extended Data Fig. 2 in Louca and Pennell (2020)

if(INCLUDE_ALROY2008_FOSSIL_DATA){
	dataset_name	= "Alroy2008_marine_invertebrates"
	lambda_tsv		= "input/Alroy2008_marine_invertebrate_fossil_rates/Alroy2008_speciation_rates.tsv"
	mu_tsv			= "input/Alroy2008_marine_invertebrate_fossil_rates/Alroy2008_extinction_rates.tsv"
	time_unit		= "Myr"
	LTT0			= 1000000 # number of tips for simulations. This does not matter much, since everything is just rescaled for Ntips
	rho0			= 1 # present-day sampling fraction to use for all calculations. The choice does not matter for our purposes, as long as we use the same value throughout.
	oldest_age		= 400 # oldest age to consider for plots and simulations

	cat2(sprintf("Examining fossil time series '%s'..\n",dataset_name));
	timeseries_output_dir=sprintf("%s/%s",output_dir,dataset_name);
	dir.create(timeseries_output_dir, showWarnings = FALSE, recursive=TRUE);

	cat2(sprintf("  Loading lambdas & mus from files..\n"))
	loaded_lambdas 	= read.table(file=lambda_tsv, sep="\t", header=FALSE, stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA","NaN"))
	loaded_mus 		= read.table(file=mu_tsv, sep="\t", header=FALSE, stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA","NaN"))
	# sort in increasing age
	loaded_lambdas 	= loaded_lambdas[order(loaded_lambdas[,1]),]
	loaded_mus		= loaded_mus[order(loaded_mus[,1]),]

	# consolidate loaded lambdas & mus onto the same age grid
	max_covered_age = min(max(loaded_lambdas[,1]), max(loaded_mus[,1]))
	loaded_ages 	= sort(c(loaded_lambdas[,1],loaded_mus[,1]))
	loaded_ages		= loaded_ages[loaded_ages<=max_covered_age] # only include ages covered by both the lambda & mu time series
	oldest_age		= min(oldest_age,tail(loaded_ages,1)) # make sure oldest_age falls within loaded_ages[]
	loaded_lambdas 	= approx(x=loaded_lambdas[,1], y=loaded_lambdas[,2], xout=loaded_ages)$y
	loaded_mus 		= approx(x=loaded_mus[,1], y=loaded_mus[,2], xout=loaded_ages)$y

	# smoothen loaded rates
	loaded_lambdas = loess_nansafe(loaded_ages, loaded_lambdas, 0.1, 2)
	loaded_mus = loess_nansafe(loaded_ages, loaded_mus, 0.1, 2)

	# plot loaded rates
	# Extended Data Fig 2a in Louca and Pennell (2020)
	plot_age_curves(file_basepath 		= sprintf("%s/loaded_lambda_mu",timeseries_output_dir),
					data_type 			= "lambdas & mus",
					case_tag 			= sprintf("%s",dataset_name),
					curves 				= list(	list(loaded_ages, loaded_lambdas),
												list(loaded_ages, loaded_mus)),
					curve_names 		= c("lambda","mu"),
					curve_colors 		= c(PLOT_COLOR_LAMBDA,PLOT_COLOR_MU),
					curve_line_types	= c(1,1),
					curve_widths 		= c(PLOT_LW_LAMBDA,PLOT_LW_MU),
					max_age 			= oldest_age,
					age_label 			= sprintf("age (%s)",time_unit),
					value_label 		= sprintf("rate (1/%s)",time_unit),
					plot_log_values 	= FALSE,
					legend_pos 			= "outside",
					plot_title 			= sprintf("%s loaded rates",dataset_name),
					data_file_comments 	= sprintf("Loaded speciation and extinction rates for timeseries '%s'",dataset_name),
					verbose 			= FALSE,
					verbose_prefix 		= "    ");


	# simulate loaded model to get PDR
	cat2(sprintf("  Simulating deterministic BD model with loaded rates..\n"))
	lsimulation = simulate_deterministic_hbd(	LTT0			= LTT0,
												oldest_age		= oldest_age,
												rho0			= rho0,
												age_grid		= loaded_ages,
												lambda			= loaded_lambdas,
												mu				= loaded_mus,
												splines_degree	= 1,
												relative_dt		= INTEGRATION_RELATIVE_DT);
	if(!lsimulation$success){
		cat2(sprintf("    ERROR: Simulation failed: %s\n",lsimulation$error))
		stop()
	}else{
		cat2(sprintf("  Plotting simulated PDR..\n"));
		plot_age_curves(file_basepath 		= sprintf("%s/PDR",timeseries_output_dir),
						data_type 			= "PDR",
						case_tag 			= sprintf("%s",dataset_name),
						curves 				= list(	list(lsimulation$ages, lsimulation$PDR)),
						curve_names 		= c("PDR"),
						curve_colors 		= c(PLOT_COLOR_PDR),
						curve_line_types	= c(1),
						curve_widths 		= c(2),
						plot_curves 		= NULL,
						max_age 			= oldest_age,
						Nages				= NULL,
						age_label 			= sprintf("age (%s)",time_unit),
						value_label 		= sprintf("PDR (1/%s)",time_unit),
						plot_log_values 	= FALSE,
						legend_pos 			= "outside",
						plot_title 			= sprintf("%s PDR, based on loaded rates",dataset_name),
						data_file_comments 	= sprintf("PDR for timeseries '%s'",dataset_name),
						verbose 			= FALSE,
						verbose_prefix 		= "    ");
	}


	# simulate alternative congruent scenario, where the trend of mu has been reversed
	# Extended Data Fig 2b in Louca and Pennell (2020)
	cat2(sprintf("  Simulating alternative congruent model (with reverse-trended mu)..\n"));
	reversed_mu = reverse_linear_trend(lsimulation$ages,lsimulation$mu)
	reversed_mu = reversed_mu - min(reversed_mu) # make sure new mu is positive
	rsimulation = simulate_deterministic_hbd(	LTT0			= LTT0,
												oldest_age		= oldest_age,
												rho0			= rho0,
												age_grid		= lsimulation$ages,
												PDR				= lsimulation$PDR,
												lambda0			= lsimulation$lambda[1],
												mu				= reversed_mu,
												splines_degree	= 1,
												relative_dt		= INTEGRATION_RELATIVE_DT);
	if(!rsimulation$success){
		cat2(sprintf("    WARNING: Simulation failed: %s\n",lsimulation$error))
		stop()
	}else{
		cat2(sprintf("  Plotting congruent model for '%s'..\n",dataset_name));
		plot_age_curves(file_basepath 		= sprintf("%s/congruent_lambda_mu_reverse_trended",timeseries_output_dir),
						data_type 			= "lambdas & mus",
						case_tag 			= sprintf("%s, reverse-trended mu",dataset_name),
						curves 				= list(	list(rsimulation$ages, rsimulation$lambda),
													list(rsimulation$ages, rsimulation$mu)),
						curve_names 		= c("lambda","mu"),
						curve_colors 		= c(PLOT_COLOR_LAMBDA,PLOT_COLOR_MU),
						curve_line_types	= c(1,1),
						curve_widths 		= c(PLOT_LW_LAMBDA,PLOT_LW_MU),
						plot_curves 		= NULL,
						max_age 			= oldest_age,
						Nages				= NULL,
						age_label 			= sprintf("age (%s)",time_unit),
						value_label 		= sprintf("rate (1/%s)",time_unit),
						plot_log_values 	= FALSE,
						legend_pos 			= "outside",
						plot_title 			= sprintf("Congruent model to %s\nreverse-trended mu",dataset_name),
						data_file_comments 	= sprintf("Speciation and extinction rates for congruent model to timeseries '%s', where mu has been reverse-trended",dataset_name),
						verbose 			= FALSE,
						verbose_prefix 		= "    ");
	}


	# simulate alternative congruent scenario, with zero mu
	# Extended Data Fig 2c in Louca and Pennell (2020)
	cat2(sprintf("  Simulating alternative congruent model (with zero mu)..\n"));
	zmsimulation = simulate_deterministic_hbd(	LTT0			= LTT0,
												oldest_age		= oldest_age,
												rho0			= rho0,
												age_grid		= lsimulation$ages,
												PDR				= lsimulation$PDR,
												lambda0			= lsimulation$lambda[1],
												mu				= 0,
												splines_degree	= 1,
												relative_dt		= INTEGRATION_RELATIVE_DT);
	if(!zmsimulation$success){
		cat2(sprintf("    WARNING: Simulation failed: %s\n",zmsimulation$error))
	}else{
		cat2(sprintf("  Plotting congruent model for '%s'..\n",dataset_name));
		plot_age_curves(file_basepath 		= sprintf("%s/congruent_lambda_mu_zero_extinction",timeseries_output_dir),
						data_type 			= "lambdas & mus",
						case_tag 			= sprintf("%s, zero mu",dataset_name),
						curves 				= list(	list(zmsimulation$ages, zmsimulation$lambda),
													list(zmsimulation$ages, zmsimulation$mu)),
						curve_names 		= c("lambda","mu"),
						curve_colors 		= c(PLOT_COLOR_LAMBDA,PLOT_COLOR_MU),
						curve_line_types	= c(1,1),
						curve_widths 		= c(PLOT_LW_LAMBDA,PLOT_LW_MU),
						plot_curves 		= NULL,
						max_age 			= oldest_age,
						Nages				= NULL,
						age_label 			= sprintf("age (%s)",time_unit),
						value_label 		= sprintf("rate (1/%s)",time_unit),
						plot_log_values 	= FALSE,
						legend_pos 			= "outside",
						plot_title 			= sprintf("Congruent model to %s\nAssuming zero mu",dataset_name),
						data_file_comments 	= sprintf("Speciation and extinction rates for congruent model to timeseries '%s', where mu is assumed to be zero",dataset_name),
						verbose 			= FALSE,
						verbose_prefix 		= "    ");
	}
}









####################################
# Analysis of seed plants tree
# obtained from Smith et al 2018

if(INCLUDE_SEED_PLANT_TREE){
	tree_name	= "Seed_Plants_Smith2018"
	tree_path	= "input/SeedPlants_Smith2018/GBMB.tre"
	time_unit	= "Myr"
	extant_N	= 422127 # number of extant seed plant species. [Govaerts (2001). How many species of seed plants are there?] estimates that there exist 422127 seed plants on Earth.
	oldest_plot_age	= 100 # oldest considered age for plotting
	oldest_fit_age	= 130 # oldest considered age for fitting

	cat2(sprintf("Examining external tree '%s'..\n",tree_name));
	tree_output_dir=sprintf("%s/%s",output_dir,tree_name);
	dir.create(tree_output_dir, showWarnings = FALSE, recursive=TRUE);

	cat2(sprintf("  Loading tree from file..\n"))
	tree		= castor::read_tree(file=tree_path)
	tree 		= castor::extend_tree_to_height(tree)$tree # make sure tree is ultrametric (correct numerical rounding errors)
	Ntips		= length(tree$tip.label)
	root_age	= castor::get_tree_span(tree)$max_distance
	rho0		 = Ntips/extant_N # present-day sampling fraction
	cat2(sprintf("    --> Loaded tree has %d tips and %d nodes, root_age=%g %s, sampling fraction=%g\n",Ntips,tree$Nnode,root_age,time_unit,rho0));

	cat2(sprintf("  Calculating the tree's LTT..\n"))
	empirical_LTT = castor::count_lineages_through_time(tree,Ntimes=10*log2(Ntips),include_slopes=TRUE);
	empirical_LTT$ages = root_age - empirical_LTT$times

	#########################
	# fit a HBD model with exponentially varying lambda and mu
	# Figure 2 in Louca and Pennell (2020)

	cat(sprintf("  Fitting exp-HBD model to tree '%s'..\n",tree_name))
	fit_dir = sprintf("%s/fitting_exp_HBD_model",tree_output_dir)
	dir.create(fit_dir, showWarnings = FALSE, recursive=TRUE);
	age_grid 	 = seq(from=0,to=oldest_fit_age,length.out=1000)
	lambda 		 = function(age,params){ pmax(0,params['alpha'] * exp(params['beta']*(age/root_age))); }
	mu 			 = function(age,params){ pmax(0,params['gamma'] * exp(params['delta']*(age/root_age))); }
	start_lambda = max(empirical_LTT$relative_slopes/rho0)
	fit = fit_hbd_model_parametric(	tree,
									age_grid			= age_grid,
									param_values		= c(alpha=NA, beta=NA, gamma=NA, delta=NA),
									param_guess			= c(start_lambda,0,0.5*start_lambda,0),
									param_min			= c(0,-5,0,-5),
									param_max			= c(100*start_lambda,5,100*start_lambda,5),
									param_scale			= 1,
									oldest_age			= oldest_fit_age,
									lambda				= lambda,
									mu					= mu,
									rho0				= function(params) rho0,
									relative_dt			= INTEGRATION_RELATIVE_DT,
									Ntrials				= DEFAULT_FITTING_NTRIALS,
									max_start_attempts	= 10,
									Nthreads			= NUMBER_OF_PARALLEL_THREADS,
									max_model_runtime	= max(1,Ntips/FITTING_TIPS_PER_RUNTIME_SECOND),
									fit_control			= list(eval.max=FITTING_NEVALUATIONS, iter.max=FITTING_NITERATIONS, rel.tol=FITTING_REL_TOLERANCE))
	if(!fit$success){
		cat2(sprintf("    ERROR: Fitting failed: %s\n",fit$error));
		stop()
	}
	save_object_to_file(fit, sprintf("%s/exp_BD_fit_results.txt",fit_dir));
	cat2(sprintf("  Simulating fitted model..\n"));
	fsimulation = simulate_deterministic_hbd(	LTT0			= Ntips,
												oldest_age		= oldest_plot_age,
												rho0			= rho0,
												age_grid		= age_grid,
												lambda			= lambda(age_grid,fit$param_fitted),
												mu				= mu(age_grid,fit$param_fitted),
												splines_degree	= 1,
												relative_dt		= INTEGRATION_RELATIVE_DT);
	if(!fsimulation$success){
		cat2(sprintf("    WARNING: Simulation failed: %s\n",fsimulation$error))
	}else{
		cat2(sprintf("  Plotting LTT of tree '%s' and fitted dLTT..\n",tree_name));
		plot_age_curves(file_basepath 		= sprintf("%s/fitted_HBD_model_LTT",fit_dir),
						data_type 			= "LTTs",
						case_tag 			= sprintf("%s, fitted dLTT",tree_name),
						curves 				= list(	list(empirical_LTT$ages, empirical_LTT$lineages),
													list(fsimulation$ages, fsimulation$LTT)),
						curve_names 		= c("tree LTT","fitted dLTT"),
						curve_colors 		= c(PLOT_COLOR_LTT,PLOT_COLOR_FIT_DLTT),
						curve_line_types	= c(1,2),
						curve_widths 		= c(PLOT_LW_LTT,PLOT_LW_FIT_DLTT),
						max_age 			= oldest_plot_age,
						age_label 			= sprintf("age (%s)",time_unit),
						value_label 		= sprintf("lineages"),
						plot_log_values 	= TRUE,
						legend_pos 			= "outside",
						plot_title 			= sprintf("LTT of tree and fitted dLTT\n%s",tree_name),
						data_file_comments 	= sprintf("LTT of external tree '%s', and fitted dLTT",tree_name),
						verbose 			= FALSE,
						verbose_prefix 		= "    ");

	}

	# Construct a congruent model with modified exponent for mu
	cat2(sprintf("  Simulating congruent model, with alternative delta..\n"));
	params = fit$param_fitted
	params['delta'] = - fit$param_fitted['delta'] # simply negate the fitted delta
	csimulation = simulate_deterministic_hbd(	LTT0			= Ntips,
												oldest_age		= oldest_plot_age,
												rho0			= rho0,
												age_grid		= fsimulation$ages,
												PDR				= fsimulation$PDR,
												lambda0 		= lambda(0,fit$param_fitted),
												mu				= mu(fsimulation$ages,params),
												splines_degree	= 1,
												relative_dt		= INTEGRATION_RELATIVE_DT);
	if(!csimulation$success){
		cat2(sprintf("      WARNING: Simulation of congruent model failed: %s\n",csimulation$error))
	}else{
		cat2(sprintf("  Plotting rates of fitted and congruent models..\n"));
		plot_age_curves(file_basepath 		= sprintf("%s/fitted_and_congruent_lambda_mu",fit_dir),
						data_type 			= "lambdas & mus",
						case_tag 			= sprintf("%s, fit & congruent",tree_name),
						curves 				= list(	list(fsimulation$ages, fsimulation$lambda),
													list(fsimulation$ages, fsimulation$mu),
													list(csimulation$ages, csimulation$lambda),
													list(csimulation$ages, csimulation$mu)),
						curve_names 		= c("lambda (fit)","mu (fit)","lambda (congruent)","mu (congruent)"),
						curve_colors 		= c(PLOT_COLOR_LAMBDA,PLOT_COLOR_MU,PLOT_COLOR_LAMBDA,PLOT_COLOR_MU),
						curve_line_types	= c(1,1,2,2),
						curve_widths 		= c(1,1,2,2),
						max_age 			= oldest_plot_age,
						age_label 			= sprintf("age (%s)",time_unit),
						value_label 		= sprintf("rate (1/%s)",time_unit),
						plot_log_values 	= FALSE,
						legend_pos 			= "outside",
						plot_title 			= sprintf("Fitted exp-HBD model & congruent model\n%s",tree_name),
						data_file_comments 	= sprintf("Speciation and extinction rates for fitted and congruent model, tree '%s'",tree_name),
						verbose 			= FALSE,
						verbose_prefix 		= "    ");
	}


	#########################
	# Fit a parametric and a grid-model to the tree
	# Extended data figure 5g-i in Louca and Pennell (2020)

	# Fit a parametric HBD model to the tree
	# Functional forms:
	#	lambda(t) = p1 * exp(-p2*t) + p3 + p4*t + p5*t^2 + p6*t^3 + p7*t^4
	#	mu(t) 	  = q1 * exp(-q2*t) + q3 + q4*t + q5*t^2 + q6*t^3 + q7*t^4
	# where t is relative age, i.e., in multiples of root_age.
	# The sampling fraction is assumed known.
	cat2(sprintf("  Fitting parametric-HBD-model to tree '%s' (this may take a while)..\n",tree_name))
	age_grid	= seq(from=0,to=oldest_fit_age,length.out=1000)
	lambda 		= function(age,params){ pmax(0,params['p1'] * exp(-params['p2']*(age/root_age)) + params['p3'] + params['p4']*(age/root_age) + params['p5']*(age/root_age)^2 + params['p6']*(age/root_age)^3 + params['p7']*(age/root_age)^4); }
	mu 			= function(age,params){ pmax(0,params['q1'] * exp(-params['q2']*(age/root_age)) + params['q3']+ params['q4']*(age/root_age) + params['q5']*(age/root_age)^2 + params['q6']*(age/root_age)^3 + params['q7']*(age/root_age)^4); }
	pfit = fit_hbd_model_parametric(tree,
									age_grid			= age_grid,
									param_values		= c(p1=NA, p2=NA, p3=NA, p4=NA, p5=NA, p6=NA, p7=NA, q1=NA, q2=NA, q3=NA, q4=NA, q5=NA, q6=NA, q7=NA),
									param_guess			= c(1,0,0,0,0,0,0,1,0,0,0,0,0,0),
									param_min			= -100,
									param_max			= +100,
									param_scale			= 1,
									oldest_age			= oldest_fit_age,
									lambda				= lambda,
									mu					= mu,
									rho0				= function(params) rho0,
									condition			= "stem",
									relative_dt			= INTEGRATION_RELATIVE_DT,
									Ntrials				= DEFAULT_FITTING_NTRIALS,
									max_start_attempts	= 10,
									Nthreads			= NUMBER_OF_PARALLEL_THREADS,
									max_model_runtime	= max(1,Ntips/FITTING_TIPS_PER_RUNTIME_SECOND),
									fit_control			= list(eval.max=FITTING_NEVALUATIONS, iter.max=FITTING_NITERATIONS, rel.tol=FITTING_REL_TOLERANCE))
	save_object_to_file(pfit, sprintf("%s/parametric_fit_results.txt",tree_output_dir))
	if(!pfit$success){
		cat2(sprintf("    ERROR: Fitting parametric model failed: %s\n",pfit$error))
		stop()
	}
	cat2(sprintf("  Simulating fitted parametric model..\n"));
	psimulation = simulate_deterministic_hbd(	LTT0			= Ntips,
												oldest_age		= oldest_plot_age,
												rho0			= rho0,
												age_grid		= age_grid,
												lambda			= lambda(age_grid,pfit$param_fitted),
												mu				= mu(age_grid,pfit$param_fitted),
												splines_degree	= 1,
												relative_dt		= INTEGRATION_RELATIVE_DT);
	if(!psimulation$success){
		cat2(sprintf("    ERROR: Simulation failed: %s\n",psimulation$error))
		stop()
	}

	# Fit an HBD model on a grid
	cat2(sprintf("  Fitting HBD-model on grid to tree '%s' (this may take a while)..\n",tree_name))
	Nfitting_ages = 8 # number of age-grid points on which lambda & mu will be fitted
	gfit = fit_hbd_model_on_grid(	tree,
									oldest_age			= oldest_fit_age,
									age_grid			= seq(from=0,to=oldest_fit_age,length.out=Nfitting_ages),
									fixed_rho0			= rho0,
									splines_degree		= 1,
									condition			= "stem",
									relative_dt			= INTEGRATION_RELATIVE_DT,
									Ntrials				= DEFAULT_FITTING_NTRIALS,
									Nthreads			= NUMBER_OF_PARALLEL_THREADS,
									max_model_runtime	= max(1,Ntips/FITTING_TIPS_PER_RUNTIME_SECOND),
									fit_control			= list(eval.max=FITTING_NEVALUATIONS, iter.max=FITTING_NITERATIONS, rel.tol=FITTING_REL_TOLERANCE))
	save_object_to_file(gfit,sprintf("%s/grid_fit_results.txt",tree_output_dir))
	if(!gfit$success){
		cat2(sprintf("      ERROR: Fitting grid model failed: %s\n",gfit$error))
		stop()
	}
	cat2(sprintf("    Simulating fitted grid model..\n"));
	gsimulation = simulate_deterministic_hbd(	LTT0			= Ntips,
												oldest_age		= oldest_plot_age,
												rho0			= rho0,
												age_grid		= gfit$age_grid,
												lambda			= gfit$fitted_lambda,
												mu				= gfit$fitted_mu,
												splines_degree	= 1,
												relative_dt		= INTEGRATION_RELATIVE_DT);
	if(!gsimulation$success){
		cat2(sprintf("      WARNING: Simulation of grid model failed: %s\n",gsimulation$error))
		stop()
	}

	# compare fitted parametric and grid-model
	# Extended Data Figs 5g-i in Louca and Pennell (2020)
	cat2(sprintf("  Plotting fitted parametric and grid models..\n"))
	plot_age_curves(file_basepath 		= sprintf("%s/fitted_parametric_and_grid_lambda_mu",tree_output_dir),
					data_type 			= "lambdas & mus",
					case_tag 			= sprintf("%s, fit parametric & grid",tree_name),
					curves 				= list(	list(gsimulation$ages, gsimulation$lambda),
												list(gsimulation$ages, gsimulation$mu),
												list(psimulation$ages, psimulation$lambda),
												list(psimulation$ages, psimulation$mu)),
					curve_names 		= c("lambda (fit grid)","mu (fit grid)","lambda (fit param)","mu (fit param)"),
					curve_colors 		= c(PLOT_COLOR_LAMBDA,PLOT_COLOR_MU,PLOT_COLOR_LAMBDA,PLOT_COLOR_MU),
					curve_line_types	= c(1,1,2,2),
					curve_widths 		= c(1,1,2,2),
					max_age 			= oldest_plot_age,
					age_label 			= sprintf("age (%s)",time_unit),
					value_label 		= sprintf("rate (1/%s)",time_unit),
					plot_log_values 	= FALSE,
					legend_pos 			= "outside",
					plot_title 			= sprintf("Fitted parametric & grid-model, lambda & mu\n%s",tree_name),
					data_file_comments 	= sprintf("Speciation and extinction rates for fitted parametric and grid-model, tree '%s'",tree_name),
					verbose 			= FALSE,
					verbose_prefix 		= "    ")
	plot_age_curves(file_basepath 		= sprintf("%s/fitted_parametric_and_grid_LTT",tree_output_dir),
					data_type 			= "LTTs",
					case_tag 			= sprintf("%s, fit parametric & grid",tree_name),
					curves 				= list(	list(empirical_LTT$ages, empirical_LTT$lineages),
												list(gsimulation$ages, gsimulation$LTT),
												list(psimulation$ages, psimulation$LTT)),
					curve_names 		= c("LTT (tree)","dLTT (fit grid)","dLTT (fit param)"),
					curve_colors 		= c(PLOT_COLOR_LTT,PLOT_COLOR_FIT_DLTT,PLOT_COLOR_FIT_DLTT),
					curve_line_types	= c(1,1,2),
					curve_widths 		= c(1,1,2),
					max_age 			= oldest_plot_age,
					age_label 			= sprintf("age (%s)",time_unit),
					value_label 		= sprintf("lineages",time_unit),
					plot_log_values 	= TRUE,
					legend_pos 			= "outside",
					plot_title 			= sprintf("Fitted parametric & grid-model, LTTs\n%s",tree_name),
					data_file_comments 	= sprintf("dLTTs for fitted parametric and grid-model, tree '%s'",tree_name),
					verbose 			= FALSE,
					verbose_prefix 		= "    ")

	# get PDR at fitting grid via finite differences, then interpolate linearly in between
	gfit$PDR = get_PDR_on_grid(gfit$age_grid, gfit$fitted_lambda, gfit$fitted_mu)
	gsimulation$PDR = approx(x=gfit$age_grid, y=gfit$PDR, xout=gsimulation$ages, method="linear", yleft=NaN, yright=NaN, rule = 1, f = 0, ties = mean)$y
	plot_age_curves(file_basepath 		= sprintf("%s/fitted_parametric_and_grid_PDR",tree_output_dir),
					data_type 			= "PDRs",
					case_tag 			= sprintf("%s, fit parametric & grid",tree_name),
					curves 				= list(	list(gsimulation$ages, gsimulation$PDR),
												list(psimulation$ages, psimulation$PDR)),
					curve_names 		= c("PDR (fit grid)","PDR (fit param)"),
					curve_colors 		= c(PLOT_COLOR_PDR,PLOT_COLOR_PDR),
					curve_line_types	= c(1,2),
					curve_widths 		= c(1,2),
					max_age 			= oldest_plot_age,
					age_label 			= sprintf("age (%s)",time_unit),
					value_label 		= sprintf("PDR (1/%s)",time_unit),
					plot_log_values 	= FALSE,
					legend_pos 			= "outside",
					plot_title 			= sprintf("Fitted parametric & grid-model, PDR\n%s",tree_name),
					data_file_comments 	= sprintf("Pulled diversification rates for fitted parametric and grid-model, tree '%s'",tree_name),
					verbose 			= FALSE,
					verbose_prefix 		= "    ")
}




#########################################
# Analysis of Cetacean tree (Extended Data Fig 3 in Louca and Pennell (2020))
# Tree taken from Steeman et al. 2009

if(INCLUDE_CETACEA_TREE){
	tree_name			= "Cetacea_Steeman2009"
	tree_path			= "input/Cetacea_Steeman2009/Cetacea_Steeman2009_timetree.tre"
	precomputed_rates	= "input/Cetacea_Steeman2009/Steeman2009lambda.tsv"
	time_unit			= "Myr"
	extant_N			= 89 # # number of extant Cetacea species. [Steeman et al. (2009). Radiation of Extant Cetaceans Driven by Restructuring of the Oceans. Systematic Biology. 58:573-585]

	cat2(sprintf("Examining external tree '%s'..\n",tree_name));
	tree_output_dir=sprintf("%s/%s",output_dir,tree_name);
	dir.create(tree_output_dir, showWarnings = FALSE, recursive=TRUE);

	cat2(sprintf("  Loading tree from file..\n"))
	tree		= castor::read_tree(file=tree_path)
	tree 		= castor::extend_tree_to_height(tree)$tree # make sure tree is ultrametric (correct numerical rounding errors)
	Ntips		= length(tree$tip.label)
	root_age	= castor::get_tree_span(tree)$max_distance
	rho0		 = Ntips/extant_N # present-day sampling fraction
	cat2(sprintf("    --> Loaded tree has %d tips and %d nodes, root_age=%g %s, sampling fraction=%g\n",Ntips,tree$Nnode,root_age,time_unit,rho0));

	cat2(sprintf("  Calculating the tree's LTT..\n"))
	empirical_LTT = castor::count_lineages_through_time(tree,Ntimes=10*log2(Ntips),include_slopes=TRUE);
	empirical_LTT$ages = root_age - empirical_LTT$times

	# load precomputed speciation (and extinction rates, if available) from file
	precomputed_rates = read.table(file=precomputed_rates, sep="\t", header=FALSE, stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA","NaN"))
	precomputed_rates = precomputed_rates[order(precomputed_rates[,1]),] # sort in ascending age
	precomputed_includes_mu = (ncol(precomputed_rates)>2)
	if(!precomputed_includes_mu) precomputed_rates = cbind(precomputed_rates,rep(0,times=nrow(precomputed_rates)))
	cat2(sprintf("  Note: Loaded previously fitted %s rates from file\n",(if(precomputed_includes_mu) "speciation & extinction" else "speciation")))
	precomputed_rates[,2] = loess_nansafe(precomputed_rates[,1], precomputed_rates[,2], 0.1, 2) # smoothen lambda
	precomputed_rates[,3] = loess_nansafe(precomputed_rates[,1], precomputed_rates[,3], 0.1, 2) # smoothen mu

	cat2(sprintf("  Simulating model with previously fitted rates..\n"))
	psimulation = simulate_deterministic_hbd(	LTT0			= Ntips,
												oldest_age		= max(precomputed_rates[,1]),
												rho0			= rho0,
												age_grid		= precomputed_rates[,1],
												lambda			= precomputed_rates[,2],
												mu				= precomputed_rates[,3],
												splines_degree	= 1,
												relative_dt		= INTEGRATION_RELATIVE_DT);
	if(!psimulation$success){
		cat2(sprintf("    WARNING: Simulation failed: %s\n",psimulation$error))
		stop()
	}
	cat2(sprintf("  Plotting tree LTT and dLTT based on loaded rates..\n"))
	plot_age_curves(file_basepath 		= sprintf("%s/fitted_HBD_model_LTT",tree_output_dir),
					data_type 			= "LTTs",
					case_tag 			= sprintf("%s, fitted dLTT",tree_name),
					curves 				= list(	list(empirical_LTT$ages, empirical_LTT$lineages),
												list(psimulation$ages, psimulation$LTT)),
					curve_names 		= c("tree LTT","fitted dLTT"),
					curve_colors 		= c(PLOT_COLOR_LTT,PLOT_COLOR_FIT_DLTT),
					curve_line_types	= c(1,2),
					curve_widths 		= c(PLOT_LW_LTT,PLOT_LW_FIT_DLTT),
					age_label 			= sprintf("age (%s)",time_unit),
					value_label 		= sprintf("lineages"),
					plot_log_values 	= TRUE,
					legend_pos 			= "outside",
					plot_title 			= sprintf("LTT of tree and dLTT of previously fitted model\n%s",tree_name),
					data_file_comments 	= sprintf("LTT of external tree '%s', and dLTT of previously fitted model",tree_name),
					verbose 			= FALSE,
					verbose_prefix 		= "    ");

	# construct congruent model
	mu_over_lambda = 0.9
	cat2(sprintf("  Simulating congruent model, by assuming that mu/lambda=%g..\n",mu_over_lambda))
	csimulation = simulate_deterministic_hbd(	LTT0			= Ntips,
												oldest_age		= tail(psimulation$ages,1),
												rho0			= rho0,
												age_grid		= psimulation$ages,
												PDR				= psimulation$PDR,
												lambda0			= psimulation$lambda[1],
												mu_over_lambda	= mu_over_lambda, # assume mu is close to the extinction rate
												splines_degree	= 1,
												relative_dt		= INTEGRATION_RELATIVE_DT);
	if(!csimulation$success){
		cat2(sprintf("    WARNING: Simulation of congruent model failed: %s\n",csimulation$error))
	}else{
		cat2(sprintf("  Plotting rates of previously fitted and congruent model..\n"));
		plot_age_curves(file_basepath 		= sprintf("%s/fitted_and_congruent_lambda_mu",tree_output_dir),
						data_type 			= "lambdas & mus",
						case_tag 			= sprintf("%s, fit & congruent",tree_name),
						curves 				= list(	list(psimulation$ages, psimulation$lambda),
													list(psimulation$ages, psimulation$mu),
													list(csimulation$ages, csimulation$lambda),
													list(csimulation$ages, csimulation$mu)),
						curve_names 		= c("lambda (fit)","mu (assumed)","lambda (congruent)","mu (congruent)"),
						curve_colors 		= c(PLOT_COLOR_LAMBDA,PLOT_COLOR_MU,PLOT_COLOR_LAMBDA,PLOT_COLOR_MU),
						curve_line_types	= c(1,1,2,2),
						curve_widths 		= c(1,1,2,2),
						age_label 			= sprintf("age (%s)",time_unit),
						value_label 		= sprintf("rate (1/%s)",time_unit),
						miny				= -0.1,
						plot_log_values 	= FALSE,
						legend_pos 			= "outside",
						plot_title 			= sprintf("Previously fitted & congruent model, lambda & mu\n%s",tree_name),
						data_file_comments 	= sprintf("Speciation and extinction rates for previously fitted and congruent model, tree '%s'",tree_name),
						verbose 			= FALSE,
						verbose_prefix 		= "    ")
		plot_age_curves(file_basepath 		= sprintf("%s/fitted_and_congruent_r",tree_output_dir),
						data_type 			= "diversification rates",
						case_tag 			= sprintf("%s, fit & congruent",tree_name),
						curves 				= list(	list(psimulation$ages, psimulation$diversification_rate),
													list(csimulation$ages, csimulation$diversification_rate)),
						curve_names 		= c("r (fit)","r (congruent)"),
						curve_colors 		= c(PLOT_COLOR_DIVERSIFICATION_RATE,PLOT_COLOR_DIVERSIFICATION_RATE),
						curve_line_types	= c(1,2),
						curve_widths 		= c(1,2),
						age_label 			= sprintf("age (%s)",time_unit),
						value_label 		= sprintf("net diversification rate (1/%s)",time_unit),
						plot_log_values 	= FALSE,
						legend_pos 			= "outside",
						plot_title 			= sprintf("Previously fitted & congruent model, net diversification rate\n%s",tree_name),
						data_file_comments 	= sprintf("Net diversification rates for previously fitted and congruent model, tree '%s'",tree_name),
						verbose 			= FALSE,
						verbose_prefix 		= "    ")
	}
	#BCO this section
  # Need to pad the rates: they only go up to 32MY, tree is 35 MY deep. So add one more interval.
	loglikelihood_hbd_steeman <- loglikelihood_hbd(
		tree=tree,
		oldest_age	= castor::get_tree_span(tree)$max_distance,
		rho0		 = Ntips/extant_N ,
		age_grid = c(precomputed_rates[,1], castor::get_tree_span(tree)$max_distance),
		lambda = c(precomputed_rates[,2], precomputed_rates[nrow(precomputed_rates),2]),
		mu = c(precomputed_rates[,3], precomputed_rates[nrow(precomputed_rates),3])
	)
	loglikelihood_hbd_congruent <- loglikelihood_hbd(
		tree=tree,
		oldest_age	= castor::get_tree_span(tree)$max_distance,
		rho0		 = Ntips/extant_N ,
		age_grid = c(csimulation$ages, castor::get_tree_span(tree)$max_distance),
		lambda = c(csimulation$lambda, csimulation$lambda[length(csimulation$lambda)]),
		mu = c(csimulation$mu, csimulation$mu[length(csimulation$mu)])
	)
	print("loglikelihood_hbd_steeman")
	print(loglikelihood_hbd_steeman)
	print("loglikelihood_hbd_congruent")
	print(loglikelihood_hbd_congruent)

# Now to optimize
library(nloptr)
# x are log lambda
negloglikelihood_hbd_for_nloptr <- function(x, ef=0.9, tree, age_grid, rho0=87/89, oldest_age=castor::get_tree_span(tree)$max_distance, badval=1e6) {
	lambda <- exp(x)
	mu <- lambda*ef
	result <- NA
	try(result <- loglikelihood_hbd(
		tree=tree,
		oldest_age	= oldest_age,
		rho0		 = rho0 ,
		age_grid = age_grid,
		lambda = lambda,
		mu= mu
	))
	if(result$success) {
		return(-result$loglikelihood)
	} else {
		return(badval)
	}
}
full_age_grid <- c(precomputed_rates[,1], castor::get_tree_span(tree)$max_distance)

full_lambda <- c(precomputed_rates[,2], precomputed_rates[nrow(precomputed_rates),2])

steeman_fixed_ef <- function(ef, full_lambda, tree, full_age_grid, rho0=87/89) {
	results <- list(
		nloptr(log(full_lambda), eval_f=negloglikelihood_hbd_for_nloptr, ef=ef, tree=tree, age_grid=full_age_grid, rho0=rho0, oldest_age=castor::get_tree_span(tree)$max_distance, badval=1e6, opts=list(print_level=0, algorithm="NLOPT_LN_NEWUOA", maxeval=1000)),
		nloptr(log((1+ef)*full_lambda), eval_f=negloglikelihood_hbd_for_nloptr, ef=ef, tree=tree, age_grid=full_age_grid, rho0=rho0, oldest_age=castor::get_tree_span(tree)$max_distance, badval=1e6, opts=list(print_level=0, algorithm="NLOPT_LN_NEWUOA", maxeval=1000)),
		nloptr(log(2*full_lambda), eval_f=negloglikelihood_hbd_for_nloptr, ef=ef, tree=tree, age_grid=full_age_grid, rho0=rho0, oldest_age=castor::get_tree_span(tree)$max_distance, badval=1e6, opts=list(print_level=0, algorithm="NLOPT_LN_NEWUOA", maxeval=1000))
	)
	best <- results[[which.min(lapply(results, "[[", 'objective'))[1]]]
	best$lambda <- exp(best$x0)
	best$mu <- ef*best$lambda
	best$netdiv <- best$lambda-best$mu
	best$turnover <- best$lambda+best$mu
	best$ef <- ef
	best$full_age_grid <- full_age_grid
	best$coarseness <- length(full_age_grid)
	return(best)
}

ef_tries <-seq(from=0, to=0.9, length.out=10)


fine_age_grid <- approx(x=full_age_grid, y=full_age_grid, n=100, rule=2)$y
fine_lambda <- approx(x=full_age_grid, y=full_lambda, xout=fine_age_grid, rule=2)$y

coarse_age_grid <- approx(x=full_age_grid, y=full_age_grid, n=20, rule=2)$y
coarse_lambda <- approx(x=full_age_grid, y=full_lambda, xout=coarse_age_grid, rule=2)$y

all_steeman_medium <- parallel::mclapply(ef_tries, steeman_fixed_ef, full_lambda=full_lambda, tree=tree, full_age_grid=full_age_grid, mc.cores=parallel::detectCores())

all_steeman_fine <- parallel::mclapply(ef_tries, steeman_fixed_ef, full_lambda=fine_lambda, tree=tree, full_age_grid=fine_age_grid, mc.cores=parallel::detectCores())

all_steeman_coarse <- parallel::mclapply(ef_tries, steeman_fixed_ef, full_lambda=coarse_lambda, tree=tree, full_age_grid=coarse_age_grid, mc.cores=parallel::detectCores())

all_steeman <- c(all_steeman_coarse, all_steeman_medium, all_steeman_fine)

summary_steeman <- data.frame(ef=unlist(lapply(all_steeman, "[[", "ef")), coarse=unlist(lapply(all_steeman, "[[", "coarseness")), neglogl=unlist(lapply(all_steeman, "[[", "objective")))
summary_steeman$deltaneglogl = summary_steeman$neglogl - min(summary_steeman$neglogl)


all.div <- unlist(lapply(all_steeman, "[[", "netdiv"))
pdf(file=sprintf("%s/%s",output_dir,"fitting.pdf"))
plot(x=rev(range(full_age_grid)), y=range(all.div), xlab='age (MYR)', ylab="net diversification rate (1/MY)", bty="n", type="n")
for (i in seq_along(all_steeman)) {
	lines(all_steeman[[i]]$full_age_grid, all_steeman[[i]]$netdiv)
}
dev.off()
#s90 <- steeman_fixed_ef(ef=.9, full_lambda=full_lambda, tree=tree, full_age_grid=full_age_grid)
#steeman_ef0.9 <- nloptr(log(full_lambda), eval_f=negloglikelihood_hbd_for_nloptr, ef=0.9, tree=tree, age_grid=full_age_grid, rho0=87/89, oldest_age=castor::get_tree_span(tree)$max_distance, badval=1e6, opts=list(print_level=1, algorithm="NLOPT_LN_NEWUOA", maxeval=1000))


#print(steeman_yule)

	# END BCO SECTION
}






###################################
# Analyze simulated trees


if(INCLUDE_SIM_MODEL_ANALYSIS){
	NM = length(SIM_MODELS);
	for(m in 1:NM){
		model = SIM_MODELS[[m]];
		cat2(sprintf("Examining simulation model '%s' (#%d)..\n",model$name,m));
		model_output_dir=sprintf("%s/%s",output_dir,model$name);
		dir.create(model_output_dir, showWarnings = FALSE, recursive=TRUE);

		# seed random number generator
		if(is.null(model$random_seed)){
			# no seed specified for this tree, so pick one randomly based on the global seed
			set.seed(global_random_seed)
			model$random_seed = sample.int(n=1000000,size=1)
		}
		set.seed(model$random_seed)
		cat2(sprintf("  Note: Seeding random generator at %d\n",model$random_seed));

		# simulate random realization of the model (generate tree)
		cat2(sprintf("  Simulating random tree..\n"));
		time_unit	= "Myr"
		rho0 		= model$rho0
		time_grid 	= seq(from=0,to=model$max_time,length.out=1000)
		parameters 	= list(birth_rate_factor=model$base_lambda, death_rate_factor=model$base_mu, rarefaction=rho0)
		tree_sim	= generate_random_tree(parameters, max_time=model$max_time, max_tips=model$max_tips, coalescent=TRUE, added_rates_times=time_grid, added_birth_rates_pc=model$added_lambda(time_grid), added_death_rates_pc=model$added_mu(time_grid), include_birth_times=TRUE, include_death_times=TRUE)
		if(!tree_sim$success) stop(sprintf("ERROR: Tree generation failed: %s",tree_sim$error))
		tree		= tree_sim$tree
		Ntips 		= length(tree$tip.label)
		root_age 	= castor::get_tree_span(tree)$max_distance
		cat2(sprintf("  Note: Tree has %d tips at final time=%g\n",Ntips,tree_sim$final_time))
		castor::write_tree(tree, file=sprintf("%s/tree.tre",model_output_dir))

		# define lambda & mu as functions of age
		age_grid = seq(from=0,to=root_age*1.1,length.out=1000)
		age2lambda = function(age){
			return(model$base_lambda+model$added_lambda(tree_sim$final_time-age))
		}
		age2mu = function(age){
			return(model$base_mu+model$added_mu(tree_sim$final_time-age))
		}
		lambda0 = age2lambda(0)

		cat2(sprintf("  Calculating empirical LTT..\n"))
		empirical_LTT 		= castor::count_lineages_through_time(tree,Ntimes=1000);
		empirical_LTT$ages 	= root_age - empirical_LTT$times
		empirical_N			= castor:::get_diversities_from_birth_and_death_events_CPP(empirical_LTT$times, birth_times=tree_sim$birth_times-tree_sim$root_time, death_times=tree_sim$death_times-tree_sim$root_time, start_diversity=1, Nsplits=2)$diversities

		cat2(sprintf("  Simulating deterministic HBD model..\n"))
		simulation = simulate_deterministic_hbd(LTT0			= Ntips,
												oldest_age		= root_age,
												rho0			= rho0,
												age_grid		= age_grid,
												lambda			= age2lambda(age_grid),
												mu				= age2mu(age_grid),
												splines_degree	= 1,
												relative_dt		= INTEGRATION_RELATIVE_DT);
		if(!simulation$success){
			stop(sprintf("ERROR: HBD simulation failed: %s\n",simulation$error))
		}

		cat2(sprintf("  Plotting tree properties..\n"))
		plot_age_curves(	file_basepath 		= sprintf("%s/tree_LTTs",model_output_dir),
							data_type 			= "tree diversity curves",
							case_tag 			= sprintf("%s - random tree",model$name),
							curves = list(	list(simulation$ages, simulation$LTT),
											list(empirical_LTT$ages, empirical_LTT$lineages),
											list(empirical_LTT$ages, empirical_N)),
							curve_names 		= c("dLTT", "tree LTT", "N"),
							curve_colors 		= c(PLOT_COLOR_LTT, PLOT_COLOR_LTT, PLOT_COLOR_TOTAL_DIVERSITY),
							curve_line_types	= c(PLOT_LINE_TYPE_DETERMINISTIC_TRUE, PLOT_LINE_TYPE_EMPIRICAL, PLOT_LINE_TYPE_EMPIRICAL),
							curve_widths 		= c(1,1,2),
							plot_curves 		= NULL,
							max_age 			= root_age,
							Nages				= PLOT_DOWNSAMPLING_RESOLUTION,
							age_label 			= "age",
							value_label 		= "diversities",
							plot_log_values 	= TRUE,
							legend_pos 			= "outside",
							plot_title 			= "Diversities of generated tree",
							data_file_comments 	= sprintf("# Diversities of randomly generated tree (%s, Ntips=%d)",model$name,Ntips),
							verbose 			= FALSE,
							verbose_prefix 		= "    ")

		plot_age_curves(	file_basepath 		= sprintf("%s/tree_diversification_rates",model_output_dir),
							data_type 			= "tree diversification rates",
							case_tag 			= sprintf("%s - random tree",model$name),
							curves = list(	list(simulation$ages, simulation$diversification_rate),
											list(simulation$ages, simulation$PDR)),
							curve_names 		= c("r", "PDR"),
							curve_colors 		= c(PLOT_COLOR_DIVERSIFICATION_RATE, PLOT_COLOR_PDR),
							curve_line_types	= c(PLOT_LINE_TYPE_DETERMINISTIC_TRUE, PLOT_LINE_TYPE_DETERMINISTIC_TRUE),
							curve_widths 		= c(1,2),
							plot_curves 		= NULL,
							max_age 			= root_age,
							Nages				= PLOT_DOWNSAMPLING_RESOLUTION,
							age_label 			= "age",
							value_label 		= "rate",
							plot_log_values 	= FALSE,
							legend_pos 			= "outside",
							plot_title 			= "Diversification rates underlying tree",
							data_file_comments 	= sprintf("# Diversification rates used to generate tree (%s, Ntips=%d)",model$name,Ntips),
							verbose 			= FALSE,
							verbose_prefix 		= "    ")

		plot_age_curves(	file_basepath 		= sprintf("%s/tree_lambdas_mus",model_output_dir),
							data_type 			= "tree speciation and extinction rates",
							case_tag 			= sprintf("%s - random tree",model$name),
							curves = list(	list(age_grid, age2lambda(age_grid)),
											list(age_grid, age2mu(age_grid))),
							curve_names 		= c("lambda", "mu"),
							curve_colors 		= c(PLOT_COLOR_LAMBDA, PLOT_COLOR_MU),
							curve_line_types	= c(PLOT_LINE_TYPE_DETERMINISTIC_TRUE, PLOT_LINE_TYPE_DETERMINISTIC_TRUE),
							curve_widths 		= c(1,2),
							plot_curves 		= NULL,
							max_age 			= root_age,
							Nages				= PLOT_DOWNSAMPLING_RESOLUTION,
							age_label 			= "age",
							value_label 		= "rate",
							plot_log_values 	= FALSE,
							legend_pos 			= "outside",
							plot_title 			= "Speciation/extinction rates underlying tree",
							data_file_comments 	= sprintf("# Speciation/extinction rates used to generate tree (%s, Ntips=%d)",model$name,Ntips),
							verbose 			= FALSE,
							verbose_prefix 		= "    ")



		# create an alternative but congruent model, with a modified extinction rate
		if(model$include_model_comparison && (!is.null(model$congruent_added_mus)) && (length(model$congruent_added_mus)>0)){
			for(ccounter in 1:length(model$congruent_added_mus)){
				cat2(sprintf("  Simulating congruent model #%d to %s..\n",ccounter,model$name))
				added_mu2 = model$congruent_added_mus[[ccounter]]
				modified_age2mu = function(age){
					return(pmax(0,model$base_mu + added_mu2(tree_sim$final_time-age)))
				}
				simulation2 = simulate_deterministic_hbd(	LTT0			= Ntips,
															oldest_age		= root_age,
															rho0			= rho0,
															age_grid		= simulation$ages,
															mu				= modified_age2mu(simulation$ages),
															PDR				= simulation$PDR,
															lambda0			= lambda0,
															splines_degree	= 1,
															relative_dt		= INTEGRATION_RELATIVE_DT);
				if(!simulation2$success){
					cat2(sprintf("    ERROR: Modified HBD simulation failed: %s\n",simulation2$error))
					next;
				}
				cat2(sprintf("  Plotting comparisons between congruent models..\n"))
				plot_age_curves(	file_basepath 		= sprintf("%s/congruent_model_%d_LTT_N",model_output_dir,ccounter),
									data_type 			= "LTTs",
									case_tag 			= sprintf("%s and congruent model #%d",model$name,ccounter),
									curves = list(	list(simulation$ages, simulation$LTT),
													list(simulation2$ages, simulation2$LTT),
													list(simulation$ages, simulation$total_diversity),
													list(simulation2$ages, simulation2$total_diversity)),
									curve_names 		= c("LTT 0", sprintf("LTT %d",ccounter), "N 0",  sprintf("N %d",ccounter)),
									curve_colors 		= c(PLOT_COLOR_LTT, PLOT_COLOR_LTT, PLOT_COLOR_TOTAL_DIVERSITY, PLOT_COLOR_TOTAL_DIVERSITY),
									curve_line_types	= c(1, 2, 1, 2),
									curve_widths 		= c(1, 1, 2, 2),
									plot_curves 		= NULL,
									max_age 			= root_age,
									Nages				= PLOT_DOWNSAMPLING_RESOLUTION,
									age_label 			= "age",
									value_label 		= "lineages",
									plot_log_values 	= TRUE,
									legend_pos 			= "outside",
									plot_title 			= sprintf("LTT and N predicted by models 0 & %d",ccounter),
									data_file_comments 	= sprintf("# Deterministic LTT and total diversity predicted by two congruent models (%s and modified #%d)",model$name,ccounter),
									verbose 			= FALSE,
									verbose_prefix 		= "    ")

				plot_age_curves(	file_basepath 		= sprintf("%s/congruent_model_%d_PSR",model_output_dir,ccounter),
									data_type 			= "PSRs",
									case_tag 			= sprintf("%s and congruent model #%d",model$name,ccounter),
									curves = list(	list(simulation$ages, simulation$PSR),
													list(simulation2$ages, simulation2$PSR)),
									curve_names 		= c("PSR 0", sprintf("PSR %d",ccounter)),
									curve_colors 		= c(PLOT_COLOR_PSR, PLOT_COLOR_PSR),
									curve_line_types	= c(1, 2),
									curve_widths 		= c(1, 2),
									plot_curves 		= NULL,
									max_age 			= root_age,
									Nages				= PLOT_DOWNSAMPLING_RESOLUTION,
									age_label 			= "age",
									value_label 		= "rate",
									plot_log_values 	= TRUE,
									legend_pos 			= "outside",
									plot_title 			= sprintf("Pulled speciation rates in models 0 & %d",ccounter),
									data_file_comments 	= sprintf("# Pulled speciation rates (PSR) in two congruent models (%s and modified #%d)",model$name,ccounter),
									verbose 			= FALSE,
									verbose_prefix 		= "    ")

				plot_age_curves(	file_basepath 		= sprintf("%s/congruent_model_%d_diversification_rates",model_output_dir,ccounter),
									data_type 			= "rates",
									case_tag 			= sprintf("%s and congruent model #%d",model$name,ccounter),
									curves = list(	list(simulation$ages, simulation$diversification_rate),
													list(simulation2$ages, simulation2$diversification_rate),
													list(simulation$ages, simulation$PDR),
													list(simulation2$ages, simulation2$PDR)),
									curve_names 		= c("r 0", sprintf("r %d",ccounter), "PDR 0", sprintf("PDR %d",ccounter)),
									curve_colors 		= c(PLOT_COLOR_DIVERSIFICATION_RATE, PLOT_COLOR_DIVERSIFICATION_RATE, PLOT_COLOR_PDR, PLOT_COLOR_PDR),
									curve_line_types	= c(1, 2, 1, 2),
									curve_widths 		= c(1, 1, 2, 2),
									plot_curves 		= NULL,
									max_age 			= root_age,
									Nages				= PLOT_DOWNSAMPLING_RESOLUTION,
									age_label 			= "age",
									value_label 		= "rate",
									plot_log_values 	= FALSE,
									legend_pos 			= "outside",
									plot_title 			= sprintf("diversification rates in models 0 & %d",ccounter),
									data_file_comments 	= sprintf("# Diversification rates in two congruent models (%s and modified #%d)",model$name,ccounter),
									verbose 			= FALSE,
									verbose_prefix 		= "    ")

				plot_age_curves(	file_basepath 		= sprintf("%s/congruent_model_%d_lambdas_mus",model_output_dir,ccounter),
									data_type 			= "rates",
									case_tag 			= sprintf("%s and congruent model #%d",model$name,ccounter),
									curves = list(	list(simulation$ages, simulation$lambda),
													list(simulation2$ages, simulation2$lambda),
													list(simulation$ages, simulation$mu),
													list(simulation2$ages, simulation2$mu)),
									curve_names 		= c("lambda 0", sprintf("lambda %d",ccounter), "mu 0", sprintf("mu %d",ccounter)),
									curve_colors 		= c(PLOT_COLOR_LAMBDA, PLOT_COLOR_LAMBDA, PLOT_COLOR_MU, PLOT_COLOR_MU),
									curve_line_types	= c(1, 2, 1, 2),
									curve_widths 		= c(1, 1, 2, 2),
									plot_curves 		= NULL,
									max_age 			= root_age,
									Nages				= PLOT_DOWNSAMPLING_RESOLUTION,
									age_label 			= "age",
									value_label 		= "rate",
									plot_log_values 	= FALSE,
									legend_pos 			= "outside",
									plot_title 			= sprintf("rates in models 0 & %d",ccounter),
									data_file_comments 	= sprintf("# Speciation/extinction rates in two congruent models (%s and modified #%d)",model$name,ccounter),
									verbose 			= FALSE,
									verbose_prefix 		= "    ");
			}
		}

		# fit PDR to tree
		if(model$include_PDR_grid_fitting){
			cat2(sprintf("  Fitting PDR to simulated timetree '%s'..\n",model$name))
			if(Ntips < model$fitting_min_lineages){
				cat2(sprintf("    WARNING: Tree size (%d) is below the threshold for HBD fitting (%d). Skipping\n",Ntips,model$fitting_min_lineages))
			}else{
				oldest_age = root_age - empirical_LTT$times[which(empirical_LTT$lineages>=model$fitting_min_lineages)[1]]
				cat2(sprintf("    Note: Max considered age = %g, root_age = %g\n", oldest_age, root_age));
				for(gs in 1:length(model$Nfitting_ages)){
					Nfitting_ages = model$Nfitting_ages[gs]
					cat2(sprintf("    Fitting with age-grid size %d..\n",Nfitting_ages));
					age_grid = seq(from=0,to=oldest_age,length.out=Nfitting_ages)
					fit = fit_hbd_pdr_on_grid(	tree,
												oldest_age			= oldest_age,
												age_grid			= age_grid,
												guess_rholambda0	= rho0*lambda0,
												guess_PDR			= approx(x=simulation$ages, y=simulation$PDR, xout=age_grid)$y,
												relative_dt			= INTEGRATION_RELATIVE_DT,
												Ntrials				= DEFAULT_FITTING_NTRIALS,
												Nthreads			= NUMBER_OF_PARALLEL_THREADS,
												max_model_runtime	= max(1,Ntips/FITTING_TIPS_PER_RUNTIME_SECOND),
												splines_degree		= 1,
												fit_control			= list(eval.max=FITTING_NEVALUATIONS, iter.max=FITTING_NITERATIONS, rel.tol=FITTING_REL_TOLERANCE));
					if(!fit$success){
						cat2(sprintf("    ERROR: Fitting failed: %s\n",fit$error));
					}else{
						save_object_to_file(fit, sprintf("%s/fitted_HBD_PDR/grid_size_%d/fit_results.txt",model_output_dir,Nfitting_ages))
						cat2(sprintf("    Simulating fitted congruence class '%s'..\n",model$name));
						fsimulation = simulate_deterministic_hbd(	LTT0			= Ntips,
																	oldest_age		= oldest_age,
																	rho0			= rho0,
																	age_grid		= fit$age_grid,
																	mu				= 0, # pick an arbitrary member of the fitted equivalence class
																	PDR				= fit$fitted_PDR,
																	lambda0			= fit$fitted_rholambda0/rho0,
																	splines_degree	= 1,
																	relative_dt		= INTEGRATION_RELATIVE_DT);
						if(!fsimulation$success) stop(sprintf("ERROR: Simulation failed: %s\n",fsimulation$error));

						cat2(sprintf("    Plotting fitted congruence class '%s'..\n",model$name));
						plot_age_curves(file_basepath 		= sprintf("%s/fitted_HBD_PDR/grid_size_%d/fitted_class_LTTs",model_output_dir,Nfitting_ages),
										data_type 			= "LTTs",
										case_tag 			= sprintf("%s - fitted class, grid size %d",model$name,Nfitting_ages),
										curves = list(	list(fsimulation$ages, fsimulation$LTT),
														list(simulation$ages, simulation$LTT),
														list(empirical_LTT$ages, empirical_LTT$lineages)),
										curve_names 		= c("fitted deterministic LTT", "true deterministic LTT", "empirical LTT"),
										curve_colors 		= c(PLOT_COLOR_LTT, PLOT_COLOR_LTT, PLOT_COLOR_LTT),
										curve_line_types	= c(PLOT_LINE_TYPE_DETERMINISTIC_FITTED, PLOT_LINE_TYPE_DETERMINISTIC_TRUE, PLOT_LINE_TYPE_EMPIRICAL),
										curve_widths 		= c(1,2,3),
										max_age 			= root_age,
										Nages				= PLOT_DOWNSAMPLING_RESOLUTION,
										age_label 			= "age",
										value_label 		= "lineages",
										plot_log_values 	= TRUE,
										legend_pos 			= "outside",
										plot_title 			= sprintf("LTT (fitted class vs tree)\ngrid size %d",Nfitting_ages),
										data_file_comments 	= sprintf("# LTTs in tree and fitted HBD-class (%s), fitting age-grid size %d",model$name,Nfitting_ages),
										verbose 			= FALSE,
										verbose_prefix 		= "      ");
						plot_age_curves(file_basepath 		= sprintf("%s/fitted_HBD_PDR/grid_size_%d/fitted_class_PDRs",model_output_dir,Nfitting_ages),
										data_type 			= "PDRs",
										case_tag 			= sprintf("%s - fitted class, grid size %d",model$name,Nfitting_ages),
										curves = list(	list(fsimulation$ages, fsimulation$PDR),
														list(simulation$ages, simulation$PDR)),
										curve_names 		= c("fitted PDR", "true PDR"),
										curve_colors 		= c(PLOT_COLOR_PDR, PLOT_COLOR_PDR),
										curve_line_types	= c(PLOT_LINE_TYPE_DETERMINISTIC_FITTED, PLOT_LINE_TYPE_DETERMINISTIC_TRUE),
										curve_widths 		= 2,
										max_age 			= root_age,
										Nages				= PLOT_DOWNSAMPLING_RESOLUTION,
										age_label 			= "age",
										value_label 		= "rate",
										plot_log_values 	= FALSE,
										legend_pos 			= "outside",
										plot_title 			= sprintf("PDRs (fitted class vs true)\ngrid size %d",Nfitting_ages),
										data_file_comments 	= sprintf("# PDR fitted to tree, model class '%s', fitting age-grid size %d",model$name,Nfitting_ages),
										verbose 			= FALSE,
										verbose_prefix 		= "      ");
					}
				}
			}
		}



		# fit PSR to tree
		if(model$include_PSR_grid_fitting){
			cat2(sprintf("  Fitting PSR to simulated timetree '%s'..\n",model$name))
			if(Ntips < model$fitting_min_lineages){
				cat2(sprintf("    WARNING: Tree size (%d) is below the threshold for HBD fitting (%d). Skipping\n",Ntips,model$fitting_min_lineages))
			}else{
				oldest_age = root_age - empirical_LTT$times[which(empirical_LTT$lineages>=model$fitting_min_lineages)[1]]
				cat2(sprintf("    Note: Max considered age = %g, root_age = %g\n", oldest_age, root_age));
				for(gs in 1:length(model$Nfitting_ages)){
					Nfitting_ages = model$Nfitting_ages[gs]
					cat2(sprintf("    Fitting with age-grid size %d..\n",Nfitting_ages));
					age_grid = seq(from=0,to=oldest_age,length.out=Nfitting_ages)
					fit = fit_hbd_psr_on_grid(	tree,
												oldest_age			= oldest_age,
												age_grid			= age_grid,
												guess_PSR			= approx(x=simulation$ages, y=simulation$PDR, xout=age_grid)$y,
												relative_dt			= INTEGRATION_RELATIVE_DT,
												Ntrials				= DEFAULT_FITTING_NTRIALS,
												Nthreads			= NUMBER_OF_PARALLEL_THREADS,
												max_model_runtime	= max(1,Ntips/FITTING_TIPS_PER_RUNTIME_SECOND),
												splines_degree		= 1,
												fit_control			= list(eval.max=FITTING_NEVALUATIONS, iter.max=FITTING_NITERATIONS, rel.tol=FITTING_REL_TOLERANCE));
					if(!fit$success){
						cat2(sprintf("    ERROR: Fitting failed: %s\n",fit$error));
					}else{
						cat2(sprintf("    --> Loglikelihood %.10g, converged %s\n",fit$loglikelihood,fit$converged));
						save_object_to_file(fit, sprintf("%s/fitted_HBD_PSR/grid_size_%d/fit_results.txt",model_output_dir,Nfitting_ages))
						cat2(sprintf("    Plotting dLTT of fitted congruence class '%s'..\n",model$name));
						refined_age_grid = seq(from=0,to=oldest_age,length.out=PLOT_DOWNSAMPLING_RESOLUTION)
						fitted_LTT = Ntips * exp(-castor:::get_antiderivative_of_splines_function(Xgrid=fit$age_grid, Xstart=0, Ygrid=fit$fitted_PSR, splines_degree=1, Xtarget=refined_age_grid));
						plot_age_curves(file_basepath 		= sprintf("%s/fitted_HBD_PSR/grid_size_%d/fitted_class_LTTs",model_output_dir,Nfitting_ages),
										data_type 			= "LTTs",
										case_tag 			= sprintf("%s - fitted class, grid size %d",model$name,Nfitting_ages),
										curves = list(	list(refined_age_grid, fitted_LTT),
														list(simulation$ages, simulation$LTT),
														list(empirical_LTT$ages, empirical_LTT$lineages)),
										curve_names 		= c("fitted deterministic LTT", "true deterministic LTT", "empirical LTT"),
										curve_colors 		= c(PLOT_COLOR_LTT, PLOT_COLOR_LTT, PLOT_COLOR_LTT),
										curve_line_types	= c(PLOT_LINE_TYPE_DETERMINISTIC_FITTED, PLOT_LINE_TYPE_DETERMINISTIC_TRUE, PLOT_LINE_TYPE_EMPIRICAL),
										curve_widths 		= c(1,2,3),
										max_age 			= root_age,
										Nages				= PLOT_DOWNSAMPLING_RESOLUTION,
										age_label 			= "age",
										value_label 		= "lineages",
										plot_log_values 	= TRUE,
										legend_pos 			= "outside",
										plot_title 			= sprintf("LTT (fitted class vs tree)\ngrid size %d",Nfitting_ages),
										data_file_comments 	= sprintf("# LTTs in tree and fitted HBD congruence class (%s), fitting age-grid size %d",model$name,Nfitting_ages),
										verbose 			= FALSE,
										verbose_prefix 		= "      ");

						cat2(sprintf("    Plotting PSR of fitted congruence class '%s'..\n",model$name));
						plot_age_curves(file_basepath 		= sprintf("%s/fitted_HBD_PSR/grid_size_%d/fitted_class_PSR",model_output_dir,Nfitting_ages),
										data_type 			= "PSRs",
										case_tag 			= sprintf("%s - fitted class, grid size %d",model$name,Nfitting_ages),
										curves = list(	list(fit$age_grid, fit$fitted_PSR),
														list(simulation$ages, simulation$PSR)),
										curve_names 		= c("fitted PSR", "true PSR"),
										curve_colors 		= c(PLOT_COLOR_PSR, PLOT_COLOR_PSR),
										curve_line_types	= c(PLOT_LINE_TYPE_DETERMINISTIC_FITTED, PLOT_LINE_TYPE_DETERMINISTIC_TRUE),
										curve_widths 		= 2,
										max_age 			= root_age,
										age_label 			= "age",
										value_label 		= "rate",
										plot_log_values 	= FALSE,
										legend_pos 			= "outside",
										plot_title 			= sprintf("PSR (fitted class vs tree)\ngrid size %d",Nfitting_ages),
										data_file_comments 	= sprintf("# PSR in tree and fitted HBD congruence class (%s), fitting age-grid size %d",model$name,Nfitting_ages),
										verbose 			= FALSE,
										verbose_prefix 		= "      ");
					}
				}
			}
		}


		# fit HBD-model to tree
		if(model$include_grid_model_fitting){
			cat2(sprintf("  Fitting HBD-model to simulated timetree '%s'..\n",model$name))
			if(Ntips < model$fitting_min_lineages){
				cat2(sprintf("    WARNING: Tree size (%d) is below the threshold for HBD fitting (%d). Skipping\n",Ntips,model$fitting_min_lineages))
			}else{
				oldest_age = root_age - empirical_LTT$times[which(empirical_LTT$lineages>=model$fitting_min_lineages)[1]]
				cat2(sprintf("    Note: Max considered age = %g, root_age = %g\n", oldest_age, root_age));
				for(gs in 1:length(model$Nfitting_ages)){
					Nfitting_ages = model$Nfitting_ages[gs]
					cat2(sprintf("    Fitting model with age-grid size %d..\n",Nfitting_ages));
					fit = fit_hbd_model_on_grid(tree,
												oldest_age			= oldest_age,
												age_grid			= seq(from=0,to=oldest_age,length.out=Nfitting_ages),
												fixed_rho0			= rho0,
												splines_degree		= 1,
												relative_dt			= INTEGRATION_RELATIVE_DT,
												Ntrials				= DEFAULT_FITTING_NTRIALS,
												Nthreads			= NUMBER_OF_PARALLEL_THREADS,
												max_model_runtime	= max(1,Ntips/FITTING_TIPS_PER_RUNTIME_SECOND),
												fit_control			= list(eval.max=FITTING_NEVALUATIONS, iter.max=FITTING_NITERATIONS, rel.tol=FITTING_REL_TOLERANCE))
					if(!fit$success){
						cat2(sprintf("    ERROR: Fitting failed: %s\n",fit$error));
					}else{
						save_object_to_file(fit, sprintf("%s/fitted_HBD_model/grid_size_%d/fit_results.txt",model_output_dir,Nfitting_ages))
						cat2(sprintf("    Simulating fitted model '%s'..\n",model$name));
						fsimulation = simulate_deterministic_hbd(	LTT0			= Ntips,
																	oldest_age		= oldest_age,
																	rho0			= rho0,
																	age_grid		= fit$age_grid,
																	lambda			= fit$fitted_lambda,
																	mu				= fit$fitted_mu,
																	splines_degree	= 1,
																	relative_dt		= INTEGRATION_RELATIVE_DT);
						if(!fsimulation$success){
							cat2(sprintf("ERROR: Simulation failed: %s\n",fsimulation$error))
						}else{
							cat2(sprintf("    Plotting fitted model '%s'..\n",model$name));
							plot_age_curves(file_basepath 		= sprintf("%s/fitted_HBD_model/grid_size_%d/fitted_dLTT",model_output_dir,Nfitting_ages),
											data_type 			= "LTTs",
											case_tag 			= sprintf("%s, fitted dLTT",model$name),
											curves 				= list(	list(empirical_LTT$ages, empirical_LTT$lineages),
																		list(simulation$ages, simulation$LTT),
																		list(fsimulation$ages, fsimulation$LTT)),
											curve_names 		= c("tree LTT","true dLTT","fitted dLTT"),
											curve_colors 		= c(PLOT_COLOR_LTT,PLOT_COLOR_LTT,PLOT_COLOR_FIT_DLTT),
											curve_line_types	= c(1,2,2),
											curve_widths 		= c(1,1,2),
											max_age 			= oldest_age,
											age_label 			= sprintf("age (%s)",time_unit),
											value_label 		= sprintf("lineages"),
											plot_log_values 	= TRUE,
											legend_pos 			= "outside",
											plot_title 			= sprintf("LTT of tree and fitted dLTT\n%s",model$name),
											data_file_comments 	= sprintf("LTT of simulated tree '%s', and fitted dLTT",model$name),
											verbose 			= FALSE,
											verbose_prefix 		= "    ")
							plot_age_curves(file_basepath 		= sprintf("%s/fitted_HBD_model/grid_size_%d/fitted_lambda_mu",model_output_dir,Nfitting_ages),
											data_type 			= "lambdas & mus",
											case_tag 			= sprintf("%s, fitted lambda & mu",model$name),
											curves 				= list(	list(simulation$ages, simulation$lambda),
																		list(simulation$ages, simulation$mu),
																		list(fsimulation$ages, fsimulation$lambda),
																		list(fsimulation$ages, fsimulation$mu)),
											curve_names 		= c("lambda (true)","mu (true)","lambda (fit)","mu (fit)"),
											curve_colors 		= c(PLOT_COLOR_LAMBDA,PLOT_COLOR_MU,PLOT_COLOR_LAMBDA,PLOT_COLOR_MU),
											curve_line_types	= c(1,1,2,2),
											curve_widths 		= c(1,1,2,2),
											max_age 			= oldest_age,
											age_label 			= sprintf("age (%s)",time_unit),
											value_label 		= sprintf("rate (1/%s)",time_unit),
											plot_log_values 	= FALSE,
											legend_pos 			= "outside",
											plot_title 			= sprintf("True and fitted lambda & mu\n%s",model$name),
											data_file_comments 	= sprintf("True and fitted lambda & mu, '%s'",model$name),
											verbose 			= FALSE,
											verbose_prefix 		= "    ")
						}
					}
				}
			}
		}

		# fit HBD-lambda to tree (assume mu is known)
		if(model$include_grid_lambda_fitting){
			cat2(sprintf("  Fitting HBD-lambda to simulated timetree '%s' (assuming known mu)..\n",model$name))
			if(Ntips < model$fitting_min_lineages){
				cat2(sprintf("    WARNING: Tree size (%d) is below the threshold for HBD fitting (%d). Skipping\n",Ntips,model$fitting_min_lineages))
			}else{
				oldest_age = root_age - empirical_LTT$times[which(empirical_LTT$lineages>=model$fitting_min_lineages)[1]]
				cat2(sprintf("    Note: Max considered age = %g, root_age = %g\n", oldest_age, root_age));
				for(gs in 1:length(model$Nfitting_ages)){
					Nfitting_ages = model$Nfitting_ages[gs]
					cat2(sprintf("    Fitting lambda with age-grid size %d..\n",Nfitting_ages));
					age_grid = seq(from=0,to=oldest_age,length.out=Nfitting_ages)
					fit = fit_hbd_model_on_grid(tree,
												oldest_age			= oldest_age,
												age_grid			= age_grid,
												fixed_rho0			= rho0,
												fixed_mu			= age2mu(age_grid),
												splines_degree		= 1,
												relative_dt			= INTEGRATION_RELATIVE_DT,
												Ntrials				= DEFAULT_FITTING_NTRIALS,
												Nthreads			= NUMBER_OF_PARALLEL_THREADS,
												max_model_runtime	= max(1,Ntips/FITTING_TIPS_PER_RUNTIME_SECOND),
												fit_control			= list(eval.max=FITTING_NEVALUATIONS, iter.max=FITTING_NITERATIONS, rel.tol=FITTING_REL_TOLERANCE))
					if(!fit$success){
						cat2(sprintf("    ERROR: Fitting failed: %s\n",fit$error));
					}else{
						save_object_to_file(fit, sprintf("%s/fitted_HBD_lambda/grid_size_%d/fit_results.txt",model_output_dir,Nfitting_ages))
						cat2(sprintf("    Simulating fitted model '%s'..\n",model$name));
						fsimulation = simulate_deterministic_hbd(	LTT0			= Ntips,
																	oldest_age		= oldest_age,
																	rho0			= rho0,
																	age_grid		= fit$age_grid,
																	lambda			= fit$fitted_lambda,
																	mu				= fit$fitted_mu,
																	splines_degree	= 1,
																	relative_dt		= INTEGRATION_RELATIVE_DT);
						if(!fsimulation$success){
							cat2(sprintf("ERROR: Simulation failed: %s\n",fsimulation$error))
						}else{
							cat2(sprintf("    Plotting fitted model '%s'..\n",model$name));
							plot_age_curves(file_basepath 		= sprintf("%s/fitted_HBD_lambda/grid_size_%d/fitted_dLTT",model_output_dir,Nfitting_ages),
											data_type 			= "LTTs",
											case_tag 			= sprintf("%s, fitted dLTT",model$name),
											curves 				= list(	list(empirical_LTT$ages, empirical_LTT$lineages),
																		list(simulation$ages, simulation$LTT),
																		list(fsimulation$ages, fsimulation$LTT)),
											curve_names 		= c("tree LTT","true dLTT","fitted dLTT"),
											curve_colors 		= c(PLOT_COLOR_LTT,PLOT_COLOR_LTT,PLOT_COLOR_FIT_DLTT),
											curve_line_types	= c(1,2,2),
											curve_widths 		= c(1,1,2),
											max_age 			= oldest_age,
											age_label 			= sprintf("age (%s)",time_unit),
											value_label 		= sprintf("lineages"),
											plot_log_values 	= TRUE,
											legend_pos 			= "outside",
											plot_title 			= sprintf("LTT of tree and fitted dLTT\n%s",model$name),
											data_file_comments 	= sprintf("LTT of simulated tree '%s', and fitted dLTT",model$name),
											verbose 			= FALSE,
											verbose_prefix 		= "    ")

							plot_age_curves(file_basepath 		= sprintf("%s/fitted_HBD_lambda/grid_size_%d/fitted_lambda_mu",model_output_dir,Nfitting_ages),
											data_type 			= "lambdas",
											case_tag 			= sprintf("%s, fitted lambda",model$name),
											curves 				= list(	list(simulation$ages, simulation$lambda),
																		list(fsimulation$ages, fsimulation$lambda)),
											curve_names 		= c("lambda (true)","lambda (fit)"),
											curve_colors 		= c(PLOT_COLOR_LAMBDA,PLOT_COLOR_LAMBDA),
											curve_line_types	= c(1,2),
											curve_widths 		= c(1,2),
											max_age 			= oldest_age,
											age_label 			= sprintf("age (%s)",time_unit),
											value_label 		= sprintf("rate (1/%s)",time_unit),
											plot_log_values 	= FALSE,
											legend_pos 			= "outside",
											plot_title 			= sprintf("True and fitted lambda, when mu is known\n%s",model$name),
											data_file_comments 	= sprintf("True and fitted lambda, when mu is known, '%s'",model$name),
											verbose 			= FALSE,
											verbose_prefix 		= "    ")
						}
					}
				}
			}
		}
	}
}
save(list=ls(), file=sprintf("%s/%s",output_dir,"everything.rda"))
cat2(sprintf("Done. All outputs were written to '%s'\n",output_dir));
