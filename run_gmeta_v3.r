# runtime
runtime <- proc.time()

run.gmeta <- function(ma.set, ma.type, data.type) {
	
	if(ma.type == 'big' & data.type == 'binary') {

		########## BIG, BINARY ##########

		#####
		# create skeleton dataframes to be filled in later
		#####

		# fixed-effects

		fixed.methods <- c("fixed-mle", "fixed-robust1")

		fixed.nrows <- length(ma.id)*length(fixed.methods)

		fixed.ma.id <- sort(rep(ma.id, length=fixed.nrows))

		fixed.effects.bb <<- data.frame(ma_id=fixed.ma.id, 
									    method=rep(fixed.methods, length=fixed.nrows),
									    mean=rep(NA, length=fixed.nrows),
									    median=rep(NA, length=fixed.nrows),
									    stdev=rep(NA, length=fixed.nrows),
									    ci_lower=rep(NA, length=fixed.nrows),
									    ci_upper=rep(NA, length=fixed.nrows),
									    ma_type=rep("big", length=fixed.nrows),
									    data_type=rep("bin", length=fixed.nrows))

		# random-effects
		random.methods <- c("random-mm", "random-reml", "random-tau2", "random-robust1")
		tau2.methods <- c("DL", "HS", "SJ", "HE", "ML", "REML", "EB")

		len.tau2.methods <- length(tau2.methods)
		random.nrows <- length(ma.id)*length(random.methods)*len.tau2.methods

		random.methods.list <- character()

		for(each.rm in random.methods) {
			random.methods.list <- c(random.methods.list, 
				                     rep(each.rm, length=len.tau2.methods))
		}

		random.ma.id <- sort(rep(ma.id, length=random.nrows))

		random.effects.bb <<- data.frame(ma_id=random.ma.id,
									     method=rep(random.methods.list, 
									 	            length=random.nrows),
									     tau2_method=rep(tau2.methods, 
									 	                 length=random.nrows),
									     mean=rep(NA, length=random.nrows),
									     median=rep(NA, length=random.nrows),
									     stdev=rep(NA, length=random.nrows),
									     ci_lower=rep(NA, length=random.nrows),
									     ci_upper=rep(NA, length=random.nrows),
									     tau2=rep(NA, length=random.nrows),
									     ma_type=rep("big", length=random.nrows),
									     data_type=rep("bin", length=random.nrows))

		#####
		# run gmeta
		#####

		for(each.ma in ma.id) {

			# grab individual one m-a at a time
			ma.raw <- ma.set[ma.set[, 'ma_id'] == each.ma, ]

			# grab 2x2 table counts from meta-analyses (ma.raw)
			ma.2by2 <- cbind(ma.raw[, "study_events1"], ma.raw[, "study_total1"], 
							 ma.raw[, "study_events2"], ma.raw[, "study_total2"])

			#####
			# compute MH odds ratios for each study in m-a using gmeta
			#####

			# apply continuity correction if study contains 0 count
			if(sum(ma.2by2 == 0) > 0) {
				row.sub <- apply(ma.2by2, 1, function(row) 0 %in% row)
				ma.2by2[row.sub, ] <- ma.2by2[row.sub, ] + 0.5
				print(paste("For ma_id", each.ma, 
					        "Haldane continuity correction was applied", sep=" "))
			}

			gmeta.MHOR <- gmeta(gmi=ma.2by2, gmi.type='2x2', method='MH')

			# grab individual study logOR means, stdevs
			study.means <- log(gmeta.MHOR$individual.means)
			study.stdevs <- gmeta.MHOR$individual.stddevs/gmeta.MHOR$individual.means
			summ.stats <- cbind(study.means, study.stdevs)

			# run gmeta fixed-effects
			for(each.fmethod in fixed.methods) {
				gmeta.fixed <- gmeta(gmi=summ.stats, 
									 gmi.type='pivot', 
									 method=each.fmethod)
				fixed.row <- fixed.effects.bb[, 'ma_id'] == each.ma & 
				             fixed.effects.bb[, 'method'] == each.fmethod
				fixed.effects.bb[fixed.row, 'mean'] <<- gmeta.fixed$combined.mean
				fixed.effects.bb[fixed.row, 'median'] <<- gmeta.fixed$combined.median
				fixed.effects.bb[fixed.row, 'stdev'] <<- gmeta.fixed$combined.sd 
				fixed.effects.bb[fixed.row, 'ci_lower'] <<- gmeta.fixed$combined.ci[1]
				fixed.effects.bb[fixed.row, 'ci_upper'] <<- gmeta.fixed$combined.ci[2]
			}

			# run gmeta random-effects

			for(each.rmethod in random.methods) {
				for(each.tmethod in tau2.methods) {
					gmeta.random <- gmeta(gmi=summ.stats, gmi.type='pivot', 
						                  method=each.rmethod, tau2=each.tmethod)
					random.row <- random.effects.bb[, 'ma_id'] == each.ma &
								  random.effects.bb[, 'method'] == each.rmethod &
								  random.effects.bb[, 'tau2_method'] == each.tmethod
					random.effects.bb[random.row, 'mean'] <<- gmeta.random$combined.mean
					random.effects.bb[random.row, 'median'] <<- gmeta.random$combined.median
					random.effects.bb[random.row, 'stdev'] <<- gmeta.random$combined.sd
					random.effects.bb[random.row, 'ci_lower'] <<- gmeta.random$combined.ci[1]
					random.effects.bb[random.row, 'ci_upper'] <<- gmeta.random$combined.ci[2]
					random.effects.bb[random.row, 'tau2'] <<- gmeta.random$tau2 
				}
			}
		}
	} else if(ma.type == "many" & data.type == "binary") {

		########## MANY, BINARY ##########

		#####
		# create skeleton dataframes to be filled in later
		#####

		# fixed-effects

		fixed.methods <- c("fixed-mle", "fixed-robust2")

		fixed.nrows <- length(ma.id)*length(fixed.methods)

		fixed.ma.id <- sort(rep(ma.id, length=fixed.nrows))

		fixed.effects.mb <<- data.frame(ma_id=fixed.ma.id, 
									    method=rep(fixed.methods, length=fixed.nrows),
									    mean=rep(NA, length=fixed.nrows),
									    median=rep(NA, length=fixed.nrows),
									    stdev=rep(NA, length=fixed.nrows),
									    ci_lower=rep(NA, length=fixed.nrows),
									    ci_upper=rep(NA, length=fixed.nrows),
									    ma_type=rep("many", length=fixed.nrows),
									    data_type=rep("bin", length=fixed.nrows))

		# random-effects
		random.methods <- c("random-mm", "random-reml", "random-tau2", "random-robust2")
		tau2.methods <- c("DL", "HS", "SJ", "HE", "ML", "REML", "EB")

		len.tau2.methods <- length(tau2.methods)
		random.nrows <- length(ma.id)*length(random.methods)*len.tau2.methods

		random.methods.list <- character()

		for(each.rm in random.methods) {
			random.methods.list <- c(random.methods.list, 
				                     rep(each.rm, length=len.tau2.methods))
		}

		random.ma.id <- sort(rep(ma.id, length=random.nrows))

		random.effects.mb <<- data.frame(ma_id=random.ma.id,
									     method=rep(random.methods.list, 
									 	            length=random.nrows),
									     tau2_method=rep(tau2.methods, 
									 	                 length=random.nrows),
									     mean=rep(NA, length=random.nrows),
									     median=rep(NA, length=random.nrows),
									     stdev=rep(NA, length=random.nrows),
									     ci_lower=rep(NA, length=random.nrows),
									     ci_upper=rep(NA, length=random.nrows),
									     tau2=rep(NA, length=random.nrows),
									     ma_type=rep("many", length=random.nrows),
									     data_type=rep("bin", length=random.nrows))

		#####
		# run gmeta
		#####

		for(each.ma in ma.id) {

			# grab individual one m-a at a time
			ma.raw <- ma.set[ma.set[, 'ma_id'] == each.ma, ]

			# grab 2x2 table counts from meta-analyses (ma.raw)
			ma.2by2 <<- cbind(ma.raw[, "study_events1"], ma.raw[, "study_total1"], 
							 ma.raw[, "study_events2"], ma.raw[, "study_total2"])

			#####
			# compute MH odds ratios for each study in m-a using gmeta
			#####

			# apply continuity correction if study contains 0 count
			if(sum(ma.2by2 == 0) > 0) {
				row.sub <- apply(ma.2by2, 1, function(row) 0 %in% row)
				ma.2by2[row.sub, ] <<- ma.2by2[row.sub, ] + 0.5
				print(paste("For ma_id", each.ma, 
					        "Haldane continuity correction was applied", sep=" "))
			}

			gmeta.MHOR <- gmeta(gmi=ma.2by2, gmi.type='2x2', method='MH')

			# grab individual study logOR means, stdevs
			study.means <- log(gmeta.MHOR$individual.means)
			study.stdevs <- gmeta.MHOR$individual.stddevs/gmeta.MHOR$individual.means
			summ.stats <- cbind(study.means, study.stdevs)

			# run gmeta fixed-effects
			for(each.fmethod in fixed.methods) {
				gmeta.fixed <- gmeta(gmi=summ.stats, 
									 gmi.type='pivot', 
									 method=each.fmethod)
				fixed.row <- fixed.effects.mb[, 'ma_id'] == each.ma & 
				             fixed.effects.mb[, 'method'] == each.fmethod
				fixed.effects.mb[fixed.row, 'mean'] <<- gmeta.fixed$combined.mean
				fixed.effects.mb[fixed.row, 'median'] <<- gmeta.fixed$combined.median
				fixed.effects.mb[fixed.row, 'stdev'] <<- gmeta.fixed$combined.sd 
				fixed.effects.mb[fixed.row, 'ci_lower'] <<- gmeta.fixed$combined.ci[1]
				fixed.effects.mb[fixed.row, 'ci_upper'] <<- gmeta.fixed$combined.ci[2]
			}

			# run gmeta random-effects

			for(each.rmethod in random.methods) {
				for(each.tmethod in tau2.methods) {
					gmeta.random <- gmeta(gmi=summ.stats, gmi.type='pivot', 
						                  method=each.rmethod, tau2=each.tmethod)
					random.row <- random.effects.mb[, 'ma_id'] == each.ma &
								  random.effects.mb[, 'method'] == each.rmethod &
								  random.effects.mb[, 'tau2_method'] == each.tmethod
					random.effects.mb[random.row, 'mean'] <<- gmeta.random$combined.mean
					random.effects.mb[random.row, 'median'] <<- gmeta.random$combined.median
					random.effects.mb[random.row, 'stdev'] <<- gmeta.random$combined.sd
					random.effects.mb[random.row, 'ci_lower'] <<- gmeta.random$combined.ci[1]
					random.effects.mb[random.row, 'ci_upper'] <<- gmeta.random$combined.ci[2]
					random.effects.mb[random.row, 'tau2'] <<- gmeta.random$tau2 
				}
			}
		}
	} else if(ma.type == "big" & data.type == "continuous") {

		########## BIG, CONTINUOUS ##########

		#####
		# create skeleton dataframes to be filled in later
		#####

		# fixed-effects

		fixed.methods <- c("fixed-mle", "fixed-robust1")

		fixed.nrows <- length(ma.id)*length(fixed.methods)

		fixed.ma.id <- sort(rep(ma.id, length=fixed.nrows))

		fixed.effects.bc <<- data.frame(ma_id=fixed.ma.id, 
									    method=rep(fixed.methods, length=fixed.nrows),
									    mean=rep(NA, length=fixed.nrows),
									    median=rep(NA, length=fixed.nrows),
									    stdev=rep(NA, length=fixed.nrows),
									    ci_lower=rep(NA, length=fixed.nrows),
									    ci_upper=rep(NA, length=fixed.nrows),
									    ma_type=rep("big", length=fixed.nrows),
									    data_type=rep("cont", length=fixed.nrows))

		# random-effects
		random.methods <- c("random-mm", "random-reml", "random-tau2", "random-robust1")
		tau2.methods <- c("DL", "HS", "SJ", "HE", "ML", "REML", "EB")

		len.tau2.methods <- length(tau2.methods)
		random.nrows <- length(ma.id)*length(random.methods)*len.tau2.methods

		random.methods.list <- character()

		for(each.rm in random.methods) {
			random.methods.list <- c(random.methods.list, 
				                     rep(each.rm, length=len.tau2.methods))
		}

		random.ma.id <- sort(rep(ma.id, length=random.nrows))

		random.effects.bc <<- data.frame(ma_id=random.ma.id,
									     method=rep(random.methods.list, 
									 	            length=random.nrows),
									     tau2_method=rep(tau2.methods, 
									 	                 length=random.nrows),
									     mean=rep(NA, length=random.nrows),
									     median=rep(NA, length=random.nrows),
									     stdev=rep(NA, length=random.nrows),
									     ci_lower=rep(NA, length=random.nrows),
									     ci_upper=rep(NA, length=random.nrows),
									     tau2=rep(NA, length=random.nrows),
									     ma_type=rep("big", length=random.nrows),
									     data_type=rep("cont", length=random.nrows))

		#####
		# run gmeta
		#####

		for(each.ma in ma.id) {

			# grab individual one m-a at a time
			ma.raw <- ma.set[ma.set[, 'ma_id'] == each.ma, ]

			summ.stats <- c(ma.raw[, 'study_means1'] - ma.raw[, 'study_means2'],
				            sqrt(ma.raw[, 'study_sd1'] + ma.raw[, 'study_sd2']))

			# run gmeta fixed-effects
			for(each.fmethod in fixed.methods) {
				gmeta.fixed <- gmeta(gmi=summ.stats, 
									 gmi.type='pivot', 
									 method=each.fmethod)
				fixed.row <- fixed.effects.bc[, 'ma_id'] == each.ma & 
				             fixed.effects.bc[, 'method'] == each.fmethod
				fixed.effects.bc[fixed.row, 'mean'] <<- gmeta.fixed$combined.mean
				fixed.effects.bc[fixed.row, 'median'] <<- gmeta.fixed$combined.median
				fixed.effects.bc[fixed.row, 'stdev'] <<- gmeta.fixed$combined.sd 
				fixed.effects.bc[fixed.row, 'ci_lower'] <<- gmeta.fixed$combined.ci[1]
				fixed.effects.bc[fixed.row, 'ci_upper'] <<- gmeta.fixed$combined.ci[2]
			}

			# run gmeta random-effects

			for(each.rmethod in random.methods) {
				for(each.tmethod in tau2.methods) {
					gmeta.random <- gmeta(gmi=summ.stats, gmi.type='pivot', 
						                  method=each.rmethod, tau2=each.tmethod)
					random.row <- random.effects.bc[, 'ma_id'] == each.ma &
								  random.effects.bc[, 'method'] == each.rmethod &
								  random.effects.bc[, 'tau2_method'] == each.tmethod
					random.effects.bc[random.row, 'mean'] <<- gmeta.random$combined.mean
					random.effects.bc[random.row, 'median'] <<- gmeta.random$combined.median
					random.effects.bc[random.row, 'stdev'] <<- gmeta.random$combined.sd
					random.effects.bc[random.row, 'ci_lower'] <<- gmeta.random$combined.ci[1]
					random.effects.bc[random.row, 'ci_upper'] <<- gmeta.random$combined.ci[2]
					random.effects.bc[random.row, 'tau2'] <<- gmeta.random$tau2 
				}
			}
		}
	} else if(ma.type == "many" & data.type == "continuous") {

		########## MANY, CONTINUOUS ##########

		#####
		# create skeleton dataframes to be filled in later
		#####

		# fixed-effects

		fixed.methods <- c("fixed-mle", "fixed-robust2")

		fixed.nrows <- length(ma.id)*length(fixed.methods)

		fixed.ma.id <- sort(rep(ma.id, length=fixed.nrows))

		fixed.effects.mc <<- data.frame(ma_id=fixed.ma.id, 
									    method=rep(fixed.methods, length=fixed.nrows),
									    mean=rep(NA, length=fixed.nrows),
									    median=rep(NA, length=fixed.nrows),
									    stdev=rep(NA, length=fixed.nrows),
									    ci_lower=rep(NA, length=fixed.nrows),
									    ci_upper=rep(NA, length=fixed.nrows),
									    ma_type=rep("many", length=fixed.nrows),
									    data_type=rep("cont", length=fixed.nrows))

		# random-effects
		random.methods <- c("random-mm", "random-reml", "random-tau2", "random-robust2")
		tau2.methods <- c("DL", "HS", "SJ", "HE", "ML", "REML", "EB")

		len.tau2.methods <- length(tau2.methods)
		random.nrows <- length(ma.id)*length(random.methods)*len.tau2.methods

		random.methods.list <- character()

		for(each.rm in random.methods) {
			random.methods.list <- c(random.methods.list, 
				                     rep(each.rm, length=len.tau2.methods))
		}

		random.ma.id <- sort(rep(ma.id, length=random.nrows))

		random.effects.mc <<- data.frame(ma_id=random.ma.id,
									     method=rep(random.methods.list, 
									 	            length=random.nrows),
									     tau2_method=rep(tau2.methods, 
									 	                 length=random.nrows),
									     mean=rep(NA, length=random.nrows),
									     median=rep(NA, length=random.nrows),
									     stdev=rep(NA, length=random.nrows),
									     ci_lower=rep(NA, length=random.nrows),
									     ci_upper=rep(NA, length=random.nrows),
									     tau2=rep(NA, length=random.nrows),
									     ma_type=rep("many", length=random.nrows),
									     data_type=rep("cont", length=random.nrows))

		#####
		# run gmeta
		#####

		for(each.ma in ma.id) {

			# grab individual one m-a at a time
			ma.raw <- ma.set[ma.set[, 'ma_id'] == each.ma, ]

			summ.stats <- c(ma.raw[, 'study_means1'] - ma.raw[, 'study_means2'],
				            sqrt(ma.raw[, 'study_sd1'] + ma.raw[, 'study_sd2']))

			# run gmeta fixed-effects
			for(each.fmethod in fixed.methods) {
				gmeta.fixed <- gmeta(gmi=summ.stats, 
									 gmi.type='pivot', 
									 method=each.fmethod)
				fixed.row <- fixed.effects.mc[, 'ma_id'] == each.ma & 
				             fixed.effects.mc[, 'method'] == each.fmethod
				fixed.effects.mc[fixed.row, 'mean'] <<- gmeta.fixed$combined.mean
				fixed.effects.mc[fixed.row, 'median'] <<- gmeta.fixed$combined.median
				fixed.effects.mc[fixed.row, 'stdev'] <<- gmeta.fixed$combined.sd 
				fixed.effects.mc[fixed.row, 'ci_lower'] <<- gmeta.fixed$combined.ci[1]
				fixed.effects.mc[fixed.row, 'ci_upper'] <<- gmeta.fixed$combined.ci[2]
			}

			# run gmeta random-effects

			for(each.rmethod in random.methods) {
				for(each.tmethod in tau2.methods) {
					gmeta.random <- gmeta(gmi=summ.stats, gmi.type='pivot', 
						                  method=each.rmethod, tau2=each.tmethod)
					random.row <- random.effects.mc[, 'ma_id'] == each.ma &
								  random.effects.mc[, 'method'] == each.rmethod &
								  random.effects.mc[, 'tau2_method'] == each.tmethod
					random.effects.mc[random.row, 'mean'] <<- gmeta.random$combined.mean
					random.effects.mc[random.row, 'median'] <<- gmeta.random$combined.median
					random.effects.mc[random.row, 'stdev'] <<- gmeta.random$combined.sd
					random.effects.mc[random.row, 'ci_lower'] <<- gmeta.random$combined.ci[1]
					random.effects.mc[random.row, 'ci_upper'] <<- gmeta.random$combined.ci[2]
					random.effects.mc[random.row, 'tau2'] <<- gmeta.random$tau2 
				}
			}
		}
	}
}

# set working directory (contains files to work with)
setwd("C:/Users/Ian/Desktop/meta/sandbox/")

# load libraries
library(gmeta)

# read in meta-analyses from database
big.studies <- read.csv("big_studies.csv") # large sample size

# grab meta-analysis IDs
ma.id <- sort(unique(big.studies[, "ma_id"]))

# run function
run.gmeta(big.studies, "big", "binary")

many.studies <- read.csv("many_studies.csv") # large # of studies
ma.id <- sort(unique(many.studies[, "ma_id"]))
run.gmeta(many.studies, "many", "binary") 



# print runtime
runtime <- proc.time() - runtime
print("Runtime:")
print(runtime)