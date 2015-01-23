# runtime
runtime <- proc.time()

# set working directory (contains files to work with)
setwd("C:/Users/Ian/Desktop/meta/sandbox/")

# read in meta-analyses from database
large.studies.bin <- read.csv("large_studies_bin.csv") # large sample size
many.studies.bin <- read.csv("many_studies_bin.csv") # large # of studies
large.studies.cont <- read.csv("large_studies_cont.csv")
many.studies.cont <- read.csv("many_studies_cont.csv")

ma.id.lbin <- sort(unique(large.studies.bin[, "ma_id"]))
ma.id.mbin <- sort(unique(many.studies.bin[, "ma_id"]))
ma.id.lcont <- sort(unique(large.studies.cont[, "ma_id"]))
ma.id.mc.BAD <- c(170)
ma.id.mcont <- sort(unique(many.studies.cont[, "ma_id"]))
ma.id.mcont <- ma.id.mcont[!(ma.id.mcont %in% ma.id.mc.BAD)]

# load libraries
library(gmeta)

run.gmeta <- function(ma.set, ma.type, data.type) {

	# large studies, binary data
	if(ma.type == "large" & data.type == "bin") {
		fixed.methods <- c("fixed-mle", "fixed-robust1")
		fixed.nrows <- length(ma.id.lbin)*length(fixed.methods)
		fixed.ma.id <- sort(rep(ma.id.lbin, length=fixed.nrows))
		fixed.effects.lbin <<- data.frame(ma_id=fixed.ma.id,
			                              method=rep(fixed.methods, length=fixed.nrows),
			                              mean=rep(NA, length=fixed.nrows),
			                              median=rep(NA, length=fixed.nrows),
			                              stdev=rep(NA, length=fixed.nrows),
			                              ci_lower=rep(NA, length=fixed.nrows),
			                              ci_upper=rep(NA, length=fixed.nrows),
			                              ma_type=rep(ma.type, length=fixed.nrows),
			                              data_type=rep(data.type, length=fixed.nrows))

		# random-effects
		random.methods <- c("random-mm", "random-reml", "random-tau2", "random-robust1")
		tau2.methods <- c("DL", "HS", "SJ", "HE", "ML", "REML", "EB")
		len.tau2.methods <- length(tau2.methods)
		random.nrows <- length(ma.id.lbin)*length(random.methods)*len.tau2.methods

		random.methods.list <- character()
		for(each.rm in random.methods) {
			random.methods.list <- c(random.methods.list, 
				                     rep(each.rm, length=len.tau2.methods))
		}

		random.ma.id <- sort(rep(ma.id.lbin, length=random.nrows))
		random.effects.lbin <<- data.frame(ma_id=random.ma.id,
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
										   ma_type=rep(ma.type, length=fixed.nrows),
				                           data_type=rep(data.type, length=fixed.nrows))

		#####
		# run gmeta
		#####

		for(each.ma in ma.id.lbin) {

			# grab individual one m-a at a time
			ma.raw <- large.studies.bin[large.studies.bin[, 'ma_id'] == each.ma, ]

			# grab 2x2 table counts from meta-analyses (ma.raw)
			ma.2by2 <- cbind(ma.raw[, "study_events1"], ma.raw[, "study_total1"], 
							 ma.raw[, "study_events2"], ma.raw[, "study_total2"])
			colnames(ma.2by2) <- c("study_events1", "study_total1", "study_events2", "study_total2")

			#####
			# compute MH odds ratios for each study in m-a using gmeta
			#####

			# apply continuity correction if study contains 0 count
			for(ma.row.num in 1:nrow(ma.2by2)) {
				tmp.ma.ckr <- cbind(ma.2by2[ma.row.num, "study_events1"] - ma.2by2[ma.row.num, "study_total1"],
				                    ma.2by2[ma.row.num, "study_events2"] - ma.2by2[ma.row.num, "study_total2"])
				ma.check1 <- 0 %in% tmp.ma.ckr
				ma.check2 <- 0 %in% ma.2by2[ma.row.num, ]
				if(ma.check1 | ma.check2) {
					ma.2by2[ma.row.num, c("study_events1", "study_events2")] <- ma.2by2[ma.row.num, c("study_events1", "study_events2")] + 0.5
					ma.2by2[ma.row.num, c("study_total1", "study_total2")] <- ma.2by2[ma.row.num, c("study_total1", "study_total2")] + 1
					print(paste("For ma_id", each.ma, 
					            "Haldane continuity correction was applied to study", 
					            ma.row.num, sep=" "))
				}
			}     							   

			gmeta.MHOR <- gmeta(gmi=ma.2by2, gmi.type='2x2', method='MH', gmo.xgrid=seq(-20, 20, by=0.001))

			# grab individual study logOR means, stdevs
			study.means <- log(gmeta.MHOR$individual.means)
			study.stdevs <- gmeta.MHOR$individual.stddevs/gmeta.MHOR$individual.means
			summ.stats <- cbind(study.means, study.stdevs)

			# run gmeta fixed-effects
			for(each.fmethod in fixed.methods) {
				gmeta.fixed <- gmeta(gmi=summ.stats, 
									 gmi.type='pivot', 
									 method=each.fmethod,
									 gmo.xgrid=seq(min(summ.stats[, 1]) - 4*max(summ.stats[, 2]),
									 	           max(summ.stats[, 1]) + 4*max(summ.stats[, 2]),
									 	           by=0.001))
				fixed.row <- fixed.effects.lbin[, 'ma_id'] == each.ma & 
				             fixed.effects.lbin[, 'method'] == each.fmethod
				fixed.effects.lbin[fixed.row, 'mean'] <<- gmeta.fixed$combined.mean
				fixed.effects.lbin[fixed.row, 'median'] <<- gmeta.fixed$combined.median
				fixed.effects.lbin[fixed.row, 'stdev'] <<- gmeta.fixed$combined.sd 
				fixed.effects.lbin[fixed.row, 'ci_lower'] <<- gmeta.fixed$combined.ci[1]
				fixed.effects.lbin[fixed.row, 'ci_upper'] <<- gmeta.fixed$combined.ci[2]
			}

			# run gmeta random-effects

			for(each.rmethod in random.methods) {
				for(each.tmethod in tau2.methods) {
					gmeta.random <- gmeta(gmi=summ.stats, gmi.type='pivot', 
						                  method=each.rmethod, tau2=each.tmethod,
									 	  gmo.xgrid=seq(min(summ.stats[, 1]) - 4*max(summ.stats[, 2]),
									 	                max(summ.stats[, 1]) + 4*max(summ.stats[, 2]),
									 	                by=0.001))
					random.row <- random.effects.lbin[, 'ma_id'] == each.ma &
								  random.effects.lbin[, 'method'] == each.rmethod &
								  random.effects.lbin[, 'tau2_method'] == each.tmethod
					random.effects.lbin[random.row, 'mean'] <<- gmeta.random$combined.mean
					random.effects.lbin[random.row, 'median'] <<- gmeta.random$combined.median
					random.effects.lbin[random.row, 'stdev'] <<- gmeta.random$combined.sd
					random.effects.lbin[random.row, 'ci_lower'] <<- gmeta.random$combined.ci[1]
					random.effects.lbin[random.row, 'ci_upper'] <<- gmeta.random$combined.ci[2]
					random.effects.lbin[random.row, 'tau2'] <<- gmeta.random$tau2 
				}
			}
		}
	}

	# many studies, binary data
	if(ma.type == "many" & data.type == "bin") {
		fixed.methods <- c("fixed-mle", "fixed-robust2")
		fixed.nrows <- length(ma.id.mbin)*length(fixed.methods)
		fixed.ma.id <- sort(rep(ma.id.mbin, length=fixed.nrows))
		fixed.effects.mbin <<- data.frame(ma_id=fixed.ma.id,
			                              method=rep(fixed.methods, length=fixed.nrows),
			                              mean=rep(NA, length=fixed.nrows),
			                              median=rep(NA, length=fixed.nrows),
			                              stdev=rep(NA, length=fixed.nrows),
			                              ci_lower=rep(NA, length=fixed.nrows),
			                              ci_upper=rep(NA, length=fixed.nrows),
			                              ma_type=rep(ma.type, length=fixed.nrows),
			                              data_type=rep(data.type, length=fixed.nrows))

		# random-effects
		random.methods <- c("random-mm", "random-reml", "random-tau2", "random-robust2")
		tau2.methods <- c("DL", "HS", "SJ", "HE", "ML", "REML", "EB")
		len.tau2.methods <- length(tau2.methods)
		random.nrows <- length(ma.id.mbin)*length(random.methods)*len.tau2.methods

		random.methods.list <- character()
		for(each.rm in random.methods) {
			random.methods.list <- c(random.methods.list, 
				                     rep(each.rm, length=len.tau2.methods))
		}

		random.ma.id <- sort(rep(ma.id.mbin, length=random.nrows))
		random.effects.mbin <<- data.frame(ma_id=random.ma.id,
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
										   ma_type=rep(ma.type, length=fixed.nrows),
				                           data_type=rep(data.type, length=fixed.nrows))

		#####
		# run gmeta
		#####

		for(each.ma in ma.id.mbin) {

			# grab individual one m-a at a time
			ma.raw <- many.studies.bin[many.studies.bin[, 'ma_id'] == each.ma, ]

			# grab 2x2 table counts from meta-analyses (ma.raw)
			ma.2by2 <- cbind(ma.raw[, "study_events1"], ma.raw[, "study_total1"], 
							 ma.raw[, "study_events2"], ma.raw[, "study_total2"])
			colnames(ma.2by2) <- c("study_events1", "study_total1", "study_events2", "study_total2")

			#####
			# compute MH odds ratios for each study in m-a using gmeta
			#####

			# apply continuity correction if study contains 0 count
			for(ma.row.num in 1:nrow(ma.2by2)) {
				tmp.ma.ckr <- cbind(ma.2by2[ma.row.num, "study_events1"] - ma.2by2[ma.row.num, "study_total1"],
				                    ma.2by2[ma.row.num, "study_events2"] - ma.2by2[ma.row.num, "study_total2"])
				ma.check1 <- 0 %in% tmp.ma.ckr
				ma.check2 <- 0 %in% ma.2by2[ma.row.num, ]
				if(ma.check1 | ma.check2) {
					ma.2by2[ma.row.num, c("study_events1", "study_events2")] <- ma.2by2[ma.row.num, c("study_events1", "study_events2")] + 0.5
					ma.2by2[ma.row.num, c("study_total1", "study_total2")] <- ma.2by2[ma.row.num, c("study_total1", "study_total2")] + 1
					print(paste("For ma_id", each.ma, 
					            "Haldane continuity correction was applied to study", 
					            ma.row.num, sep=" "))
				}
			}     

			gmeta.MHOR <- gmeta(gmi=ma.2by2, gmi.type='2x2', method='MH', gmo.xgrid=seq(-20, 20, by=0.001))

			# grab individual study logOR means, stdevs
			study.means <- log(gmeta.MHOR$individual.means)
			study.stdevs <- gmeta.MHOR$individual.stddevs/gmeta.MHOR$individual.means
			summ.stats <- cbind(study.means, study.stdevs)

			# run gmeta fixed-effects
			for(each.fmethod in fixed.methods) {
				gmeta.fixed <- gmeta(gmi=summ.stats, 
									 gmi.type='pivot', 
									 method=each.fmethod,
									 gmo.xgrid=seq(min(summ.stats[, 1]) - 4*max(summ.stats[, 2]),
									 	           max(summ.stats[, 1]) + 4*max(summ.stats[, 2]),
									 	           by=0.001))
				fixed.row <- fixed.effects.mbin[, 'ma_id'] == each.ma & 
				             fixed.effects.mbin[, 'method'] == each.fmethod
				fixed.effects.mbin[fixed.row, 'mean'] <<- gmeta.fixed$combined.mean
				fixed.effects.mbin[fixed.row, 'median'] <<- gmeta.fixed$combined.median
				fixed.effects.mbin[fixed.row, 'stdev'] <<- gmeta.fixed$combined.sd 
				fixed.effects.mbin[fixed.row, 'ci_lower'] <<- gmeta.fixed$combined.ci[1]
				fixed.effects.mbin[fixed.row, 'ci_upper'] <<- gmeta.fixed$combined.ci[2]
			}

			# run gmeta random-effects

			for(each.rmethod in random.methods) {
				for(each.tmethod in tau2.methods) {
					gmeta.random <- gmeta(gmi=summ.stats, gmi.type='pivot', 
						                  method=each.rmethod, tau2=each.tmethod,
									 	  gmo.xgrid=seq(min(summ.stats[, 1]) - 4*max(summ.stats[, 2]),
									 	                max(summ.stats[, 1]) + 4*max(summ.stats[, 2]),
									 	                by=0.001))
					random.row <- random.effects.mbin[, 'ma_id'] == each.ma &
								  random.effects.mbin[, 'method'] == each.rmethod &
								  random.effects.mbin[, 'tau2_method'] == each.tmethod
					random.effects.mbin[random.row, 'mean'] <<- gmeta.random$combined.mean
					random.effects.mbin[random.row, 'median'] <<- gmeta.random$combined.median
					random.effects.mbin[random.row, 'stdev'] <<- gmeta.random$combined.sd
					random.effects.mbin[random.row, 'ci_lower'] <<- gmeta.random$combined.ci[1]
					random.effects.mbin[random.row, 'ci_upper'] <<- gmeta.random$combined.ci[2]
					random.effects.mbin[random.row, 'tau2'] <<- gmeta.random$tau2 
				}
			}
		}
	}

	# large studies, continuous data
	if(ma.type == "large" & data.type == "cont") {
		fixed.methods <- c("fixed-mle", "fixed-robust1")
		fixed.nrows <- length(ma.id.lcont)*length(fixed.methods)
		fixed.ma.id <- sort(rep(ma.id.lcont, length=fixed.nrows))
		fixed.effects.lcont <<- data.frame(ma_id=fixed.ma.id,
			                               method=rep(fixed.methods, length=fixed.nrows),
			                               mean=rep(NA, length=fixed.nrows),
			                               median=rep(NA, length=fixed.nrows),
			                               stdev=rep(NA, length=fixed.nrows),
			                               ci_lower=rep(NA, length=fixed.nrows),
			                               ci_upper=rep(NA, length=fixed.nrows),
			                               ma_type=rep(ma.type, length=fixed.nrows),
			                               data_type=rep(data.type, length=fixed.nrows))

		# random-effects
		random.methods <- c("random-mm", "random-reml", "random-tau2", "random-robust1")
		tau2.methods <- c("DL", "HS", "SJ", "HE", "ML", "REML", "EB")
		len.tau2.methods <- length(tau2.methods)
		random.nrows <- length(ma.id.lcont)*length(random.methods)*len.tau2.methods

		random.methods.list <- character()
		for(each.rm in random.methods) {
			random.methods.list <- c(random.methods.list, 
				                     rep(each.rm, length=len.tau2.methods))
		}

		random.ma.id <- sort(rep(ma.id.lcont, length=random.nrows))
		random.effects.lcont <<- data.frame(ma_id=random.ma.id,
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
										    ma_type=rep(ma.type, length=fixed.nrows),
				                            data_type=rep(data.type, length=fixed.nrows))

		#####
		# run gmeta
		#####

		for(each.ma in ma.id.lcont) {

			# grab individual one m-a at a time
			ma.raw <- large.studies.cont[large.studies.cont[, 'ma_id'] == each.ma, ]

			# grab summary statistics
			summ.stats <<- cbind(ma.raw[, "study_mean1"] - ma.raw[, "study_mean2"], 
							     sqrt((ma.raw[, "study_sd1"])^2 + (ma.raw[, "study_sd2"])^2))


			# run gmeta fixed-effects
			for(each.fmethod in fixed.methods) {
				gmeta.fixed <- gmeta(gmi=summ.stats, 
									 gmi.type='pivot', 
									 method=each.fmethod,
									 gmo.xgrid=seq(min(summ.stats[, 1]) - 4*max(summ.stats[, 2]),
									 	           max(summ.stats[, 1]) + 4*max(summ.stats[, 2]),
									 	           by=0.001))
				fixed.row <- fixed.effects.lcont[, 'ma_id'] == each.ma & 
				             fixed.effects.lcont[, 'method'] == each.fmethod
				fixed.effects.lcont[fixed.row, 'mean'] <<- gmeta.fixed$combined.mean
				fixed.effects.lcont[fixed.row, 'median'] <<- gmeta.fixed$combined.median
				fixed.effects.lcont[fixed.row, 'stdev'] <<- gmeta.fixed$combined.sd 
				fixed.effects.lcont[fixed.row, 'ci_lower'] <<- gmeta.fixed$combined.ci[1]
				fixed.effects.lcont[fixed.row, 'ci_upper'] <<- gmeta.fixed$combined.ci[2]
			}

			# run gmeta random-effects

			for(each.rmethod in random.methods) {
				for(each.tmethod in tau2.methods) {
					gmeta.random <- gmeta(gmi=summ.stats, gmi.type='pivot', 
						                  method=each.rmethod, tau2=each.tmethod,
									 	  gmo.xgrid=seq(min(summ.stats[, 1]) - 4*max(summ.stats[, 2]),
									 	           	    max(summ.stats[, 1]) + 4*max(summ.stats[, 2]),
									 	           	    by=0.001))
					random.row <- random.effects.lcont[, 'ma_id'] == each.ma &
								  random.effects.lcont[, 'method'] == each.rmethod &
								  random.effects.lcont[, 'tau2_method'] == each.tmethod
					random.effects.lcont[random.row, 'mean'] <<- gmeta.random$combined.mean
					random.effects.lcont[random.row, 'median'] <<- gmeta.random$combined.median
					random.effects.lcont[random.row, 'stdev'] <<- gmeta.random$combined.sd
					random.effects.lcont[random.row, 'ci_lower'] <<- gmeta.random$combined.ci[1]
					random.effects.lcont[random.row, 'ci_upper'] <<- gmeta.random$combined.ci[2]
					random.effects.lcont[random.row, 'tau2'] <<- gmeta.random$tau2 
				}
			}
		}
	}

	# many studies, continuous data
	if(ma.type == "many" & data.type == "cont") {
		fixed.methods <- c("fixed-mle", "fixed-robust2")
		fixed.nrows <- length(ma.id.mcont)*length(fixed.methods)
		fixed.ma.id <- sort(rep(ma.id.mcont, length=fixed.nrows))
		fixed.effects.mcont <<- data.frame(ma_id=fixed.ma.id,
			                               method=rep(fixed.methods, length=fixed.nrows),
			                               mean=rep(NA, length=fixed.nrows),
			                               median=rep(NA, length=fixed.nrows),
			                               stdev=rep(NA, length=fixed.nrows),
			                               ci_lower=rep(NA, length=fixed.nrows),
			                               ci_upper=rep(NA, length=fixed.nrows),
			                               ma_type=rep(ma.type, length=fixed.nrows),
			                               data_type=rep(data.type, length=fixed.nrows))

		# random-effects
		random.methods <- c("random-mm", "random-reml", "random-tau2", "random-robust2")
		tau2.methods <- c("DL", "HS", "SJ", "HE", "ML", "REML", "EB")
		len.tau2.methods <- length(tau2.methods)
		random.nrows <- length(ma.id.mcont)*length(random.methods)*len.tau2.methods

		random.methods.list <- character()
		for(each.rm in random.methods) {
			random.methods.list <- c(random.methods.list, 
				                     rep(each.rm, length=len.tau2.methods))
		}

		random.ma.id <- sort(rep(ma.id.mcont, length=random.nrows))
		random.effects.mcont <<- data.frame(ma_id=random.ma.id,
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
										    ma_type=rep(ma.type, length=fixed.nrows),
				                            data_type=rep(data.type, length=fixed.nrows))

		#####
		# run gmeta
		#####

		for(each.ma in ma.id.mcont) {

			# grab individual one m-a at a time
			ma.raw <<- many.studies.cont[many.studies.cont[, 'ma_id'] == each.ma, ]

			# grab summary statistics
			summ.stats <<- cbind(ma.raw[, "study_mean1"] - ma.raw[, "study_mean2"], 
							     sqrt((ma.raw[, "study_sd1"])^2 + (ma.raw[, "study_sd2"])^2))

			# run gmeta fixed-effects
			for(each.fmethod in fixed.methods) {
				gmeta.fixed <<- gmeta(gmi=summ.stats, 
									  gmi.type='pivot', 
									  method=each.fmethod,
									  gmo.xgrid=seq(min(summ.stats[, 1]) - 4*max(summ.stats[, 2]),
									 	            max(summ.stats[, 1]) + 4*max(summ.stats[, 2]),
									 	            by=0.001))
				fixed.row <- fixed.effects.mcont[, 'ma_id'] == each.ma & 
				             fixed.effects.mcont[, 'method'] == each.fmethod
				fixed.effects.mcont[fixed.row, 'mean'] <<- gmeta.fixed$combined.mean
				fixed.effects.mcont[fixed.row, 'median'] <<- gmeta.fixed$combined.median
				fixed.effects.mcont[fixed.row, 'stdev'] <<- gmeta.fixed$combined.sd 
				fixed.effects.mcont[fixed.row, 'ci_lower'] <<- gmeta.fixed$combined.ci[1]
				fixed.effects.mcont[fixed.row, 'ci_upper'] <<- gmeta.fixed$combined.ci[2]
			}

			# run gmeta random-effects

			for(each.rmethod in random.methods) {
				for(each.tmethod in tau2.methods) {
					gmeta.random <- gmeta(gmi=summ.stats, gmi.type='pivot', 
						                  method=each.rmethod, tau2=each.tmethod,
										  gmo.xgrid=seq(min(summ.stats[, 1]) - 4*max(summ.stats[, 2]),
										 	            max(summ.stats[, 1]) + 4*max(summ.stats[, 2]),
										 	            by=0.001))
					random.row <- random.effects.mcont[, 'ma_id'] == each.ma &
								  random.effects.mcont[, 'method'] == each.rmethod &
								  random.effects.mcont[, 'tau2_method'] == each.tmethod
					random.effects.mcont[random.row, 'mean'] <<- gmeta.random$combined.mean
					random.effects.mcont[random.row, 'median'] <<- gmeta.random$combined.median
					random.effects.mcont[random.row, 'stdev'] <<- gmeta.random$combined.sd
					random.effects.mcont[random.row, 'ci_lower'] <<- gmeta.random$combined.ci[1]
					random.effects.mcont[random.row, 'ci_upper'] <<- gmeta.random$combined.ci[2]
					random.effects.mcont[random.row, 'tau2'] <<- gmeta.random$tau2 
				}
			}
		}
	}

}

# run the function

run.gmeta(large.studies.bin, "large", "bin")
run.gmeta(many.studies.bin, "many", "bin")
run.gmeta(large.studies.cont, "large", "cont")
run.gmeta(many.studies.cont, "many", "cont")

# print runtime
runtime <- proc.time() - runtime
print("Runtime:")
print(runtime)