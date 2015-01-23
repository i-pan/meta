# runtime
runtime <- proc.time()

# set working directory (contains files to work with)
setwd("C:/Users/Ian/Desktop/meta/sandbox/")

# read in meta-analyses from database
large.studies <- read.csv("large_studies_bin.csv") # large sample size

# grab meta-analysis IDs
ma.id <- sort(unique(large.studies[, "ma_id"]))

# load libraries
library(gmeta)

#####
# create skeleton dataframes to be filled in later
#####

# fixed-effects

fixed.methods <- c("fixed-mle", "fixed-robust1")

fixed.nrows <- length(ma.id)*length(fixed.methods)

fixed.ma.id <- sort(rep(ma.id, length=fixed.nrows))

fixed.effects <- data.frame(ma_id=fixed.ma.id, 
							method=rep(fixed.methods, length=fixed.nrows),
							mean=rep(NA, length=fixed.nrows),
							median=rep(NA, length=fixed.nrows),
							stdev=rep(NA, length=fixed.nrows),
							ci_lower=rep(NA, length=fixed.nrows),
							ci_upper=rep(NA, length=fixed.nrows))

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

random.effects <- data.frame(ma_id=random.ma.id,
							 method=rep(random.methods.list, 
							 	        length=random.nrows),
							 tau2_method=rep(tau2.methods, 
							 	             length=random.nrows),
							 mean=rep(NA, length=random.nrows),
							 median=rep(NA, length=random.nrows),
							 stdev=rep(NA, length=random.nrows),
							 ci_lower=rep(NA, length=random.nrows),
							 ci_upper=rep(NA, length=random.nrows),
							 tau2=rep(NA, length=random.nrows))

#####
# run gmeta
#####

for(each.ma in ma.id) {

	# grab individual one m-a at a time
	ma.raw <- large.studies[large.studies[, 'ma_id'] == each.ma, ]

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
		fixed.row <- fixed.effects[, 'ma_id'] == each.ma & 
		             fixed.effects[, 'method'] == each.fmethod
		fixed.effects[fixed.row, 'mean'] <- gmeta.fixed$combined.mean
		fixed.effects[fixed.row, 'median'] <- gmeta.fixed$combined.median
		fixed.effects[fixed.row, 'stdev'] <- gmeta.fixed$combined.sd 
		fixed.effects[fixed.row, 'ci_lower'] <- gmeta.fixed$combined.ci[1]
		fixed.effects[fixed.row, 'ci_upper'] <- gmeta.fixed$combined.ci[2]
	}

	# run gmeta random-effects

	for(each.rmethod in random.methods) {
		for(each.tmethod in tau2.methods) {
			gmeta.random <- gmeta(gmi=summ.stats, gmi.type='pivot', 
				                  method=each.rmethod, tau2=each.tmethod)
			random.row <- random.effects[, 'ma_id'] == each.ma &
						  random.effects[, 'method'] == each.rmethod &
						  random.effects[, 'tau2_method'] == each.tmethod
			random.effects[random.row, 'mean'] <- gmeta.random$combined.mean
			random.effects[random.row, 'median'] <- gmeta.random$combined.median
			random.effects[random.row, 'stdev'] <- gmeta.random$combined.sd
			random.effects[random.row, 'ci_lower'] <- gmeta.random$combined.ci[1]
			random.effects[random.row, 'ci_upper'] <- gmeta.random$combined.ci[2]
			random.effects[random.row, 'tau2'] <- gmeta.random$tau2 
		}
	}
}

# print runtime
runtime <- proc.time() - runtime
print("Runtime:")
print(runtime)