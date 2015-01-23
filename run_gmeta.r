# runtime
runtime <- proc.time()

# set working directory
setwd("C:/Users/Ian/Desktop/meta/sandbox/")

# read in database files
big.studies <- read.csv("big_studies.csv")

# get M-A IDs
ma.id <- unique(big.studies[, 'ma_id'])

# for each M-A, do the following
##	1. calculate log odds ratio and std dev
## 	2. logit fixed conventional
##	3. logit random conventional
## 		a. each method for estimating tau-squared
## 	4. logit fixed robust
## 	5. logit random robust
## 		a. each method for estimating tau-squared

library(gmeta)

meta.CD.fixed.results <- matrix(nrow=length(ma.id), ncol=3)
colnames(meta.CD.fixed.results) <- c("mean", "median", "stdev")

meta.CD.exactOR.results <- matrix(nrow=length(ma.id), ncol=5)
colnames(meta.CD.exactOR.results) <- c("mean", "median", "stdev", "ci.LL", 
									                     "ci.UL")

meta.CD.MH.results <- matrix(nrow=length(ma.id), ncol=3)
colnames(meta.CD.MH.results) <- c("mean", "median", "stdev")

meta.CD.rand.list <- list(length=6)

run.gmeta <- function(ma, or.method, ma.method, count, tau2.method) {
	meta.2by2 <<- cbind(ma[, 'study_events1'], ma[, 'study_total1'], 
			                ma[, 'study_events2'], ma[, 'study_total2'])

	if(or.method == 'MH') {
		gmeta.MH <- gmeta(gmi = meta.2by2, gmi.type = '2x2', method = or.method,
			                gmo.xgrid = seq(-1, 1, by = 0.001))
		meta.CD.MH.results[count, "mean"] <<- log(gmeta.MH$combined.mean)
		meta.CD.MH.results[count, "median"] <<- log(gmeta.MH$combined.median)
		stdev.MH <- gmeta.MH$combined.sd/gmeta.MH$combined.mean
		meta.CD.MH.results[count, "stdev"] <<- stdev.MH
		if(sum(meta.2by2==0) > 0) {
			# tmp.nrow <- nrow(meta.2by2)
			# meta.2by2[which(meta.2by2==0)] <- NA
			# meta.2by2 <- meta.2by2[complete.cases(meta.2by2), ]
			# tmp.nrow <- tmp.nrow - nrow(meta.2by2)
			# if(tmp.nrow==1) {
			# 	print(paste("For ma_id", each.ma, tmp.nrow, "study was removed", 
			# 		          sep=" "))
			# } else {
			# 	print(paste("For ma_id", each.ma, tmp.nrow, "studies were removed", 
			# 				      sep=" "))
			#	}
			row_sub <<- apply(meta.2by2, 1, function(row) 0 %in% row)
			meta.2by2[row_sub, ] <<- meta.2by2[row_sub, ] + 0.5
			print(paste("For ma_id", each.ma, 
				          "Haldane continuity correction was applied"), sep=" ")	
			}
		}

	gmeta.logOR <<- gmeta(gmi = meta.2by2, gmi.type = '2x2', method = or.method,
						            gmo.xgrid = seq(-1, 1, by = 0.001))

	if(or.method == 'MH') {
		gmeta.means <- log(gmeta.logOR$individual.means)
		gmeta.sd <- gmeta.logOR$individual.stddevs/gmeta.logOR$individual.means
	} else {
		gmeta.means <- gmeta.logOR$individual.means
		gmeta.sd <- gmeta.logOR$individual.stddevs
	}

	if(or.method != 'exact1') {
		summ.stats <<- cbind(gmeta.means, gmeta.sd)
		gmeta.CD <<- gmeta(gmi = summ.stats, gmi.type = 'pivot', 
				               method = ma.method, tau2 = tau2.method,
				               gmo.xgrid = seq(-1, 1, by = 0.001))
		if(!grepl('random', ma.method)) {
			meta.CD.fixed.results[count, "mean"] <<- gmeta.CD$combined.mean
			meta.CD.fixed.results[count, "median"] <<- gmeta.CD$combined.median
			meta.CD.fixed.results[count, "stdev"] <<- gmeta.CD$combined.sd
		} else if(grepl('random', ma.method)) {
			meta.CD.random.results[count, "mean"] <<- gmeta.CD$combined.mean
			meta.CD.random.results[count, "median"] <<- gmeta.CD$combined.median
			meta.CD.random.results[count, "stdev"] <<- gmeta.CD$combined.sd
			meta.CD.random.results[count, "tau2"] <<- gmeta.CD$tau2
		}
	} else if(or.method == 'exact1') {
		meta.CD.exactOR.results[count, "mean"] <<- gmeta.logOR$combined.mean
		meta.CD.exactOR.results[count, "median"] <<- gmeta.logOR$combined.median
		meta.CD.exactOR.results[count, "stdev"] <<- gmeta.logOR$combined.sd
		meta.CD.exactOR.results[count, "ci.LL"] <<- gmeta.logOR$combined.ci[1]
		meta.CD.exactOR.results[count, "ci.UL"] <<- gmeta.logOR$combined.ci[2]
	}
}

## asymptotic normality assumption
or.method <- 'MH'
ma.method <- 'fixed-robust1'
count <- 1
for(each.ma in ma.id) {
	temp <- big.studies[big.studies[, 'ma_id'] == each.ma, ]
	run.gmeta(temp, or.method, ma.method, count, NULL)
	count <- count + 1
}

ma.method <- 'random-reml'
tau2.methods <- c('DL', 'HS', 'SJ', 'HE', 'ML', 'REML', 'EB')
list.count <- 1
for(each.tau2 in tau2.methods) {
	count <- 1
	meta.CD.random.results <- matrix(nrow=length(ma.id), ncol=4)
	colnames(meta.CD.random.results) <- c("mean", "median", "stdev", "tau2")
	for(each.ma in ma.id) {
		temp <- big.studies[big.studies[, 'ma_id'] == each.ma, ]
		run.gmeta(temp, or.method, ma.method, count, each.tau2)
		count <- count + 1
	}
	meta.CD.rand.list[[list.count]] <- meta.CD.random.results
	list.count <- list.count + 1
}

names(meta.CD.rand.list) <- tau2.methods

## exact LOR
or.method <- 'exact1'
count <- 1
for(each.ma in ma.id) {
	temp <- big.studies[big.studies[, 'ma_id'] == each.ma, ]
	run.gmeta(temp, or.method, NULL, count, NULL)
	count <- count + 1
}

colnames(meta.CD.exactOR.results)[4] <- paste(gmeta.logOR$ci.level*100, 
									                         colnames(meta.CD.exactOR.results)[4], 
										                       sep=".")
colnames(meta.CD.exactOR.results)[5] <- paste(gmeta.logOR$ci.level*100, 
                                           colnames(meta.CD.exactOR.results)[5], 
										                       sep=".")               

# runtime
runtime <- proc.time() - runtime
print("Runtime:")
print(runtime)