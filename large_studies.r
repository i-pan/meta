# set working directory to folder with database file
setwd("C:/Users/Ian/Desktop/meta/sandbox/")

extract.large.meta <- function(database, data.type) {

	# read database file 
	if(grepl("dta", database)) {
		library(foreign)
		all.meta <- read.dta(database)
	} else if(grepl("csv", database)) {
		all.meta <- read.csv(database)
	}

	# extract all meta-analyses with >= 8 studies
	large.meta <- all.meta[all.meta['ma_n_studies'] >= 8, ]

	if(data.type == "bin") {
	# extract all studies with >= 300 per arm
		large.meta.N <- large.meta[large.meta['study_total1'] >= 400 &
								large.meta['study_total2'] >= 400, ]
	} else if(data.type == "cont") {
	# extract all studies with >= 250 total 
		large.meta.N <- large.meta[large.meta['study_totaln'] >= 400, ]
	}

	# count how many studies met criteria in each M-A
	ma.id.df <- data.frame(table(large.meta.N['ma_id']))

	# extract the total number of studies in each M-A

	for(i in 1:length(ma.id.df[, 1])) {
		ma.id.df[i, 3] <- unique(
							    large.meta.N[which(large.meta.N[,'ma_id'] == 
												  ma.id.df[i, 1]), 
							    				  'ma_n_studies']
							    )
	}

	# get IDs of M-As where at least 75% of studies met criteria
	frac <- ma.id.df[, 2]/ma.id.df[, 3]
	ma.id.df <- cbind(ma.id.df, frac)
	colnames(ma.id.df) <- c('ma_id', 'geq300', 'total_n_studies', 'fracgeq300')
	large.studies.id <- ma.id.df[which(ma.id.df[ ,4] >= 0.75), 1]
	 
	large.studies.id <- as.numeric(levels(large.studies.id))[large.studies.id]

	# extract all data for those studies
	tmp.index <- which(all.meta['ma_id'][ ,1] %in% large.studies.id) 
	large.studies <<- all.meta[tmp.index, ]

	file.name <- paste("large_studies_", data.type, ".csv", sep="")

	write.csv(large.studies, file.name, row.names=F)

}

extract.large.meta("2012_Q1_database_binary.dta", "bin")
print(paste("number of meta-analyses present (binary):", length(unique(large.studies$ma_id))))
extract.large.meta("2012_Q1_continuous_data.csv", "cont")
print(paste("number of meta-analyses present (continuous):", length(unique(large.studies$ma_id))))