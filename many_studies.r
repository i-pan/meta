# set working directory to folder with database file
setwd("C:/Users/Ian/Desktop/meta/sandbox/")

extract.many.meta <- function(database, data.type) {

	# read database file 
		if(grepl("dta", database)) {
			library(foreign)
			all.meta <- read.dta(database)
		} else if(grepl("csv", database)) {
			all.meta <- read.csv(database)
		}

	# extract all M-As >= 25 studies
	all.meta <- all.meta[all.meta[, 'estimable'] == 'YES', ]
	many.studies <<- all.meta[all.meta[, 'ma_n_studies'] >= 25, ]

	file.name <- paste("many_studies_", data.type, ".csv", sep="")

	write.csv(many.studies, file.name, row.names=F)

}

extract.many.meta("2012_Q1_database_binary.dta", "bin")
print(paste("number of meta-analyses present (binary):", length(unique(many.studies$ma_id))))
extract.many.meta("2012_Q1_continuous_data.csv", "cont")
print(paste("number of meta-analyses present (continuous):", length(unique(many.studies$ma_id))))