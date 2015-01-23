# load library
library(foreign)

# load database
cont.ma <- read.dta("C:/Users/Ian/Desktop/meta/sandbox/2012_Q1_continuous_data.dta")
cont.ma <- cont.ma[cont.ma[, 'estimable' == 'YES', ]]
cont.ma <- data.frame(cont.ma, ma_id = rep(NA, length=nrow(cont.ma)), ma_n_studies = rep(NA, length=nrow(cont.ma)))

# separate based on out_pool
cont.ma.Y <- cont.ma[cont.ma[, 'out_pool'] == 'Y', ]
cont.ma.S <- cont.ma[cont.ma[, 'out_pool'] == 'S', ]

#####
# start assigning in cont.ma.Y 
#####

# grab review_id
review.id.Y <- unique(cont.ma.Y[, 'review_id'])
idnum <- 1

for(each.rid.Y in review.id.Y) {
	same.rid.Y <- cont.ma.Y[cont.ma.Y[, 'review_id'] == each.rid.Y, ]
	comp.id.Y <- unique(same.rid.Y[, 'comp_id'])
	for(each.compid.Y in comp.id.Y) {
		same.compid.Y <- same.rid.Y[same.rid.Y[, 'comp_id'] == each.compid.Y, ]
		out.id.Y <- unique(same.compid.Y[, 'out_id'])
		for(each.outid.Y in out.id.Y) {
			index <- cont.ma[, 'out_pool'] == 'Y' &
					 cont.ma[, 'review_id'] == each.rid.Y &
					 cont.ma[, 'comp_id'] == each.compid.Y &
					 cont.ma[, 'out_id'] == each.outid.Y
			cont.ma[index, 'ma_id'] <- idnum
			idnum <- idnum + 1
		}
	}
}

review.id.S <- unique(cont.ma.S[, 'review_id'])

for(each.rid.S in review.id.S) {
	same.rid.S <- cont.ma.S[cont.ma.S[, 'review_id'] == each.rid.S, ]
	comp.id.S <- unique(same.rid.S[, 'comp_id'])
	for(each.compid.S in comp.id.S) {
		same.compid.S <- same.rid.S[same.rid.S[, 'comp_id'] == each.compid.S, ]
		out.id.S <- unique(same.compid.S[, 'out_id'])
		for(each.outid.S in out.id.S) {
			same.outid.S <- same.compid.S[same.compid.S[, 'out_id'] == each.outid.S, ]
			sub.id.S <- unique(same.outid.S[, 'subgroup_id'])
			for(each.subid.S in sub.id.S) {
				index <- cont.ma[, 'out_pool'] == 'S' &
						 cont.ma[, 'review_id'] == each.rid.S &
						 cont.ma[, 'comp_id'] == each.compid.S &
						 cont.ma[, 'out_id'] == each.outid.S &
						 cont.ma[, 'subgroup_id'] == each.subid.S
				cont.ma[index, 'ma_id'] <- idnum
				idnum <- idnum + 1
			}
		}
	}
}


cont.ma <- cont.ma[cont.ma[, 'out_pool'] != 'N', ]

ma.id.df <- data.frame(table(cont.ma[, 'ma_id']))

counter <- 1
for(each.id in ma.id.df[, 1]) {
	row.num <- cont.ma[, 'ma_id'] == each.id
	cont.ma[row.num, 'ma_n_studies'] <- ma.id.df[counter, 2]
	counter <- counter + 1
}

cont.ma <- cont.ma[with(cont.ma, order(ma_id)), ]
cont.ma <- data.frame(cont.ma, study_totaln=rep(NA, length=nrow(cont.ma)))
cont.ma[, "study_totaln"] <- cont.ma[, "study_total1"] + cont.ma[, "study_total2"]
write.csv(cont.ma, "2012_Q1_continuous_data.csv", row.names=F)
