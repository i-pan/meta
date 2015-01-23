nf1 = ma.raw[, "study_total1"]
nf2 = ma.raw[, "study_total2"]
mf1 = ma.raw[, "study_mean1"]
mf2 = ma.raw[, "study_mean2"]
var1 = ma.raw[, "study_sd1"]^2
var2 = ma.raw[, "study_sd2"]^2
s.pool = sqrt(((nf1 - 1)*var1 + (nf2 - 1)*var2)/(nf1 + nf2 - 2))
J = 1 - 3/(4*(nf1 + nf2) - 9)
hedges.g = J*(mf1 - mf2)/s.pool
var.g = 1/nf1 + 1/nf2 + hedges.g^2/(2*(nf1 + nf2))

summ.stats = cbind(hedges.g, var.g)
