#### this function computes the inter-correlation matrix for the l-correlations. These computations take a LONG time, so I'm computing them once, outputting to a 
#### file, and loading into the monte carlo so I don't have to keep recomputing them

clear()
cors = c(0, .1, .3, .5)
skew = c(0, .3, .5, .7)
params = expand.grid(cors=cors, skew=skew)
inter.correlations = list()
source("research/RPackages/lcorrelations/R/Headrick.R")
i=1
for (i in 1:nrow(params)){
	skew = params$skew[i]
	cors = params$cors[i]
	
	sig = matrix(c(1, cors, cors, 1), nrow=2)
	inter = compute.inter(sig, c(skew, skew), c(0,0))	
	inter.correlations[i] = inter[1,2]; names(inter.correlations)[i]= paste0("skew=", skew, "; cors=", cors)
	print(i)
}

save(inter.correlations, file="research/caseiii/data/lcorrelations.Rdat")
