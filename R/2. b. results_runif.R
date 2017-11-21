set.seed(11212017)
		#### this file takes the results from the monte carlo and sees which parameters are predictive of bias
d = read.csv("research/interactions/data/mc_runif.csv")


			##### perform random forest, only using the first 1000 rows (otherwise it takes a LONG time)
thresh = rfThresh(corrected~a+b+c+cor+skew+p.missing+n+mu.z+mu.x, data=d, nruns=10, importance="gini")
interp = rfInterp(thresh, nruns=10, importance="gini")
plot(interp)
			##### skewness, proportion missing, cor, and b alone are important
