set.seed(11212017)
		#### this file takes the results from the monte carlo and sees which parameters are predictive of bias
d = read.csv("research/interactions/data/mc_runif.csv")
head(d)

require(tidyverse)
p = ggplot(data=d, aes(x=c, y=CaseIII))
p + geom_point(alpha=.15) + geom_smooth() + scale_y_continuous(limits = c(-.1, .1)) + geom_hline(yintercept=0, col="red")


getwd()

nrow(d)
source("research/RPackages/fifer/R/rf.interp.r")

			##### perform random forest, only using the first 1000 rows (otherwise it takes a LONG time)
thresh = rfThresh(corrected~ a+b+c+cor+p.missing+n+mu.z+mu.x, data=d, nruns=20, importance="gini")
interp = rfInterp(thresh, nruns=20, importance="gini")
			### p.missing, c, a, cor, and mu.z are related (albeit barely)
save(interp, file= "research/interactions/data/rf_results.Rdat")
