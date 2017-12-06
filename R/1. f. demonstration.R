set.seed(11212017)
		#### this file takes the results from the monte carlo and sees which parameters are predictive of bias
d = read.csv("research/interactions/data/mc_runif.csv")
head(d)

require(tidyverse)
p = ggplot(data=d, aes(x=c, y=CaseIII))
p + geom_point(alpha=.15) + geom_smooth() + scale_y_continuous(limits = c(-.1, .1)) + geom_hline(yintercept=0, col="red")

coplot(CaseIII~c|p.missing + a, data=d, panel=panel.smooth, ylim=c(-.15, .15), pch=16, cex=.5, lwd=3)
#### a/c/p.missing have an effect together


coplot(CaseIII~c|p.missing + cor, data=d, panel=panel.smooth, ylim=c(-.15, .15), pch=16, cex=.5, lwd=3)
#### cor/c/p.missing have an effect together		

coplot(CaseIII~c|p.missing + mu.z, data=d, panel=function(x,y,...){panel.smooth(x,y);abline(h=0)}, ylim=c(-.15, .15), pch=16, cex=.5, lwd=3)
#### mu.z/c/p.missing have an effect together		

source("research/RPackages/fifer/R/rf.interp.r")

			##### perform random forest, only using the first 1000 rows (otherwise it takes a LONG time)
thresh = rfThresh(case~ a+b+c+cor+p.missing+n+mu.z+mu.x, data=d, nruns=20, importance="gini")
interp = rfInterp(thresh, nruns=20, importance="gini")
			### p.missing, c, a, cor, and mu.z are related (albeit barely)
save(interp, file= "research/interactions/data/rf_results.Rdat")



#high missingness, low c, mean less than 20

d.messed = subset(d, p.missing>.6 & c>.5)
coplot(CaseIII~mu.z|p.missing + a, data=d, panel=panel.smooth, ylim=c(-.25, .25), pch=16, cex=.5, lwd=3)

coplot(CaseIII~cor|p.missing + b, data=d, panel=panel.smooth, ylim=c(-.25, .25), pch=16, cex=.5, lwd=3)



coplot(CaseIII~c|p.missing + mu.x, data=d, panel=panel.smooth, ylim=c(-.25, .25), pch=16, cex=.5, lwd=3)




coplot(CaseIII~mu.z|p.missing + a, data=d, panel=panel.smooth, ylim=c(-.25, .25), pch=16, cex=.5, lwd=3)
coplot(CaseIII~mu.z|p.missing + a, data=d, panel=panel.smooth, ylim=c(-.25, .25), pch=16, cex=.5, lwd=3)





nrow(d)
#60% missing, c = -.5, b = .1, mu.z = 3, mu.x = 0 (standardized scale)



coplot(CaseIII~c|mu.z + mu.x, data=subset(d, p.missing>.6), panel=panel.smooth, ylim=c(-.25, .25), pch=16, cex=.5, lwd=3)
