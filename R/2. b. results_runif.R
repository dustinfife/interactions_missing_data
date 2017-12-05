set.seed(11212017)
		#### this file takes the results from the monte carlo and sees which parameters are predictive of bias
d = read.csv("research/interactions/data/mc_runif.csv")
source("research/RPackages/fifer/R/rf.interp.r")

			##### perform random forest, only using the first 1000 rows (otherwise it takes a LONG time)
thresh = rfThresh(corrected~a+b+c+cor+skew+p.missing+n+mu.z+mu.x, data=d, nruns=20, importance="gini")
require(randomForest)
interp = rfInterp(thresh, nruns=20, importance="gini")
plot(interp)

interp$err.interp


labs = c(	expression(mu[X]), 
			expression(beta[y.xz]),
			expression(p[missing]),
			expression(beta[y.x]))

require(ggplot2)
require(ggrepel)
interp2 = data.frame(err.interp = interp$err.interp, varselect.interp = interp$varselect.interp)
p = ggplot(data=interp,
	mapping = aes(x=reorder(varselect.interp, err.interp, decreasing=T),
					y = err.interp))
p + geom_point() + labs(x="Variable", y="Stepwise OOB Error") + 
cartesian_coord()

scale_x_discrete(labels=c(labs[1], labs[2], labs[3], labs[4], "skew"))
			##### skewness, proportion missing, cor, and b alone are important

		### let's visualize each of these
p = ggplot(data=d,
	mapping = aes(x=mu.x,
					y = corrected))
p + geom_point() + geom_smooth() + geom_hline(yintercept=0, size=2, alpha=.2)

p = ggplot(data=d,
	mapping = aes(x=skew,
					y = corrected))
p + geom_point() + geom_smooth() + geom_hline(yintercept=0, size=2, alpha=.2)


p = ggplot(data=d,
	mapping = aes(x=b,
					y = corrected))
p + geom_point() + geom_smooth() + geom_hline(yintercept=0, size=2, alpha=.2) + labs(x=expression(beta[y.x]))

p = ggplot(data=d,
	mapping = aes(x=p.missing,
					y = corrected))
p + geom_point() + geom_smooth() + geom_hline(yintercept=0, size=2, alpha=.2) + labs(x=expression(p[missing]))


p = ggplot(data=d,
	mapping = aes(x=c,
					y = corrected))
p + geom_point() + geom_smooth() + geom_hline(yintercept=0, size=2, alpha=.2) + labs(x=expression(b[y.xz]))

p = ggplot(data=d,
	mapping = aes(x=cor,
					y = corrected))
p + geom_point() + geom_smooth() + geom_hline(yintercept=0, size=2, alpha=.2) + labs(x=expression(r[xz]))

colMeans(d)

head(d)