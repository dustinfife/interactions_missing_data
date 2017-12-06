clear()
load("research/interactions/data/rf_results.Rdat")

labs = c(	expression(mu[Z]), 
			expression(r[xz]),
			expression(beta[y.z]),			
			expression(beta[y.xz]),
			expression(p[missing]))

interp2 = data.frame(err.interp = interp$err.interp, varselect.interp = interp$varselect.interp)
p = ggplot(data=interp2,
	mapping = aes(x=reorder(varselect.interp, err.interp, decreasing=T),
					y = err.interp))
p + geom_point() + labs(x="Variable", y="Stepwise OOB Error") + scale_x_discrete(labels=labs) + labs(x="")


		### let's visualize each of these
d = read.csv("research/interactions/data/mc_runif.csv")		
p = ggplot(data=d,
	mapping = aes(x=mu.z,
					y = corrected))
p + geom_point() + geom_smooth() + geom_hline(yintercept=0, size=2, alpha=.2) + scale_y_continuous(limits=c(-.05, .05))
		#### the realtionship is really tiny (if it exists at all)
		
		
p = ggplot(data=d,
	mapping = aes(x=cor,
					y = corrected))
p + geom_point() + geom_smooth() + geom_hline(yintercept=0, size=2, alpha=.2) + scale_y_continuous(limits=c(-.01, .01))
		#### also a very tiny relationship
		
		
		
p = ggplot(data=d,
	mapping = aes(x=a,
					y = corrected))
p + geom_point() + geom_smooth() + geom_hline(yintercept=0, size=2, alpha=.2) + scale_y_continuous(limits=c(-.01, .01))
		#### also a very tiny relationship		


p = ggplot(data=d,
	mapping = aes(x=c,
					y = corrected))
p + geom_point() + geom_smooth() + geom_hline(yintercept=0, size=2, alpha=.2) + scale_y_continuous(limits=c(-.01, .01)) + facet_wrap(~ cut(p.missing, 6))
head(d)


coplot(corrected~c|p.missing + a, data=d, panel=panel.smooth, ylim=c(-.1, .1), pch=16, cex=.5, lwd=3)
		#### a/c/p.missing have an effect together
		
		
coplot(corrected~c|p.missing + cor, data=d, panel=panel.smooth, ylim=c(-.1, .1), pch=16, cex=.5, lwd=3)
		#### cor/c/p.missing have an effect together		

coplot(corrected~mu.z|p.missing + cor, data=d, panel=panel.smooth, ylim=c(-.1, .1), pch=16, cex=.5, lwd=3)
		#### mu.z/c/p.missing have an effect together		
		
		
		#### the final parameters to be estimated are mu.z, c, p.missing, cor, and skew
