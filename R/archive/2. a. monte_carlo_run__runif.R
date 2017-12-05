
#### this simulation randomly varies a bunch of parameters (slopes between each of predictors, means, correlations, etc). This is a preliminary simulation to determine which 
#### parameters actually make a difference in bias (or standard errors). 

iterations = 1000	
		#### preallocated
bias.mat = data.frame(iteration=1:iterations, a=NA, b=NA, c=NA, cor=NA, skew=NA, p.missing=NA, n=1000, mu.z=NA, mu.x=NA)
i=1
for (i in 1:nrow(bias.mat)){

		#### randomly sample parameters
	a = runif(1, 0, .3); bias.mat$a[i]=a
	b = runif(1, 0, .3); bias.mat$b[i]=b
	c = runif(1, 0, .3); bias.mat$c[i]=c
		##### randomly select means
	mu.z = runif(1, 0, 100); bias.mat$mu.z[i]=mu.z
	mu.x = runif(1, 0, 100); bias.mat$mu.x[i]=mu.x
	mu.y = runif(1, 0, 100); bias.mat$mu.y[i]=mu.y
		##### randomly select sds
	sd.z = runif(1, 1, 20); bias.mat$sd.z[i]=sd.z
	sd.x = runif(1, 1, 20); bias.mat$sd.x[i]=sd.x
	sd.y = runif(1, 1, 20); bias.mat$sd.y[i]=sd.y	
	p.missing = runif(1, 0, .9); bias.mat$p.missing[i] = p.missing			
	skew=runif(1, -100, 100); bias.mat$skew[i] = skew
	cor = runif(1, 0, .5); 	bias.mat$cor[i] = cor
	n = round(runif(1, 50, 500))
	
	
			##### create a correlated X and Z
	sig = matrix(c(1, cor, cor, 1), nrow=2)

		#### create new weights to retain proper metrix (a' = a*sx/sy)	
	ap = a*(sd.y/sd.x)
	bp = b*(sd.y/sd.z)
	cov = cor*sd.x*sd.z
	
	#### generate skewed data
	require(sn)
	d = data.frame(rmsn(n=bias.mat$n[i], xi=c(mu.z, mu.x), Omega = matrix(c(sd.z^2,cov, cov, sd.x^2), nrow=2), alpha=c(skew, skew)))
	#d = data.frame(mvrnorm(n, mu=c(mu.z, mu.x), Sigma=matrix(c(sd.z^2,cov, cov, sd.x^2), nrow=2)))
	names(d) = c("z", "x")	

		##### create interaction term
	d$zx = d$z*d$x

		##### compute expected covariance and variance of interaction term (must be done empirically because skewness throws things off)
	covx.xz =sd.x^2*mu.z + cov*mu.x
	covz.xz =sd.z^2*mu.x + cov*mu.z
	var.xz = sd.z^2*mu.x^2 + sd.x^2*mu.z^2 + 2*cov*mu.x*mu.z + sd.x^2*sd.z^2 + cov^2		

			##### compute correlations (standardized) with interaction term
	corx.xz = covx.xz/(sd.x*sqrt(var.xz))	
	corz.xz = covz.xz/(sd.z*sqrt(var.xz))		
	var.xz.standardized = 1 + cor^2		#### variance of interaction if x/z are standardized

			##### finish recreating weights	
	cp = c*(sd.y/sqrt(var.xz))
	mu.xz = mu.x*mu.z + cov


			##### compute expected population correlation  (computed empirically because var.xz and covz.xz, etc. assume skewness)
	pop = (ap*sd.x^2 + bp*cov + cp*covx.xz)
	pop.cov = pop
	
	pop = pop/(sd.x*sd.y)	

			##### compute explained variance of y IF ALL VARIABLES ARE STANDARDIZES
	exp.y.st = (a^2 + b^2 + c^2 + 2*a*b*cor + 2*a*c*corx.xz + 2*b*c*corz.xz)

			##### compute explained varince in y (unstandardized)
	exp.y.unst = (ap^2*sd.x^2 + bp^2*sd.z^2 + cp^2*var.xz + 2*ap*bp*cov + 2*ap*cp*covx.xz + 2*bp*cp*covz.xz)
	
			##### compute total variance of y (if unstandardized)
	total.y.var = exp.y.unst/exp.y.st

			##### compute intercept of y
	b0 = mu.y - (a* mu.x + b* mu.z + c*mu.xz)

			##### create y
	d$y = b0 + ap*d$x + bp*d$z + cp*d$zx + rnorm(length(d$x), 0, sqrt(sd.y^2-exp.y.unst))	



			##### now select on z and create the datasets
	n.selected = which(d$z<quantile(d$z, bias.mat$p.missing[i]))		
	full = d
	rest = full; rest[n.selected,c("x", "zx", "y")] = NA
	sample = full[sample(1:nrow(full), size=nrow(na.omit(rest))),]

			##### correct using Case III
	require(selection)			
	bias.mat$CaseIII[i] = caseIII(rest, y=4)
	bias.mat$Sample[i] = cor(sample$x, sample$y)

			#### now do the correction
			#### 1. do PL on the 2x2 matrix
	pl = mv.correction(na.omit(rest[,1:2]), 1, matrix(var(d$z), nrow=1))			

			#### 2. Use West and Aiken's formulas for the remainder of the matrix
	covx.xz = pl[2,2]*mean(d$z) + pl[1,2]*mean(d$x)
	covz.xz = pl[1,1]*mean(d$x) + pl[1,2]*mean(d$z)
	var.xz = round(pl[1,1]*mean(d$x)^2 + pl[2,2]*mean(d$z)^2 + 2*pl[1,2]*mean(d$x)*mean(d$z) + pl[2,2]*pl[1,1] + pl[1,2]^2, digits=2)
	
			##### 3. create sigma prime
	corrected = rbind(pl, c(covz.xz, covx.xz))	
	corrected = cbind(corrected, c(covz.xz, covx.xz, var.xz)) 

			##### 4. Use PL again to correct the matrix
	final.corrected =cov2cor(mv.correction(data.matrix(cov(na.omit(rest))), p=3, v.pp=data.matrix(corrected)))
	bias.mat$corrected[i] = final.corrected[2,4]
	
	
			##### compute slope differences (to see if they're affected)
	mod.full = lm(y~x+z+zx, data=full)		
	mod.rest = lm(y~x+z+zx, data=rest)

		


		# #### track progress
	if (i/100 == round(i/100)){
		cat(paste0("Iteration ", i, " of ", nrow(bias.mat), "\n"))
		#write.csv(bias.mat, paste0("research/caseiii/data/monte_carlo_caseiii_", num, ".csv"), row.names=F)		
	}
	pop
}
head(bias.mat)
d = bias.mat

write.csv(bias.mat, "research/interactions/data/mc_runif.csv", row.names=F)		