#### this simulation randomly varies a bunch of parameters (slopes between each of predictors, means, correlations, etc). This is a preliminary simulation to determine which 
#### parameters actually make a difference in bias (or standard errors). 

iterations = 10000

		#### preallocated
bias.mat = data.frame(iteration=1:iterations, a=NA, b=NA, c=NA, cor=NA, skew=NA, p.missing=NA, n=1000, mu.z=NA, mu.x=NA)

for (i in 1:nrow(bias.mat)){

		#### randomly sample parameters
	a = runif(1, 0, .3); bias.mat$a[i]=a
	b = runif(1, 0, .3); bias.mat$b[i]=b
	c = runif(1, 0, .3); bias.mat$c[i]=c
	mu.z = runif(1, 0, 100); bias.mat$mu.z[i]=mu.z
	mu.x = runif(1, 0, 100); bias.mat$mu.x[i]=mu.x
	sd.y = 1;sd.x=1; sd.z=1; mu.y=0
	p.missing = runif(1, 0, .9); bias.mat$p.missing[i] = p.missing			
	skew=runif(1, -100, 100); bias.mat$skew[i] = skew
	cor = runif(1, 0, .5); 	bias.mat$cor[i] = cor
	n = round(runif(1, 50, 500))
	
			##### create a correlated X and Z
	sig = matrix(c(1, cor, cor, 1), nrow=2)

	#### generate skewed data
	require(sn)
	d = data.frame(rmsn(n=bias.mat$n[i], xi=c(mu.z, mu.x), Omega = matrix(c(1, cor, cor, 1), nrow=2), alpha=c(skew, skew)))
	names(d) = c("z", "x")

	#### standardize them
	d$z = d$z-mean(d$z)
	d$x = d$x-mean(d$x)
	mu.xn = 0; mu.zn = 0
	

		##### create interaction term
	d$zx = d$z*d$x
	
		#### create new weights to retain proper metrix (a' = a*sx/sy)	
	ap = a*(sd.y/sd(d$x))
	bp = b*(sd.y/sd(d$z))
	cov = cor(d)[1,2]

			##### compute expected covariance and variance of interaction term (must be done empirically because skewness throws things off)
	covx.xz =sd(d$x)^2*mu.zn + cov*mu.xn
	covz.xz =sd(d$z)^2*mu.xn + cov*mu.zn
	var.xz = sd(d$z)^2*mu.xn^2 + sd(d$x)^2*mu.zn^2 + 2*cov*mu.xn*mu.zn + sd(d$x)^2*sd(d$z)^2 + cov^2	


			##### finish recreating weights	
	cp = c*(sd.y/sqrt(var.xz))
	mu.xz = mu.xn*mu.zn + cov


			##### compute expected population correlation  (computed empirically because var.xz and covz.xz, etc. assume skewness)
	pop = (ap*sd(d$x)^2 + bp*cov + cp*cov(d$z, d$zx) )/(sd(d$x)*sd.y)

			##### compute variance of y (this will be under a variance of one when there's skewness)
	vy = sd.y^2-(a^2 + b^2 + c^2*var(d$zx) + 2*a*b*cov + 2*a*c*cov(d$z,d$zx) + 2*b*c*cov(d$z, d$zx))

			##### compute intercept of y
	b0 = mu.y - (a* mu.xn + b* mu.zn + c*mu.xz)


			##### create y
	d$y = b0 + ap*d$x + bp*d$z + cp*d$zx + rnorm(length(d$x), 0, sqrt(vy))	


			##### now select on z and create the datasets
	n.selected = which(d$z<quantile(d$z, bias.mat$p.missing[i]))		
	full = d
	rest = full; rest[n.selected,c("x", "zx", "y")] = NA
	sample = full[sample(1:nrow(full), size=nrow(na.omit(rest))),]

			##### correct using Case III
	bias.mat$CaseIII[i] = caseIII(rest, y=4)-pop
	bias.mat$Sample[i] = cor(sample$x, sample$y)-pop

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


	bias.mat$corrected[i] = final.corrected[2,4]-pop
		


		# #### track progress
	if (i/100 == round(i/100)){
		cat(paste0("Iteration ", i, " of ", nrow(bias.mat), "\n"))
		#write.csv(bias.mat, paste0("research/caseiii/data/monte_carlo_caseiii_", num, ".csv"), row.names=F)		
	}
	pop
}

d = bias.mat

write.csv(bias.mat, "research/interactions/data/mc_runif.csv", row.names=F)		