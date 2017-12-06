#clear()

#### determine which datasets have been run already
load("research/interactions/data/last_done.Rdat")
num = mc_caseiii[mc_caseiii[,2]==1,1]
num = num[length(num)] + 1


		##### load next dataset
bias.mat = read.csv(paste0("research/interactions/data/monte_carlo_caseiii_", num, ".csv"))		
head(bias.mat)
		##### note that I've started another iteration
mc_caseiii = rbind(mc_caseiii, c(num, 0))
save(mc_caseiii, file="research/interactions/data/last_done.Rdat")

		##### read in global params
sd.z = sd.x = sd.y = 1; mu.y = 0
for (i in 1:nrow(bias.mat)){

		#### extract name objects
	a = bias.mat$a[i]
	b = bias.mat$b[i]
	c = bias.mat$c[i]
	cor = bias.mat$cor[i]
	n = bias.mat$n[i]
	mu.z = bias.mat$mu.z[i]
	mu.x = bias.mat$mu.x[i]
	skew = bias.mat$skew[i]
						
	sig = matrix(c(sd.z^2, sd.x*sd.z*cor, sd.x*sd.z*cor, sd.x^2), nrow=2)

	#### generate skewed data
	require(sn)
	d = data.frame(rmsn(n=n, xi=c(mu.z, mu.x), Omega = matrix(c(1, cor, cor, 1), nrow=2), alpha=c(skew, skew)))
	names(d) = c("z", "x")	

			##### create interaction term
	d$zx = d$z*d$x
	
			#### create new weights to retain proper metrix (a' = a*sx/sy)	
	ap = a*(sd.y/sd.x)
	bp = b*(sd.y/sd.z)
	cov = cor*(sd.x*sd.z)

			##### compute expected covariance and variance of interaction term
	covx.xz =sd.x^2*mu.z + cov*mu.x
	covz.xz =sd.z^2*mu.x + cov*mu.z
	var.xz = sd.z^2*mu.x^2 + sd.x^2*mu.z^2 + 2*cov*mu.x*mu.z + sd.x^2*sd.z^2 + cov^2	

			##### finish recreating weights	
	cp = c*(sd.y/sqrt(var.xz))
	mu.xz = mu.x*mu.z + cov


			##### compute expected population correlation
	pop = (ap*sd.x^2 + bp*cov + cp*(covx.xz) )/(sd.x*sd.y)


			##### compute variance of y
	vy = sd.y^2-(ap^2*sd.x^2 + bp^2*sd.z^2 + cp^2*var(d$zx) + 2*ap*bp*cov + 2*ap*cp*cov(d$x, d$zx) + 2*bp*cp*cov(d$z, d$zx))

			##### compute intercept of y
	b0 = mu.y - (ap* mu.x + bp* mu.z + cp*mu.xz)


			##### create y
	d$y = b0 + ap*d$x + bp*d$z + cp*d$zx + rnorm(length(d$x), 0, sqrt(vy))	

			##### now select on z and create the datasets
	n.selected = which(d$z<quantile(d$z, bias.mat$p.missing[i]))		
	full = d
	rest = full; rest[n.selected,c("x", "zx", "y")] = NA
	sample = full[sample(1:nrow(full), size=nrow(na.omit(rest))),]

			##### correct using Case III
	bias.mat$CaseIII[i] = caseIII(rest[,c(1,2,4)])
	bias.mat$Sample[i] = cor(sample$x, sample$y)
	bias.mat$uncorrected[i] = cor(rest$x, rest$y, use="complete.obs")

			#### now do the correction
			#### 1. do PL on the 2x2 matrix
	pl = mv.correction(na.omit(rest[,1:2]), 1, matrix(var(d$z), nrow=1))			
			#### 2. Use West and Aiken's formulas for the remainder of the matrix
	covx.xz = pl[2,2]*mean(d$z) + pl[1,2]*mean(d$x)
	covz.xz = pl[1,1]*mean(d$x) + pl[1,2]*mean(d$z)
	var.xz = round(pl[1,1]*mean(d$x)^2 + pl[2,2]*mean(d$z)^2 + 2*pl[1,2]*mean(d$x)*mean(d$z) + pl[2,2]*pl[1,1] + pl[1,2]^2, digits=2)
	corrected = rbind(pl, c(covz.xz, covx.xz))	
	corrected = cbind(corrected, c(covz.xz, covx.xz, var.xz)) 
			##### 3. Use PL again to correct the matrix
	final.corrected =cov2cor(mv.correction(data.matrix(cov(na.omit(rest))), p=3, v.pp=data.matrix(corrected)))


	bias.mat$corrected[i] = final.corrected[2,4]
		


		#### track progress
	# if (i/1000 == round(i/1000)){
		# cat(paste0("Iteration ", i, " of ", nrow(bias.mat), "\n"))
		# 
	# }
	
}

write.csv(bias.mat, paste0("research/interactions/data/monte_carlo_caseiii_", num, ".csv"), row.names=F)		

	#### if the monte carlo has been completed, update the tally
if (i==nrow(bias.mat)){
	mc_caseiii[nrow(mc_caseiii), 2] = 1
	save(mc_caseiii, file="research/interactions/data/last_done.Rdat")
}	
