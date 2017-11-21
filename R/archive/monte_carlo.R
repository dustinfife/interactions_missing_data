clear()


		##### read in MC parameters from knitr file
params = read.csv("research/caseiii/data/mc_parameters.csv")
fix.func = function(x){as.numeric(unlist(strsplit(as.character(x), ", ")))}
a = fix.func(params[1,2])
b= fix.func(params[2,2])
c= fix.func(params[3,2])
cor= fix.func(params[4,2])
n= fix.func(params[7,2])
iterations=100
sd.z = 1; mu.z = fix.func(params[6,2])
sd.x = 1; mu.x = fix.func(params[5,2])
sd.y = 1; mu.y = 0


		##### read in packages
require(selection)
require(fifer)

		##### preallocate
bias.mat = expand.grid(a=a,b=b,c=c,cor=cor,n=n,mu.z=mu.z, mu.x=mu.x, p.missing=fix.func(params[8,2]), i=1:iterations, Sample=NA, caseIII=NA, corrected=NA)

		##### randomly sort bias.mat so we get equal probability of each combination
bias.mat = bias.mat[sample(1:nrow(bias.mat), nrow(bias.mat)), ]
i=1
for (i in 1:nrow(bias.mat)){

		#### extract name objects
	a = bias.mat$a[i]
	b = bias.mat$b[i]
	c = bias.mat$c[i]
	cor = bias.mat$cor[i]
	n = bias.mat$n[i]
	mu.z = bias.mat$mu.z[i]
	mu.x = bias.mat$mu.x[i]
							
			##### create a correlated X and Z
	sig = matrix(c(sd.z^2, sd.x*sd.z*cor, sd.x*sd.z*cor, sd.x^2), nrow=2)
	d = data.frame(mvrnorm(n=n, mu=c(mu.z, mu.x), Sigma=sig))
	names(d) = c("z", "x")
	
	#### standardize them
	d$z = d$z-mean(d$z)
	d$x = d$x-mean(d$x)
	mu.xn = 0; mu.zn = 0

			##### create interaction term
	d$zx = d$z*d$x
	
			#### create new weights to retain proper metrix (a' = a*sx/sy)	
	ap = a*(sd.y/sd.x)
	bp = b*(sd.y/sd.z)
	cov = cor*(sd.x*sd.z)

			##### compute expected covariance and variance of interaction term
	covx.xz =sd.x^2*mu.zn + cov*mu.xn
	covz.xz =sd.z^2*mu.xn + cov*mu.zn
	var.xz = sd.z^2*mu.xn^2 + sd.x^2*mu.zn^2 + 2*cov*mu.xn*mu.zn + sd.x^2*sd.z^2 + cov^2	

			##### finish recreating weights	
	cp = c*(sd.y/sqrt(var.xz))
	mu.xz = mu.xn*mu.zn + cov


			##### compute expected population correlation
	pop = (ap*sd.x^2 + bp*cov + cp*(covx.xz) )/(sd.x*sd.y)


			##### compute variance of y
	vy = sd.y^2-(ap^2*sd.x^2 + bp^2*sd.z^2 + cp^2*var(d$zx) + 2*ap*bp*cov + 2*ap*cp*cov(d$x, d$zx) + 2*bp*cp*cov(d$z, d$zx))

			##### compute intercept of y
	b0 = mu.y - (ap* mu.xn + bp* mu.zn + cp*mu.xz)


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
	corrected = rbind(pl, c(covz.xz, covx.xz))	
	corrected = cbind(corrected, c(covz.xz, covx.xz, var.xz)) 
			##### 3. Use PL again to correct the matrix
	final.corrected =cov2cor(mv.correction(data.matrix(cov(na.omit(rest))), p=3, v.pp=data.matrix(corrected)))


	bias.mat$corrected[i] = final.corrected[2,4]
		


		#### track progress
	if (i/100 == round(i/100)){
		cat(paste0("Iteration ", i, " of ", iterations, "\n"))
		write.csv(bias.mat, "research/caseIII/data/monte_carlo_caseIII.csv", row.names=F)		
	}
	
}
i