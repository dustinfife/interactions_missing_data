clear()
require(MASS)
require(fifer)
require(xtable)
a =.3;
b=.3;
c=.3;
cor=.3
n=20000

			##### create a correlated X and Z
	sig = matrix(c(1, cor, cor, 1), nrow=2)
	d = data.frame(mvrnorm(n=n, mu=c(0,0), Sigma=sig))
	names(d) = c("z", "x")
	
			##### create interaction term
	d$zx = d$z*d$x

	
			##### create a y as a function of x,z, and x/z
	vy = 1-(a^2 + b^2 + c^2*var(d$zx) + 2*a*b*cor + 2*a*c*cov(d$x, d$zx) + 2*b*c*cov(d$z, d$zx))
	d$y = a*d$x + b*d$z + c*d$zx + rnorm(length(d$x), 0, sqrt(vy))

			# ##### convert all to sensible metrics
	# d$z = rescale(d$z, 500, 100)
	# d$x = rescale(d$x, 3, .4)
	# d$y = rescale(d$y, 3, .4)
	# d$zx = d$x*d$x

	#### convert all to mean deviant form (rescale back)
	n.selected = which(d$z<median(d$z))		
	full = d
	rest = full; rest[n.selected,c("y", "x", "zx")] = NA
	sample = full[sample(1:nrow(full), size=nrow(na.omit(rest))),]	
	


	
	
	require(selection)
	em.results = em(data.matrix(rest[,1:2]))
	pl =(mv.correction(na.omit(rest[,1:2]), 1, matrix(var(d$z), nrow=1)))
	
	covx.xz =var(d$x)*mean(d$z) + cov(d$z, d$x)*mean(d$x)
	covz.xz = var(d$z)*mean(d$x) + cov(d$z, d$x)*mean(d$z)
	
	##### estimate var(xz) from West and Aiken, page 179
	var.xz = round(var(d$z)*mean(d$x)^2 + var(d$x)*mean(d$z)^2 + 2*cov(d$z, d$x)*mean(d$x)*mean(d$z) + var(d$x)*var(d$z) + cov(d$x, d$z)^2, digits=2)

	corrected = rbind(pl, c(covz.xz, covx.xz))	
	corrected = cbind(corrected, c(covz.xz, covx.xz, var.xz)) 
	
	corrected
	cov(d)[1:3, 1:3]
	
	final.corrected =(mv.correction(data.matrix(cov(na.omit(rest))), p=3, v.pp=data.matrix(corrected)))
cov(na.omit(rest))[1:3,1:3]	

final.corrected[1:3, 1:3]

	n
	f = final.corrected
	final.corrected = data.frame(final.corrected)
	names(final.corrected) = c("SAT", "HSGPA", "SAT $\\times$ HSGPA", "First Year GPA")
	row.names(final.corrected) = names(final.corrected)
	corrected.cor = cov2cor(f)[2,4]
cor(na.omit(rest))[2,4]



	# names(rest) = c("SAT", "HSGPA", "SAT $\\times$ HSGPA", "First Year GPA")

# xx = xtable(rest, caption="Simulated Dataset of SAT, High School GPA, and First Year GPA Scores. Missing Cells Indicate Those Scores that Fall Below the Median of SAT.",
	# label="tab:combs", align="lcccc", round=c(1,2,2,2))
# print(xx, latex.environments=NULL, caption.placement="top", sanitize.text.function=function(x){x}, include.rownames=F)	
# cov(d)


# names(pl) = c("SAT", "GPA")
# row.names(pl) = names(pl)
# xx = xtable(pl, caption="Variance/Covariance Matrix of SAT and High School GPA, After Correcting for Range Restriction on SAT.",
	# label="tab:pl1", align="lcc", round=2)

