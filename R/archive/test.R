	clear()
	
	a =.3;
	b=.3;
	c=.3;
	cor=.3
	n=100000
	
	
	bias.mat = matrix(nrow=50, ncol=9)
	for (i in 1:nrow(bias.mat)){
	
				##### create a correlated X and Z
		sig =(matrix(c(1, cor, cor, 1), nrow=2))
		d = data.frame(mvrnorm(n=n, mu=c(0, 0), Sigma=sig))
		names(d) = c("z", "x")
		
				##### create interaction term
		d$zx = d$z*d$x
		
				##### create a y as a function of x,z, and x/z
		vy = 1-(a^2 + b^2 + c^2 + 2*a*b*cor)
		# d$y = a*10*d$x + b*15*d$z + c*12*d$x*d$z + rnorm(length(d$x), 0, sqrt(27254468))
		# var(d$y)
		
				##### now convert standardized variables to non-standard variables
		d$y = a*d$x + b*d$z + c*d$x*d$z + rnorm(length(d$x), 0, vy)		
		d$y = 2.8 + .25*d$y
		d$x = 50 + 10*d$x
		d$z = 100 + 15*d$z
		cor(d)
		
				##### now select on z and create the datasets
		n.selected = which(d$z<median(d$z))		
		full = d
		rest = full; rest[n.selected,c("y", "x", "zx")] = NA
		
				##### use the PL (correcting for all three. This will work)
		require(selection)
		pl = mv.correction(na.omit(rest), p=3, v.pp=cov(full[,1:3]))
		pl
		cov(full)
		
		
				#### test whether expected value comes out right
		mod = lm(x~z, data=rest)
		ve = summary(mod)$sigma^2
		ez = mean(full$z)
		ez2 = mean(full$z^2)
		ez3 = mean(full$z^3)
		b0 = coef(mod)[1]; b1 = coef(mod)[2]
		exxz = ve*ez + b0^2*ez + 2*b0*b1*ez2 + b1^2*ez3
		covxxz = exxz - (b0 +b1*ez)*mean(full$x*full$z)
			###this works
		
		
			#### now see if we can get the other estimate.
		covxzz = var(full$z)*((b0 +b1*ez))+ez*cov(full)[1,2]
			#### this seems to work
			
			#### now get variance of xy
		varxz = mean(full$x)^2*var(full$z) + (mean(full$z)^2 + var(full$z))*var(full$x)


			#### now put them all in one matrix
		xz = data.matrix(mv.correction(na.omit(rest[,1:2]), 1, var(full[,1])))
		xz = cbind(xz,c(covxzz, covxxz))
		xz = rbind(xz, c(covxzz, covxxz, varxz))

			### do pl correction
		pl2 = mv.correction(na.omit(rest), p=3, v.pp=xz)	
		
				### do MI
		# require(mice)
		# mice.impute.prod = function(y, ry, x){return(x["x"]*x["z"])}
		# mi = mice(rest, m=5, method=c("norm", "norm", "prod", "norm"))
		# fit = with(data=mi, exp=cor(cbind(z,x,zx,y))[2,4])
		# mi = mean(unlist(fit$analyses))
	
		
		bias.mat[i,1] = covxzz 
		bias.mat[i,2] = cov(full)[1,3]
		bias.mat[i,3] = covxxz
		bias.mat[i,4] = cov(full)[2,3]
		bias.mat[i,5] = varxz
		bias.mat[i,6] = var(full$zx)
		bias.mat[i,7] = cov2cor(cov(full))[2,4]
		bias.mat[i,8] = cov2cor(pl2)[2,4]
		bias.mat[i,9] = cor(rest, use="pairwise.complete.obs")[2,4]				
	}
	
	options(scipen=1111)
	sqrt(colMeans(bias.mat))

(b0+b1*ez)*(mean(full$x*full$z))
mean(full$z*full$x*full$x)
exxz
cov(full)
mean(full$x)









		###### use Jorge's correction
restricted.data=rest; x="x"; z="z"
jorge.correction = function(restricted.data, x=x, z=z){
		#### 1. create linear model
	f = make.formula(x, z)		
	xzmod = lm(f, data=restricted.data); b0 = coef(xzmod)[1]; b1 = coef(xzmod)[2]

		### 2. compute e(x), e(x^2), e(x^3)
	ez = mean(restricted.data$z); ez2 = mean(restricted.data$z^2); ez3 = mean(restricted.data$z^3)		
	vare = summary(xzmod)$sigma^2	
	
	#### correct z,zx
	covxxz = (vare + b0^2)*ez + 2*b0*b1*ez2 + b1^2*ez3	
	cov(full)
	head(full)
	
	##### correct x,zx
	
	
	
	
cov(d)	
covxxz
cov(full)



		#### now apply case II to correct x/z correlation
	covxz.r = cor(restricted.data[,z], restricted.data[,x], use="pairwise.complete.obs")	
	sdz = sd(restricted.data[,z]); sdzr = sd(na.omit(restricted.data)[,z])
	covxz = caseII(covxz.r, sdz, sdzr)

		#### put into a matrix and input into mv correction
	corrected.mat = c(sdz^2, covxz, covxxz,
					covxz, NA, )

}






		##### apply PL sequentially, starting with z/x
pl = mv.correction(na.omit(rest[,1:2]), p=1, v.pp=1)
pl
cor(rest, use="pairwise.complete.obs")
cor(full)


lm(zx~z + x, data=rest)

var(full$y, na.rm=T)*sqrt(1-summary(lm(y~x*z, data=full))$r.squared)
summary(lm(y~x*z, data=full))$sigma


var(rest$y, na.rm=T)*sqrt(1-summary(lm(y~x*z, data=rest))$r.squared)
summary(lm(y~x*z, data=rest))$sigma





10.16/15 + .5*(1-10.16/15)


summary(lm(y~x*z, data=full))$sigma

		##### now correct
require(selection)
pop = cor(full)[2,4]
c3 = caseIII(data=rest, x=2, y=4, z=1)
em.est = em(data.matrix(rest))[2,4]
pl = cov2cor(mv.correction(na.omit(rest), p=3, v.pp = cov(full)[1:3,1:3]))[2,4]
require(mice)
mice.impute.prod = function(y, ry, x){return(x["x"]*x["z"])}
mi = mice(rest, m=5, method=c("norm", "norm", "prod", "norm"))
fit = with(data=mi, exp=cor(cbind(z,x,zx,y))[2,4])
mi = mean(unlist(fit$analyses))
pop
c3
em.est
pl
mi

		#### see if equation 4 in fife, mendoza, terry actually works
sigq = sd(rest[,"zx"], na.rm=T)
rzq = cor(rest, use="pairwise.complete.obs")[1,3]
varz = var(full$z); varzr = var(rest$z[!is.na(rest$x)])
correction = sigq*(1-rzq^2 + rzq^2*(varz/varzr))
sd(full$zx)
correction

head(rest)






hist(full$y)


			##### when I make x:z missing, I'm deleting the cause of missingness (since z was a part of the selection process). Even though x:z wasn't directly selected on, because x is missing, it too contains missing information not inherent in the dataset
			##### But, if we can impute x, then create a copy of z:x, we can turn it into a MAR situation. Right?
			##### OR, we can correct the r(z,zx) correlation for direct range restriction and subsequently use Case III.
			##### If we get rid of X OR X/Z, we move into an NMAR situation because the cause of missingness in zx is NOT accounted for. The product creates a unique variable that is not a direct function of the missing variable.
			##### So here's what we do: we impute the X's, then we multiply the imputed X's by Z to create xz, then we impute Y. 
			##### OR, we just compute the covariance matrix from the regression weights (since the regression weights are unbiased). 
			##### Put differently, all missing data MUST be explainable from the observed data. If xz values are missing, some missingness is attributable to them because they have a unique relationship to y, above and beyond just z (that is, afterall, what an interaction means). 
			##### Therefore, if xz are missing, we cannot adequetly explain the missingness in Y. 
			##### Actually, can't do from regression weights; too few equations for too many unknowns (only have four beta parameters and 10 elements in v/covar matrix). 
			


mice.impute.norm