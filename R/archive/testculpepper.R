clear()

a =.3;
b=.3;
c=.3;
cor=.3
n=200
iterations=1000


x.2=1:2
harvey76=function(dat,x, x.2, y, tol=.00001){
	Xmat = na.omit(data.matrix(cbind(1, dat[,x])))
	Yvec = na.omit(dat[,y])
	Zmat = na.omit(data.matrix(cbind(1, dat[,x.2])))
	b0 = solve(t(Xmat)%*%Xmat)%*%(t(Xmat)%*%Yvec)
	g0 = c(1,rep(0,ncol(Zmat)-1))
	mdiff=1
	while( mdiff>tol ){
		ezmat = exp(-Zmat%*%g0)
		eXmat = diag(ezmat)*Xmat
		b1 = b0 + (solve(t(eXmat)%*%Xmat))%*%(t(eXmat)%*%(Yvec-Xmat%*%b0))
		g1 = g0 + solve(t(Zmat)%*%Zmat)%*%(t(Zmat)%*%(ezmat*(Yvec-Xmat%*%b1)^2 -1)   )
		mdiff = max(abs(c(b1-b0,g1-g0)))
		b0=b1
		g0=g1
	}
	VCb = solve(t(eXmat)%*%Xmat)
	VCg = 2*solve(t(Zmat)%*%Zmat)
	LL =-2*(nrow(Xmat)/2*log(2*pi)-.5*sum(Zmat%*%g0)-.5*sum(exp(-Zmat%*%g0)*(Yvec-Xmat%*%b0)^2))
	Spop.emp = cov(Xmat)  
	BtSB = t(b0)%*%Spop.emp%*%b0
	#3. Computing Average Prediction Error
	mean_sigma2hat = mean(exp((g0[1]+Zmat*g0[2])))
	#4. Computing Ratio of Predicted to Total Variance
	Vyc = BtSB + mean_sigma2hat
	R2c = BtSB/Vyc
	return(R2c)
}



	i=1
bias.mat = data.frame(rmsr.mi=1:iterations, covxy.mi=NA, covzy.mi=NA, covxzy.mi=NA,	
					rmsr.aiken = NA, covxy.aiken=NA, covzy.aiken=NA, covxzy.aiken=NA,
					rmsr.caseiii = NA, covxy.caseiii=NA, covzy.caseiii=NA, covxzy.caseiii=NA, culpepper=NA, rsq=NA)
for (i in 1:nrow(bias.mat)){

			##### create a correlated X and Z
	sig = matrix(c(1, cor, cor, 1), nrow=2)
	d = data.frame(mvrnorm(n=n, mu=c(0,0), Sigma=sig))
	names(d) = c("z", "x")

			##### create interaction term
	d$zx = d$z*d$x
	
			##### create a y as a function of x,z, and x/z (see proof_y.)
	vy = 1-(a^2 + b^2 + c^2*var(d$zx) + 2*a*b*cor)

	##### now convert standardized variables to non-standard variables
	d$y = a*d$x + b*d$z + c*d$zx + rnorm(length(d$x), 0, sqrt(vy))		

	
			##### now select on z and create the datasets
	n.selected = which(d$z<median(d$z))		
	full = d
	rest = full; rest[n.selected,c("y", "x", "zx")] = NA

			##### use the PL (correcting for all three. This will work)
	require(selection)
	pl = mv.correction(na.omit(rest), p=3, v.pp=cov(full[,1:3]))
	pl
	
	
	#### estimate cov(x,xz) and cov(z,xz) from Cohen, Cohen, West, and Aiken, p. 264
	covx.xz = var(d$x)*mean(d$z) + cov(d$z, d$x)*mean(d$x)
	covz.xz = var(d$z)*mean(d$x) + cov(d$z, d$x)*mean(d$z)
	
	
	##### estimate var(xz) from West and Aiken, page 179
	var.xz = var(d$z)*mean(d$x)^2 + var(d$x)*mean(d$z)^2 + 2*cov(d$z, d$x)*mean(d$x)*mean(d$z) + var(d$x)*var(d$z) + cov(d$x, d$z)^2
	
	##### estimate var/covar of xz (without interaction)
	caseiii = mv.correction(na.omit(rest)[,(1:2)], p=1, v.pp=var(d$z))


	#### put into a variance/covariance matrix (partially corrected)
	corrected = rbind(caseiii, c(covz.xz, covx.xz))
	corrected = cbind(corrected, c(covz.xz, covx.xz, var.xz))
	corrected 

	##### do pearson-lawley
	final.corrected = mv.correction(na.omit(rest), p=3, v.pp=corrected)
	head(rest)
	
		###### do culpepper correction
		two = rest; rest[,2] = rest[,1]^2
	bias.mat[i,"culpepper"] = harvey76(rest, 1:2, 1, 4)
	bias.mat[i,"rsq"] = summary(lm(y~x+z+x*z, data=full))$r.squared

			### do MI	
		#### do my own friggin imputation
		#### 1. predict x from z with error
		#### 2. multiply x.n by z
		#### 3. use pl to fill in y
	f = function(z){x = .3*z + rnorm(length(z), 0, sqrt(1-.3^2))}	
	x = .3*d$z + rnorm(length(d$z), 0, sqrt(1-.3^2))
	xz = x*d$z
	pl.imputed = mv.correction(na.omit(rest), p=3, v.pp=cov(cbind(d$z, x, xz)))
	bias = pl.imputed-cov(d)
	bias.mat[i,1] = sqrt(sum((cov(d) - pl.imputed)^2)/16)
	bias.mat[i,2] = bias[2,4]
	bias.mat[i,3] = bias[1,4]
	bias.mat[i,4] = bias[3,4]
	
	###### record aiken	
	bias.mat[i,5] = sqrt(sum((cov(d) - final.corrected)^2)/16)
	bias = final.corrected-cov(d)
	bias.mat[i,6] = bias[2,4]
	bias.mat[i,7] = bias[1,4]
	bias.mat[i,8] = bias[3,4]	
	
	###### compare to mv correction
	bias.mat[i,9] = sqrt(sum((cov(d) - pl)^2)/16)
	bias = pl-cov(d)
	bias.mat[i,10] = bias[2,4]
	bias.mat[i,11] = bias[1,4]
	bias.mat[i,12] = bias[3,4]				
}
bias.mat = bias.mat[,order(names(bias.mat))]
colMeans(bias.mat)

boxplot(bias.mat$culpepper, bias.mat$rsq)


means = colMeans(bias.mat)[1:9]
sds = apply(bias.mat, 2, sd)[1:9]
par2()
par(mar=c(5,3.5,1,1))
plotCI(1:9, means, sds, xaxt="n", ylab="Bias", xlab="", ylim=c(-.5, .5))
staxlab(1, 1:9, labels=names(bias.mat)[1:9], srt=45)
abline(h=0, col="red")


