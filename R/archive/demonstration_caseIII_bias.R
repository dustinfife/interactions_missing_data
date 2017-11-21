clear()
a =.1;
b=.1;
c=.7;
cor=.1
n=100
iterations=1000
sd.z = 100; mu.z = 500
sd.x = .4; mu.x = 3
sd.y = .4; mu.y = 3

bias.mat = data.frame(iteration=1:iterations, CaseIII=NA, EM = NA, Sample=NA)
for (i in 1:nrow(bias.mat)){

			##### create a correlated X and Z
	sig = matrix(c(sd.z^2, sd.x*sd.z*cor, sd.x*sd.z*cor, sd.x^2), nrow=2)
	d = data.frame(mvrnorm(n=n, mu=c(mu.z, mu.x), Sigma=sig))
	names(d) = c("z", "x")

			##### create interaction term
	d$zx = d$z*d$x
	cov(d)

			#### create new weights to retain proper metrix (a' = a*sx/sy)
	cov(d)			
	ap = a*(sd.y/sd.x)
	bp = b*(sd.y/sd.z)
	cp = c*(sd.y/sd(d$zx))
	cov = cor*(sd.x*sd.z)

			##### compute variance of y
	vy = sd.y^2-(ap^2*sd.x^2 + bp^2*sd.z^2 + cp^2*var(d$zx) + 2*ap*bp*cov + 2*ap*cp*cov(d$x, d$zx) + 2*bp*cp*cov(d$z, d$zx))
	
			##### compute intercept of y
	b0 = mu.y - (ap*mu.x + bp*mu.z + cp*mean(d$zx))
	
			##### create y
	d$y = b0 + ap*d$x + bp*d$z + cp*d$zx + rnorm(length(d$x), 0, sqrt(vy))	
	mean(d$y)
	sd(d$y)



cov(d)
#plot(d$x, d$y)
			##### now select on z and create the datasets
	n.selected = which(d$z<median(d$z))		
	full = d
	rest = full; rest[n.selected,c("y", "x", "zx")] = NA
	sample = full[sample(1:nrow(full), size=nrow(na.omit(rest))),]
	
			##### correct using Case III
	require(selection)
	pl = data.frame(mv.correction(cov(rest[,1:2], use="pairwise.complete.obs"), p=1, v.pp=matrix(var(d$z), nrow=1)))
	names(pl) = c("SAT", "GPA")
	row.names(pl) = names(pl)


	
	bias.mat$CaseIII[i] = caseIII(rest, y=4)
	bias.mat$Sample[i] = cor(sample$x, sample$y)
	sink(tempfile())
	bias.mat$EM[i] = em(data.matrix(rest))[2,4]
	sink()

		#### track progress
	if (i/100 == round(i/100)){
		cat(paste0("Iteration ", i, " of ", iterations, "\n"))
	}
	
}

head(bias.mat)
mns = apply(bias.mat, 2, mean)
mns
mns[2] - mns[4]
pdf("research/caseIII/writing/plots/demonstration.pdf")
par1()
boxplot(bias.mat$CaseIII, bias.mat$EM, bias.mat$Sample, xaxt="n", ylab="Estimated Correlation Between X and Y")
axis(1, at=1:3, labels=c("Case III Corrected", "EM Corrected", "Random Sample"))
dev.off()



a =.3;
b=.3;
c=.3;
cor=.3
n=2000
iterations = 1
i=1

			##### create a correlated X and Z
	sig = matrix(c(1, cor, cor, 1), nrow=2)
	d = data.frame(mvrnorm(n=n, mu=c(0,0), Sigma=sig))
	names(d) = c("z", "x")
	
			##### create interaction term
	d$zx = d$z*d$x

	
			##### create a y as a function of x,z, and x/z
	vy = 1-(a^2 + b^2 + c^2*var(d$zx) + 2*a*b*cor + 2*a*c*cov(d$x, d$zx) + 2*b*c*cov(d$z, d$zx))
	d$y = a*d$x + b*d$z + c*d$zx + rnorm(length(d$x), 0, sqrt(vy))

			##### convert all to sensible metrics
	d$z = rescale(d$z, 500, 100)
	d$x = rescale(d$x, 3, .4)
	d$y = rescale(d$y, 3, .4)
	d$zx = d$x*d$x		
			##### now select on z and create the datasets
	n.selected = which(d$z<median(d$z))		
	full = d
	rest = full; rest[n.selected,c("y", "x", "zx")] = NA
	sample = full[sample(1:nrow(full), size=nrow(na.omit(rest))),]
	
		#### plot that sucka
pdf("research/caseIII/writing/plots/explanation.pdf")
	plot(d$x, d$y, col=string.to.colors(is.na(rest$y), colors=c(rgb(0, 53, 149, 40, maxColorValue=255), rgb(255,140,0,40, maxColorValue=255))), pch=16, xlab="X", ylab="Y", cex=.8)		
	abline(lm(y~x, data=d[n.selected,]), col=rgb(0,53,149, maxColorValue=255))
	abline(lm(y~x, data=d[which(d$z>median(d$z)),]), col=rgb(255, 140, 0, maxColorValue=255))
	lines(lowess(d$x, d$y), col='black', lwd=2)
dev.off()

string.to.colors(is.na(rest$y))[2]