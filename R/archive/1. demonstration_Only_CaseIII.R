clear()
source('research/interactions/R/0.parameters_setup.R')
a =demo.slopes;
b=demo.slopes;
c=demo.slopes;
cor=demo.cor
n=demo.n
iterations=demo.iterations

sd.z = demo.sd[2]; mu.z = demo.mean[2]
sd.x = demo.sd[1]; mu.x = demo.mean[1]
sd.y = demo.sd[3]; mu.y = demo.mean[3]

bias.mat = data.frame(iteration=1:iterations, standardize = rep(c(T,F), times=iterations/2), CaseIII=NA, EM = NA, Sample=NA)
i=2
for (i in 1:nrow(bias.mat)){

			##### create a correlated X and Z

	sig = matrix(c(sd.z^2, sd.x*sd.z*cor, sd.x*sd.z*cor, sd.x^2), nrow=2)
	d = data.frame(mvrnorm(n=n, mu=c(mu.z, mu.x), Sigma=sig))
	names(d) = c("z", "x")
	
			##### standardize?
	if (bias.mat$standardize[i]){
		d$z = d$z-mean(d$z)
		d$x = d$x-mean(d$x)
		mu.xn = 0; mu.zn = 0
	} else {
		mu.xn = mu.x; mu.zn=mu.z
	}
	

			##### create interaction term
	d$zx = d$z*d$x


		
			#### create new weights to retain proper metrix (a' = a*sx/sy)	
	ap = a*(sd.y/sd.x)
	bp = b*(sd.y/sd.z)
	cov = cor*(sd.x*sd.z)		
	var.xz = sd.z^2*mu.xn^2 + sd.x^2*mu.zn^2 + 2*cov*mu.xn*mu.zn + sd.x^2*sd.z^2 + cov^2	
	cp = c*(sd.y/sqrt(var.xz))



	covx.xz =sd.x^2*mu.zn + cov*mu.xn
	covz.xz =sd.z^2*mu.xn + cov*mu.zn
	mu.xz = mu.xn*mu.zn + cov



	pop = (ap*sd.x^2 + bp*cov + cp*(covx.xz) )/(sd.x*sd.y)


			##### compute variance of y
	vy = sd.y^2-(ap^2*sd.x^2 + bp^2*sd.z^2 + cp^2*var(d$zx) + 2*ap*bp*cov + 2*ap*cp*cov(d$x, d$zx) + 2*bp*cp*cov(d$z, d$zx))

			##### compute intercept of y
	b0 = mu.y - (ap* mu.xn + bp* mu.zn + cp*mu.xz)


			##### create y
	d$y = b0 + ap*d$x + bp*d$z + cp*d$zx + rnorm(length(d$x), 0, sqrt(vy))	

			##### now select on z and create the datasets
	n.selected = which(d$z<mu.zn)		
	full = d
	rest = full; rest[n.selected,c("x", "zx", "y")] = NA
	sample = full[sample(1:nrow(full), size=nrow(na.omit(rest))),]
	
			##### correct using Case III
	require(selection)
	pl = data.frame(mv.correction(cov(rest[,1:2], use="pairwise.complete.obs"), p=1, v.pp=matrix(var(d$z), nrow=1)))
	names(pl) = c("SAT", "GPA")
	row.names(pl) = names(pl)

	pop.cov = matrix(c(sd.z^2, cov, covz.xz,	
						cov, sd.x^2, covx.xz,
						covz.xz, covx.xz, var.xz), nrow=3)
	pl = mv.correction(cov(rest, use="complete.obs"), p=3, v.pp=pop.cov)						
	bias.mat$CaseIII[i] = caseIII(rest, y=4)-pop
	bias.mat$Sample[i] = cor(sample$x, sample$y)-pop

	sink(tempfile())
	bias.mat$EM[i] = em(data.matrix(rest))[2,4]-pop
	sink()

		#### track progress
	if (i/100 == round(i/100)){
		cat(paste0("Iteration ", i, " of ", iterations, "\n"))
	}
	
}

write.csv(bias.mat, "research/interactions/data/demonstration_results.csv", row.names=F)
aggregate(CaseIII~standardize, FUN=mean, data=bias.mat)
d = reshape(bias.mat, varying=c("CaseIII", "EM", "Sample"), v.names="Estimate", times=c("CaseIII", "EM", "Sample"), timevar = ("Method"),direction="long")		
	
pdf("research/interactions/writing/plots/demonstration.pdf")		
par1()		
#par(mfrow=c(3,2), oma=c(0,0,0,1.5))

boxplot(Estimate~Method, data=d[d$standardize,], at=c(1,4,7), xlim=c(0,9), xaxt="n", ylim=range(d$Estimate), ylab="Bias")
boxplot(Estimate~Method, data=d[!(d$standardize),], at=c(2, 5, 8), add=T, xaxt="n", col="lightgray")
axis(1, c(1.5, 4.5, 7.8), labels=c("Case III", "EM", "Random Sample"))
abline(h=0)

dev.off()

mns = colMeans(bias.mat[!(bias.mat$standardize),])
mns[3] - mns[5]

options(scipen=10)
mns = apply(bias.mat[!(bias.mat$standardize),], 2, sd)
mns = apply(bias.mat[(bias.mat$standardize),], 2, sd)
mns



# # pdf("research/interactions/writing/plots/demonstration.pdf")
# par1()
# boxplot(bias.mat$CaseIII, bias.mat$EM, bias.mat$Sample, xaxt="n", ylab="Estimated Correlation Between X and Y")
# axis(1, at=1:3, labels=c("Case III Corrected", "EM Corrected", "Random Sample"))
# dev.off()
