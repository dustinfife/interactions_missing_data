clear()

		#### read in parameters for the monte carlo
source('research/interactions/R/0.parameters_setup.R')
require(selection)

		#### set up global parameters
a =demo.slopes[1];
b=demo.slopes[2];
c=demo.int;
cor=demo.cor
n=demo.n
iterations=demo.iterations
sd.z = demo.sd[2]; mu.z = demo.mean[2]
sd.x = demo.sd[1]; mu.x = demo.mean[1]
sd.y = demo.sd[3]; mu.y = demo.mean[3]

i=2
		##### pre-populate matrix
bias.mat = data.frame(iteration=1:iterations, standardize = rep(c(T,F), times=iterations/2), CaseIII=NA, EM = NA, Sample=NA)
for (i in 1:nrow(bias.mat)){

			##### create a correlated X and Z
	sig = matrix(c(sd.z^2, sd.x*sd.z*cor, sd.x*sd.z*cor, sd.x^2), nrow=2)
	d = data.frame(mvrnorm(n=n, mu=c(mu.z, mu.x), Sigma=sig))
	names(d) = c("z", "x")
	mu.xn = mu.x; mu.zn=mu.z
	
			##### create interaction term
	d$zx = d$z*d$x
		
			#### create new weights to retain proper metrix (a' = a*sx/sy)	
	ap = a*(sd.y/sd.x)
	bp = b*(sd.y/sd.z)
	cov = cor*(sd.x*sd.z)		
	var.xz = sd.z^2*mu.xn^2 + sd.x^2*mu.zn^2 + 2*cov*mu.xn*mu.zn + sd.x^2*sd.z^2 + cov^2	
	cp = c*(sd.y/sqrt(var.xz))
	
			##### estimate covariance between X/Z and the interaction term (accordig to Aiken and West)
	covx.xz =sd.x^2*mu.zn + cov*mu.xn
	covz.xz =sd.z^2*mu.xn + cov*mu.zn
	mu.xz = mu.xn*mu.zn + cov

			##### compute expected population correlation
	pop = (ap*sd.x^2 + bp*cov + cp*(covx.xz) )/(sd.x*sd.y)


			##### compute expected variance of y
	vy = sd.y^2-(ap^2*sd.x^2 + bp^2*sd.z^2 + cp^2*var.xz + 2*ap*bp*cov + 2*ap*cp*covx.xz + 2*bp*cp*covz.xz)

			##### compute expected intercept of y
	b0 = mu.y - (ap* mu.xn + bp* mu.zn + cp*mu.xz)

			##### create y
	d$y = b0 + ap*d$x + bp*d$z + cp*d$x*d$z + rnorm(length(d$x), 0, sqrt(vy))	


			##### now select on z and create the datasets
	n.selected = which(d$z<quantile(d$z, .7))		
	full = d
	rest = full; rest[n.selected,c("x", "zx", "y")] = NA
	sample = full[sample(1:nrow(full), size=nrow(na.omit(rest))),]
	
			##### center the variables
	if (bias.mat$standardize[i]){
		d$z = d$z-mean(d$z)
		d$x = d$x-mean(d$x)
		
		rest$z = rest$z-mean(rest$z, na.rm=T)
		rest$x = rest$x-mean(rest$x, na.rm=T)
				
		sample$z = sample$z-mean(sample$z)
		sample$x = sample$x-mean(sample$x)
	}

			#### create the interaction term (based on the new centered variables)
			#### note: it makes a big different in the results of the simulation if we center before generating Y. Centering biases things tremendously. 
	d$xz = d$x*d$z	
	
			##### correct using Pearson Lawley modified correction (I don't report this in the paper)
	pop.cov = matrix(c(sd.z^2, cov, covz.xz,	
						cov, sd.x^2, covx.xz,
						covz.xz, covx.xz, var.xz), nrow=3)
	pl = cov2cor(mv.correction(cov(rest, use="complete.obs"), p=3, v.pp=pop.cov))[2,4]
	bias.mat$dustin[i] = pl-pop

			#### correct using case III/EM
	bias.mat$CaseIII[i] = caseIII(rest, y=4)-pop
	sink(tempfile())
	bias.mat$EM[i] = em(data.matrix(rest))[2,4]-pop
	sink()
	
			#### output the sample estimate
	bias.mat$Sample[i] = cor(sample$x, sample$y)-pop	


			#### see whether interaction terms are unbiased (not reported)
	mod.full = lm(y~x + z + zx, data=d)
	mod.rest = lm(y~x + z + zx, data=rest)	
	bias.mat$interaction.coef[i] = coef(mod.full)[4]-coef(mod.rest)[4]
	
	#### track progress
	if (i/100 == round(i/100)){
		cat(paste0("Iteration ", i, " of ", iterations, "\n"))
	}
	
}

		##### export results
write.csv(bias.mat, "research/interactions/data/demonstration_results.csv", row.names=F)

		#### load tidyverse stuff
require(tidyverse)

		###### create Figure 1 in paper
d = reshape(bias.mat, varying=c("CaseIII", "EM", "Sample"), v.names="Estimate", times=c("CaseIII", "EM", "Sample"), timevar = ("Method"),direction="long")		

		###### group the outliers
d2 <-
  d %>%
  group_by(Method) %>%
  mutate(outlier = Estimate > quantile(Estimate, .75) + IQR(Estimate) * 1.5 | Estimate < quantile(Estimate, .25) - IQR(Estimate) * 1.5) %>%
  ungroup
require(ggplot2)  
		###### produce the plot
theme_set(theme_bw(base_size=14,base_family='Times New Roman'))
p = ggplot(data = d2, mapping=aes(x=Method, y=Estimate))
p + geom_boxplot(outlier.color="lightgray", outlier.shape = NA, width=.5) +  # NO OUTLIERS
  geom_jitter(data = function(x) dplyr::filter_(x, ~ outlier), width=.05, alpha=.15) +
  theme(text=element_text(family="Times")) + geom_hline(yintercept=0, col="lightgray") + labs(x="") + scale_y_continuous(labels=seq(-.3, .5, .1), minor_breaks = seq(-.3, .5, .1))
head(d2)
head(d)
prism.plots(Estimate~Method, data=d)


ggsave("research/interactions/writing/plots/demonstration.pdf")
