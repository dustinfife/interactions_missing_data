clear()
set.seed(732017)
		###the final parameters to be estimated are mu.z, c, p.missing, cor, and skew
		##### read in MC parameters from knitr file
params = read.csv("research/interactions/data/mc_parameters.csv")
fix.func = function(x){as.numeric(unlist(strsplit(as.character(x), ", ")))}
a = fix.func(params[1,2])
b= fix.func(params[2,2])
c= fix.func(params[3,2])
cor= fix.func(params[4,2])
n= fix.func(params[7,2])
skew = fix.func(params[9,2])
iterations=500
sd.z = 1; mu.z = fix.func(params[6,2])
sd.x = 1; mu.x = fix.func(params[5,2])
sd.y = 1; mu.y = 0

		##### preallocate
bias.mat = expand.grid(a=a,b=b,c=c,cor=cor,n=n,mu.z=mu.z, mu.x=mu.x, skew=skew, p.missing=fix.func(params[8,2]), i=1:iterations, Sample=NA, caseIII=NA, corrected=NA, uncorrected=NA)
nrow(bias.mat)

		##### break up into several chunks (of 10,000 rows each)
breakups = seq(from=1, to=nrow(bias.mat), by=(nrow(bias.mat)/iterations)*10)		

		##### loop through and create a separate dataset for each
for (i in 1:length(breakups)){
	if (i<length(breakups)){
		new.d = bias.mat[breakups[i]:(breakups[i+1]-1),]
	} else {
		new.d = bias.mat[breakups[i]:nrow(bias.mat),]
	}
	write.csv(new.d, paste0("research/interactions/data/monte_carlo_caseiii_", i, ".csv"), row.names=F)
}

mc_caseiii = data.frame(last.done=0, completed=1)
save(mc_caseiii, file="research/interactions/data/last_done.Rdat")
