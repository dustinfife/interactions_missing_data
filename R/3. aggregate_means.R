clear()
files = list.files("research/caseiii/data/", pattern="monte_carlo_caseiii_", full.names=T)
files = gsub("//", "/", files)

		#### read one initially to get things started
d = read.csv(files[1])

			#### find out which variables are varying and put into formulas
varying = which(apply(d, 2, sd)>0)
varying = names(varying[-which(names(varying) %in% c("i", "Sample", "corrected", "uncorrected", "CaseIII"))])

sample = make.formula("Sample", c(varying))
caseiii = make.formula("CaseIII", c(varying))
corrected = make.formula("corrected", c(varying))
uncorrected = make.formula("uncorrected", c(varying))

			##### get size of first 2 dimensions
res = aggregate(sample, FUN=mean, data=d)

		#### preallocate three arrays
results_sample = array(dim=c(dim(res), length(files)))
results_caseiii = results_sample
results_corrected = results_sample
results_uncorrected = results_sample
sd_sample =sd_caseiii=sd_corrected=sd_uncorrected = results_sample

for (i in 1:length(files)){
	d = read.csv(files[i])
	b = data.matrix(aggregate(sample, FUN=mean, data=d, na.rm=F))
	
	### record means
	results_sample[,,i] = b
	results_caseiii[,,i] = data.matrix(aggregate(caseiii, FUN=mean, data=d, na.rm=T))
	results_corrected[,,i] = data.matrix(aggregate(corrected, FUN=mean, data=d, na.rm=T))
	results_uncorrected[,,i] = data.matrix(aggregate(uncorrected, FUN=mean, data=d, na.rm=T))	
	
	### record standard errors
	sd_sample[,,i] = data.matrix(aggregate(sample, FUN=sd, data=d, na.rm=T))
	sd_caseiii[,,i] = data.matrix(aggregate(caseiii, FUN=sd, data=d, na.rm=T))
	sd_corrected[,,i] = data.matrix(aggregate(corrected, FUN=sd, data=d, na.rm=T))		
	sd_uncorrected[,,i] = data.matrix(aggregate(uncorrected, FUN=sd, data=d, na.rm=T))			
	print(i)
}


		#### aggregate across 3rd dimension
sample_results = data.frame(apply(results_sample, c(1,2), FUN=mean))
caseiii_results = data.frame(apply(results_caseiii, c(1,2), FUN=mean))
results_corrected = data.frame(apply(results_corrected, c(1,2), FUN=mean))
results_uncorrected = data.frame(apply(results_uncorrected, c(1,2), FUN=mean))
names(sample_results) = names(caseiii_results) = names(results_corrected) = names(results_uncorrected) = c(varying, "Estimate")
names(results_corrected)[5] = "p"
names(caseiii_results)[5] = "p"


sd_sample = data.frame(apply(sd_sample, c(1,2), FUN=mean))
sd_caseiii = data.frame(apply(sd_caseiii, c(1,2), FUN=mean))
sd_corrected = data.frame(apply(sd_corrected, c(1,2), FUN=mean))
sd_uncorrected = data.frame(apply(sd_uncorrected, c(1,2), FUN=mean))
names(sd_sample) = names(sd_caseiii) = names(sd_corrected) = names(sd_uncorrected) = c(varying, "SD")
names(results_corrected)[5] = "p"
names(caseiii_results)[5] = "p"

		#### export the results
save(sample_results, caseiii_results, results_corrected, results_uncorrected,
	sd_sample, sd_caseiii, sd_corrected, sd_uncorrected,
	file="research/caseiii/data/MC_final_results.Rdat")


