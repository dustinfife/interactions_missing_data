clear()

#### see what iteration we're at
load("research/caseiii/data/last_done.Rdat")
num = mc_caseiii[mc_caseiii[,2]==1,1]
num = num[length(num)] 

#### see how many are left
max.num = subsetString(list.files("research/caseiii/data", pattern="monte_carlo_caseiii_"), sep="_", position=4)
max.num = max(as.numeric(gsub(".csv", "", max.num)))

require(selection)
##### now iterate
for (l in num:(max.num)){
	source("research/caseiii/R/2b. monte_carlo_run_c3.R")
	cat(paste0("Iteration ", l, " of ", max.num, "\n"))
}