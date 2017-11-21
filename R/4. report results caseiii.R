clear()
load("research/caseiii/data/MC_final_results.Rdat")
		#### subtract caseiii and corrected from random sample
caseiii_results$Estimate = caseiii_results$Estimate-sample_results$Estimate
results_corrected$Estimate = results_corrected$Estimate-sample_results$Estimate
results_uncorrected$Estimate = results_uncorrected$Estimate-sample_results$Estimate

source("research/caseiii/R/multi_way_interaction_plot_gen.R")


names(results_corrected)[2] = "r"
names(results_corrected)[4] = "s"
names(caseiii_results) = names(results_uncorrected) = names(results_corrected)
pdf("research/caseiii/writing/knitr/bias_correction.pdf")
results_corrected
m = mw.interaction.plot(y="Estimate", sort.vars=c("s","r", "b"), smallest.to.biggest=T, FUN=mean, d=results_corrected, add=F, col=rgb(0,0,0,.25), ylim=c(-.2, .3), line.loc=c(0,3.25), label.points=T, sort.by.effects=F, pch=16)
abline(h=0, col="lightgray")
m = mw.interaction.plot(y="Estimate", sort.vars=c("s","r", "b"), smallest.to.biggest=T, FUN=mean, d=caseiii_results, add=T, label.points=F, sort.by.effects=F, pch=1)
#m = mw.interaction.plot(y="Estimate", sort.vars=c("s","r", "b"), FUN=mean, d=results_uncorrected, add=T, label.points=F, sort.by.effects=F, pch=2, col="lightgray")
legend("bottomleft", legend=c(expression(italic('r'[c3])), expression(italic('r'[pl]))), col=c("black", rgb(0,0,0,.25)), pch=c(1,16), bty="n")
dev.off()


names(sd_corrected)[5] = names(sd_caseiii)[5] = names(sd_sample)[5] = "p"
names(sd_corrected)[4] = names(sd_caseiii)[4] = names(sd_sample)[4] = "s"
pdf("research/caseiii/writing/knitr/se_correction.pdf")

m = mw.interaction.plot(y="SD", sort.vars=c("p","n", "s"), FUN=mean, d=sd_corrected, add=F, col=rgb(0,0,0,.25),ylim=c(0, .6), line.loc=c(0,3.25), label.points=F, sort.by.effects=F, pch=16)
text(1:5, m$SD[1:5]+c(-.02, 0, 0, -.02,.01), paste0("s=", unique(m$s)))
m = mw.interaction.plot(y="SD", sort.vars=c("p","n", "s"), FUN=mean, d=sd_caseiii,  add=T, label.points=F, sort.by.effects=F, pch=1)

abline(h=0, col="lightgray")

legend("topright", legend=c(expression(italic('r'[c3])), expression(italic('r'[pl]))), col=c("black", rgb(0,0,0,.25)), pch=c(1,16), bty="n")
dev.off()




require(Hmisc)
colored.table(results_corrected, "research/caseiii/writing/knitr/corrected.Rnw", dep.var="Estimate", row.factors=c("skew","cor"), col.factors=c("b"), landscape=T)	
colored.table(caseiii_results, "research/caseiii/writing/knitr/caseiii.Rnw", dep.var="Estimate", row.factors=c("skew","cor"), col.factors=c("b"), landscape=T)	

