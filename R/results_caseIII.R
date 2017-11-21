files = list.files("research/caseiii/data/", pattern="monte_carlo_caseiii_")






clear()
require(party)
require(edarf)
d = read.csv("research/caseiii/data/runif_simulation_caseiii.csv")
fit <- cforest(corrected~a+b+c+cor+skew+p.missing+n+mu.z+mu.x, controls = cforest_unbiased(mtry = 2), data=d)
imp = varimp(fit)
options(scipen=10)
sort(imp)
names(fit)

# pd_int <- partial_dependence(fit, c("skew", "cor", "b"), n = c(10, 25), interaction = TRUE)
# plot_pd(pd_int)






?cforest
imp <- variable_importance(fit, nperm = 10)
plot_imp(imp)






thresh = rfThresh(corrected~a+b+c+cor+skew+p.missing+n+mu.z+mu.x, data=d[1:1000,], nruns=25, importance="gini")

require(edarf)
mod = thresh$model
pd <- partial_dependence(mod, vars = "n", n = c(10, 25), data=d)
plot_pd(pd)


require(party)






vi = varimp(fit)
vi
pd <- partial_dependence(fit, vars = "b", n = c(10, 25))
plot_pd(pd)
? partial_dependence

thresh
plot(thresh)
rf = VSURF(, data=d)