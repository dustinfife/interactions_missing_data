#Function for Estimating Parameters of Harvey (1976)
  #Xmat = Design matrix for modeling conditional mean
  #Zmat = Design matrix for modeling conditional error variance
  #Yvec = Criterion Variable
  #tol  = tolerance level for convergence 
harvey76=function(Xmat,Zmat,Yvec,tol=.00001){
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
  
  list(b0,g0,LL)
}

#R Code for Generating Empirical Example
set.seed(20150215)
N=5000 #sample size
dat = matrix(,N,2)
dat[,1] = rnorm(N)
cmax = max(dat[,1])
cmin = min(dat[,1])
b1 = .5
Gs = c(0,b1,-b1/(2*cmax))
gs = c(0,-.25)
dat[,2] = rnorm(N,mean=Gs[1]+dat[,1]*Gs[2]+dat[,1]^2*Gs[3],
                sd=exp((gs[1]+dat[,1]*gs[2])/2) )
colnames(dat) = c('X','Y')
head(dat)
summary(lm(Y~X, data=data.frame(dat)))
#introducing truncation at x = q(1-.8)
pis = which(dat[,1]>quantile(dat[,1],probs=.80))
dats = dat[pis,] 
#True population correlation
rho2 = summary(lm(dat[,'Y']~poly(dat[,'X'],2,raw=T)))$r.square

dat_supp = dat
dat_supp = as.data.frame(cbind(dat,matrix(,N,1)))
head(dat_supp)
dat_supp[pis,3] = dat[pis,2]
colnames(dat_supp)=c('X','Y','Ys')
EmpDat = dat_supp
#Standard Correction for Range Restriction
Vc = sd(EmpDat[,1])/sd(EmpDat[!is.na(EmpDat[,3]),1])
rs = cor(EmpDat[,1],EmpDat[,3],use="complete.obs")
R2linc = (Vc*rs/sqrt(1+rs^2*(Vc^2-1)) )^2
sqrt(R2linc) #traditional correction based upon Equation 1

#1. Computing Unrestricted Predictor Variance-Covariance Matrix 
Xpop = cbind(1,EmpDat[,1],EmpDat[,1]^2)
head(Xpop)
Spop.emp = cov(Xpop)
#2. Estimating Harvey (1976) Parameters
fit_h = harvey76(Xmat=cbind(1,EmpDat[!is.na(EmpDat[,3]),1],EmpDat[!is.na(EmpDat[,3]),1]^2),		### non-missing design matrix data
                 Zmat=cbind(1,EmpDat[!is.na(EmpDat[,3]),1]),								### non-missing design matrix for all but x^2
                 Yvec=EmpDat[!is.na(EmpDat[,3]),3],tol=.00001)								### outcome variable
BtSB = t(fit_h[[1]])%*%Spop.emp%*%fit_h[[1]] 
#3. Computing Average Prediction Error
head(EmpDat)
mean_sigma2hat = mean(exp((fit_h[[2]][1]+EmpDat[,1]*fit_h[[2]][2])))
#4. Computing Ratio of Predicted to Total Variance
Vyc = BtSB + mean_sigma2hat
R2c = BtSB/Vyc 
(R2c)
head(EmpDat)


x = rnorm(5000)
y = .5*x + .2*x^2 + rnorm(length(x))


x = 1
dat[,3] = dat[,1]^2

	### dat = a dataset containing all the data (including nonlinear terms such as x*x)
	### x = the indices (for dat) that indicate which columns are included as predictors (should include x*x)
	### y = the index for the y column
harvey76=function(dat,x, y, tol=.00001){
	Xmat = cbind(1, dat[,x])
	Yvec = dat[,y]
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
  
  list(b0,g0,LL)
}































#