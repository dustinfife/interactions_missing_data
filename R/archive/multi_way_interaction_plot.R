clear()
d = read.csv("research/johnson/data/mvDataPostSC2.csv", stringsAsFactors=F)

y = "Bias"
sort.vars = c("c", "rzs", "rzt")
FUN = mean
d = d
highest.first=F
sort.by.effects=T		#### if false, the order inputted in sort.vars is used
add = F

mw.interaction.plot = function(y, sort.vars, FUN, d, highest.first=F, sort.by.effects=T, ...){

			# 1. sort according to size of the effect (using linear models)
	if (sort.by.effects){
		f = make.formula(y, c(sort.vars, "-1"))		## -1 = remove intercept
		effects = lm(f, data=d)
		res = summary(effects)$coefficients[,3]
		ord = names(res)[order(res, decreasing=!highest.first)]
	} else {
		ord = sort.vars
	}

	
			## 2. sort dataset
	ind = which(d[,lines]==line.types[i])
	m = d[ind,]
	m = m[with(m, order(-rzs, -rzt, -c)), ]
	m$sequence = 1:nrow(m)


	par1()
	par(mar=c(5,3.5,1,1), tck=-.01, mgp=c(2, .5, 0))
	plot(1:nrow(m),(m$Bias), type="o", xaxt="n", xlab="", ylab="Bias", cex=.4, col="white", ...)
	#### fastest sequence: c
	r1 = unique(m$rzt)
	r2 = unique(m$rzs)
	a=1;b=1
	for (a in 1:length(r1)){
	for (b in 1:length(r2)){	
		l = subset(m, rzt==r1[a] & rzs == r2[b])
		
		lines(l$sequence, l$Bias, type="o", cex=.8, pch=16, lty=i)	
		if (a==1 & b==1){
			text(l$sequence, l$Bias, labels = paste0("c = ", l$c), cex=0.8, adj=c(-.25, .91))
		}
	}}
	
	#### divide by number of fastest changing one
	rzt.seq = m$rzt[seq(from=1, to=nrow(m), by=length(unique(m$c)))]
	rzt.loc = seq(from=mean(c(1, length(unique(m$c)))), to=nrow(m), by=length(unique(m$c)))
	axis(1, at=rzt.loc, labels=paste0("rzt = ", rzt.seq), las=1, cex.axis=.8)
	
	rzs.seq = m$rzs[seq(from=1, to=nrow(m), by=length(unique(m$c))*length(unique(m$rzt)))]
	rzs.loc = seq(from=mean(c(1, length(unique(m$c))*length(unique(m$rzt)))), to=nrow(m), by=length(unique(m$c))*length(unique(m$rzt)))
	axis(1, at=rzs.loc, labels=paste0("rzs = ", rzs.seq), las=1, cex.axis=.8, line=2, lwd=0, lwd.tick=1)
	
	width.grp.line = rzt.loc[2]-rzt.loc[1]
		#### hack the labels
	for (a in 1:length(rzs.loc)){
		axis(1, at=c(rzs.loc[a]-width.grp.line, rzs.loc[a]+width.grp.line), labels=rep("", 2), las=1, cex.axis=.8, line=2, lwd=1, lwd.tick=0)	
	}

lapply(c(.2, .4, .6), FUN=function(x){expression(paste0(r=x))})	
	
	}
}


x <- sample(1:100, 10, replace = T) # just 10 random numbers
y <- sample(1:100, 10, replace = T) # 10 more random numbers
par(mar = c(10, 5, 5, 5)) 
    # increasing the 1st number to 10 makes 10 lines below axis 1
plot(x~y) # normal plot
axis(1, at = c(20, 40, 60, 80), labels = c("1", "2", "3", "4"), line = 5, col = 4) 



		## 2. re-sort table according to order
require(dplyr)		
res = d[,c(sort.vars, lines, y)]
res %>% arrange_(lines,ord)

		## 3. make an empty plot
plot(1:nrow(res), )		



a <- (1:5)
b <- (6:10)
c <- (11:15)
d <- (16:20)
df <- as.data.frame(cbind(a,b,c,d))
x <- names(df)[1:2]
df %>%
   arrange_(x)

require(dplyr)



effects[with(effects, order(ord)), ]


order(effects, order(abs(res), decreasing=T))

?order
effects = effects[order(effects[,order]),]
		## 2. make empty plot (as long as the number of conditions)
length = unique(d[,sort.vars])
