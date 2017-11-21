# clear()
# d = read.csv("research/johnson/data/mvDataPostSC2.csv", stringsAsFactors=F)
# d=conditions
# y="seh"
# sort.vars=c("condom","preg","Gender","Condition")
# FUN = mean
# highest.first=F
# sort.by.effects=T		#### if false, the order inputted in sort.vars is used
# add = F
# smallest.to.biggest=F
# label.pos=c(-.25, .91)
# col.lines="black"
# line.loc=seq(from=0, to=length(sort.vars)-1)
# line.loc=c(0,6,9)
# mar=c(12,3.5,1,1)
mw.interaction.plot = function(y, sort.vars, FUN, d, highest.first=F, sort.by.effects=T, smallest.to.biggest=F, add=F, return.table=T, line.loc=seq(from=0, to=length(sort.vars)-1),
	col.lines="black", label.points=T, label.pos=c(-.25, .91), mar=c(7,3.5,1,1), spec.label=NULL, ...){

	
	if (length(sort.vars)<=2){
		quit("sort.vars must have at least three elements. If you're only using two, try using a two way interaction plot instead.")
	}

		#### reduce dataset to aggregated
	f = make.formula(y, c(sort.vars, "-1"))		## -1 = remove intercept
	d = aggregate(f, data=d, FUN=FUN)
			# 1. sort according to size of the effect (using linear models)
	if (sort.by.effects){
		effects = lm(f, data=d)
		res = anova(effects)[,"F value"]
		ord = na.omit(sort.vars[order(res, decreasing=!highest.first)])
	} else {
		ord = sort.vars
	}

	
			## 2. sort dataset			
	m = d[do.call("order", c(d[ord], decreasing=!(smallest.to.biggest))), ]
	m$sequence = 1:nrow(m)

	if (!add){

		#### make empty plot
		par1()
		par(mar=mar, tck=-.01, mgp=c(2, .5, 0))
		plot(1:nrow(m),(m[,y]), type="o", xaxt="n", xlab="", ylab=y, cex=.4, col="white",...)

	}
		#### loop through all to draw lines
		total.lines = seq(from=1, to=nrow(m), by=length(unique(m[,ord[length(ord)]])))		
		for (i in 1:length(total.lines)){
			if (i==length(total.lines)){
				nums = total.lines[i]:nrow(m)
			} else {
				nums = total.lines[i]:(total.lines[i+1]-1)
			}
			lines(nums, m[nums,y], type="o", cex=.8, col=col.lines,...)
			
			### if the first iteration, label the points
			if (i==1 & label.points){
				text(nums, m[nums,y], labels = paste0(ord[length(ord)], " = ", unique(m[,ord[length(ord)]])), cex=0.8, adj=label.pos)
			}
			
		}


		#### now label the axes
		remaining = rev(ord[-length(ord)])
		current = length(unique(m[,ord[length(ord)]]))
		if (!add){
		for (i in 1:length(remaining)){
			
			if (i>1){	
				width.grp.line = new.loc[2]-new.loc[1] + .75*new.loc[2]-new.loc[1]		### this goes here before I replace new.loc (because I need it LAST time, not this time)
			}
			tick = ifelse(i==1, .1,0)
			new.seq = m[,remaining[i]][seq(from=1, to=nrow(m), by=current)]
			new.loc = seq(from=mean(c(1, current)), to=nrow(m), by=current)
			current2 = length(unique(m[,remaining[i]]))
			current = current*current2
			labs = paste0(remaining[i], " = ", new.seq)
			axis(1, at=new.loc, labels=labs, las=1, cex.axis=.8, line=line.loc[i], lwd=0, lwd.tick=tick, las=2)

						#### figure out how much "padding" there should be between axis labels
			#max_nchar <- max(nchar(labs)) *2.5

		
				#### now hack the lines
			if (i>1){
				for (a in 1:length(new.loc)){
					axis(1, at=c(new.loc[a]-width.grp.line, new.loc[a]+width.grp.line), labels=rep("", 2), las=1, cex.axis=.8, line=line.loc[i], lwd=1, lwd.tick=0)	
				}			
			}
		}}


	if (return.table){
		return(m)
	}
	
}


