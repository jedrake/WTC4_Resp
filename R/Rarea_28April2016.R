#- read in the leaf respiration data recorded at 15 deg C for WTC4
source("R/load.R")


#----------------------------------------------------------------
#- get the Rleaf15 data
Rleaf <- returnRleaf()

#- update the trendlogs based on the Rleaf date measurements.
#  Note that this can take awhile.
updateTrendlogsRleaf(measuredates = unique(dat.m$Date))
                   
#-- loop over each observation in the Rleaf dataset, calculate the mean temperature of the preceding X dates (3)
ndays <- 3
Rleaf$Tpre <- NA
for (i in 1:nrow(Rleaf)){
  searchdate <- Rleaf$Date[i]
  mindate <- searchdate-ndays
  inds <- which(TrendlogChDF$Date >=mindate & TrendlogChDF$Date <=searchdate & TrendlogChDF$chamber == Rleaf$chamber[i])
  Rleaf$Tpre[i] <- mean(TrendlogChDF[inds,"Tair_al"],na.rm=T)
  
}  
  
#----------------------------------------------------------------



#----------------------------------------------------------------
#- average across chambers
dat.m <- summaryBy(.~chamber+Date+T_treatment,data=Rleaf,FUN=mean,keep.names=T)

windows(20,40);par(mar=c(5,6,1,1),cex.lab=2)
boxplot(Rarea~T_treatment,data=dat.m,ylim=c(0,1.5),ylab="Rarea (umol m-2 s-1)",col=c("grey","red"))

#----------------------------------------------------------------


#----------------------------------------------------------------
#- get Mike's data from WTC3. Takes a few moments to read and interpolate all the met data
mja <- returnMJAdata(ndays=3)
#----------------------------------------------------------------



#----------------------------------------------------------------
#- plot treatment averages and SE's of Rleaf15 vs. thermal history
Rleaf.m <- summaryBy(Rarea+Tpre~T_treatment+Date,data=dat.m,FUN=c(mean,standard.error))
mja.m <- summaryBy(Rarea25+Rarea15+Tpre~T_treatment+Date,data=mja,FUN=c(mean,standard.error))


windows(25,35);par(mar=c(8,6,1,1),cex.lab=1.7,cex.axis=1.1)
size=1.2
xlims=c(5,30)
ylims=c(0,1.5)
#- plot Mike's data
plotBy(Rarea15.mean~Tpre.mean|T_treatment,data=mja.m,legend=F,pch=1,xaxt="n",yaxt="n",type="p",
       cex=size,ylim=ylims,xlim=xlims,ylab="",xlab="",
       panel.first=adderrorbars(x=mja.m$Tpre.mean,y=mja.m$Rarea15.mean,SE=mja.m$Rarea15.standard.error,direction="updown"))
magaxis(side=c(1,2,3,4),labels=c(1,1,0,0),las=1)

#- add my data
plotBy(Rarea.mean~Tpre.mean|T_treatment,data=Rleaf.m,legend=F,pch=16,xaxt="n",yaxt="n",type="p",add=T,
       cex=size,ylim=ylims,xlim=xlims,ylab="",xlab="",
       panel.first=adderrorbars(x=Rleaf.m$Tpre.mean,y=Rleaf.m$Rarea.mean,SE=Rleaf.m$Rarea.standard.error,direction="updown"))
title(xlab=expression(atop(Thermal~history,
                             T[air]~of~preceding~3~days~(degree*C))),line=5)
title(ylab=expression(R[leaf]~at~15~degree*C~(mu*mol~CO[2]~m^-2~s^-1)))

legend("topright",c("A-WTC3","W-WTC4","A-WTC4","W-WTC4"),ncol=2,col=c("black","red","black","red"),
       pch=c(1,1,16,16))
#----------------------------------------------------------------


