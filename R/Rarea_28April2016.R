#- read in the leaf respiration data recorded at 15 deg C for WTC4
source("R/load.R")


#----------------------------------------------------------------
#- get the Rleaf15 data
Rleaf <- returnRleaf()
#----------------------------------------------------------------



#----------------------------------------------------------------
#- average across chambers
dat.m <- summaryBy(.~chamber+Date+T_treatment,data=Rleaf,FUN=mean,keep.names=T)

windows(20,40);par(mar=c(5,6,1,1),cex.lab=2)
boxplot(Rarea25~T_treatment,data=dat.m,ylim=c(0,3),ylab="Rarea (umol m-2 s-1)",col=c("grey","red"))

#----------------------------------------------------------------


#----------------------------------------------------------------
#- get Mike's data from WTC3. Takes a few moments to read and interpolate all the met data
mja <- returnMJAdata(ndays=3)
#----------------------------------------------------------------
