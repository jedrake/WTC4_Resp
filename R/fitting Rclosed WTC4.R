#- read in the minutely whole-canopy respiration data for WTC4
#  measured on the night of 28 April 2016 to 29 April 2016
setwd("Work/HFE/WTC4/")
library(doBy)
library(HIEv)
library(plotBy)
library(RColorBrewer)
library(DEoptim)
library(calibrate)
setToken(tokenfile="C://Repos/wtc3_flux/HIEv_token.txt")
palette(c("black",rev(brewer.pal(11,"Spectral"))))









#- download the data
downloadHIEv(searchHIEv("WTC_AUTO_ALL_REF01MIN_R_20160430.dat"),topath="Data/WCR/")
downloadHIEv(searchHIEv("WTC_AUTO_ALL_CH01MIN_R_20160430.dat"),topath="Data/WCR/")

#- read in the data
refdat <- readTOA5("Data/WCR/WTC_AUTO_ALL_REF01MIN_R_20160430.dat")
chdat <- readTOA5("Data/WCR/WTC_AUTO_ALL_CH01MIN_R_20160430.dat")

#- subset to just my measurements
starttime <- as.POSIXct("2016-04-28 19:00:00", tz="UTC")
starttime2 <- as.POSIXct("2016-04-28 22:30:00", tz="UTC")
starttime3 <- as.POSIXct("2016-04-28 22:48:00", tz="UTC")
endtime <- as.POSIXct("2016-04-29 01:05:00",tz="UTC")

refdat <- subset(refdat,DateTime>starttime & DateTime < endtime)
refdat$DateTime <- nearestTimeStep(refdat$DateTime,nminutes=1,align="floor")
chdat <- subset(chdat,DateTime>starttime & DateTime < endtime)

#- just subset the first set of measuremens
#----------------------------------------------------------------------------------------------
#model the flux and the leak rate, return the squared residual sum for minimiation by DEoptim()
closedWTC.mod <- function(mu, V, Ca, Ci,fit=1,breakpoint){
  
  R <- mu[1]
  theta <- mu[2]
  
  
  resid <- rep(NA,length(Ca))
  resid[1] <- 0
  pred <- rep(NA,length(Ca))
  pred[1] <- Ci[1]
  
  #- fit the first bit (before the break point)
  for (i in 2:breakpoint){
    dCO2 <- (-theta*(Ci[i-1]-Ca[i-1])+R)/V
    pred[i] <- pred[i-1] + dCO2
    
    resid[i] <- (Ci[i] - pred[i])^2
  }
  
  #- define the breakpoint values so the second loop will run
  pred[breakpoint+1] <- Ci[breakpoint+1]
  resid[breakpoint+1] <- 0
  
  #- fit the second bit (after the break point)
  for (i in (breakpoint+2):length(Ca)){
    dCO2 <- (-theta*(Ci[i-1]-Ca[i-1])+R)/V
    pred[i] <- pred[i-1] + dCO2
    
    resid[i] <- (Ci[i] - pred[i])^2
  }
  resid.sum <- sum(resid)
  if (fit==1) return(resid.sum)
  if (fit==0) return(data.frame(Ci=Ci,Ca=Ca,pred=pred))
}
#----------------------------------------------------------------------------------------------

chdat1 <- subset(chdat,DateTime<starttime2)
chdat2 <- subset(chdat,DateTime>starttime3)

plotBy(CO2L~DateTime|chamber,data=chdat2,type="o",legend=F)


chdat1.l <- split(chdat1,chdat1$chamber)
chdat2.l <- split(chdat2,chdat2$chamber)

lower <- c(0,0)
upper <- c(150,1)

fits <- list()
preds <- list()
for(i in 1:length(chdat1.l)){
  tofit <- rbind(chdat1.l[[i]],chdat2.l[[i]])
  
  #- round to nearest minute
  tofit$DateTime <- nearestTimeStep(tofit$DateTime,nminutes=1,align="floor")
  tofit <- tofit[,c("chamber","DateTime","CO2L")]
  
  #- merge in the reference data
  refdat2 <- refdat[,c("DateTime","REFCO2")]
  tofit <- merge(tofit,refdat2,by="DateTime")

  #- fit the model, extract the best parameters, and rerun the model to get the predicted values  
  DEfit <- DEoptim(fn=closedWTC.mod,lower=lower,upper=upper,V=53,Ca=tofit$REFCO2,Ci=tofit$CO2L,
                   breakpoint=210,
                   DEoptim.control(NP = 50))
  DEfit.best <- unname(DEfit$optim$bestmem)
  DEpred <- closedWTC.mod(mu=DEfit.best,V=53,Ca=tofit$REFCO2,Ci=tofit$CO2L,fit=0,breakpoint=210)
  DEpred$chamber <- tofit$chamber[1]
  
  #- extract the parameters and predictions
  fits[[i]] <- data.frame(chamber=tofit$chamber[1],Rd = DEfit.best[1], theta = DEfit.best[2])
  preds[[i]] <- DEpred
  
}

Rfits <- do.call(rbind,fits)
Rpreds <- do.call(rbind,preds)
Rfits$T_treatment <- factor(ifelse(Rfits$chamber %%2 == 1,"ambient","warmed"))


#---- figure out the molar density of air
Ta <- mean(refdat$TAREF)+273.15      # air temperature in K
vapor_pressure <- 6.11*10^((7.5*mean(refdat$DEWPNTC)/(237.3+mean(refdat$DEWPNTC))))/10*1000 #Pa

Dry_air_density <- 0.0035*mean(refdat$AIRPRESS)*1000/Ta  # in kg m-3
Wet_air_density <- 0.0022*vapor_pressure/Ta              # in kg m-3
#Rgas <- 287.05 # J kg-1 K-1

#- calculate the density of air, given the observed moisture. Recall that wet air is LESS dense than dry air.
#    see http://www.engineeringtoolbox.com/density-air-d_680.html
#    This uses the outside reference data. This is not really appropriate, as the chambers are at 15 and maybe a different humidity.
#    I need to download and use the real met data.
#    It would be WAAAY easier to assume 44.6 mol m-3, as I did previously.
# the following had units of kg m-3
Mix_air_density <- Dry_air_density*(1+Wet_air_density/Dry_air_density)/(1+1.609*Wet_air_density/Dry_air_density)
Mix_air_density_mol <- Mix_air_density/28.97*1000  # convert to mol air m-3. See http://www.engineeringtoolbox.com/molecular-mass-air-d_679.html

Rfits$Rd_umol <- Rfits$Rd*Mix_air_density_mol/60  # convert to umol s-1. Recall that the data are recorded minutely, but I want per second

boxplot(Rd_umol~T_treatment,data=Rfits,col=c("blue","red"))
plotBy(Rd_umol~theta|chamber,data=Rfits,legend=F,pch=15)



#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------
#- merge in the tree volume data?
initd <- getwd()
setwd("C://Repos/wtc4sheets")

# Load packages, code.
source("R/load.R")

# Map the HIE-Data2 network drive. Store drive letter here.
hiedrive <- "W:"

# Flag to email the figure or not.
email_figures <- FALSE

# Directory with height / diameter measurements
hddir <- file.path(hiedrive, "WORKING_DATA/WTC/WTC4/Share/RAW_GROWTH_DATA/WTC4 tree measure data")

# Make full dataset
hddata <- make_hddata(hddir)

# Calculate stem volume (for planted trees only)
vol <- return_volume(hddata)

vol2 <- subset(vol,Date>=as.Date("2016-04-13"))
vol2$inc <- c(NA,diff(vol2$vol)) #- get the increment
vol3 <- subset(vol2,Date==as.Date("2016-04-27"))


vol3$Rd <- Rfits$Rd_umol
vol3$theta <- Rfits$theta
vol3$Rd_v <- with(vol3,Rd/vol*10000)
vol3$Rd_inc <- with(vol3,Rd/inc*10000)

windows(90,50);par(mfrow=c(1,3),mar=c(6,7,2,1),cex.lab=2,cex.axis=1.5)
plotBy(CO2L~DateTime|chamber,data=chdat,type="o",pch=16,ylim=c(410,1200),legend=F)
legend("topleft",levels(as.factor(chdat1$chamber)),ncol=6,pch=15,col=palette()[1:12])
lines(REFCO2~DateTime,data=refdat,col="black",lty=2)
legend("topright",lty=2,"Atmosphere")
# R per unit stem volume
plotBy(Rd~vol|T_treatment,pch=15,col=c("blue","red"),cex=1.5,data=vol3,
       xlim=c(400,2500),
       xlab="Stem volume (cm3)",ylab="Rcanopy raw (umol s-1)")
textxy(X=vol3$vol,Y=vol3$Rd,labs=vol3$chamber,cex=0.7)

# R per unit stem volume INCREMENT
plotBy(Rd~inc|T_treatment,pch=15,col=c("blue","red"),cex=1.5,data=vol3,
       xlim=c(000,700),
       xlab="Stem volume growth (cm3)",ylab="Rcanopy raw (umol s-1)")
textxy(X=vol3$inc,Y=vol3$Rd,labs=vol3$chamber,cex=0.7)
#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------


