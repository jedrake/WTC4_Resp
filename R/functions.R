#------------------------------------------------------------------------------------------------------------------
#- Functions that do useful things for the leaf and canopy respiration project in WTC4
#------------------------------------------------------------------------------------------------------------------



#------------------------------------------------------------------------------------------------------------------
#- function to read and download MJA's WTC3 respiraiton data (at 25 deg C). Also return the thermal history
#   (i.e., the mean temperature of a specified number of days preceding the measurements)
#   Returns a dataframe.
#------------------------------------------------------------------------------------------------------------------
returnMJAdata <- function(ndays=3){
  
  #----------------------------------------------------------------
  #- get Mike's data
  #downloadHIEv(searchHIEv("GX-Rdark25"),topath="Data")
  mja <- read.csv("Data/WTC_TEMP_CM_GX-Rdark25_20130617-20140402_L2.csv")
  mja$Rarea25 <- mja$Photo * -1
  mja$Date <- as.Date(mja$date)
  
  mja.m <-summaryBy(Rarea25~chamber+Date+T_treatment,data=subset(mja,Water_treatment =="control"),FUN=mean,keep.names=T)
  
  
  #- download the met data from WTC3 and read them in
  #downloadHIEv(searchHIEv("WTC_TEMP_CM_WTCMET"),topath="Data")
  
  metfiles <- list.files(path="Data",pattern="WTC_TEMP_CM_WTCMET",full.names=T)
  
  print("Reading met data")
  #- read in the data
  met.l <- list()
  for(i in 1:length(metfiles)){
    met.l[[i]] <- read.csv(metfiles[i])
  }
  met <- do.call(rbind,met.l)
  met$Date <- as.Date(met$DateTime)
  
  print("Interpolating met data")
    #-- loop over each observation in Mike's dataset, calculate the mean temperature of the preceding X dates (3)
  ndays <- 3
  mja.m$Tpre <- NA
  for (i in 1:nrow(mja.m)){
    searchdate <- mja.m$Date[i]
    mindate <- searchdate-ndays
    inds <- which(met$Date >=mindate & met$Date <=searchdate & met$chamber == mja.m$chamber[i])
    mja.m$Tpre[i] <- mean(met[inds,"Tair_al"],na.rm=T)
    
  }  
  print("Done")

  
  return(mja.m)
}
#------------------------------------------------------------------------------------------------------------------








#------------------------------------------------------------------------------------------------------------------
#- Return the WTC4 leaf respiration data, measured at 15 deg C
#------------------------------------------------------------------------------------------------------------------
returnRleaf <- function(){

  #- read in the data
  #- find the right csv files and extract the date
  files <- list.files(path="Data",pattern="GX-Rleaf",full.names=T)
  files <- files[grep(".csv",files)]
  dates <- as.Date(substr(gsub("[^0-9]","",files),start=1,stop=8),format="%Y%m%d")
  
  #- read in the data
  dat.l <- list()
  for(i in 1:length(files)){
    dat.l[[i]] <- read.csv(files[i])
    dat.l[[i]]$Date <- dates[i]
  }
  dat <- do.call(rbind,dat.l)
  dat$chamber <- factor(dat$chamber,levels=c("C01","C02","C03","C04","C05","C06","C07","C08",'C09',
                                             "C10","C11","C12"))
  linkdf <- data.frame(chamber=levels(dat$chamber),T_treatment=rep(c("ambient","elevated"),6))
  dat.m <- merge(dat.m,linkdf)
  
  dat.m$Rarea <- -1*dat.m$Photo
  dat.m$Rarea25 <- with(dat.m,Rarea*2) # assumes Q10 = 2
  
  return(dat.m)
}
#------------------------------------------------------------------------------------------------------------------
