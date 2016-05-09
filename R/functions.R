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
  mja.m$Rarea15 <- mja.m$Rarea25*2.1^((15-25)/10) # assumes Q10 = 2.1 (justified by data)
  
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
  dat.m <- merge(dat,linkdf)
  
  dat.m$Rarea <- -1*dat.m$Photo
  dat.m$Rarea25 <- with(dat.m,Rarea*2) # assumes Q10 = 2
  
  return(dat.m)
}
#------------------------------------------------------------------------------------------------------------------






#------------------------------------------------------------------------------------------------------------------
# Download and convert WTC trendlog data
#------------------------------------------------------------------------------------------------------------------
getWTCtrendlog <- function(Date=Sys.Date()-1, timestep=15){
  
  
  # ODBC stuff only works on 32 bit
  if(.Platform$r_arch == "x64")
    stop("This function only works on 32bit Windows.")
  
  # Download and unzip a zipped trendlog file
  wtctrend <- searchHIEv("WTC_Trendlogs_", quiet=TRUE)
  
  # List all dates available.
  date <- as.Date(ymd(str_extract(wtctrend$filename, "[0-9]{8}")))
  
  #
  Date <- as.Date(Date)
  if(is.na(Date))stop("Date malformed, try again (YYYY-MM-DD)")
  if(!Date %in% date)stop("Date not found in WTC trendlogs.")
  i <- match(Date,date)
  
  # Download record
  message("Downloading data...", appendLF=FALSE)
  d <- downloadHIEv(wtctrend, i, quiet=TRUE)
  wtcdir <- "wtctrendlogs"
  unlink(wtcdir)
  message("done.")
  
  # Unzip
  dir.create(wtcdir, showWarnings = FALSE)
  u <- unzip(d, exdir=wtcdir)
  
  # A file with the trendlog key. 
  # This is saved as loglist.RData. Redo this if needed.
  #     loglist <- read.csv("C:/Repos/WTCautoscript/trendloglist.csv", stringsAsFactors=FALSE)
  #     loglist <- subset(loglist, trendlog != "xxx")
  #     loglist$chamber <- as.factor(loglist$chamber)
  #     lv <- levels(loglist$chamber)
  #     levels(loglist$chamber) <- gsub("ch","C",lv)
  #     loglist$trendlog <- as.numeric(loglist$trendlog)
  #     save(loglist, file="data/loglist.RData")
  
  # The list of trendlogs we are interested in
  trendlog_nrs_list <- suppressWarnings(as.numeric(as.character(loglist$trendlog)))
  
  
  # Function to extract a table from an mdb, add chamber and variable fields.
  getData <- function(mdbfile){
    
    on.exit(if(exists("odb"))odbcClose(odb))
    
    odb <- try(odbcConnectAccess(mdbfile))
    if(inherits(odb, "try-error"))return(NULL)
    
    # Find variable name from loglist
    lognr <- as.numeric(str_extract(mdbfile, "[0-9]{10}"))
    varname <- loglist[trendlog_nrs_list ==lognr,"varname"]
    if(length(varname) == 0)return(NULL)
    if(length(varname) > 1){
      stop("Fatal error: some duplicate trendlog IDs in loglist.")
    }
    cham <- loglist[trendlog_nrs_list ==lognr,"chamber"]
    
    # Fetch table from mdb
    tabs <- sqlTables(odb)$TABLE_NAME
    if(!any(grepl("tblTrendlog_",tabs)))return(NULL)
    data <- sqlFetch(odb,  tabs[grep("tblTrendlog_",tabs)])
    
    # Reread Date - to ensure it is in UTC timezone.
    data$TimeOfSample <- ymd_hms(as.character(data$TimeOfSample))
    if(all(is.na(data$TimeOfSample)))return(NULL)
    
    # Nearest timestep (forward searching)
    data$DateTime <- nearestTimeStep(data$TimeOfSample, nminutes=timestep)
    
    # Keep only 15minutely averaged data.
    data <- aggregate(SampleValue ~ DateTime, FUN=mean,data=data)
    data$variable <- varname
    data$chamber <- cham
    
    return(data)
  }
  
  # Read all mdbs. Discard NULL ones. Row-bind them all.
  message("Converting data...", appendLF=FALSE)
  dats <- lapply(u,getData)
  dats <- dats[! sapply(dats, is.null)]
  dats <- suppressWarnings(as.data.frame(rbind_all(dats)))
  if(nrow(dats) == 0){
    warning("No data for ", Date)
    return(NULL)
  }
  message("done.")
  
  # clean up
  unlink(d)
  unlink(u)
  unlink(wtcdir)
  
  
  message("Reshaping and aggregating...", appendLF=FALSE)
  # Reshape. 
  logdatach <- reshape(dats[which(dats$chamber!='ref'),], direction="wide",
                       idvar=c("DateTime","chamber"),
                       timevar="variable")
  
  
  logdataref <- reshape(dats[which(dats$chamber=='ref'),], direction="wide",
                        idvar=c("DateTime","chamber"),
                        timevar="variable")
  
  
  # I don't know how to avoid the SampleValue.xxx names in the df; clean up.
  nm <- names(logdatach)
  names(logdatach) <- gsub("SampleValue.","",nm)
  nm <- names(logdataref)
  names(logdataref) <- gsub("SampleValue.","",nm)
  # tidy dataframe
  #  logdata$chamber <- as.factor(logdata$chamber)
  #  treatdf <- data.frame(chamber=c("ref", paste0("C", sprintf("%02.0f", 1:12))),
  #                        T_treatment=c("reference", rep(c("ambient","elevated"),6)))
  
  #  logdata <- merge(logdata, treatdf, by="chamber")
  logdatach <- logdatach[order(logdatach$DateTime),]
  logdataref <- logdataref[order(logdataref$DateTime),]
  logdataref$chamber<-NULL #drop redudant chamber column from ref data
  message("done.")
  
  return(list(logdatach,logdataref))
}
#------------------------------------------------------------------------------------------------------------------







#------------------------------------------------------------------------------------------------------------------
#- get the trendlogs for the measurement dates, and the 3 dates preceding the measurement dates
#------------------------------------------------------------------------------------------------------------------
updateTrendlogsRleaf <- function(measuredates){
  #- load the trendlog data frames
  load("Data/TrendlogChDF.RData")
  load("Data/TrendlogRefDF.RData")
  
  for(i in 1:length(measuredates)){
    focaldate <- measuredates[i]
    print(paste("getting trendlogs for measurement date",focaldate,sep=" "))
    
    datestoget <- c(focaldate-3,focaldate-2,focaldate-1,focaldate) # define the dates to get (3 days prior to measurement date)
    
    
    #----
    #- download the 15-minutely trendlogs for each date in turn, add to the dataframe.
    for (j in 1:length(datestoget)){
      
      datesintrendlogs <- unique(TrendlogChDF$Date)
      
      #- if the focal date is not already in teh trendlogs, get it.
      if(datestoget[[j]] %in% datesintrendlogs == F){
        print(paste("getting",datestoget[j],sep=" "))
        
        dat  <- getWTCtrendlog(datestoget[j],timestep=15)
        dat[[1]]$Date <- as.Date(dat[[1]]$DateTime)
        dat[[2]]$Date <- as.Date(dat[[2]]$DateTime)
        
        #- add to trendlogs
        TrendlogChDF<-rbind(TrendlogChDF,dat[[1]])
        TrendlogRefDF<-rbind(TrendlogRefDF,dat[[2]])
      }
      
    }
    #----
    
  }
  TrendlogChDF$chamber <- factor(TrendlogChDF$chamber) # get rid of the "ref" level of this factor
  save(TrendlogChDF,file="Data/TrendlogChDF.RData")
  save(TrendlogRefDF,file="Data/TrendlogRefDF.RData")
  
  print("done")
}
#------------------------------------------------------------------------------------------------------------------





#----------------------------------------------------------------------------------------------------------------
standard.error <- function(dat,na.rm=F,...){
  if(na.rm==T){
    dat <- subset(dat,is.na(dat)==F)
  }
  std <- sd(dat)
  n <- length(dat)
  se <- std/sqrt(n)
  return(se)
}
#----------------------------------------------------------------------------------------------------------------





#----------------------------------------------------------------------------------------------------------------
# Adds error bars to a plot
adderrorbars <- function(x,y,SE,direction,barlen=0.04,...){
  
  if(length(direction)>1)stop("direction must be of length one.")
  #if(direction == "updown")
  #  direction <- c("up","down")
  if(direction == "rightleft" | direction == "leftright")direction <- c("left","right")
  
  if("up" %in% direction)
    arrows(x0=x, x1=x, y0=y, y1=y+SE, code=3, angle=90, length=barlen,...)
  if("down" %in% direction) 
    arrows(x0=x, x1=x, y0=y, y1=y-SE, code=3, angle=90, length=barlen,...)
  if("updown" %in% direction) 
    arrows(x0=x, x1=x, y0=y+SE, y1=y-SE, code=3, angle=90, length=barlen,...)
  if("left" %in% direction) 
    arrows(x0=x, x1=x-SE, y0=y, y1=y, code=3, angle=90, length=barlen,...)
  if("right" %in% direction)
    arrows(x0=x, x1=x+SE, y0=y, y1=y, code=3, angle=90, length=barlen,...)  
  
}
#----------------------------------------------------------------------------------------------------------------

