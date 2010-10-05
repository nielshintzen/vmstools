## DCFIndicator5.r
## by Fabrizio Manco, 22/09/2010
## calculates the DCF5 indicator : for each month, the total area where fishing activity occurs

DCFIndicator5 <- function ( tacsat,
                            minThreshold=10,            # if time interval has been calculated (and named SI_INTV), it's a minimal nb of minutes, otherwise, it's minimal number of points
                            cellresX=0.05,
                            cellresY=0.05,
                            plotMapTF=FALSE,
                            exportTableName=""
                            )
  {

  
  #remove rows with NA
  tacsat<-tacsat[complete.cases(tacsat),]
  
  if (!"SI_DATIM" %in% colnames(tacsat)) {tacsat$SI_DATIM  <- as.POSIXct(paste(tacsat$SI_DATE, tacsat$SI_TIME,sep = " "), tz = "GMT", format = "%d/%m/%Y  %H:%M:%S")}
  
  ltMonth<-unique(as.numeric(format(tacsat$SI_DATIM, format="%m")))
  ltMonth<-sort(ltMonth)
  tableResultDCF5=data.frame(month=ltMonth, DCF5=rep(0, length(ltMonth)))
  
    for (x in 1:length(ltMonth)){
      currMonth<-ltMonth[x]
      monthlyTacsat<-subset(tacsat, as.numeric(format(tacsat$SI_DATIM, format="%m"))==currMonth)
    
      if ("SI_INTV" %in% colnames(tacsat)) { nameVarToSum="SI_INTV"} else {nameVarToSum=""}
      
      if (plotMapTF) {windows(record=TRUE)}
      monthlyVmsGrid<-vmsGridCreate(monthlyTacsat, nameLon = "SI_LONG", nameLat = "SI_LATI", cellsizeX=cellresX, cellsizeY=cellresY, nameVarToSum, plotMap=plotMapTF, plotPoints = FALSE)
      
      # calculate the area of each cell in square km
      monthlyVmsGrid<-calcAreaOfCells(monthlyVmsGrid)

      filteredFishingValues<-subset(monthlyVmsGrid@data, monthlyVmsGrid@data$fishing>minThreshold)
      tableResultDCF5[x,2]<-sum(filteredFishingValues$fishing)
      }
  if (exportTableName!="") {write.csv(tableResultDCF5, exportTableName)}
  return(tableResultDCF5)
}