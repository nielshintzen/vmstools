## DCFIndicator6.r
## by Fabrizio Manco, 22/09/2010
## calculates the DCF6 indicator : for each month, the area where 90% of the pings were obtained

DCFIndicator6 <- function ( tacsat,
                            pctThreshold=90,              # percentage of points to include
                            cellresX=0.05, 
                            cellresY=0.05,
                            calcAreaMethod="Trapezoid",   # "Trapezoid" (fast and less accurate, good for small cellsizes) or "UTM" (accurate but slow, good for huge cellsizes)
                            plotMapTF = FALSE,
                            exportTableName=""
                            )
  {

  #remove rows with NA
  tacsat<-tacsat[complete.cases(tacsat),]
  
  if (!"SI_DATIM" %in% colnames(tacsat)) {tacsat$SI_DATIM  <- as.POSIXct(paste(tacsat$SI_DATE, tacsat$SI_TIME,sep = " "), tz = "GMT", format = "%d/%m/%Y  %H:%M:%S")}
  ltMonth<-unique(as.numeric(format(tacsat$SI_DATIM, format="%m")))
  ltMonth<-sort(ltMonth)
  tableResultDCF6=data.frame(month=ltMonth, DCF6=rep(0, length(ltMonth)))

  for (x in 1:length(ltMonth)){
    currMonth<-ltMonth[x]
    monthlyTacsat<-subset(tacsat, as.numeric(format(tacsat$SI_DATIM, format="%m"))==currMonth)
    
    # flag the pings inside the MCP according to the pctThreshold
    monthlyTacsat<-tacsatMinimumConvexPolygon(monthlyTacsat, pctThreshold)
    
    # grid the vms pings inside the MCP
    monthlyTacsat<-subset(monthlyTacsat, monthlyTacsat$INMCP!=0)
    if (plotMapTF) {windows(record=TRUE)}
    monthlyVmsGrid<-vmsGridCreate(monthlyTacsat, nameLon = "SI_LONG", nameLat = "SI_LATI", cellsizeX=cellresX, cellsizeY=cellresY, nameVarToSum = "INMCP", plotMap=plotMapTF, plotTitle=paste("Month ", currMonth), plotPoints = FALSE)
    if (plotMapTF) {plot_mcp(plotnew=FALSE, plotpoints=FALSE)}
    
      # calculate the area of each cell in square km
      monthlyVmsGrid<-surface(monthlyVmsGrid, method=calcAreaMethod, includeNA=FALSE)
            
    tableResultDCF6[x,2]<-sum(monthlyVmsGrid@data$cellArea[!is.na(monthlyVmsGrid@data$fishing) & monthlyVmsGrid@data$fishing!=0])
    }

  if (exportTableName!="") {write.csv(tableResultDCF6, exportTableName)}
  return(tableResultDCF6)
}