countPings <- function(formula,tacsat,grid=NULL){

                  origTacsat <- tacsat
                  #Turn selected variabels into element list
                  form <- formula
                  if (form[[1]] != "~")
                      stop("Error: Formula must be one-sided.")
                  formc <- as.character(form[2])
                  formc <- gsub(" ", "", formc)
                  if (!is.element(substring(formc, 1, 1), c("+", "-")))
                      formc <- paste("+", formc, sep = "")
                  vars <- unlist(strsplit(formc, "[\\+\\-]"))
                  vars <- vars[vars != ""]
                  signs <- formc
                  for (i in 1:length(vars)) {
                      signs <- gsub(vars[i], "", signs)
                  }
                  signs <- unlist(strsplit(signs, "")) #Currently we do not use signs

                  #Define which variables selected are column names, time variables or spatial variables
                  Vars      <- vars[which(!vars %in% c("day","week","month","quarter","year","gridcell","icesrectangle","icesarea"))]
                  timeVars  <- vars[which(vars %in% c("day","week","month","quarter","year"))]
                  spatVars  <- vars[which(vars %in% c("gridcell","icesrectangle","icesarea"))]

                  #Add time notation if you want this as output
                  if(length(timeVars)>0){
                    if(!length(grep("SI_DATIM",colnames(tacsat)))>0) tacsat$SI_DATIM <- as.POSIXct(paste(tacsat$SI_DATE,  tacsat$SI_TIME,   sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M")
                    if("day" %in% timeVars){       tacsat$SI_DAY   <- an(format(tacsat$SI_DATIM,format="%j"))   ; timeVars[which(timeVars=="day")]      <- "SI_DAY"}
                    if("week" %in% timeVars){      tacsat$SI_WEEK  <- an(format(tacsat$SI_DATIM,format="%W"))   ; timeVars[which(timeVars=="week")]     <- "SI_WEEK" }
                    if("month" %in% timeVars){     tacsat$SI_MONTH <- an(format(tacsat$SI_DATIM,format="%m"))   ; timeVars[which(timeVars=="month")]    <- "SI_MONTH"}
                    if("quarter" %in% timeVars){   tacsat$SI_QUART <- an(substr(quarters(tacsat$SI_DATIM),2,2)) ; timeVars[which(timeVars=="quarter")]  <- "SI_QUART"}
                    if("year" %in% timeVars){      tacsat$SI_YEAR  <- an(format(tacsat$SI_DATIM,format="%Y"))   ; timeVars[which(timeVars=="year")]     <- "SI_YEAR" }
                  }
                  #Add spatial notation if you want this as output
                  if(length(spatVars)>0){
                    if("gridcell" %in% spatVars & is.null(grid) == T) stop("Grid needs to be specified to use the 'gridcell' option")
                    if("gridcell" %in% spatVars & is.null(grid) == F){
                      #Create coordinates of tacsat data
                      coords                    <- cbind(x=tacsat$SI_LONG,y=tacsat$SI_LATI)
                      sPDF                      <- SpatialPointsDataFrame(coords,data=tacsat)
                      #Turn grid into a spatial pixel dataframe
                      grid                      <- as(grid,"SpatialPixels");
                      grid                      <- as(grid,"SpatialPixelsDataFrame")
                      #Overlay the two spatial frameworks to see to which gridcell each tacsat coordinate belongs
                      gridCellIndex             <- overlay( grid, sPDF )
                      newCoords                 <- sPDF@coords[gridCellIndex,]
                      
                      tacsat$GR_LONG            <- newCoords[,1]
                      tacsat$GR_LATI            <- newCoords[,2]
                      spatVars[which(spatVars=="gridcell")] <- "GR_LONG"; spatVars <- c(spatVars,"GR_LATI")
                    }
                    if("icesrectangle" %in% spatVars){   tacsat$LE_RECT <- ICESrectangle(tacsat); spatVars[which(spatVars=="icesrectangle")] <- "LE_RECT" }
                    if("icesarea" %in% spatVars){        tacsat$LE_AREA <- ICESarea(tacsat)     ; spatVars[which(spatVars=="icesarea")] <- "LE_AREA"}
                  }

                  tacsat$SUM  <- 1
                  buildList   <- list()
                  totVars     <- c(Vars,timeVars,spatVars)
                  for(iVars in 1:length(totVars)) buildList[[iVars]] <- tacsat[,totVars[iVars]]

                  #Do the counting of pings
                  res           <- aggregate(tacsat$SUM,buildList,sum)
                  colnames(res) <- c(totVars,"pings")

              return(res)}

