#' Count the number of VMS pings in a selection
#' 
#' Count the number of VMS pings in any selection made in time and spatial
#' frame.
#' 
#' Formula has form: ~Variable+timeVariable+anotherTimeVariable+spatialVariable
#' \cr
#' 
#' options in formula for Variable: any column name of tacsat file \cr options
#' in formula for timeVariable: day, week, month, quarter, year \cr options in
#' formula for spatialVariable: icesrectangle, icesarea, gridcell \cr
#' 
#' @param formula specify the elements you want to use as axis of the count
#' @param tacsat tacsat data file, possibly with additional columns
#' @param grid if in formulate 'gridcell' is chosen, a SpatialGrid must be
#' provided
#' @return Returns the matrix with counted pings by each specified variable
#' @note if Tacsat is a big file, the overlay function might fail resulting in
#' terminating the function
#' @author Niels T. Hintzen
#' @seealso \code{\link{createGrid}}, \code{\link{vmsGridCreate}}
#' @references Hintzen et al. 2010 Improved estimation of trawling tracks using
#' cubic Hermite spline interpolation of position registration data, EU lot 2
#' project
#' @examples
#' 
#' data(tacsat)
#' 
#' #make the tacsat file a bit smaller
#' tacsat  <- tacsat[1:10000,]
#' 
#' grid          <- createGrid(range(tacsat$SI_LONG,na.rm=TRUE),
#'                   range(tacsat$SI_LATI,na.rm=TRUE),0.5,0.5,type="SpatialGrid")
#' 
#' result        <- countPings(~VE_REF+year+gridcell,tacsat,grid=grid)
#' result        <- countPings(~VE_REF+week+year+icesrectangle,tacsat)
#' 
#' @export countPings
countPings <- function(formula,tacsat,grid=NULL,by=NULL){
                  require(sf)
                  require(data.table)
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
                    if("day" %in% timeVars & !"SI_DAY" %in% colnames(tacsat)){       tacsat$SI_DAY   <- an(format(tacsat$SI_DATIM,format="%j"))};      if("day" %in% timeVars){   ; timeVars[which(timeVars=="day")]      <- "SI_DAY"}
                    if("week" %in% timeVars & !"SI_WEEK" %in% colnames(tacsat)){      tacsat$SI_WEEK  <- an(format(tacsat$SI_DATIM,format="%W"))};     if("week" %in% timeVars){   ; timeVars[which(timeVars=="week")]     <- "SI_WEEK" }
                    if("month" %in% timeVars & !"SI_MONTH" %in% colnames(tacsat)){     tacsat$SI_MONTH <- an(format(tacsat$SI_DATIM,format="%m"))};    if("month" %in% timeVars){   ; timeVars[which(timeVars=="month")]    <- "SI_MONTH"}
                    if("quarter" %in% timeVars & !"SI_QUART" %in% colnames(tacsat)){   tacsat$SI_QUART <- an(substr(quarters(tacsat$SI_DATIM),2,2))};  if("quarter" %in% timeVars){ ; timeVars[which(timeVars=="quarter")]  <- "SI_QUART"}
                    if("year" %in% timeVars & !"SI_YEAR" %in% colnames(tacsat)){      tacsat$SI_YEAR  <- an(format(tacsat$SI_DATIM,format="%Y"))};     if("year" %in% timeVars){   ; timeVars[which(timeVars=="year")]     <- "SI_YEAR" }
                  }
                  #Add spatial notation if you want this as output
                  if(length(spatVars)>0){
                    if("gridcell" %in% spatVars & is.null(grid) == TRUE) stop("Grid needs to be specified to use the 'gridcell' option")
                    if("gridcell" %in% spatVars & is.null(grid) == FALSE){
                      #Create coordinates of tacsat data
                      coords                    <- cbind(x=tacsat$SI_LONG,y=tacsat$SI_LATI)
                      sPDF                      <- st_as_sf(tacsat,coords=c("SI_LONG","SI_LATI"))
                      idx                       <- sapply(st_intersects(sPDF,grid), function(z) if (length(z)==0) NA_integer_ else z[1])
                      newCoords                 <- st_coordinates(st_centroid(grid))[idx,]

                      tacsat$GR_LONG            <- newCoords[,1]
                      tacsat$GR_LATI            <- newCoords[,2]
                      spatVars[which(spatVars=="gridcell")] <- "GR_LONG"; spatVars <- c(spatVars,"GR_LATI")
                    }
                    if("icesrectangle" %in% spatVars){
                      if(!"LE_RECT" %in% colnames(tacsat))
                        tacsat$LE_RECT <- ICESrectangle(tacsat)
                      spatVars[which(spatVars=="icesrectangle")] <- "LE_RECT"
                    }
                    if("icesarea" %in% spatVars){
                      if(!"LE_AREA" %in% colnames(tacsat))
                        tacsat$LE_AREA <- ICESarea(tacsat)
                      spatVars[which(spatVars=="icesarea")] <- "LE_AREA"
                    }
                  }

                  if(is.null(by)){
                    tacsat$SUM      <- 1
                  } else {
                    tacsat$SUM      <- tacsat[,by]
                    }
                  totVars         <- c(Vars,timeVars,spatVars)

                  #Do the counting of pings
                  for(iVars in 1:length(totVars)) tacsat[,totVars[iVars]] <- af(ac(tacsat[,totVars[iVars]]))
                  DT              <- data.table(tacsat)
                  eq              <- c.listquote(totVars)

                  res             <- DT[,sum(SUM),by=eval(eq)]
                  setnames(res,colnames(res),c(totVars,"pings"))
                  #colnames(res)   <- c(totVars,"pings")

              return(data.frame(res))}
