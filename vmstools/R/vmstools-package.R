

#' Vessel logbook dataset
#' 
#' A dataset consisting of logbook (landings / values) registrations of
#' disguised origin for 2 consecutice years.
#' 
#' This example dataset if for example purposes only. Without prior permission
#' of the authors it is not allowed to use this data other than for example
#' non-publishable purposes.
#' 
#' @name eflalo
#' @docType data
#' @format A data frame with 4539 observations on the following 189 variables.
#' \describe{ \item{list("VE_REF")}{Vessel reference / ID}
#' \item{list("VE_FLT")}{Fleet reference (DCF regulation)}
#' \item{list("VE_COU")}{Flag nation of vessel} \item{list("VE_LEN")}{Vessel
#' length} \item{list("VE_KW")}{Vessel power (in kW)}
#' \item{list("VE_TON")}{Vessel Tonnage (GT)} \item{list("FT_REF")}{Fishing
#' trip reference number} \item{list("FT_DCOU")}{Departure country}
#' \item{list("FT_DHAR")}{Departure harbour} \item{list("FT_DDAT")}{Departure
#' date} \item{list("FT_DTIME")}{Departure time} \item{list("FT_LCOU")}{Landing
#' country} \item{list("FT_LHAR")}{Landing harbour}
#' \item{list("FT_LDAT")}{Landing date} \item{list("FT_LTIME")}{Landing time}
#' \item{list("LE_ID")}{Log event ID} \item{list("LE_CDAT")}{Catch date}
#' \item{list("LE_STIME")}{Log event start time} \item{list("LE_ETIME")}{Log
#' event end time} \item{list("LE_SLAT")}{Log event start position latitude}
#' \item{list("LE_SLON")}{Log event start position longitude}
#' \item{list("LE_ELAT")}{Log event end position latitude}
#' \item{list("LE_ELON")}{Log event end position longitude}
#' \item{list("LE_GEAR")}{Gear} \item{list("LE_MSZ")}{Mesh size}
#' \item{list("LE_RECT")}{ICES rectangle} \item{list("LE_DIV")}{ICES division}
#' \item{list("LE_MET_level6")}{Fishing activity}
#' \item{list("LE_KG_SPEC")}{Landing weight estimate of species (FAO species
#' codes)} \item{list("LE_EURO_SPEC")}{Landing value of species (FAO species
#' codes)} }
#' @references EU Lot 2 project, Niels T. Hintzen (niels.hintzen@wur.nl)
#' @source NA
#' @keywords datasets
#' @examples
#' 
#' data(eflalo)
#' eflalo <- formatEflalo(eflalo)
#' str(eflalo)
#' 
NULL





#' Dataframe containing > 3500 harbours primarily Europe
#' 
#' A list of harbour positions primarily located in Europe, including -name,
#' -latitude, -longitude and a certain -range (in km) the harbour spans
#' (radius).
#' 
#' 
#' @name euharbours
#' @aliases euharbours harbours
#' @docType data
#' @references EU Lot 2 project, Niels T. Hintzen (niels.hintzen@wur.nl)
#' @keywords datasets
#' @examples
#' 
#' data(tacsat)
#' data(euharbours); euharbours <- harbours
#' 
#' #modify the range of the harbours to e.g. 2km
#' euharbours$range <- 2
#' 
#' pointInHarbour(tacsat$SI_LONG,tacsat$SI_LATI,euharbours)
#' 
NULL





#' Polygon dataset of Europe
#' 
#' Polygon dataset with the coastlines of Europe.
#' 
#' 
#' @name europa
#' @docType data
#' @format A data frame with 83878 observations on the following 5 variables.
#' \describe{ \item{list("PID")}{a numeric vector} \item{list("SID")}{a numeric
#' vector} \item{list("POS")}{a numeric vector} \item{list("X")}{Longitude
#' values} \item{list("Y")}{Latitude values} }
#' @source NOAA coastline datafile
#' @keywords datasets
#' @examples
#' 
#' data(europa)
#' require(PBSmapping)
#' 
#' eurPols     <- lonLat2SpatialPolygons(lst=lapply(as.list(sort(unique(europa$SID))),
#'                         function(x){data.frame(
#'                           SI_LONG=subset(europa,SID==x)$X,
#'                           SI_LATI=subset(europa,SID==x)$Y)}))
#' 
#' plot(eurPols,col="green",xlim=c(-4,10),ylim=c(48,62))
#' 
NULL





#' SpatialPolygons dataset of ICESareas
#' 
#' A pre-defined SpatialPolygons dataset with the ICES areas.
#' 
#' 
#' @name ICESareas
#' @docType data
#' @source www.ices.dk
#' @keywords datasets
#' @examples
#' 
#' require(sp)
#' data(ICESareas)
#' plot(ICESareas)
#' 
NULL





#' Merge TACSAT formatted data to EFLALO2 formatted data
#' 
#' Method to assign the EFLALO2 tripnumber to the accompanying TACSAT records,
#' which were registered between the departure date and arrival date of a trip.
#' 
#' Method only assignes a tripnumber to the tacsat data if there are vms pings
#' that occur between the departure date and arrival date of a specific vessel.
#' If no related trip can be found, tripnumber assigned to tacsat equals NA.
#' 
#' @name mergeEflalo2Tacsat
#' @docType data
#' @param eflalo2 Trip and landings data in the Eflalo2 format
#' @param tacsat Vms data in the Tacsat format
#' @author Niels T. Hintzen, Francois Bastardie
#' @seealso \code{\link{findEndTacsat}}, \code{\link{sortTacsat}},
#' \code{\link{filterTacsat}}
#' @references EU lot 2 project
#' @examples
#' 
#'   #Load the eflalo2 and tacsat data
#' data(eflalo)
#' data(tacsat)
#' 
#' #-Remove duplicated records from tacsat
#' myf       <- paste(tacsat$VE_REF,tacsat$SI_LATI,tacsat$SI_LONG,
#'                    tacsat$SI_DATE,tacsat$SI_TIME);
#' tacsat    <- tacsat[!duplicated(myf),];
#' 
#' #-Merge the tacsat and eflalo2 data
#' #computation time approx 30 seconds
#' result <- mergeEflalo2Tacsat(eflalo,tacsat)
#' 
NULL





#' Useful functions for the multivariate analysis of logbooks data for
#' identifying metiers.
#' 
#' This function contains several functions needed for the multivariate
#' analysis of logbooks data for identifying metiers.
#' 
#' 
#' @param transformation_proportion Transform quantities to percentage values
#' (between 0 and 100) of each species in the logevent total catch.
#' @param table_variables Transpose the dataset (change variables into
#' individuals)
#' @param scree Implementation of "scree-test"
#' @param select_species Remove the cluster with the smallest mean of capture
#' @param building_tab_pca Build the table with the main species
#' @param test.values Compute the test-value for each species by cluster
#' @param targetspecies Determine the species with a test-value > 1.96 by
#' cluster
#' @param withinVar Calculate the cluster's within-variance
#' @note A number of libraries are initially called for the whole metier
#' analyses and must be installed :
#' (FactoMineR),(cluster),(SOAR),(amap),(MASS),(mda)
#' @author Nicolas Deporte, Sebastien Demaneche, Stephanie Mahevas (IFREMER,
#' France), Clara Ulrich, Francois Bastardie (DTU Aqua, Denmark)
#' @references Development of tools for logbook and VMS data analysis. Studies
#' for carrying out the common fisheries policy No MARE/2008/10 Lot 2
NULL





#' Vessel Monitoring by Satellite system dataset
#' 
#' A dataset consisting of VMS registrations of disguised origin for 2
#' consecutice years.
#' 
#' This example dataset if for example purposes only. Without prior permission
#' of the authors it is not allowed to use this data other than for example
#' non-publishable purposes.
#' 
#' @name tacsat
#' @docType data
#' @format A data frame with 97015 observations on the following 8 variables.
#' \describe{ \item{list("VE_COU")}{Flag nation of vessel}
#' \item{list("VE_REF")}{Vessel reference / ID}
#' \item{list("SI_LATI")}{Latitudinal position (in decimal degrees)}
#' \item{list("SI_LONG")}{Longitudinal position (in decimal degrees)}
#' \item{list("SI_DATE")}{Date} \item{list("SI_TIME")}{Time (24 hours clock)}
#' \item{list("SI_SP")}{Speed (in knots per hour)} \item{list("SI_HE")}{Heading
#' (in degrees)} }
#' @references EU Lot 2 project, Niels T. Hintzen (niels.hintzen@wur.nl)
#' @source NA
#' @keywords datasets
#' @examples
#' 
#' data(tacsat)
#' tacsat <- formatTacsat(tacsat)
#' str(tacsat)
#' 
NULL





#' Small high frequency test VMS dataset CANNOT BE DISTRIBUTED WITHOUT PRIOR
#' PERMISSION -> SEE Author
#' 
#' A short VMS dataset with high frequency interval rate (every 6 minutes) used
#' to test spatial and VMS functions. Dataframe contains vessel ID, decimal
#' longitude and decimal latitude position, speed and heading at ping and date
#' of ping.
#' 
#' 
#' @name VMShf
#' @docType data
#' @format A data frame with 702 observations on the following 6 variables.
#' \describe{ \item{list("ship")}{ID of ship} \item{list("declat")}{Decimal
#' notation of longitude position} \item{list("declon")}{Decimal notation of
#' latitude postion} \item{list("speed")}{Speed at ping}
#' \item{list("heading")}{Heading at ping} \item{list("date")}{Date and time of
#' ping} }
#' @references niels.hintzen@wur.nl
#' @source niels.hintzen@wur.nl
#' @examples
#' 
#' data(VMShf)
#' str(VMShf)
#' 
#' 
NULL





#' development of R scripts to setup the workflow of vms and logbook data
#' 
#' -Design of routines to import VMS and Logbook data\cr -Design of routines to
#' cluster metiers\cr -Design of routines to link VMS and Logbook data\cr
#' -Design of routines to define vessel activity\cr -Design of routines to
#' interpolate VMS data\cr -Design of routines to compute pressure
#' indicators\cr
#' 
#' \tabular{ll}{ Package: \tab VMStools\cr Type: \tab Package\cr Version: \tab
#' 0.58\cr Date: \tab 2011-12-30\cr License: \tab \cr LazyLoad: \tab yes\cr }
#' 
#' @name VMStools-package
#' @aliases VMStools-package VMStools
#' @docType package
#' @author Niels T. Hintzen, Doug Beare, Francois Bastardie, Nicolas Deporte,
#' Andy South, Neil Campbell, Fabrizio Manco, Clara Ulrich Rescan, Stephanie
#' Mahevas, Sebastien Demaneche, Josefine Egekvist, Hans Gerritsen, Matthew
#' Parker-Humphreys
#' 
#' Maintainer: Niels T. Hintzen <niels.hintzen@@wur.nl>, Francois Basterdie
#' <fba@@aqua.dtu.dk>
#' @references EU Lot 2 project
NULL



