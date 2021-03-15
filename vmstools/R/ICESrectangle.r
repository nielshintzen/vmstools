#' Get ICES rectangle from coordinates
#' 
#' Get the ICES rectangle from any lon,lat position, given that this position
#' is within the ICES region.
#' 
#' 
#' @param dF dataframe given that they have 'SI_LONG' and 'SI_LATI' columns
#' @return Returns the rectangles as a vector
#' @author Neil Campbell
#' @seealso \code{\link{ICESarea}}, \code{\link{ICESrectangle2LonLat}}
#' @references EU Lot 2 project
#' @examples
#' 
#' res <- data.frame(SI_LONG = c(1,2,2,4,2),
#'                   SI_LATI = c(53,53.2,54,56.7,55.2))
#' ICESrectangle(res)
#' 
#' 
#' @export ICESrectangle
ICESrectangle <- function(dF){
                    rectChar1n2 <- sprintf("%02i",as.integer(2 * (dF[, "SI_LATI"] - 35.5)))
                    rectChar3 <- ifelse(dF[, "SI_LONG"] > -50 & dF[, "SI_LONG"]<= -40, "A",
                                 ifelse(dF[, "SI_LONG"] > -40 & dF[, "SI_LONG"]<= -30, "B",
                                 ifelse(dF[, "SI_LONG"] > -30 & dF[, "SI_LONG"]<= -20, "C",
                                 ifelse(dF[, "SI_LONG"] > -20 & dF[, "SI_LONG"]<= -10, "D",
                                 ifelse(dF[, "SI_LONG"] > -10 & dF[, "SI_LONG"]<    0, "E",
                                 #-Note that at 0 meridian the allocation of points at the meridian switch
                                 ifelse(dF[, "SI_LONG"] >=  0 & dF[, "SI_LONG"]<  10, "F",
                                 ifelse(dF[, "SI_LONG"] >= 10 & dF[, "SI_LONG"]<  20, "G",
                                 ifelse(dF[, "SI_LONG"] >= 20 & dF[, "SI_LONG"]<  30, "H", "J"))))))))
                    rectChar4  <- rep(NA,nrow(dF))
                    idxlowzero <- which(dF[,"SI_LONG"] <  0)
                    idxabozero <- which(dF[,"SI_LONG"] >= 0)
                    if(length(idxlowzero)>0) rectChar4[idxlowzero] <- ceiling(dF[idxlowzero,"SI_LONG"] %% 10 -1 + 10)%%10
                    if(length(idxabozero)>0) rectChar4[idxabozero] <- floor(dF[idxabozero,"SI_LONG"] %% 10)
                    rectID <- paste(rectChar1n2, rectChar3, rectChar4, sep = "")
                 return(rectID)}


#dF <- as.data.frame(cbind(SI_LONG = seq(-20,-9,0.5),SI_LATI = rep(54,length(seq(-20,-9,0.5)))))
#dF <- as.data.frame(cbind(SI_LONG = seq(0,12,0.5),SI_LATI = rep(54,length(seq(0,12,0.5)))))
#dF$ICES <- ICESrectangle(dF)
