#' Convert the output (a list) from interpolateTacsat to tacsat
#' 
#' This function takes the list of tracks output by interpolateTacsat and
#' converts them back to tacsat format. This means that all the functions that
#' work on tacsat can work with this new, larger interpolated dataset. For
#' example you could now distribute the landings among the new pings.  The
#' npoints argument is the optional number of pings between each real or actual
#' ping.
#' 
#' The distances between the new vessel positions are assumed to be equal.
#' 
#' @param interpolation A list of tracks between actual pings output by the
#' interpolateTacsat function.
#' @param tacsat Tacsat data used when the interpolateTacsat function was
#' called.
#' @param npoints The number of pings or positions required between each real
#' or actual vessel position or ping
#' @param equalDist Logical. Whether the number of positions returned should be
#' equally spaced or not
#' @return A data frame in tacsat format is returned with estimated locations,
#' speeds and headings between each real or actual vessel position. The column
#' added HL_ID refers to the ID given to each individual interpolation.
#' @author Niels Hintzen, Doug Beare
#' @seealso \code{\link{interpolateTacsat}}
#' @references EU Lot 2 project
#' @examples
#' 
#' 
#' data(tacsat)
#' # Take out first 100 rows
#' input   <- tacsat[1:100,]
#' # set resolution for number of locations or pings required
#' #   between each real or actual vessel location.
#' res1    <- 100
#' 
#' #Do the interpolation
#' interpolation <- interpolateTacsat(tacsat=input,interval=120,
#'                                    margin=12,res=res1,method="cHs",
#'                                    params=list(fm=0.4,distscale=20,sigline=0.2,
#'                                    st=c(4,8)),headingAdjustment=0)
#'                                           
#' #Convert the interpolated data back to tacsat with 10 locations or
#' # pings between each real or actual vessel location.
#' tacsat        <- interpolation2Tacsat(interpolation=interpolation,
#'                                       tacsat=input,npoints=10,equalDist=TRUE)
#'   
#' 
#' @export interpolation2Tacsat
interpolation2Tacsat <- function (interpolation, tacsat, npoints = 10, equalDist = TRUE)
{
    tacsat <- sortTacsat(tacsat)
    if (!"HL_ID" %in% colnames(tacsat))
        tacsat$HL_ID <- 1:nrow(tacsat)
    if (!"SI_DATIM" %in% colnames(tacsat))
        tacsat$SI_DATIM <- as.POSIXct(paste(tacsat$SI_DATE, tacsat$SI_TIME,
            sep = " "), tz = "GMT", format = "%d/%m/%Y  %H:%M")

    #- Convert the interpolation to equally distance spaced points
    if (equalDist) {
        interpolationEQ <- equalDistance(interpolation, npoints)
    } else {
        interpolationEQ <- lapply(interpolation, function(x) {
            idx <- round(seq(2, nrow(x), length.out = npoints))
            return(x[c(1, idx), ])
        })
    }
    #- Get the interpolated longitudes and latitudes out
    int <- lapply(interpolationEQ, function(x) {
        x[-1, ]
    })
    lenint <- lapply(int, length)
    idx <- which(unlist(lenint) == (npoints * 2))
    int <- lapply(interpolationEQ[idx], function(x) {
        x[-1, ]
    })
    lenEQ <- length(interpolationEQ[idx])
    intmin2 <- lapply(int, function(x) {
        x <- x[-c(1, nrow(x)), ]
        colnames(x) <- c("SI_LONG", "SI_LATI")
        return(x)
    })
    intidx <- lapply(interpolationEQ[idx], function(x) {
        x[1, ]
    })

    #- Make a unique ID per interpolation
    INT_ID <- as.data.frame(do.call(rbind,intidx))
    INT_ID$INT_ID <- 1:nrow(INT_ID)
    tacsat$INT_ID <- NA
    for(iRow in 1:nrow(INT_ID))
      tacsat[INT_ID[iRow,1]:INT_ID[iRow,2],"INT_ID"] <- INT_ID$INT_ID[iRow]

    #- Aggregate tacsat information for all pings within an interpolation
    idxNoIntTacsat  <- which(is.na(tacsat$INT_ID))
    if(length(idxNoIntTacsat)>0){
      noIntTacsat     <- tacsat[idxNoIntTacsat,]
      intTacsat       <- tacsat[-idxNoIntTacsat,]

      matchClnms  <- c("INTV","LE_SURF","LE_SUBSURF",colnames(tacsat)[kgeur(colnames(tacsat))])[which(c("INTV","LE_SURF","LE_SUBSURF",colnames(tacsat)[kgeur(colnames(tacsat))]) %in% colnames(tacsat))]
      aggInttacsat<- aggregate(tacsat[-idxNoIntTacsat,matchClnms],by=list(tacsat[-idxNoIntTacsat,"INT_ID"]),FUN=sum,na.rm=T)
      colnames(aggInttacsat)[1] <- "INT_ID"
    } else {
      noIntTacsat <- numeric()
      intTacsat   <- tacsat
      matchClnms  <- c("INTV","LE_SURF","LE_SUBSURF",colnames(tacsat)[kgeur(colnames(tacsat))])[which(c("INTV","LE_SURF","LE_SUBSURF",colnames(tacsat)[kgeur(colnames(tacsat))]) %in% colnames(tacsat))]
      aggInttacsat<- aggregate(tacsat[,matchClnms],by=list(tacsat[,"INT_ID"]),FUN=sum,na.rm=T)
      colnames(aggInttacsat)[1] <- "INT_ID"
    }

    #- get the column names that cannot be interpolated
    clnames <- colnames(tacsat)
    b <- clnames[!clnames %in% c("SI_LONG", "SI_LATI",
            "SI_HE", "SI_SP", "SI_DATE", "SI_TIME",
            "SI_DATIM")]

    #- Replicate values that cannot be interpolated
    bvals <- lapply(as.list(1:lenEQ), function(x) {
        matrix(unlist(tacsat[intidx[[x]][1], b]), nrow = npoints -
            2, ncol = length(b), byrow = TRUE, dimnames = list(round(seq(intidx[[x]][1],
            intidx[[x]][2], length.out = npoints - 2), 3), b))
    })
    #- Apply a simple equal speading to the data values (Ideally we would make this dependent on SI_SP but then the equalDist should be modified too)
    bvals <- lapply(as.list(1:lenEQ),function(x){
              dat2replace <- subset(aggInttacsat,INT_ID == bvals[[x]][1,"INT_ID"])
              if(length(kgeur(b))>0){
                for(iSpec in b[kgeur(b)])
                  bvals[[x]][,iSpec] <- c(dat2replace[,iSpec]/(npoints-1))
                mode(bvals[[x]][,kgeur(b)]) <- "numeric"
              }
              if("INTV" %in% b){
                bvals[[x]][,"INTV"] <- c(dat2replace[,"INTV"]/(npoints-1))
                mode(bvals[[x]][,"INTV"]) <- "numeric"
              }
              if("LE_SURF" %in% b){
                bvals[[x]][,"LE_SURF"] <- c(dat2replace[,"LE_SURF"]/(npoints-1))
                mode(bvals[[x]][,"LE_SURF"]) <- "numeric"
              }
              if("LE_SUBSURF" %in% b){
                bvals[[x]][,"LE_SUBSURF"] <- c(dat2replace[,"LE_SUBSURF"]/(npoints-1))
                mode(bvals[[x]][,"LE_SUBSURF"]) <- "numeric"
              }
              return(bvals[[x]])
    })

    SI_DATIMs <- data.frame(from = tacsat$SI_DATIM[do.call(rbind,
        intidx)[, 1]], to = tacsat$SI_DATIM[do.call(rbind, intidx)[,
        2]])
    SI_DATIMs <- lapply(as.list(1:lenEQ), function(x) {
        seq(SI_DATIMs[x, 1], SI_DATIMs[x, 2], length.out = npoints)[2:(npoints -
            1)]
    })
    SI_DATE <- lapply(SI_DATIMs, function(x) {
        format(x, format = "%d/%m/%Y")
    })
    timeNotation <- ifelse(length(unlist(strsplit(tacsat$SI_TIME[1],
        ":"))) > 2, "secs", "mins")
    if (timeNotation == "secs")
        SI_TIME <- lapply(SI_DATIMs, function(x) {
            format(x, format = "%H:%M:%S")
        })
    if (timeNotation == "mins")
        SI_TIME <- lapply(SI_DATIMs, function(x) {
            format(x, format = "%H:%M")
        })
    SI_SPs <- as.matrix(data.frame(from = tacsat$SI_SP[do.call(rbind,
        intidx)[, 1]], to = tacsat$SI_SP[do.call(rbind, intidx)[,
        2]]))
    SI_SP <- mapply(seq, from = SI_SPs[, 1], to = SI_SPs[, 2],
        length.out = npoints - 2)
    SI_HE <- lapply(as.list(1:lenEQ), function(x) {
        y <- int[[x]]
        return(bearing(y[2:(nrow(y) - 1), 1], y[2:(nrow(y) -
            1), 2], y[3:nrow(y), 1], y[3:nrow(y), 2]))
    })
    ret <- lapply(as.list(1:lenEQ), function(x) {
        data.frame(bvals[[x]], intmin2[[x]], SI_DATIM = SI_DATIMs[[x]],
            SI_DATE = SI_DATE[[x]], SI_TIME = SI_TIME[[x]], SI_SP = SI_SP[,
                x], SI_HE = SI_HE[[x]], stringsAsFactors = F)
    })

    #- Need to apply the same equal spreading to the original tacsat data
    matchClnms  <- c("INTV","LE_SURF","LE_SUBSURF",
                     colnames(intTacsat)[kgeur(colnames(intTacsat))])[which(c("INTV","LE_SURF","LE_SUBSURF",colnames(intTacsat)[kgeur(colnames(intTacsat))]) %in% colnames(intTacsat))]
    intTacsat[,matchClnms] <- intTacsat[,matchClnms] / (npoints-1)

    #- Combine all together again
    interpolationTot <- do.call(rbind, ret)
    interpolationTot <- formatTacsat(interpolationTot)
    if(length(noIntTacsat)>0){
      tacsatInt <- rbindTacsat(rbindTacsat(intTacsat,noIntTacsat),interpolationTot)
    } else {
      tacsatInt <- rbindTacsat(intTacsat,interpolationTot)
    }
    for(iCol in c("INTV","LE_SURF","LE_SUBSURF",colnames(tacsatInt)[kgeur(colnames(tacsatInt))]))
      tacsatInt[,iCol] <- anf(tacsatInt[,iCol])
    tacsatInt <- sortTacsat(tacsatInt)

    #- Do consistency check
    if(length(kgeur(colnames(tacsatInt)))>0){
      ratio     <- sum(tacsatInt[,kgeur(colnames(tacsatInt))],na.rm=T)/sum(tacsat[,kgeur(colnames(tacsatInt))],na.rm=T)
      if(ratio <= 0.99 | ratio >= 1.01 | is.na(ratio))
        warnings("Some kilo/euros/surface area swept is misplaced")
    }

    return(tacsatInt)
}
