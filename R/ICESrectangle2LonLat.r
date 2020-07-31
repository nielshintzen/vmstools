ICESrectangle2LonLat <-
function (statsq,midpoint=FALSE)
{
    #Split the code into its respective parts
    latpart <- substr(statsq, 1, 2)
    lonpart <- substr(statsq, 3, 4)

    #Handle the latitude first, as its easier
    #The grid increments continuously from
    #south to north in 0.5 degree intervals
    latlabels <- sprintf("%02i",1:99)
    lat.mids <- seq(36, 85, 0.5) + 0.25
    lat.idx <- match(latpart, latlabels)
    lat <- lat.mids[lat.idx]

    #The longtitudinal structure is not so easy
    #There are three main traps:
    #  - A4-A9 do not exist
    #  - There are no I squares
    #  - Grid ends at M8
    lonlabels <- paste(rep(LETTERS[c(2:8,10:13)], each=10), rep(0:9,
        7), sep = "")
    lonlabels <- c("A0","A1","A2","A3",lonlabels[-length(lonlabels)])
    lon.mids <- -44:68 + 0.5
    lon.idx <- match(lonpart, lonlabels)
    lon <- lon.mids[lon.idx]

    #Check whether it worked
    #If any part of the code is not recognised, both
    #lon and lat should be NA
    failed.codes <- is.na(lat) | is.na(lon)
    if (any(failed.codes)) {
        warning("Some stat squares are not valid. Please check the help files for ICESrectangle2LonLat() for more information about the formal definition of valid ICES rectangles.")
        lat[failed.codes] <- NA
        lon[failed.codes] <- NA
    }

    #Correct for midpoints
    if(midpoint == FALSE){

                        lat <- lat - 0.25
                        lon <- lon - 0.5

                        }
    #Done
    return(data.frame(SI_LATI=lat,SI_LONG=lon))
}
