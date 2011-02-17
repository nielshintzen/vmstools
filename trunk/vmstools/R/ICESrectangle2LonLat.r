ICESrectangle2LonLat <-
function (statsq,midpoint=F)
{
    part1 <- substr(statsq, 1, 2)
    part2 <- substr(statsq, 3, 4)
    labels <- 0:90
    latlabels <- ifelse(labels < 10, paste("0", labels, sep = ""), 
        as.character(labels))
    latvalues <- seq(35.5, 80.5, 0.5) + 0.25
    lonlabels <- paste(rep(LETTERS[2:8], rep(10, 7)), rep(0:9, 
        7), sep = "")
    lonvalues <- (-40:29) + 0.5
    indx <- match(part1, latlabels)
    lat <- latvalues[indx]
    indx <- match(part2, lonlabels)
    lon <- lonvalues[indx]
    if (any(is.na(lat)) | any(is.na(lon))) 
        warning("Some stat squares have not been recognised.")
        
    if(midpoint == F){

                        lat <- lat - 0.25
                        lon <- lon - 0.5

                        }
    return(data.frame(SI_LATI=lat,SI_LONG=lon))
}
