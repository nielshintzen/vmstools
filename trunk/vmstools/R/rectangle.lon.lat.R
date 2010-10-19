rectangle.lon.lat <- function(Rectangle, midpoint = TRUE){

                    lat  <- (as.numeric(substr(Rectangle, 1, 2)) / 2) +35.5

                    lon  <- seq(-40, 40, by=10)[match(substr(Rectangle, 3, 3), 
                               c("A", "B", "C", "D", "E", "F", "G", "H", "I"))] +
                               (as.numeric(substr(Rectangle, 4, 4)) %% 10)
                    
                    if(midpoint == TRUE){
                    
                        lat <- lat + 0.25
                        lon <- lon + 0.5
                        
                        }

                    return(cbind(lat, lon))
                    
                    }   
