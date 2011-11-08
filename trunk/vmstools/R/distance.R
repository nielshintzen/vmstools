`distance` <-
function(lon,lat,lonRef,latRef){

                    pd <- pi/180

                    a1<- sin(((latRef-lat)*pd)/2)
                    a2<- cos(lat*pd)
                    a3<- cos(latRef*pd)
                    a4<- sin(((lonRef-lon)*pd)/2)
                    a <- a1*a1+a2*a3*a4*a4

                                      c <- 2*atan2(sqrt(a),sqrt(1-a));
                    return(6371*c)}