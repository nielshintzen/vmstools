`distance` <-
function(lon,lat,lonRef,latRef){
                    x1 <- lon
                    y1 <- lat
                    x2 <- lonRef
                    y2 <- latRef
                    
                    a <- sin(((y2-y1)*pi/180)/2)*sin(((y2-y1)*pi/180)/2) + cos(y1*pi/180)*cos(y2*pi/180)*
                         sin(((x2-x1)*pi/180)/2)*sin(((x2-x1)*pi/180)/2);
                                      c <- 2*atan2(sqrt(a),sqrt(1-a));
                                      R <- 6371;
                                      dx1 <- R*c
                    res <- dx1
                    return(res)}