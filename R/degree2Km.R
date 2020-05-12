`degree2Km` <-
function(lon,lat,degree){
                      x1 <- lon
                      y1 <- lat
                      
                      a <- cos(y1*pi/180)*cos(y1*pi/180)*sin((1*pi/180)/2)*sin((1*pi/180)/2);
                        c <- 2*atan2(sqrt(a),sqrt(1-a));
                        R <- 6371;
                        dx1 <- R*c
                        
              return(dx1 * degree)}

