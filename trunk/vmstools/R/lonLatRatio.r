`lonLatRatio` <-
function(lon,lat){
                                  x1 <- lon[1]
                                  x2 <- lon[2]
                                  y1 <- lat[1]
                                  y2 <- lat[2]
                                  
                                    #Based on the Haversine formula  
                                      #At the position, the y-position remains the same, hence, cos(y1)*cos(y1) instead of cos(y1) * cos(y2)
                                  a <- cos(y1*pi/180)*cos(y1*pi/180)*sin((0.1*pi/180)/2)*sin((0.1*pi/180)/2); 
                                      c <- 2*atan2(sqrt(a),sqrt(1-a)); 
                                      R <- 6371; 
                                      dx1 <- R*c
                                  a <- cos(y2*pi/180)*cos(y2*pi/180)*sin((0.1*pi/180)/2)*sin((0.1*pi/180)/2); 
                                      c <- 2*atan2(sqrt(a),sqrt(1-a)); 
                                      R <- 6371; 
                                      dx2 <- R*c 
                                      
                                  ratio1 <- dx1/11.12
                                  ratio2 <- dx2/11.12
                              return(c(ratio1,ratio2))}

