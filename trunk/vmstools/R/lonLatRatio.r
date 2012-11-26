`lonLatRatio` <-
    function(x1,y1){
      #Based on the Haversine formula
      #At the position, the y-position remains the same, hence, cos(y1)*cos(y1) instead of cos(y1) * cos(y2)
      a <- cos(y1*pi/180)*cos(y1*pi/180)*sin((0.1*pi/180)/2)*sin((0.1*pi/180)/2);
      c <- 2*atan2(sqrt(a),sqrt(1-a));
      R <- 6371;
      dx1 <- R*c

    return(c(dx1/11.12))}

