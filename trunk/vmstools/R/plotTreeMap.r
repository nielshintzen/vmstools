
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
# F.Bastardie
# using stolen code on the net
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!##

plotTreeMap <- 
    function(x,  gridcell=c(0.1,0.1), gear="OTB", 
                xlim= c(-1,17), ylim= c(52,62),  acolors = rainbow(7),
                  species.to.keep= c("LE_KG_COD","LE_KG_NEP", "LE_KG_PLE", "LE_KG_SOL") ){  
    
   
   chop <- function(x) rev( rev(x)[-1] )

   simple.hook <- function(z, xl, yl, xu, yu) {
     #rect(xl, yl, xu, yu, lwd=3, border="blue")
     #text((xl+xu)/2, (yl+yu)/2, z$one, cex=2.5, col="green")
     rect(xl, yl, xu, yu,  col=as.character(z$one), border=NA)
     }

   squarified.treemap <- function(z, x=0, y=0, w=1, h=1, hook) {
     cz <- cumsum(z$size) / sum(z$size)
     n <- which.min(abs(log(max(w/h, h/w) * sum(z$size) * cz^2/z$size)))
     more <- n < length(z$size)
     a <- c(0, cz[1:n]) / cz[n]
     if (h > w) {
       hook(z[1:n, ], x+w*chop(a), rep(y,n), x+w*a[-1], rep(y+h*cz[n],n))
       if (more) Recall(z[-(1:n), ], x, y+h*cz[n], w, h*(1-cz[n]), hook)
     } else {
       hook(z[1:n, ], rep(x,n), y+h*chop(a), rep(x+w*cz[n],n), y+h*a[-1])
       if (more) Recall(z[-(1:n), ], x+w*cz[n], y, w*(1-cz[n]), h, hook)
     }    
    }
   ## test squarified.treemap()
   if(FALSE){
        z <- data.frame(size=c(5,8,3,12,2,6,7), one=colors()[100:106])
        z <- z[order(-z$size), ]
        plot(0:2, 0:2, type="n", axes=FALSE, xlab="", ylab="")
      squarified.treemap(z, x=1, y=1, hook=simple.hook)
   }


    x<-  x[x$SI_STATE==1,]  # keep only fishing points
  
    
    x$SI_LONG             <- anf(x$SI_LONG)
    x$SI_LATI             <- anf(x$SI_LATI)
    x <- x[!is.na(x$SI_LATI),]

    
    idxx                  <- which(x$SI_LONG >= xlim[1] & x$SI_LONG <= xlim[2])
    idxy                  <- which(x$SI_LATI >= ylim[1] & x$SI_LATI <= ylim[2])
    x                     <- x[idxx[which(idxx %in% idxy)],]
   
    xx <- x[x$LE_GEAR %in% gear,] # keep only a given (set of) gear(s)
    
    grids                 <- createGrid(xrange=range(xx$SI_LONG),yrange=range(xx$SI_LATI), gridcell[1],gridcell[2],type="SpatialPixelsDataFrame")
    coords                <- SpatialPointsDataFrame(cbind(x=an(ac(xx$SI_LONG)),y=an(ac(xx$SI_LATI))),data=xx)
    coords$dens           <- overlay(grids, coords)
      
    #- Sum by gridcell
    DT                    <- data.table(data.frame(coords))
    DT$x                  <- af(ac(grids@coords[DT$dens,1]))
    DT$y                  <- af(ac(grids@coords[DT$dens,2]))
      
    idx.col               <- grep('KG', names(coords))
    eq1                   <- c.listquote(paste("sum(",names(coords[,idx.col]),",na.rm=TRUE)",sep=""))
    eq2                   <- c.listquote(c("x","y"))
        
    byRect                <- data.frame(DT[,eval(eq1),by=eval(eq2)])
    colnames(byRect)      <- c("SI_LONG","SI_LATI", names(coords)[idx.col])
    byRect$SI_LONG        <- signif(anf(byRect$SI_LONG))
    byRect$SI_LATI        <- signif(anf(byRect$SI_LATI))
    idx.col               <- grep('KG', colnames(byRect))
    rangeRect             <- range(apply(byRect[idx.col],1,sum,na.rm=T))
    rangeRect             <- c(0,rangeRect[2])
    
    A.sum <- apply(byRect[,idx.col], 1, sum, na.rm=TRUE) 
    A.sum2 <- apply(byRect[,idx.col], 2, sum, na.rm=TRUE) 
   
 
    # create OTH
    species.to.merge <- names(A.sum2)[!names(A.sum2) %in% species.to.keep]
    byRect$LE_KG_OTH <- apply(byRect[,species.to.merge], 1, sum, na.rm=TRUE)
    
    # compute the proportion
    byRect   <-byRect[, !names(byRect) %in%  species.to.merge] 
    idx.col               <- grep('KG', names(byRect))
    byRect[,idx.col] <- sweep(byRect[,idx.col], 1, A.sum, FUN="/")
   
    
    X11(15,15)
    library(mapdata)
    map("worldHires",res=1,xlim=range(x$SI_LONG),ylim=range(x$SI_LATI),fill=T,col="darkgreen");map.axes();box()
   
 
    for(iRect in 1:nrow(byRect)){
  
      x1 <- an(ac(byRect[iRect,"SI_LONG"]))
      y1 <- an(ac(byRect[iRect,"SI_LATI"]))

      size <- an(ac(byRect[iRect, idx.col]))
      size <- replace(size, is.na(size) | size==0,1)
  
      z <- data.frame(size=size, one=acolors[1:((1+length(idx.col))-1)])
      z <- z[order(-z$size), ]
      squarified.treemap(z, x=x1, y=y1, w=gridcell[1], h=gridcell[2], hook=simple.hook)

      }
    map("worldHires",add=TRUE,res=1,xlim=range(x$SI_LONG),ylim=range(x$SI_LATI),fill=T,col="darkgreen");map.axes();box()
    legend("topright", legend=names(byRect[,idx.col]), fill=acolors[1:((1+length(idx.col))-1)])
  
 
 
    if(FALSE){
    color <- rainbow(9)
    library(mapdata)
    map("worldHires",res=1,xlim=range(x$SI_LONG),ylim=range(x$SI_LATI),fill=T,col="darkgreen");map.axes();box()
    for(iRect in 1:nrow(byRect)){
         i       <- round((sum(byRect[iRect,idx.col],na.rm=T)-ifelse(rangeRect[1]==0,0,rangeRect[1])) /
                       (rangeRect[2] - ifelse(rangeRect[1]==0,0,rangeRect[1]))     *(length(color)-1)) +1
 
    polygon(x=c(an(ac(byRect[iRect,"SI_LONG"])),
             an(ac(byRect[iRect,"SI_LONG"]))+gridcell[1],
             an(ac(byRect[iRect,"SI_LONG"]))+gridcell[1],
             an(ac(byRect[iRect,"SI_LONG"]))
             ),
          y=c(an(ac(byRect[iRect,"SI_LATI"])),
              an(ac(byRect[iRect,"SI_LATI"])),
              an(ac(byRect[iRect,"SI_LATI"]))+gridcell[2],
              an(ac(byRect[iRect,"SI_LATI"]))+gridcell[2]
             ),
          col=color[i],lwd=1,border=NA)         
     }
    } # end FALSE  

return()
}


##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##


   # load("C:\\merging\\redo3merged2010\\all_merged_weight_2010.RData")
   # load("C:\\merging\\redo3merged2010\\all_merged_value_2010.RData")
   # library(vmstools)
   # plotTreeMap (all.merged,  gridcell=c(0.1,0.1), gear="OTB", 
   #             xlim= c(-1,17), ylim= c(52,62), acolors=rainbow(7),
   #               species.to.keep= c("LE_KG_COD","LE_KG_NEP", "LE_KG_PLE", "LE_KG_SOL") ) 
 



 