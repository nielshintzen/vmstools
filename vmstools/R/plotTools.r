#' Plot eflalo or tacsat files
#' 
#' Plot eflalo or tacsat given a grid and column name as a map of intensity
#' 
#' 
#' @param x Eflalo or tacsat dataframe
#' @param level Aggregating level: "ICESrectangle" or specified "gridcell"
#' between xlim and ylim in steps gricell for tacsat or "ICESrectangle" only
#' for eflalo
#' @param xlim two element numeric vector giving a range of longitudes,
#' expressed in degrees, to which drawing should be restricted. Longitude is
#' measured in degrees east of Greenwich, so that, in particular, locations in
#' the USA have negative longitude. If fill = TRUE, polygons selected by region
#' must be entirely inside the xlim range. The default value of this argument
#' spans the entire longitude range of the database.
#' @param ylim two element numeric vector giving a range of latitudes,
#' expressed in degrees, to which drawing should be restricted. Latitude is
#' measured in degrees north of the equator, so that, in particular, locations
#' in the USA have positive latitude. If fill = TRUE, polygons selected by
#' region must be entirely inside the ylim range. The default value of this
#' argument spans the entire latitude range of the database.
#' @param zlim the minimum and maximum z values for which colors should be
#' plotted, defaulting to the range of the finite values of z. Each of the
#' given colors will be used to color an equispaced interval of this range. The
#' midpoints of the intervals cover the range, so that values just outside the
#' range will be plotted.
#' @param log whether values to plot need to be logged, TRUE or FALSE
#' @param gridcell two element numeric vector giving the steps in longitudinal
#' degrees and steps in latitudinal degrees for the grid to plot on. Only
#' needed when level = "gridcell"
#' @param color colors range. default = brewer.pal(9,"YlOrRd")
#' @param control.tacsat list with attribute "clm": column names to display
#' aggregated results over
#' @param control.eflalo list with attribute "clm": column names to display
#' aggregated results over
#' @param returnRange Logical: return range of plotted values (which can be
#' used to define own legend
#' @param las Direction of axis legends. See ?par for more info
#' @author Niels T. Hintzen
#' @seealso \code{\link{plotTreeMap}}, \code{\link{Grid2KLM}},
#' \code{\link{landingsMap2GIFanim}}, \code{\link{pings2EffortMaps}},
#' \code{\link{pings2LandingsMaps}}
#' @examples
#' 
#' data(tacsat)
#' data(eflalo)
#' 
#' plotTools(tacsat,level="ICESrectangle",xlim=c(-5,10),ylim=c(48,62),zlim=NULL,
#'           log=FALSE,gridcell=c(0.1,0.05),color=NULL,control.tacsat=list(clm=NULL))
#' x11()
#' plotTools(eflalo,level="ICESrectangle",xlim=c(-5,10),ylim=c(48,62),zlim=NULL,
#'           log=FALSE,gridcell=c(0.1,0.05),color=NULL,control.tacsat=list(clm=NULL))
#' x11()
#' plotTools(tacsat,level="gridcell",xlim=c(-5,10),ylim=c(48,62),zlim=NULL,
#'           log=FALSE,gridcell=c(0.1,0.05),color=NULL,control.tacsat=list(clm=NULL))
#' 
#' x11()
#' plotTools(eflalo,level="ICESrectangle",xlim=c(-5,10),ylim=c(48,62),zlim=NULL,
#'           log=FALSE,gridcell=c(0.1,0.05),color=NULL,control.tacsat=list(clm=NULL),
#'           control.eflalo=list(clm=c("LE_KG_COD","LE_KG_PLE")))
#' 
#' 
#' @export plotTools
plotTools <- function(x,level="ICESrectangle",xlim,ylim,zlim=NULL,log=FALSE,gridcell=c(0.1,0.05),color=NULL,control.tacsat=list(clm=NULL),control.eflalo=list(clm=NULL),returnRange=FALSE,las=1){
  if(is.null(color)==TRUE) color <- rev(heat.colors(9))
  require(data.table)
  require(sf)
    #TACSAT
  if(all(c("SI_LATI","SI_LONG") %in% colnames(x))){
      #- Get left lower point from ICES rectangle in gps format
    x$LE_RECT             <- ICESrectangle(x)
    
      #- Limit the spatial boundaries of tacsat
    idxx                  <- which(x$SI_LONG >= xlim[1] & x$SI_LONG <= xlim[2])
    idxy                  <- which(x$SI_LATI >= ylim[1] & x$SI_LATI <= ylim[2])
    x                     <- x[idxx[which(idxx %in% idxy)],]
    x$SI_LONGF            <- af(ac(x$SI_LONG))
    x$SI_LATIF            <- af(ac(x$SI_LATI))
    if(is.null(control.tacsat$clm)==TRUE){
        control.tacsat$clm  <- "idx"
        x$idx               <- rep(1,nrow(x))
    }
      
    #---------------------------------------------------------------------------
    #- Sum by rectangle
    #---------------------------------------------------------------------------
    if(level == "ICESrectangle"){
      x$SI_LONGF            <- af(ICESrectangle2LonLat(ac(x$LE_RECT))[,2])
      x$SI_LATIF            <- af(ICESrectangle2LonLat(ac(x$LE_RECT))[,1])

      DT                    <- data.table(x)
      eq1                   <- c.listquote(paste("sum(",control.tacsat$clm,",na.rm=TRUE)",sep=""))
      eq2                   <- c.listquote(c("SI_LONGF","SI_LATIF"))
      
      byRect                <- data.frame(DT[,eval(eq1),by=eval(eq2)]); colnames(byRect) <- c("SI_LONG","SI_LATI",control.tacsat$clm)
      byRect                <- byRect[which(is.na(byRect$SI_LONG)==FALSE & is.na(byRect$SI_LATI) == FALSE),]
      if(length(control.tacsat$clm)>1)  rangeRect             <- range(apply(byRect[,control.tacsat$clm],1,sum,na.rm=TRUE))
      if(length(control.tacsat$clm)==1) rangeRect             <- range(byRect[,control.tacsat$clm],na.rm=TRUE)
      rangeRect             <- c(0,rangeRect[2])
    }
    #---------------------------------------------------------------------------
    #- Sum by grid cell
    #---------------------------------------------------------------------------
    if(level == "gridcell"){
      grids                 <- createGrid(xlim,ylim,gridcell[1],gridcell[2],type="GridDF")
      coords                <- st_as_sf(x,coords=c("SI_LONG","SI_LATI"))
      coords$dens           <- vmstools:::st_over(coords, st_geometry(grids))
      
        #- Sum by gridcell
      DT                    <- data.table(data.frame(coords))
      DT$x                  <- af(ac(st_coordinates(st_centroid(grids))[DT$dens,1]))
      DT$y                  <- af(ac(st_coordinates(st_centroid(grids))[DT$dens,2]))
      
      eq1                   <- c.listquote(paste("sum(",control.tacsat$clm,",na.rm=TRUE)",sep=""))
      eq2                   <- c.listquote(c("x","y"))
      
      byRect                <- data.frame(DT[,eval(eq1),by=eval(eq2)]); colnames(byRect) <- c("SI_LONG","SI_LATI",control.tacsat$clm)
      byRect$SI_LONG        <- signif(anf(byRect$SI_LONG))
      byRect$SI_LATI        <- signif(anf(byRect$SI_LATI))
      if(length(control.tacsat$clm)>1)  rangeRect             <- range(apply(byRect[,control.tacsat$clm],1,sum,na.rm=TRUE))
      if(length(control.tacsat$clm)==1) rangeRect             <- range(byRect[,control.tacsat$clm],na.rm=TRUE)
      rangeRect             <- c(0,rangeRect[2])
    }
    ctrl                    <- control.tacsat
  }


    #EFLALO
  if(all(c("VE_FLT","VE_KW") %in% colnames(x))){
    if(level != "ICESrectangle") stop("Data supplied is 'eflalo' which only has ICESrectangle")
    x$SI_LONG             <- ICESrectangle2LonLat(ac(x$LE_RECT))[,2]
    x$SI_LATI             <- ICESrectangle2LonLat(ac(x$LE_RECT))[,1]
    
      #- Limit the spatial boundaries of eflalo
    idxx                  <- which(x$SI_LONG >= xlim[1] & x$SI_LONG <= xlim[2])
    idxy                  <- which(x$SI_LATI >= ylim[1] & x$SI_LATI <= ylim[2])
    x                     <- x[idxx[which(idxx %in% idxy)],]
    
      #- Sum by rectangle
    if(is.null(control.eflalo$clm)==TRUE) control.eflalo$clm <- colnames(x[,kgeur(colnames(x))])
    DT                    <- data.table(x)
    eq1                   <- c.listquote(paste("sum(",control.eflalo$clm,",na.rm=TRUE)",sep=""))
    eq2                   <- c.listquote(c("SI_LONG","SI_LATI"))
    DT$SI_LONG            <- af(ac(DT$SI_LONG)); DT$SI_LATI <- af(ac(DT$SI_LATI))
    
    byRect                <- data.frame(DT[,eval(eq1),by=eval(eq2)]); colnames(byRect) <- c("SI_LONG","SI_LATI",control.eflalo$clm)
    byRect                <- byRect[which(is.na(byRect$SI_LONG)==FALSE & is.na(byRect$SI_LATI) == FALSE),]
    if(length(control.eflalo$clm)>1)  rangeRect             <- range(apply(byRect[,control.eflalo$clm],1,sum,na.rm=TRUE))
    if(length(control.eflalo$clm)==1) rangeRect             <- range(byRect[,control.eflalo$clm],na.rm=TRUE)
    rangeRect             <- c(0,rangeRect[2])
    ctrl                  <- control.eflalo
  }

  data(europa)
  plt <- list()
  cols <- rep(NA,nrow(byRect))
  for(iRect in 1:nrow(byRect)){
    if(log){
      if(is.null(zlim)==TRUE){ i         <- round((log(sum(byRect[iRect,ctrl$clm],na.rm=TRUE))-ifelse(rangeRect[1]==0,0,log(rangeRect[1])))
                                               /(log(rangeRect[2]) - ifelse(rangeRect[1]==0,0,log(rangeRect[1])))*(length(color)-1)) +1
      } else {
          i                           <- round((log(sum(byRect[iRect,ctrl$clm],na.rm=TRUE))-ifelse(zlim[1]==0,0,log(zlim[1])))
                                               /(log(zlim[2]) -      ifelse(zlim[1]==0,0,log(zlim[1])))          *(length(color)-1)) +1
        }
    } else {
        if(is.null(zlim)==TRUE){ i       <- round((sum(byRect[iRect,ctrl$clm],na.rm=TRUE)-ifelse(rangeRect[1]==0,0,rangeRect[1]))
                                               /(rangeRect[2] -      ifelse(rangeRect[1]==0,0,rangeRect[1]))     *(length(color)-1)) +1
        } else {
            i                         <- round((sum(byRect[iRect,ctrl$clm],na.rm=TRUE)-ifelse(zlim[1]==0,0,zlim[1]))
                                               /(zlim[2] -           ifelse(zlim[1]==0,0,zlim[1]))               *(length(color)-1)) +1
          }
      }
    cols[iRect]                       <- color[i]
    if(level == "ICESrectangle") plt[[iRect]] <-  lonLat2SFPolygons(SI_LONG=c(an(ac(byRect[iRect,"SI_LONG"])),an(ac(byRect[iRect,"SI_LONG"]))+1,an(ac(byRect[iRect,"SI_LONG"]))+1,an(ac(byRect[iRect,"SI_LONG"])),an(ac(byRect[iRect,"SI_LONG"]))),
                                                                               SI_LATI=c(rep(an(ac(byRect[iRect,"SI_LATI"])),2),rep(an(ac(byRect[iRect,"SI_LATI"]))+0.5,2),an(ac(byRect[iRect,"SI_LATI"]))))#,col=color[i],lwd=1,border=NA)
    if(level == "gridcell") plt[[iRect]]      <-  lonLat2SFPolygons(SI_LONG=c(an(ac(byRect[iRect,"SI_LONG"]))-gridcell[1]/2,an(ac(byRect[iRect,"SI_LONG"]))+gridcell[1]/2,an(ac(byRect[iRect,"SI_LONG"]))+gridcell[1]/2,an(ac(byRect[iRect,"SI_LONG"]))-gridcell[1]/2,an(ac(byRect[iRect,"SI_LONG"]))-gridcell[1]/2),
                                                                               SI_LATI=c(rep(an(ac(byRect[iRect,"SI_LATI"])),2)-gridcell[2]/2,rep(an(ac(byRect[iRect,"SI_LATI"]))+gridcell[2]/2,2),rep(an(ac(byRect[iRect,"SI_LATI"])),1)-gridcell[2]/2))#,col=color[i],lwd=1,border=NA)
  }

  pols <- st_as_sfc(plt)
  st_crs(pols) <- 4326

  require(ggplot2)
  require(ggthemes)
  theme_plot <- 
  (theme_foundation(base_size=12)
    + theme(plot.title       = element_text(face = "bold",size = rel(1.0), hjust = 0),
            plot.margin      = unit(c(5,2,2,2),"mm"),
            plot.background  = element_rect(colour = NA),
            text             = element_text(),
            axis.title       = element_text(face = "bold",size = rel(1)),
            axis.title.y     = element_text(angle=90,vjust =2),
            axis.title.x     = element_text(vjust = -0.2),
            axis.text        = element_text(), 
            axis.line        = element_line(colour="black"),
            axis.ticks       = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            panel.border     = element_rect(colour="black" , linewidth=0.1),
            panel.background = element_rect(colour = NA),
            strip.background = element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text       = element_text(face="bold"),
            legend.key       = element_rect(colour = NA),
            legend.position  = "bottom",
            legend.direction = "horizontal",
            legend.key.size  = unit(0.2, "cm"),
            legend.spacing   = unit(0, "cm"),  # updated from legend.margin which is deprecated
            legend.title     = element_text(face="italic", size=rel(0.8))
    ))
  print(ggplot() +
    geom_sf(data=europa) + theme_plot + 
    geom_sf(data=pols,colour=NA,fill=cols) +
    geom_sf(data=europa) +
    scale_x_continuous(limits = xlim) + 
    scale_y_continuous(limits = ylim))

  if(returnRange) return(rangeRect)
}
