
## Francois Bastardie - fba@aqua.dtu.dk
## Jan 2012

## -------! paths and data inputs (e.g. eflalo and tacsat) to be adapted by the user !------


library(vmstools)  # version 0.58
library(help=vmstools)  # getting help...
memory.limit(4000)

# set your own path here:
main.path <- "C:\\merging"

  ###--------------------------------###
  ###--------------------------------###
  ###--------------------------------###
  ### LOAD TACSAT & EFLALO DATA SETS ###
  ###--------------------------------###
  ###--------------------------------###
  ###--------------------------------###
  y <- "2010"
  load(paste(main.path,"\\VMStools_DTUAqua_course_Jan12\\EflaloAndTacsat\\tacsat_10vessels_",y,".RData",sep=''))
  load(paste(main.path,"\\VMStools_DTUAqua_course_Jan12\\EflaloAndTacsat\\eflalo_10vessels_",y,".RData",sep=''))
 
  
  ###--------------------------------###
  ###--------------------------------###
  ###--------------------------------###
  ### CREATE A PORT DATA.FRAME       ###
  ###--------------------------------###
  ###--------------------------------###
  ###--------------------------------###
  EUports <- read.table(file=file.path(main.path,
                     "EU_ports_location.csv"), header=TRUE,sep=",") 
                     # from http://ec.europa.eu/fisheries/cfp/control/codes/index_en.htm
  plot(EUports$Longitude, EUports$Latitude, xlim=c(-10,20), ylim=c(50,70))
  library(mapdata)
  map("worldHires", add=TRUE)
  
  if(FALSE){
  ## add missing harbour found from suspicious lat/long
  missing.ports <- EUports[1:41,] ; missing.ports[]<- NA  #init
 
  # long/lat
  coord.missing <- data.frame(NA, ncol=3, nrow=nrow(missing.ports))
  coord.missing[1,1] <- 10.3087 ; coord.missing[1,2] <- 56.9913    ; coord.missing[1,3] <- "DNK"
  coord.missing[2,1] <-  10.6720 ; coord.missing[2,2] <- 54.7507      ; coord.missing[2,3] <-  "DNK"
  coord.missing[3,1] <-  9.7120 ; coord.missing[3,2] <- 55.2627      ; coord.missing[3,3] <- "DNK"
  coord.missing[4,1] <- 6.2040 ; coord.missing[4,2] <- 53.4080      ; coord.missing[4,3] <-  "NLD"
  coord.missing[5,1] <-  11.9293 ; coord.missing[5,2] <- 54.5727         ; coord.missing[5,3] <- "DNK"
  coord.missing[6,1] <-   12.0420 ; coord.missing[6,2] <- 54.8920      ; coord.missing[6,3] <-  "DNK"
  coord.missing[7,1] <- 11.198 ; coord.missing[7,2] <-  55.2113     ; coord.missing[7,3] <- "DNK"
  coord.missing[8,1] <-  -1.1440 ; coord.missing[8,2] <- 60.1707      ; coord.missing[8,3] <-  "DNK"
  coord.missing[9,1] <-10.7127  ; coord.missing[9,2] <- 56.5333      ; coord.missing[9,3] <- "DNK"
  coord.missing[10,1] <-  10.9253 ; coord.missing[10,2] <-57.2967   ; coord.missing[10,3] <- "DNK"
  coord.missing[11,1] <- 10.9253 ; coord.missing[11,2] <-57.2960    ; coord.missing[11,3] <- "DNK"
  coord.missing[12,1] <- 10.6713 ; coord.missing[12,2] <- 54.752   ; coord.missing[12,3] <-  "DNK"
  coord.missing[13,1] <-  10.2407 ; coord.missing[13,2] <-55.0940    ; coord.missing[13,3] <-  "DNK"
  coord.missing[14,1] <- 10.6713 ; coord.missing[14,2] <- 54.7507    ; coord.missing[14,3] <-  "DNK"
  coord.missing[15,1] <-  10.6713 ; coord.missing[15,2] <- 54.7513    ; coord.missing[15,3] <-  "DNK"
  coord.missing[16,1] <- 9.9880  ; coord.missing[16,2] <- 54.6853    ; coord.missing[16,3] <-   "DNK"
  coord.missing[17,1] <- -4.3267 ; coord.missing[17,2] <- 48.0993      ; coord.missing[17,3] <- "GBR"
  coord.missing[18,1] <-  9.965 ; coord.missing[18,2] <-57.594        ; coord.missing[18,3] <-  "DNK"
  coord.missing[19,1] <-8.304  ; coord.missing[19,2] <- 56.55      ; coord.missing[19,3] <-    "DNK"
  coord.missing[20,1] <- 9.63  ; coord.missing[20,2] <- 55.513      ; coord.missing[20,3] <- "DNK"
  coord.missing[21,1] <- 10.052 ; coord.missing[21,2] <-55.822          ; coord.missing[21,3] <- "DNK"
  coord.missing[22,1] <- 9.783 ; coord.missing[22,2] <- 54.913      ; coord.missing[22,3] <- "DNK"
  coord.missing[23,1] <- 15.137 ; coord.missing[23,2] <-  55.058         ; coord.missing[23,3] <-  "DNK"
  coord.missing[24,1] <- -0.066 ; coord.missing[24,2] <-53.581      ; coord.missing[24,3] <-  "DNK"
  coord.missing[25,1] <-  8.565 ; coord.missing[25,2] <- 55.088         ; coord.missing[25,3] <- "DNK"
  coord.missing[26,1] <-  10.588 ; coord.missing[26,2] <-57.718        ; coord.missing[26,3] <- "DNK"
  coord.missing[27,1] <-10.2585  ; coord.missing[27,2] <-54.942           ; coord.missing[27,3] <-"DNK"
  coord.missing[28,1] <-  12.467 ; coord.missing[28,2] <- 54.953        ; coord.missing[28,3] <- "DNK"
  coord.missing[29,1] <-  8.565 ; coord.missing[29,2] <- 55.088           ; coord.missing[29,3] <- "DNK"
  coord.missing[30,1] <-  12.375 ; coord.missing[30,2] <-55.253          ; coord.missing[30,3] <- "DNK"
  coord.missing[31,1] <- 8.223  ; coord.missing[31,2] <- 56.703        ; coord.missing[31,3] <-"DNK"
  coord.missing[32,1] <-  8.219 ; coord.missing[32,2] <- 56.698      ; coord.missing[32,3] <-  "DNK"
  coord.missing[33,1] <-   8.222; coord.missing[33,2] <- 56.697       ; coord.missing[33,3] <-  "DNK"
  coord.missing[34,1] <- 10.502 ; coord.missing[34,2] <- 57.494        ; coord.missing[34,3] <- "DNK"
  coord.missing[35,1] <- 14.692 ; coord.missing[35,2] <- 55.094       ; coord.missing[35,3] <-   "DNK"
  coord.missing[36,1] <- 8.5485 ; coord.missing[36,2] <- 56.583       ; coord.missing[36,3] <-  "DNK"
  coord.missing[37,1] <- 9.8885 ; coord.missing[37,2] <-55.271       ; coord.missing[37,3] <-  "DNK"
  coord.missing[38,1] <- 11.709 ; coord.missing[38,2] <-55.72            ; coord.missing[38,3] <- "DNK"
  coord.missing[39,1] <- 11.35 ; coord.missing[39,2] <-55.88          ; coord.missing[39,3] <-  "DNK"
  coord.missing[40,1] <- 14.835 ; coord.missing[40,2] <-55.249           ; coord.missing[40,3] <- "DNK"
  coord.missing[41,1] <-10.6645 ; coord.missing[41,2] <-55.449           ; coord.missing[41,3] <- "DNK"
 
  for (i in 1 : nrow(coord.missing)){
      missing.ports[i,"Longitude"] <-coord.missing[i,1];  missing.ports[i,"Latitude"] <-coord.missing[i,2] ; 
      missing.ports [i,"ISO.3.Country.Code"] <- coord.missing[i,3]  ;  missing.ports [i,"Description"] <- paste("NA",i,sep='')
      points(missing.ports$Longitude[i], missing.ports$Latitude[i], xlim=c(-1,20), ylim=c(50,65), pch=16, cex=0.5, col=2)
      ## browser()
  }

  EUports <- rbind ( EUports, missing.ports)
  } # end FALSE
  
  colnames(EUports) [colnames(EUports) %in% "Latitude"]  <- "lat"
  colnames(EUports) [colnames(EUports) %in% "Longitude"] <- "lon"
  
  EUports$range <- 3 # in km
  
  
  ###--------------------------------###
  ###--------------------------------###
  ###--------------------------------###
  ### DETECT POSITIONS IN HARBOUR    ###
  ###--------------------------------###
  ###--------------------------------###
  ###--------------------------------###
  tacsat$SI_HARB <- NA
  EUports$harbour  <- EUports$Description
  tacsat$SI_HARB   <- pointInHarbour(lon=tacsat$SI_LONG, lat=tacsat$SI_LATI, harbours=EUports, rowSize=30, returnNames=TRUE)
  inHarb <- tacsat$SI_HARB
  inHarb <- replace(inHarb, !is.na(inHarb), 1)
  inHarb <- replace(inHarb, is.na(inHarb), 0)
  inHarb <- as.numeric(inHarb)

  # assign a trip identifier
  tacsat$SI_FT <- 1 # init
  idx <- which(inHarb==0)
  tacsat[idx,"SI_FT"] <- cumsum(inHarb) [idx] # add a SI_FT index

  # keep 'out of harbour' points only
  # (but keep the departure point and the arrival point lying in the harbour)
  startTrip <- c(diff(tacsat[,"SI_FT"]), 0)
  endTrip   <- c(0, diff(tacsat[,"SI_FT"]))
  tacsat[which(startTrip>0),"SI_FT"]  <-  tacsat[which(startTrip>0)+1,"SI_FT"] # tricky here
  tacsat[which(endTrip<0),"SI_FT"]    <-  tacsat[which(endTrip<0)-1,"SI_FT"] # tricky here
  tacsat <- tacsat[which(inHarb==0 |  startTrip>0 |  endTrip<0),]
    
  ###--------------------------------###
  ###--------------------------------###
  ###--------------------------------###
  ### ASSIGN A FISHING STATE         ###
  ###--------------------------------###
  ###--------------------------------###
  ###--------------------------------###
  tacsat$SI_STATE <- 2 # init (1: fishing; 2: steaming)
  tacsat$SI_STATE [(tacsat$SI_SP>4 & tacsat$SI_SP<8)] <-1
    ## => fill in SI_STATE with a fake rule just to initiate
 
 
  ###--------------------------------###
  ###--------------------------------###
  ###--------------------------------###
  ### FEW ADJUSTEMENTS...            ###
  ###--------------------------------###
  ###--------------------------------###
  ###--------------------------------###
  eflalo$LE_MET_level6 <- eflalo$LE_MET
  eflalo <- eflalo[eflalo$LE_MET!="No_logbook6",]
  eflalo$VE_FLT<-"fleet1"



  ###--------------------------------###
  ###--------------------------------###
  ###--------------------------------###
  ### VESSEL BY VESSEL MERGING       ###
  ###--------------------------------###
  ###--------------------------------###
  ###--------------------------------###
  mergeEflalo2Pings (eflalo=eflalo, tacsat=tacsat, 
          general=list(output.path=file.path(main.path, "VMStools_DTUAqua_course_Jan12", paste('merged','2010',sep='')), 
                visual.check=TRUE, detectFishing=TRUE, speed="segment", what.speed="calculated", conserve.all=TRUE))

  ###--------------------------------###
  ###--------------------------------###
  ###--------------------------------###
  ### LOOK AT ONE GIVEN VESSEL       ###
  ###--------------------------------###
  ###--------------------------------###
  ###--------------------------------###
  load(paste(main.path, "\\VMStools_DTUAqua_course_Jan12\\merged2010\\merged_DNK000003612_2010.RData",sep=''))
  head(merged[,c("VE_REF", "FT_REF", "VE_FLT", "VE_KW", "LE_MET_level6", "LE_GEAR",
                 "SI_LATI", "SI_LONG", "SI_SP", "SI_HE", "SI_HARB","SI_STATE","SI_RECT","LE_EFF_VMS","flag", 
                 "SI_DATE","SI_TIME",
                 "LE_KG_COD","LE_KG_DAB")], 40)


  
  ###--------------------------------###
  ###--------------------------------###
  ###--------------------------------###
  ### BUILD ONE BIG DATASET          ###
  ###--------------------------------###
  ###--------------------------------###
  ###--------------------------------###
  
  # ...selecting relevant species:
  spp <- c('COD', 'CSH', 'DAB', 'ELE', 'FLE', 'HAD', 'HER', 'HKE', 'HOM',
         'LEM', 'MAC', 'MON', 'MUS', 'NEP', 'NOP', 'OYF', 'PLE', 'POK', 'PRA', 'SAN',
           'SOL', 'SPR', 'TUR', 'WHB', 'WIT')

                          
  load(paste(main.path,"\\VMStools_DTUAqua_course_Jan12\\EflaloAndTacsat\\tacsat_10vessels_",y,".RData",sep=''))
  vnames <- as.character(unique(tacsat$VE_REF))
   

  tmp <- bindAllMergedTables (vessels=vnames, a.year=y, species.to.keep=spp, 
             folder = file.path(main.path, "VMStools_DTUAqua_course_Jan12", paste('merged','2010',sep='')), all.in.one.table=FALSE)
   ##=> all_merged_value_2010.RData and all_merged_weight_2010.RData are now created in the output folder...

  
  
  ###--------------------------------###
  ###--------------------------------###
  ###--------------------------------###
  ### PLOTTING ON A GRID             ###
  ###--------------------------------###
  ###--------------------------------###
  ###--------------------------------###

  a.year <- '2010'

  # 1:  fishing effort in hours
  load(file.path(main.path, "VMStools_DTUAqua_course_Jan12",
     paste("merged",a.year,sep=''), paste("all_merged_weight_",a.year,".RData",sep='')))
  all.merged$SI_LONG <- anf(all.merged$SI_LONG ) ; all.merged$SI_LATI <- anf(all.merged$SI_LATI )
  all.merged$LE_MET <- all.merged$LE_MET_level6

  # assign a trip identifier
  all.merged$tripnum <- 0 # init
  all.merged[all.merged$SI_HARB!="NA", "tripnum"] <- 1
  all.merged[, "tripnum"]  <- cumsum( all.merged[, "tripnum"]  )
  # keep only fishing positions (and the first point in the harbour)
  all.merged <- all.merged[all.merged$SI_STATE==1 | !all.merged$SI_HARB=="NA",]
  all.merged$SI_LONG <- anf(all.merged$SI_LONG)
  all.merged$SI_LATI <- anf(all.merged$SI_LATI)
  all.merged[,"LE_EFF_VMS"] <- as.numeric(as.character(all.merged[,"LE_EFF_VMS"] ))
  all.merged[, "LE_EFF_VMS"] <- replace(all.merged[,"LE_EFF_VMS"], is.na(all.merged[, "LE_EFF_VMS"]) | all.merged[, "LE_EFF_VMS"] < 0, 0)
  all.merged$LE_EFF_VMS <-   an(all.merged$LE_EFF_VMS)/ 60 # in hours
  dnk_area1 <- vmsGridCreate(all.merged, nameVarToSum="LE_EFF_VMS",  numCats=10,  plotPoints =FALSE, legendtitle="fishing (hours)",
          colLand="darkolivegreen4",  addICESgrid=TRUE,
              nameLon="SI_LONG", nameLat="SI_LATI", cellsizeX =0.05, cellsizeY =0.05, we = -3.5,ea = 3.5, so = 57, no = 61,
               breaks0=c(0,0.5,1,2,4,6,12,24,32, 100000), legendncol=2,
                outGridFile= paste("dnk",a.year,"_area1.asc",sep='') )
  rect(1.005,59.32,1.20,59.37, border=4, lwd=1.7)
   #59.32 N to 59.37N and 1.00 E to 1.20E
   # save
  #savePlot(filename = paste("dnk",a.year,"_area1",sep=''), type ="jpeg")

  # 2: A ZOOM
  dnk_area1_zoom <- vmsGridCreate(all.merged, nameVarToSum="LE_EFF_VMS",  numCats=10,  plotPoints =FALSE, legendtitle="fishing (hours)",
          colLand="darkolivegreen4",  addICESgrid=TRUE,
              nameLon="SI_LONG", nameLat="SI_LATI", cellsizeX =0.05, cellsizeY =0.05, we = -2,ea = 2, so = 58, no = 60.5,
               breaks0=c(0,2,4,6,12,24,32, 100000), legendncol=2,
                outGridFile= paste("dnk",a.year,"_area1_zoom.asc",sep='') )
  rect(1.005,59.32,1.20,59.37, border=4, lwd=1.7)
  #savePlot(filename =  paste("dnk",a.year,"_area1_zoom",sep=''), type ="jpeg")

  ###--------------------------------###
  ###--------------------------------###
  ###--------------------------------###
  ### PLOTTING THE VMS TRACKS        ###
  ###--------------------------------###
  ###--------------------------------###
  ###--------------------------------###

  # 3: PLOT THE TRACKS
  load(file.path(main.path, "VMStools_DTUAqua_course_Jan12",
     paste("merged",a.year,sep=''), paste("all_merged_weight_",a.year,".RData",sep='')))
  all.merged$SI_LONG <- anf(all.merged$SI_LONG ) ; all.merged$SI_LATI <- anf(all.merged$SI_LATI )
   # assign a trip identifier
  all.merged$tripnum <- 0 # init
  all.merged[all.merged$SI_HARB!="NA", "tripnum"] <- 1
  all.merged[, "tripnum"]  <- cumsum( all.merged[, "tripnum"]  )

  windows(6,6)
  plot(0,0, type="n", xlim=c(-1.5,1.5),ylim=c(58,60.5), ylab="Degree North", xlab="Degree West", axes=FALSE)
  idx <- point.in.polygon(point.x = anf(all.merged$SI_LONG), point.y = anf(all.merged$SI_LATI),
        pol.x = c(-1.5, 1.5, 1.5, -1.5), pol.y = c(58, 58, 60.5, 60.5)) > 0
  focus <- all.merged[idx,]
  lst <- split(focus, focus$tripnum, drop=TRUE)
  for(i in 1:length(lst)){
    x <- lst[[i]]
    x$SI_STATE <- replace (x$SI_STATE, x$SI_STATE==2, -1) # disable
    if(nrow(x)>0) segments(x$SI_LONG, x$SI_LATI,
                                   c(x$SI_LONG [-1], x$SI_LONG [length(x$SI_LONG)]) ,
                                       c(x$SI_LATI[-1],x$SI_LATI[length(x$SI_LATI)]), col=an(x$SI_STATE))
  }
  rect(1.005,59.32,1.20,59.37, border=4, lwd=2)
  axis(1)
  axis(2, las=2)
  box()
  map("worldHires", add = TRUE, col = "darkolivegreen4", fill = TRUE,
        bg = "white",  xlim=c(-1.5,1.5),ylim=c(58,60), regions = c("uk",
            "ireland", "france", "germany", "netherlands", "norway",
            "belgium", "spain", "luxembourg", "denmark", "sweden",
            "iceland", "portugal", "italy", "sicily", "ussr",
            "sardinia", "albania", "monaco", "turkey", "austria",
            "switzerland", "czechoslovakia", "finland", "libya",
            "hungary", "yugoslavia", "poland", "greece", "romania",
            "bulgaria", "slovakia", "morocco", "tunisia", "algeria",
            "egypt"))
  #savePlot(filename =  paste("dnk",a.year,"_area1_zoom_with_tracks",sep=''), type ="jpeg")




  ###--------------------------------###
  ###--------------------------------###
  ###--------------------------------###
  ### PLOTTING IN GOOGLE MAP         ###
  ###--------------------------------###
  ###--------------------------------###
  ###--------------------------------###

  
  Grid2KML <- function(output.mat=output.mat, what.quantity = 'effort', kmlfile='vms.kml',imagefile='vms.png')  {
     #Takes the output from vmsGridCreate ie. when plotMap is set to FALSE.
      #output.mat[["fishing"]] <- log(output.mat[["fishing"]])
  dd <- output.mat@grid
  d1 <- dd@cells.dim[1]
  d2 <- dd@cells.dim[2]
  fishing <- output.mat@data$fishing
  mat <- matrix((fishing),byrow=T,ncol=d1,nrow=d2)
  mat <- t(mat[d2:1,])
  bbox <- output.mat@bbox
  xxx <- seq(bbox[1,1],bbox[1,2],length=d1)
  yyy <- seq(bbox[2,1],bbox[2,2],length=d2)
  rr <- range(mat[mat!='-Inf'],na.rm=T)
  labs<-seq(rr[1],rr[2],length=9)
  image(xxx,yyy,mat,zlim=c(rr[1],rr[2]),xlab="",ylab="",col=rainbow(9))
  gd <- list(x=xxx,y=yyy,z=mat)
  gd.1 <- as.SpatialGridDataFrame.im(as.im(gd))
  proj4string(gd.1) <- CRS("+proj=longlat +datum=WGS84")
  vms.kml <- GE_SpatialGrid(gd.1)
  tf <- tempfile()
  png(file=paste(tf,imagefile,sep=''), width=vms.kml$width, height=vms.kml$height, bg="transparent",res=576)
  par(mar=c(0,0,0,0), xaxs="i", yaxs="i",cex=.25)
  image(as.image.SpatialGridDataFrame(gd.1[1]), col=heat.colors(9),xlim=vms.kml$xlim, ylim=vms.kml$ylim)
  dev.off()
  kmlOverlay(vms.kml, kmlfile=paste(tf,kmlfile,sep=''), imagefile=paste(tf,imagefile,sep=''), name=what.quantity)

  legend(x='topleft', legend=as.character(labs), pch = 22,pt.bg=heat.colors(length(labs)),
  title=what.quantity, ncol=2,bg="transparent",pt.cex=1.5 )
  }

  # call it...
  dnk_area1_zoom <- vmsGridCreate(all.merged, plotMap=FALSE, nameVarToSum="LE_EFF_VMS",  numCats=10,  plotPoints =FALSE, legendtitle="fishing (hours)",
          colLand="darkolivegreen4",  addICESgrid=TRUE,
              nameLon="SI_LONG", nameLat="SI_LATI", cellsizeX =0.05, cellsizeY =0.05, we = -2,ea = 2, so = 58, no = 60.5,
               breaks0=c(0,2,4,6,12,24,32, 100000), legendncol=2,
                outGridFile= paste("dnk",a.year,"_area1_zoom.asc",sep='') )

  Grid2KML (dnk_area1_zoom)  
  
  
  ###--------------------------------###
  ###--------------------------------###
  ###--------------------------------###
  ### PLOTTING THE ORIGIN OF LANDINGS###
  ###--------------------------------###
  ###--------------------------------###
  ###--------------------------------###
  sp <- "LE_EURO_COD"
  what <- "value"
  a.unit <- "(EURO)"
  cellsizeX = 0.05
  cellsizeY = 0.05 
  we <- 9.8; ea <- 12.7; no <- 55.2; so <- 54;
  breaks0 <- c(0, 100, 100 * (2^1), 100 * (2^2), 100 * (2^3), 100 * (2^4), 100 * 
            (2^5), 100 * (2^6), 100 * (2^7), 100 * (2^8), 100 * 
            (2^9), 1e+07) 
                                  

  load(file.path(main.path, "VMStools_DTUAqua_course_Jan12",
  paste("merged",a.year,sep=''), paste("all_merged_value_",a.year,".RData",sep='')))
   
  df1 <- all.merged[all.merged$SI_STATE == 1, c("SI_LATI", "SI_LONG", sp)]
  df1$SI_LATI <- anf(df1$SI_LATI)
  df1$SI_LONG <- anf(df1$SI_LONG)
  df1[, sp] <- replace(df1[, sp], is.na(df1[, sp]) | df1[, sp] < 0, 0)
  vmsGridCreate(df1, nameVarToSum = sp, numCats = 10, plotPoints = FALSE, 
        legendtitle = paste("landings", what, a.unit, sep = " "), 
        colLand = "darkolivegreen4", addICESgrid = TRUE, nameLon = "SI_LONG", 
        nameLat = "SI_LATI", cellsizeX = cellsizeX, cellsizeY = cellsizeY, 
        we = we, ea = ea, no = no, so = so, breaks0 = breaks0, 
        legendncol = 2)
    title(sp)



  ###--------------------------------###
  ###--------------------------------###
  ###--------------------------------###
  ### PLOTTING SQUARIFIED TREE MAP   ###
  ###--------------------------------###
  ###--------------------------------###
  ###--------------------------------###
  
plotTreeMap <-
function (x, gridcell = c(0.1, 0.1), gear = "OTB", xlim = c(-1,
    17), ylim = c(52, 62), acolors = rainbow(7), species.to.keep = c("LE_KG_COD",
    "LE_KG_NEP", "LE_KG_PLE", "LE_KG_SOL"))
{
    chop <- function(x) rev(rev(x)[-1])
    simple.hook <- function(z, xl, yl, xu, yu) {
        rect(xl, yl, xu, yu, col = as.character(z$one), border = NA)
    }
    squarified.treemap <- function(z, x = 0, y = 0, w = 1, h = 1,
        hook) {
        cz <- cumsum(z$size)/sum(z$size)
        n <- which.min(abs(log(max(w/h, h/w) * sum(z$size) *
            cz^2/z$size)))
        more <- n < length(z$size)
        a <- c(0, cz[1:n])/cz[n]
        if (h > w) {
            hook(z[1:n, ], x + w * chop(a), rep(y, n), x + w *
                a[-1], rep(y + h * cz[n], n))
            if (more)
                Recall(z[-(1:n), ], x, y + h * cz[n], w, h *
                  (1 - cz[n]), hook)
        }
        else {
            hook(z[1:n, ], rep(x, n), y + h * chop(a), rep(x +
                w * cz[n], n), y + h * a[-1])
            if (more)
                Recall(z[-(1:n), ], x + w * cz[n], y, w * (1 -
                  cz[n]), h, hook)
        }
    }
    x <- x[x$SI_STATE == 1, ]
    x$SI_LONG <- anf(x$SI_LONG)
    x$SI_LATI <- anf(x$SI_LATI)
    x <- x[!is.na(x$SI_LATI), ]
    idxx <- which(x$SI_LONG >= xlim[1] & x$SI_LONG <= xlim[2])
    idxy <- which(x$SI_LATI >= ylim[1] & x$SI_LATI <= ylim[2])
    x <- x[idxx[which(idxx %in% idxy)], ]
    xx <- x[x$LE_GEAR %in% gear, ]
    grids <- createGrid(xrange = xlim, yrange = ylim,
        gridcell[1], gridcell[2], type = "SpatialPixelsDataFrame")
    coords <- SpatialPointsDataFrame(cbind(x = an(ac(xx$SI_LONG)),
        y = an(ac(xx$SI_LATI))), data = xx)
    coords$dens <- overlay(grids, coords)
    DT <- data.table(data.frame(coords))
    DT$x <- af(ac(grids@coords[DT$dens, 1]))
    DT$y <- af(ac(grids@coords[DT$dens, 2]))
    idx.col <- grep("KG", names(coords))
    eq1 <- c.listquote(paste("sum(", names(coords[, idx.col]),
        ",na.rm=TRUE)", sep = ""))
    eq2 <- c.listquote(c("x", "y"))
    byRect <- data.frame(DT[, eval(eq1), by = eval(eq2)])
    colnames(byRect) <- c("SI_LONG", "SI_LATI", names(coords)[idx.col])
    byRect$SI_LONG <- signif(anf(byRect$SI_LONG))
    byRect$SI_LATI <- signif(anf(byRect$SI_LATI))
    idx.col <- grep("KG", colnames(byRect))
    rangeRect <- range(apply(byRect[idx.col], 1, sum, na.rm = T))
    rangeRect <- c(0, rangeRect[2])
    A.sum <- apply(byRect[, idx.col], 1, sum, na.rm = TRUE)
    A.sum2 <- apply(byRect[, idx.col], 2, sum, na.rm = TRUE)
    species.to.merge <- names(A.sum2)[!names(A.sum2) %in% species.to.keep]
    byRect$LE_KG_OTH <- apply(byRect[, species.to.merge], 1,
        sum, na.rm = TRUE)
    byRect <- byRect[, !names(byRect) %in% species.to.merge]
    idx.col <- grep("KG", names(byRect))
    byRect[, idx.col] <- sweep(byRect[, idx.col], 1, A.sum, FUN = "/")
    X11(7, 7)
    library(mapdata)
    map("worldHires", res = 1, xlim = xlim, ylim = ylim,
        fill = T, col = "darkgreen")
    map.axes()
    box()
    for (iRect in 1:nrow(byRect)) {
        x1 <- an(ac(byRect[iRect, "SI_LONG"]))
        y1 <- an(ac(byRect[iRect, "SI_LATI"]))
        size <- an(ac(byRect[iRect, idx.col]))
        size <- replace(size, is.na(size) | size <= 0, 1e-04) ## DEBUG: NO NEGATIVE size ALLOWED
        z <- data.frame(size = size, one = acolors[1:((1 + length(idx.col)) -
            1)])
        z <- z[order(-z$size), ]
        print(z)
        squarified.treemap(z, x = x1, y = y1, w = gridcell[1],
            h = gridcell[2], hook = simple.hook)
    }
    for (i in seq(xlim[1], xlim[2], by = gridcell[1])) abline(v = i,
        col = grey(0.9))
    for (i in seq(ylim[1], ylim[2], by = gridcell[2])) abline(h = i,
        col = grey(0.9))
    map("worldHires", add = TRUE, res = 1, xlim = xlim,
        ylim = ylim, fill = T, col = "darkgreen")
    map.axes()
    box()
    legend("topright", legend = gsub('LE_KG_','', names(byRect[, idx.col])), fill = acolors[1:((1 +
        length(idx.col)) - 1)], bg = "white")

   
  
    return()
}

 ############
 ############
 ############
 # load the data
  load(file.path(main.path, "VMStools_DTUAqua_course_Jan12",
     paste("merged",a.year,sep=''), paste("all_merged_weight_",a.year,".RData",sep='')))
  all.merged$SI_LONG <- anf(all.merged$SI_LONG ) ; all.merged$SI_LATI <- anf(all.merged$SI_LATI )
  all.merged$LE_MET  <- all.merged$LE_MET_level6


 # who is the main metier landing plaice?
 dd <- round(tapply(all.merged$LE_KG_PLE, all.merged$LE_MET_level6, sum, na.rm=TRUE))
  dd[order(dd, decreasing=TRUE)]
 
 ## add an area
 all.merged$SI_LONG  <- anf(all.merged$SI_LONG ) ; all.merged$SI_LATI <- anf(all.merged$SI_LATI )
 all.merged$ICESarea <- ICESarea(  all.merged ) # in vmstools
 all.merged$LE_MET   <- all.merged$LE_MET_level6
 all.merged <- all.merged[,sort(colnames(all.merged))]


 # choice of species
 apply(all.merged[all.merged$ICESarea=="3an" & all.merged$LE_GEAR=="OTB",c(7:31)],2,sum,na.rm=T)

 # install.packages("RColorBrewer")
 library (RColorBrewer)

 # call to plotTreeMap()
 a.metier <- "OTB_DEF_90-119_0_0"
 a.gear   <-"OTB"
  plotTreeMap (all.merged[all.merged$LE_MET_level6 %in% a.metier & all.merged$SI_STATE==1,],  gridcell=c(0.1,0.05), gear=a.gear,
                xlim= c(6.5,12), ylim= c(56.5,59.5), acolors=brewer.pal(7,"Accent"),
                  species.to.keep= c("LE_KG_COD","LE_KG_PLE", "LE_KG_NEP","LE_KG_POK","LE_KG_PRA") )
 mtext("Latitude", 2, 3) ;  mtext("Longitude",1, 2)



 # ...and load the merged output table for all vessels
  load(file.path("C:","output",paste("all_merged__",a.year,".RData",sep='')))


  ###--------------------------------###
  ###--------------------------------###
  ###--------------------------------###
  ### AUTOMATIC HIERARCHY OF PLOTS   ###
  ###--------------------------------###
  ###--------------------------------###
  ###--------------------------------###

 # load the data
  load(file.path(main.path, "VMStools_DTUAqua_course_Jan12",
     paste("merged",a.year,sep=''), paste("all_merged_weight_",a.year,".RData",sep='')))
  all.merged$SI_LONG <- anf(all.merged$SI_LONG ) ; all.merged$SI_LATI <- anf(all.merged$SI_LATI )

  graphics.off()
  
  # generate the effort maps (jpeg files) and store in a hierarchy of folders
  pings2EffortMaps (all.merged=all.merged,  output=
    file.path("C:","output"),
    cellsizeX =0.05, cellsizeY =0.05, we=9.8, ea=12.7, no=55.2, so=54.0,
     breaks0=c(0,25, 50,100,200,400,800,1600, 3200,6400,12800, 100000))

  # generate the Landings maps (jpeg files) and store in a hierarchy of folders
  pings2LandingsMaps (all.merged=all.merged,  output=
    file.path("C:","output"),
    cellsizeX =0.05, cellsizeY =0.05, we=9.8, ea=12.7, no=55.2, so=54.0,
     breaks0=c(0,25, 50,100,200,400,800,1600, 3200,6400,12800, 100000))
            
            
            
                      