
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
# F.Bastardie
# using stolen code on the net
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!##

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
    coords@data$dens      <- over(as(coords,"SpatialPoints"), as(grids,"SpatialPixels"))
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
    rangeRect <- range(apply(byRect[idx.col], 1, sum, na.rm = TRUE))
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
    require(mapdata)
    map("worldHires", resolution = 1, xlim = xlim, ylim = ylim,
        fill = TRUE, col = "darkgreen")
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
    map("worldHires", add = TRUE, resolution = 1, xlim = xlim,
        ylim = ylim, fill = TRUE, col = "darkgreen")
    map.axes()
    box()
    legend("topright", legend = gsub('LE_KG_','', names(byRect[, idx.col])), fill = acolors[1:((1 +
        length(idx.col)) - 1)], bg = "white")

   
  
    return()
}


##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##


  ## do mergeEflalo2Pings()....
  # ...and load the merged output table for all vessels
  #load(file.path("C:","output","all_merged__1800.RData"))


  #graphics.off()
  #df1 <- all.merged[, c("LE_MET_level6", "LE_GEAR","SI_STATE", "SI_LATI","SI_LONG","LE_KG_COD","LE_KG_PLE")]
  #df1$SI_LONG <- anf(df1$SI_LONG)
  #df1$SI_LATI <- anf(df1$SI_LATI)
  #df1 <-   df1[ !is.na(df1$SI_LATI),]
  #df1 <-   df1[ !is.na(df1$SI_LONG),]
  #df1$LE_MET  <- df1$LE_MET_level6


  ## call to plotTreeMap()
  #a.metier <- "TBB_DEF_70-99_0_0"
  #a.gear   <-"TBB"
  #plotTreeMap (df1[df1$LE_MET_level6 %in% a.metier & df1$SI_STATE==1,],  gridcell=c(0.1,0.05), gear=a.gear,
  #              xlim= c(3,6), ylim= c(50,54), acolors=rainbow(7),
  #                species.to.keep= c("LE_KG_COD","LE_KG_PLE") )
  #mtext("Latitude", 2, 3) ;  mtext("Longitude",1, 2)



 