#This is taken from gpclib package
mcp.area <- function(xy, id, percent = seq(20, 100, by = 5),
                     unin = c("m","km"), unout = c("ha", "km2", "m2"), plotit = TRUE){
    unin <- match.arg(unin)
    unout <- match.arg(unout)
    if (length(id) != nrow(xy))
        stop("xy and id should be of the same length")
    xy <- xy[!is.na(xy[, 1]), ]
    xy <- xy[!is.na(xy[, 2]), ]
    id <- id[!is.na(xy[, 1])]
    id <- id[!is.na(xy[, 2])]
    lev <- percent
    res <- list()
    ar <- matrix(0, nrow = length(lev), ncol = nlevels(factor(id)))
    lixy <- split(xy, id)
    le <- names(lixy)
    for (i in 1:length(lev)) {
        ar[i, ] <- unlist(lapply(lixy, function(z) {
            res <- mcp(z, rep(1, nrow(z)), percent = lev[i])
            class(res) <- "data.frame"
            return(area.poly(as(res[, 2:3], "gpc.poly")))
        }))
    }
    ar <- as.data.frame(ar)
    names(ar) <- le
    if (unin == "m") {
        if (unout == "ha")
            ar <- ar/10000
        if (unout == "km2")
            ar <- ar/1e+06
    }
    if (unin == "km") {
        if (unout == "ha")
            ar <- ar * 100
        if (unout == "m2")
            ar <- ar * 1e+06
    }
    row.names(ar) <- lev
    class(ar) <- c("hrsize", "data.frame")
    attr(ar, "units") <- unout
    if (plotit)
        plot(ar)
    return(ar)
}
