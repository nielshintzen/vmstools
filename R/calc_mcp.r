calc_mcp <- function (id = 1, points = NULL, filename = "MCP_Output.txt",
    verbose = FALSE, pct = 100)
{
    require(adehabitat)
    
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


    errorcode <- 1000
    if ((pct > 100) || (pct < 0)) {
        errorcode <- 100
        if (verbose) {
            cat("\n\nWARNING: Supplied percentage must be between 0 and 100 (inclusive).")
            cat("\nERROR CODE: ", errorcode, "\n\n", sep = "")
        }
        return("ERROR")
    }
    if (length(dim(points)) != 2) {
        errorcode <- 71
        if (verbose) {
            cat("\n\nWARNING: Provided points input matrix has fewer than 2 columns.")
            cat("\nERROR CODE: ", errorcode, "\n\n", sep = "")
        }
        return("ERROR")
    }
    if (dim(points)[2] != 2) {
        errorcode <- 70
        if (verbose) {
            cat("\n\nWARNING: Provided points input matrix has too many columns, only 2 are allowed.")
            cat("\nERROR CODE: ", errorcode, "\n\n", sep = "")
        }
        return("ERROR")
    }
    else {
        temp <- as.data.frame(cbind(1, points))
        temp[, 1] <- as.factor(temp[, 1])
        MCP <- (mcp(temp[, 2:3], temp[, 1], percent = pct))
        area <- mcp.area(temp[, 2:3], temp[, 1], percent = pct,
            unin = "m", unout = "km2", plotit = FALSE)
    }
    coordsMCP <- cbind(id, MCP)
    tmp <- coordsMCP[, 2:4]
    outtabMCP <- cbind(id, tmp)
    write.table(outtabMCP, sep = ",", append = TRUE, file = filename,
        col.names = FALSE)
    r.MCP <- list(MCP = MCP, points = points, id = id, MCP.area = area,
        MCP.pct = pct)
    assign("r.MCP", r.MCP, pos = 1)
    mcp.result <- list(id = id, MCP.area = area, MCP.pct = pct,
        MCP.coords = MCP)
    return(mcp.result)
}
