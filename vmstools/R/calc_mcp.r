#' Computing the Minimum Convex Polygon (MCP)
#' 
#' This function computes the Minimum Convex Polygon (MCP) from a set of
#' points. The MCP is the minimum area polygon containing a set of point
#' locations.
#' 
#' This function is most powerful when used repetitively within a loop to
#' compute the MCP for subsets of points stored in a large data table.
#' 
#' @param id Provide a unique integer to identify an MCP from others that you
#' may construct with other data points
#' @param points Two-column matrix or data frame of point coordinates
#' @param filename A character name for an ASCII output file
#' @param verbose Boolean: set to TRUE if extended processing feedback is
#' wanted
#' @param pct Integer 0 <= pct <=100, the percentage of the MCP for which area
#' is provided
#' @return The returned result is a list: %% If it is a LIST, use
#' \item{MCP.area }{The area of the MCP in square kilometers} \item{MCP.pct
#' }{The desired percentage of the MCP for which the area is computed}
#' \item{MPC.coords}{A matrix containing MCP vertices. Each row represents a
#' unique point, the first column contains x-coordinates, and the second,
#' y-coordinates }
#' @note Results are stored in the r.MCP object (required for graphical
#' visualization using plot_mcp). This function can be used on its own (once)
#' or repetitively in a loop to process grouped point data stored in a larger
#' table. When used repetitively, be sure to increment the id parameter to
#' ensure that each MCP has a unique identifier. The output ASCII coordinate
#' file can be further processed using the makeshapes function to generate an
#' ESRI Shapefile for MCP polygons.
#' @author Randy Bui, Ron N. Buliung, Tarmo K. Remmel
#' @references Builds upon MCP functions available in the adehabitat package
#' @examples
#' 
#' data(tacsat)
#' calc_mcp(id=1, points = tacsat[1:10,c("SI_LONG","SI_LATI")], filename="MCP_Output.txt",
#' verbose = FALSE, pct = 100)
#' 
#' @export calc_mcp
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
