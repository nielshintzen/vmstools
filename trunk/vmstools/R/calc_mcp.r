calc_mcp <- function (id = 1, points = NULL, filename = "MCP_Output.txt",
    verbose = FALSE, pct = 100)
{
    require(adehabitat)
    require(gpclib)
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
