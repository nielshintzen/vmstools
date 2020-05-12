plot_mcp <- function (plotnew = TRUE, plotpoints = TRUE, points.col = "black",
    points.pch = 1, titletxt = "Title", xaxis = "Easting (m)",
    yaxis = "Northing (m)", mcp.col = "black", mcp.lwd = 2, fill.col = NA,
    jpeg = FALSE, ...)
{
    par(...)
    if (jpeg) {
        jpeg(filename = paste("MCP", r.MCP$id, ".jpg", sep = ""),
            width = 600, height = 600, pointsize = 12, quality = 90,
            bg = "white", res = NA)
    }
    if (plotnew) {
        plot(r.MCP$MCP, xlab = xaxis, ylab = yaxis, colpol = fill.col,
            col = mcp.col, lwd = mcp.lwd)
        title(paste(titletxt, sep = ""))
    }
    else {
        default.parameters <- par(no.readonly = TRUE)
        xlim.max <- c(default.parameters$usr[1], default.parameters$usr[2])
        ylim.max <- c(default.parameters$usr[3], default.parameters$usr[4])
        plot(r.MCP$MCP, xlab = xaxis, ylab = yaxis, colpol = fill.col,
            col = mcp.col, lwd = mcp.lwd, xlim = xlim.max, ylim = ylim.max,
            add = TRUE)
        title(paste(titletxt, sep = ""))
    }
    if (plotpoints) {
        points(r.MCP$points, col = points.col, pch = points.pch)
    }
    if (jpeg) {
        dev.off()
    }
}
