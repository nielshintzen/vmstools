#' Plot the Minimum Convex Polygon
#' 
#' This function plots the MCP as a polygon, which covers the geographical
#' extent of a set of points on a Cartesian plane.
#' 
#' The r.MCP object (generated in calc_mcp function) is required to plot the
#' MCP.
#' 
#' @param plotnew Boolean: Set to TRUE to create a new plot. Set to FALSE to
#' overlay current plot.
#' @param plotpoints Boolean: Set to TRUE if the point observations are to be
#' plotted
#' @param points.col Specify a colour for the point observations
#' @param points.pch Specify a plotting symbol for the point observations
#' @param titletxt A string to use as the title on the plot
#' @param xaxis A string to label the x-axis of the plot
#' @param yaxis A string to label the y-axis of the plot
#' @param mcp.col Specify the line colour for the MCP
#' @param mcp.lwd Specify the line width for the MCP
#' @param fill.col Specify a fill colour for the MCP
#' @param jpeg Boolean: Set to TRUE if the plot should be saved in JPEG format
#' @param \dots Arguments to be passed to graphical parameters
#' @author Randy Bui, Ron N. Buliung, Tarmo K. Remmel
#' @examples
#' 
#' data(tacsat)
#' 
#' calc_mcp(id=1, points = tacsat[1:10,c("SI_LONG","SI_LATI")], filename="MCP_Output.txt",
#' verbose = FALSE, pct = 100)
#' plot_mcp(plotnew=TRUE, plotpoints=TRUE, titletxt="Title",
#' xaxis= "Easting (m)", yaxis="Northing (m)")
#' 
#' @export plot_mcp
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
