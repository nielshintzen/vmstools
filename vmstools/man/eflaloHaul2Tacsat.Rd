\name{eflaloHaul2Tacsat}
\alias{eflaloHaul2Tacsat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
eflaloHaul2Tacsat(x, npoints = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{npoints}{
%%     ~~Describe \code{npoints} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, npoints = NULL) 
{
    if (is.null(npoints) == T) 
        npoints <- 10
    if (length(npoints) != nrow(x)) 
        npoints <- rep(npoints, nrow(x))
    x$LE_SDATIM <- as.POSIXct(paste(x$LE_CDAT, x$LE_STIME), format = "\%d/\%m/\%Y \%H:\%M")
    x$LE_EDATIM <- as.POSIXct(paste(x$LE_CDAT, x$LE_ETIME), format = "\%d/\%m/\%Y \%H:\%M")
    y <- data.frame(do.call(cbind, as.list(rep(NA, 9))))
    colnames(y) <- c("VE_COU", "VE_REF", "SI_LATI", "SI_LONG", 
        "SI_DATE", "SI_TIME", "SI_SP", "SI_HE", "FT_REF")
    taLst <- lapply(as.list(1:nrow(x)), function(z) {
        ta <- y[rep(1, npoints[z]), ]
        ta$VE_COU[] <- x$VE_COU[z]
        ta$VE_REF[] <- x$VE_REF[z]
        ta$SI_LATI[] <- seq(x$LE_SLAT[z], x$LE_ELAT[z], length.out = npoints[z])
        ta$SI_LONG[] <- seq(x$LE_SLON[z], x$LE_ELON[z], length.out = npoints[z])
        ta$SI_DATE[] <- format(seq(x$LE_SDATIM[z], x$LE_EDATIM[z], 
            length.out = npoints[z]), "\%d/\%m/\%Y")
        ta$SI_TIME[] <- format(seq(x$LE_SDATIM[z], x$LE_EDATIM[z], 
            length.out = npoints[z]), "\%H:\%M")
        ta$SI_SP[] <- c(distance(x$LE_SLON[z], x$LE_SLAT[z], 
            x$LE_ELON[z], x$LE_ELAT[z])/1.852)/c(difftime(x$LE_EDATIM[z], 
            x$LE_SDATIM[z], units = "hours"))
        ta$SI_HE[] <- bearing(x$LE_SLON[z], x$LE_SLAT[z], x$LE_ELON[z], 
            x$LE_ELAT[z])
        ta$FT_REF[] <- x$FT_REF[z]
        return(ta)
    })
    ta <- do.call(rbind, taLst)
    return(ta)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
