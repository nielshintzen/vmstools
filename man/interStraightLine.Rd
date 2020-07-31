\name{interStraightLine}
\alias{interStraightLine}
\title{
Interpolate using Straight line
}
\description{
This is an internal function of interpolateTacsat.
Function to interpolate VMS data using a straight line.
}
\usage{
interStraightLine(spltx, spltCon, res)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spltx}{
Tacsat dataset split up by vessel using 'split'.
}
  \item{spltCon}{
List of VMS records that are connected
}
  \item{res}{
Number of additional points to use for interpolation
}
}
\author{
Niels T. Hintzen
}
