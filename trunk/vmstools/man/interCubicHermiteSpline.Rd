\name{interCubicHermiteSpline}
\alias{interCubicHermiteSpline}
\title{Interpolate using Cubic Hermite Spline
}
\description{
This is an internal function of interpolateTacsat.
Function to interpolate VMS data using a cubic hermite spline function.
}
\usage{
interCubicHermiteSpline(spltx, spltCon, res, params, headingAdjustment)
}
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
  \item{params}{
Parameters to use for interpolation
}
  \item{headingAdjustment}{
Parameter to adjust the choice of heading depending on its own or previous point
}
}
\author{
Niels T. Hintzen
}
