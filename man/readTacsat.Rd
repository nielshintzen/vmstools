\name{readTacsat}
\alias{readTacsat}
\title{
Read Tacsat data into R
}
\description{
Reading tacsat data from delimited file into R, checking for obligatory columns and formatting all data columns.
}
\usage{
readTacsat(file, sep = ",", dec = ".")
}
\arguments{
  \item{file}{file path + file name}
  \item{sep}{delimiter used in file (default to ',')}
  \item{dec}{decimal notation used in file (default to '.')}
}
\value{
Returns the formatted Tacsat dataset}
\references{EU lot 2 project}
\author{Niels T. Hintzen}
\seealso{\code{readEflalo()}}
\examples{
data(tacsat)
dir.create("C:/tmpTacsat")

#temporarily write tacsat file to disk to thereafter read it back in again
tacsat$SI_HE[which(tacsat$SI_HE>360 | tacsat$SI_HE < 0)] <- NA
tacsat$SI_LONG[which(tacsat$SI_LONG < -180 | tacsat$SI_LONG > 180)] <- NA
tacsat$SI_LATI[which(tacsat$SI_LATI < -90 | tacsat$SI_LATI > 90)] <- NA
write.table(tacsat, file = "C:/tmpTacsat/tacsat.csv", quote = TRUE, sep = ",",
        eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = TRUE)

#Read in tacsat file
tacsat <- readTacsat("C:/tmpTacsat/tacsat.csv")


}
