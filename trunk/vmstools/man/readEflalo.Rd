\name{readEflalo}
\alias{readEflalo}
\title{
Read Eflalo data into R}
\description{
Reading eflalo data from delimited file into R, checking for obligatory columns and formatting all data columns.
}
\usage{
readEflalo(file, sep = ",", dec = ".")
}
\arguments{
  \item{file}{file path + file name}
  \item{sep}{delimiter used in file (default to ',')}
  \item{dec}{decimal notation used in file (default to '.')}
}

\value{
Returns the formatted Eflalo dataset}
\references{EU lot 2 project}
\author{Niels T. Hintzen}
\seealso{\code{\link{readTacsat}}}
\examples{
data(eflalo)
dir.create("C:/tmpEflalo")

#temporarily write tacsat file to disk to thereafter read it back in again
write.table(eflalo, file = "C:/tmpEflalo/eflalo.csv", quote = TRUE, sep = ",",
    eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = TRUE)

#Read in tacsat file
eflalo <- readEflalo("C:/tmpEflalo/eflalo.csv")

}

