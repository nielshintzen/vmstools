\name{summarizeTacsat}
\alias{summarizeTacsat}
\title{
Summarize primary aspects of tacsat file}
\description{
Summary nr vessels, countries, spatial range, temporal range, effort from tacsat}
\usage{
summarizeTacsat(tacsat)
}
\arguments{
  \item{tacsat}{
tacsat dataframe
}
}
\author{
Niels T. Hintzen}
\seealso{
\code{\link{summarizeEflalo}},\code{\link{summary}}
}
\examples{
data(tacsat)

summarizeTacsat(tacsat)
}