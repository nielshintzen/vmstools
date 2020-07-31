\name{summarizeEflalo}
\alias{summarizeEflalo}
\title{
Summarize primary aspects of eflalo file}
\description{
Summary nr vessels, countries, temporal range, effort, top 3 catches and values}
\usage{
summarizeEflalo(eflalo)
}
\arguments{
  \item{eflalo}{
eflalo dataframe}
}
\author{
Niels T. Hintzen}
\seealso{
\code{\link{summarizeTacsat}},\code{\link{summary}}
}
\examples{
data(eflalo)

summarizeEflalo(eflalo)
}