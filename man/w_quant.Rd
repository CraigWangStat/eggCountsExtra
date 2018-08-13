\name{w_quant}
\alias{w_quant}
\title{
Compute outlier weights using quantiles
}
\description{
Compute outlier weights using quantiles.
}
\usage{
w_quant(postFEC)
}
\arguments{
  \item{postFEC}{vector of post-treatment faecal egg counts}
}
\value{
The returned value is a list consists of
\item{w}{a vector of outlier weights}
\item{wmo}{weighted mean of outliers}
\item{postmean}{mean of post-treatment counts}
}
\author{
Anja Fallegger \cr
Craig Wang
}
\examples{
data(epgs)
w_quant(epgs$after)
}
\seealso{
\code{\link[eggCounts]{fecr_stanExtra}} for using the model specification.
}
