\name{w_ratio}
\alias{w_ratio}
\title{
Compute outlier weights using ratios
}
\description{
Compute outlier weights using pre- and post-treatment count ratios.
}
\usage{
w_ratio(preFEC, postFEC)
}
\arguments{
  \item{postFEC}{vector of post-treatment faecal egg counts}
  \item{preFEC}{vector of pre-treatment faecal egg counts}
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
w_ratio(epgs$before, epgs$after)
}
\seealso{
\code{\link[eggCounts]{fecr_stanExtra}} for using the model specification.
}
