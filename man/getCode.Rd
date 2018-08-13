\name{getCode}
\alias{getCode}
\title{
Get stan model code for egg count data analysis dealing with outliers
}
\description{
Obtain stan model specifications with handling of ponential outliers that can be used in the \strong{eggCounts} package.
}
\usage{
getCode(modelName)
}
\arguments{
  \item{modelName}{string. One of four availale models ("Po", "UPo", "ZIPo", "ZIUPo").}
}
\value{
A character string containing the stan model specification.
}
\author{
Craig Wang \email{craig.wang@uzh.ch}
}
\note{
The models are developed by Anja Falleger.
}
\examples{
model_one <- getCode("Po")
writeLines(model_one)
}
\seealso{
\code{\link[eggCounts]{fecr_stanExtra}} for using the model specification.
}
