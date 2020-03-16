\name{methods.list.rma}
\alias{methods.list.rma}
\alias{as.data.frame.list.rma}
\alias{as.matrix.list.rma}
\alias{[.list.rma}
\alias{head.list.rma}
\alias{tail.list.rma}
\alias{$<-.list.rma}
\title{Methods for 'list.rma' Objects}
\description{
   Methods for objects of class \code{"list.rma"}.
}
\usage{
\method{as.data.frame}{list.rma}(x, \dots)
\method{as.matrix}{list.rma}(x, \dots)
\method{[}{list.rma}(x, i, \dots)
\method{head}{list.rma}(x, n=6L, \dots)
\method{tail}{list.rma}(x, n=6L, \dots)
\method{$}{list.rma}(x, name) <- value
}
\arguments{
   \item{x}{an object of class \code{"list.rma"}.}
   \item{\dots}{other arguments.}
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{http://www.metafor-project.org/}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{https://www.jstatsoft.org/v036/i03}.
}
\keyword{internal}
