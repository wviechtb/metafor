\name{hettest}
\alias{hettest}
\title{Blah}
\description{
   Blah. \loadmathjax
}
\usage{
hettest(x, vi, sei, subset, data, method="REML",
        test="score", boot=TRUE, progbar=TRUE, digits, \dots)
}
\arguments{
   \item{x}{a vector with the observed effect sizes or outcomes or an object of class \code{"rma"}.}
   \item{vi}{vector with the corresponding sampling variances (ignored if \code{x} is an object of class \code{"rma"}).}
   \item{sei}{vector with the corresponding standard errors (note: only one of the two, \code{vi} or \code{sei}, needs to be specified).}
   \item{subset}{optional (logical or numeric) vector to specify the subset of studies that should be used for the calculation (ignored if \code{x} is an object of class \code{"rma"}).}
   \item{data}{optional data frame containing the variables given to the arguments above.}
   \item{method}{character string.}
   \item{test}{optional character string. See \sQuote{Details}.}
   \item{boot}{logical.}
   \item{progbar}{logical.}
   \item{digits}{optional integer to specify the number of decimal places to which the printed results should be rounded.}
   \item{\dots}{other arguments.}
}
\details{

   Blah.

}
\value{
   An object of class \code{"hettest"}. The object is a list containing the following components:
   \item{x2}{blah.}
   \item{df}{blah.}
   \item{pval}{the p-value.}
   \item{\dots}{some additional elements/values.}
}
\author{
   Wolfgang Viechtbauer (\email{wvb@metafor-project.org}, \url{https://www.metafor-project.org}).
}
\keyword{htest}
