\name{dfround}
\alias{dfround}
\title{Round Variables in a Data Frame}
\description{
   Function to round the numeric variables in a data frame.
}
\usage{
dfround(x, digits, drop0=TRUE)
}
\arguments{
   \item{x}{a data frame.}
   \item{digits}{either a single integer or a numeric vector of the same length as there are columns in \code{x}.}
   \item{drop0}{logical (or a vector thereof) to specify whether trailing zeros after the decimal mark should be removed (the default is \code{TRUE}).}
}
\details{
   A simple convenience function to round the numeric variables in a data frame, possibly to different numbers of digits. Hence, \code{digits} can either be a single integer (which will then be used to round all numeric variables to the specified number of digits) or a numeric vector (of the same length as there are columns in \code{x}) to specify the number of digits to which each variable should be rounded.

   Non-numeric variables are skipped. If \code{digits} is a vector, some arbitrary value (or \code{NA}) can be specified for those variables.

   Note: When \code{drop0=FALSE}, then \code{\link{formatC}} is used to format the numbers, which turns them into character variables.
}
\value{
   Returns the data frame with variables rounded as specified.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\examples{
dat <- dat.bcg
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat)
res <- rma(yi, vi, mods = ~ ablat + year, data=dat)
coef(summary(res))
dfround(coef(summary(res)), digits=c(2,3,2,3,2,2))
}
\keyword{manip}
