\name{print.matreg}
\alias{print.matreg}
\alias{summary.matreg}
\alias{print.summary.matreg}
\title{Print and Summary Methods for 'matreg' Objects}
\description{
   Functions to print objects of class \code{"matreg"} and \code{"summary.matreg"}. \loadmathjax
}
\usage{
\method{print}{matreg}(x, digits, signif.stars=getOption("show.signif.stars"),
      signif.legend=signif.stars, \dots)

\method{summary}{matreg}(object, digits, \dots)

\method{print}{summary.matreg}(x, digits, signif.stars=getOption("show.signif.stars"),
      signif.legend=signif.stars, \dots)
}
\arguments{
   \item{x}{an object of class \code{"matreg"} or \code{"summary.matreg"} (for \code{print}).}
   \item{object}{an object of class \code{"matreg"} (for \code{summary}).}
   \item{digits}{integer to specify the number of decimal places to which the printed results should be rounded. If unspecified, the default is to take the value from the object.}
   \item{signif.stars}{logical to specify whether p-values should be encoded visually with \sQuote{significance stars}. Defaults to the \code{show.signif.stars} slot of \code{\link{options}}.}
   \item{signif.legend}{logical to specify whether the legend for the \sQuote{significance stars} should be printed. Defaults to the value for \code{signif.stars}.}
   \item{\dots}{other arguments.}
}
\details{
   The output is a table with the estimated coefficients, corresponding standard errors, test statistics, p-values, and confidence interval bounds. When using \code{summary}, the output includes additional statistics, including \mjseqn{R^2} and the omnibus test of the model coefficients (either an F- or a chi-square test).
}
\value{
   The function does not return an object.
}
\author{
   Wolfgang Viechtbauer (\email{wvb@metafor-project.org}, \url{https://www.metafor-project.org}).
}
\seealso{
   \code{\link{matreg}} for the function to create \code{matreg} objects.
}
\keyword{print}
