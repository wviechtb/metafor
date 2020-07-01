\name{print.permutest.rma.uni}
\alias{print.permutest.rma.uni}
\title{Print Method for 'permutest.rma.uni' Objects}
\description{
   Print method for objects of class \code{"permutest.rma.uni"}.
}
\usage{
\method{print}{permutest.rma.uni}(x, digits=x$digits, signif.stars=getOption("show.signif.stars"),
      signif.legend=signif.stars, \dots)
}
\arguments{
   \item{x}{an object of class \code{"permutest.rma.uni"}.}
   \item{digits}{integer specifying the number of decimal places to which the printed results should be rounded (the default is to take the value from the object).}
   \item{signif.stars}{logical indicating whether p-values should be encoded visually with \sQuote{significance stars}. Defaults to the \code{show.signif.stars} slot of \code{\link{options}}.}
   \item{signif.legend}{logical indicating whether the legend for the \sQuote{significance stars} should be printed. Defaults to the value for \code{signif.stars}.}
   \item{\dots}{other arguments.}
}
\details{
   The output includes:
   \itemize{
   \item the results of the omnibus test of moderators. Suppressed if the model includes only one coefficient (e.g., only an intercept, like in the fixed- and random-effects model). The p-value based on the permutation test is indicated by \code{p-val*}.
   \item a table with the estimated coefficients, corresponding standard errors, test statistics, p-values, and confidence interval bounds. The p-values based on the permutation tests are indicated by \code{pval*}. When permutation-based CIs have been obtained, then the CI bounds are indicated with \code{ci.lb*} and \code{ci.ub*}.
   }
}
\value{
   The function does not return an object.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{http://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{permutest.rma.uni}}
}
\keyword{print}
