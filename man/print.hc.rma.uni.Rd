\name{print.hc.rma.uni}
\alias{print.hc.rma.uni}
\title{Print Method for 'hc.rma.uni' Objects}
\description{
   Print method for objects of class \code{"hc.rma.uni"}. \loadmathjax
}
\usage{
\method{print}{hc.rma.uni}(x, digits=x$digits, \dots)
}
\arguments{
   \item{x}{an object of class \code{"hc.rma.uni"} obtained with \code{\link[=hc.rma.uni]{hc}}.}
   \item{digits}{integer to specify the number of decimal places to which the printed results should be rounded (the default is to take the value from the object).}
   \item{\dots}{other arguments.}
}
\details{
   The output is a data frame with two rows, the first (labeled \code{rma}) corresponding to the results based on the usual estimation method, the second (labeled \code{hc}) corresponding to the results based on the method by Henmi and Copas (2010). The data frame includes the following variables:
   \itemize{
   \item the method used to estimate \mjseqn{\tau^2} (always \code{DL} for \code{hc})
   \item the estimated amount of heterogeneity
   \item the estimated average true outcome
   \item the corresponding standard error (\code{NA} when \code{transf} argument has been used)
   \item the lower and upper confidence interval bounds
   }
}
\value{
   The function returns the data frame invisibly.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{hc.rma.uni}}
}
\keyword{print}
