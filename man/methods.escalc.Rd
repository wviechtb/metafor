\name{methods.escalc}
\alias{methods.escalc}
\alias{[.escalc}
\alias{$<-.escalc}
\alias{cbind.escalc}
\alias{rbind.escalc}
\title{Methods for 'escalc' Objects}
\description{
   Methods for objects of class \code{"escalc"}.
}
\usage{
\method{[}{escalc}(x, i, \dots)
\method{$}{escalc}(x, name) <- value
\method{cbind}{escalc}(\dots, deparse.level=1)
\method{rbind}{escalc}(\dots, deparse.level=1)
}
\arguments{
   \item{x}{an object of class \code{"escalc"}.}
   \item{\dots}{other arguments.}
}
\note{
   For the \code{`[`} method, any variables specified as part of the \code{i} argument will be searched for within object \code{x} first (see \sQuote{Examples}).
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### select rows where variable 'alloc' is equal to 'random'
dat[dat$alloc == "random",]

### variables specified are automatically searched for within the object itself
dat[alloc == "random",]

### note: this behavior is specific to 'escalc' objects; this doesn't work for regular data frames
}
\keyword{internal}
