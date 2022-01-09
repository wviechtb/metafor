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
### copy data into 'dat' and examine data
dat <- dat.viechtbauer2021

### calculate log odds ratios and corresponding sampling variances
dat <- escalc(measure="OR", ai=xTi, n1i=nTi, ci=xCi, n2i=nCi, add=1/2, to="all", data=dat)

### fit mixed-effects meta-regression model
res <- rma(yi, vi, mods = ~ dose, data=dat)

### get studentized residuals
sav <- rstudent(res)
sav

### studies with studentized residuals larger than +-1.96
sav[abs(sav$z) > 1.96,]

### variables specified are automatically searched for within the object itself
sav[abs(z) > 1.96,]

### note: this behavior is specific to 'rma.list' objects; this doesn't work for regular data frames
}
\keyword{internal}
