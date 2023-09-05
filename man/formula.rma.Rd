\name{formula.rma}
\alias{formula}
\alias{formula.rma}
\title{Extract the Model Formula from 'rma' Objects}
\description{
   Function to extract the model formula from objects of class \code{"rma"}.
}
\usage{
\method{formula}{rma}(x, type="mods", \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma"}.}
   \item{type}{the formula which should be returned; either \code{"mods"} (default), \code{"yi"} (in case argument \code{yi} was used to specify a formula), or \code{"scale"} (only for location-scale models).}
   \item{\dots}{other arguments.}
}
\value{
   The requested formula.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.glmm}}, and \code{\link{rma.mv}} for functions to fit models for which a model formula can be extracted.
}
\examples{
### copy BCG vaccine data into 'dat'
dat <- dat.bcg

### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat,
              slab=paste(author, ", ", year, sep=""))

### mixed-effects meta-regression model
res <- rma(yi, vi, mods = ~ ablat + alloc, data=dat)
formula(res, type="mods")

### specify moderators via 'yi' argument
res <- rma(yi ~ ablat + alloc, vi, data=dat)
formula(res, type="yi")
}
\keyword{models}
