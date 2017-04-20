\name{model.matrix.rma}
\alias{model.matrix.rma}
\title{Model Matrix for 'rma' Objects}
\description{
   The function extracts the model matrix for objects of class \code{"rma"}.
}
\usage{
\method{model.matrix}{rma}(object, \dots)
}
\arguments{
   \item{object}{an object of class \code{"rma"}.}
   \item{\dots}{other arguments.}
}
\value{
   The model matrix.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \cr
   package website: \url{http://www.metafor-project.org/} \cr
   author homepage: \url{http://www.wvbauer.com/}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
   \code{\link{fitted.rma}}
}
\examples{
### meta-analysis of the log risk ratios using a mixed-effects meta-regression model
### with multiple moderators (absolute latitude, publication year, and allocation method)
res <- rma(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg,
           mods = ~ ablat + year + alloc, data=dat.bcg)
model.matrix(res)
}
\keyword{models}
