\name{simulate.rma}
\alias{simulate}
\alias{simulate.rma}
\title{Simulate Method for 'rma' Objects}
\description{
   Function to simulate effect sizes or outcomes based on \code{"rma"} model objects.
}
\usage{
\method{simulate}{rma}(object, nsim=1, seed=NULL, olim, \dots)
}
\arguments{
   \item{object}{an object of class \code{"rma"}.}
   \item{nsim}{number of response vectors to simulate (defaults to 1).}
   \item{seed}{an object to specify if and how the random number generator should be initialized (\sQuote{seeded}). Either \code{NULL} or an integer that will be used in a call to \code{set.seed} before simulating the response vectors. If set, the value is saved as the \code{"seed"} attribute of the returned value. The default, \code{NULL} will not change the random generator state, and return \code{\link{.Random.seed}} as the \code{"seed"} attribute; see \sQuote{Value}.}
   \item{olim}{optional argument to specify observation/outcome limits for the simulated values. If unspecified, no limits are used.}
   \item{\dots}{other arguments.}
}
\details{
   The model specified via \code{object} must be a model fitted with either the \code{\link{rma.uni}} or \code{\link{rma.mv}} functions.
}
\value{
   A data frame with \code{nsim} columns with the simulated effect sizes or outcomes.

   The data frame comes with an attribute \code{"seed"}. If argument \code{seed} is \code{NULL}, the attribute is the value of \code{\link{.Random.seed}} before the simulation was started; otherwise it is the value of the \code{seed} argument with a \code{"kind"} attribute with value \code{as.list(RNGkind())}.
}
\note{
   If the outcome measure used for the analysis is bounded (e.g., correlations are bounded between -1 and +1, proportions are bounded between 0 and 1), one can use the \code{olim} argument to enforce those observation/outcome limits when simulating values (simulated values cannot exceed those bounds then).
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{rma.uni}} and \code{\link{rma.mv}} for functions to fit models for which simulated effect sizes or outcomes can be generated.
}
\examples{
### copy BCG vaccine data into 'dat'
dat <- dat.bcg

### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat)
dat

### fit random-effects model
res <- rma(yi, vi, data=dat)
res

### simulate 5 sets of new outcomes based on the fitted model
newdat <- simulate(res, nsim=5, seed=1234)
newdat
}
\keyword{datagen}
