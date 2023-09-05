\name{coef.matreg}
\alias{coef.matreg}
\alias{vcov.matreg}
\title{Extract the Model Coefficients and Variance-Covariance Matrix from 'matreg' Objects}
\description{
   Methods for objects of class \code{"matreg"}.
}
\usage{
\method{coef}{matreg}(object, \dots)
\method{vcov}{matreg}(object, \dots)
}
\arguments{
   \item{object}{an object of class \code{"matreg"}.}
   \item{\dots}{other arguments.}
}
\details{
   The \code{coef} function extracts the estimated model coefficients from objects of class \code{"matreg"}. The \code{vcov} function extracts the corresponding variance-covariance matrix.
}
\value{
   Either a vector with the estimated model coefficients or a variance-covariance matrix.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{matreg}} for the function to create \code{matreg} objects.
}
\examples{
### fit a regression model with lm() to the 'mtcars' dataset
res <- lm(mpg ~ hp + wt + am, data=mtcars)
coef(res)
vcov(res)

### covariance matrix of the dataset
S <- cov(mtcars)

### fit the same regression model using matreg()
res <- matreg(y="mpg", x=c("hp","wt","am"), R=S, cov=TRUE,
              means=colMeans(mtcars), n=nrow(mtcars))
coef(res)
vcov(res)
}
\keyword{models}
