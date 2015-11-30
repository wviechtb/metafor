\name{hatvalues.rma.mv}
\alias{hatvalues.rma.mv}
\title{Hat values for 'rma.mv' Objects}
\description{
   The function computes the hat values for objects of class \code{"rma.mv"}.
}
\usage{
\method{hatvalues}{rma.mv}(model, type="diagonal", \dots)
}
\arguments{
   \item{model}{an object of class \code{"rma.mv"}.}
   \item{type}{character string indicating whether to return only the diagonal of the hat matrix (\code{"diagonal"}) or the entire hat matrix (\code{"matrix"}).}
   \item{\dots}{other arguments.}
}
\value{
   Either a vector with the diagonal elements of the hat matrix or the entire hat matrix.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \cr
   package website: \url{http://www.metafor-project.org/} \cr
   author homepage: \url{http://www.wvbauer.com/}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.

   Viechtbauer, W., & Cheung, M. W.-L. (2010). Outlier and influence diagnostics for meta-analysis. \emph{Research Synthesis Methods}, \bold{1}, 112--125.
}
\keyword{models}
