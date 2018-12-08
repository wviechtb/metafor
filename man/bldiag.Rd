\name{bldiag}
\alias{bldiag}
\title{Construct Block Diagonal Matrix}
\description{
   Function to construct a block diagonal matrix from (a list of) matrices.
}
\usage{
bldiag(\dots)
}
\arguments{
   \item{\dots}{individual matrices or a list of matrices.}
}
\author{
   Posted to R-help by Berton Gunter (2 Sep 2005).
}
\seealso{
   \code{\link{rma.mv}}
}
\examples{
### copy data into 'dat'
dat <- dat.berkey1998

### construct list of the variance-covariance matrices of the observed outcomes for the studies
V <- lapply(split(dat[,c("v1i", "v2i")], dat$trial), as.matrix)
V

### construct block diagonal matrix
V <- bldiag(V)
V
}
\keyword{manip}
