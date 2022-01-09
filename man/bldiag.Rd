\name{bldiag}
\alias{bldiag}
\title{Construct Block Diagonal Matrix}
\description{
   Function to construct a block diagonal matrix from (a list of) matrices.
}
\usage{
bldiag(\dots, order)
}
\arguments{
   \item{\dots}{individual matrices or a list of matrices.}
   \item{order}{optional argument to specify a variable based on which a square block diagonal matrix should be ordered.}
}
\author{
   Posted to R-help by Berton Gunter (2 Sep 2005) with some further adjustments by Wolfgang Viechtbauer
}
\seealso{
   \code{\link{rma.mv}} for the model fitting function that can take such a block diagonal matrix as input (for the \code{V} argument).

   \code{\link{blsplit}} for a function that can split a block diagonal matrix into a list of sub-matrices.
}
\examples{
### copy data into 'dat'
dat <- dat.berkey1998
dat

### construct list with the variance-covariance matrices of the observed outcomes for the studies
V <- lapply(split(dat[c("v1i","v2i")], dat$trial), as.matrix)
V

### construct block diagonal matrix
V <- bldiag(V)
V

### if we split based on 'author', the list elements in V are in a different order than tha data
V <- lapply(split(dat[c("v1i","v2i")], dat$author), as.matrix)
V

### can use 'order' argument to reorder the block-diagonal matrix into the correct order
V <- bldiag(V, order=dat$author)
V
}
\keyword{manip}
