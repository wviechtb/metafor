\name{blsplit}
\alias{blsplit}
\title{Split Block Diagonal Matrix}
\description{
   Function to split a block diagonal matrix into a list of sub-matrices.
}
\usage{
blsplit(x, cluster, fun, args, sort=FALSE)
}
\arguments{
   \item{x}{a block diagonal matrix.}
   \item{cluster}{vector to specify the clustering variable to use for splitting.}
   \item{fun}{optional argument to specify a function to apply to each sub-matrix.}
   \item{args}{optional argument to specify any additional argument(s) for the function specified via \code{fun}.}
   \item{sort}{logical to indicate whether to sort the list by the unique cluster values (the default is \code{FALSE}).}
}
\value{
   A list of one or more sub-matrices.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\seealso{
   \code{\link{bldiag}} for a function to create a block diagonal matrix based on sub-matrices.

   \code{\link{vcalc}} for a function to construct a variance-covariance matrix of dependent effect sizes or outcomes, which often has a block diagonal structure.
}
\examples{
### copy data into 'dat'
dat <- dat.assink2016

### assume that the effect sizes within studies are correlated with rho=0.6
V <- vcalc(vi, cluster=study, obs=esid, data=dat, rho=0.6)

### split V matrix into list of sub-matrices
Vs <- blsplit(V, cluster=dat$study)
Vs[1:2]
lapply(Vs[1:2], cov2cor)

### illustrate the use of the fun and args arguments
blsplit(V, cluster=dat$study, cov2cor)[1:2]
blsplit(V, cluster=dat$study, round, 3)[1:2]
}
\keyword{manip}
