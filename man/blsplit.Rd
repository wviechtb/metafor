\name{blsplit}
\alias{blsplit}
\title{Split Block Diagonal Matrix}
\description{
   Split a block diagonal matrix into a list of matrices.
}
\usage{
blsplit(x, cluster, sort=FALSE)
}
\arguments{
   \item{x}{a block diagonal matrix.}
   \item{cluster}{vector to specify the clustering variable to use for splitting.}
   \item{sort}{logical to indicate whether to sort the list by the unique cluster values (the default is \code{FALSE}).}
}
\value{
   A list of one or more matrices.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\seealso{
   \code{\link{bldiag}}, \code{\link{vcalc}}
}
\examples{
### copy data into 'dat'
dat <- dat.assink2016

### assume that the effect sizes within studies are correlated with rho=0.6
V <- vcalc(vi, cluster=study, obs=esid, data=dat, rho=0.6)

### split V matrix into list of matrices
Vs <- blsplit(V, cluster=dat$study)
Vs[1:2]
lapply(Vs[1:2], cov2cor)
}
\keyword{manip}
