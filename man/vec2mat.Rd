\name{vec2mat}
\alias{vec2mat}
\title{Convert a Vector into a Square Matrix}
\description{
   Function to convert a vector into a square matrix by filling up the lower triangular part of the matrix.
}
\usage{
vec2mat(x, diag=FALSE, corr=!diag, dimnames)
}
\arguments{
   \item{x}{a vector of the correct length.}
   \item{diag}{logical to specify whether the vector also contains the diagonal values of the lower triangular part of the matrix (the default is \code{FALSE}).}
   \item{corr}{logical to specify whether the diagonal of the matrix should be replaced with 1's (the default is to do this when \code{diag=FALSE}).}
   \item{dimnames}{optional vector of the correct length with the dimension names of the matrix.}
}
\details{
   The values in \code{x} are filled into the lower triangular part of a square matrix with the appropriate dimensions (which are determined based on the length of \code{x}). If \code{diag=TRUE}, then \code{x} is assumed to also contain the diagonal values of the lower triangular part of the matrix. If \code{corr=TRUE}, then the diagonal of the matrix is replaced with 1's.
}
\value{
   A matrix.
}
\author{
   Wolfgang Viechtbauer (\email{wvb@metafor-project.org}, \url{https://www.metafor-project.org}).
}
\examples{
vec2mat(1:6, corr=FALSE)
vec2mat(seq(0.2, 0.7, by=0.1), corr=TRUE)
vec2mat(1:10, diag=TRUE)
vec2mat(1:6, corr=FALSE, dimnames=c("A","B","C","D"))
}
\keyword{manip}
