\name{replmiss}
\alias{replmiss}
\title{Replace Missing Values in a Vector}
\description{
   Function to replace missing (\code{NA}) values in a vector.
}
\usage{
replmiss(x, y)
}
\arguments{
   \item{x}{vector that may include one or more missing values.}
   \item{y}{either a scalar or a vector of the same length as \code{x} with the value(s) to replace missing values with.}
}
\value{
   Vector \code{x} with the missing values replaced based on the scalar or vector \code{y}.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{http://www.metafor-project.org/}
}
\examples{
x <- c(4,2,7,NA,1,NA,5)
x <- replmiss(x,0)
x

x <- c(4,2,7,NA,1,NA,5)
y <- c(2,3,6,5,8,1,2)
x <- replmiss(x,y)
x
}
\keyword{manip}
