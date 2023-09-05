\name{contrmat}
\alias{contrmat}
\title{Construct Contrast Matrix for Two-Group Comparisons}
\description{
   Function to construct a matrix that indicates which two groups have been contrasted against each other in each row of a dataset.
}
\usage{
contrmat(data, grp1, grp2, last, shorten=FALSE, minlen=2, check=TRUE, append=TRUE)
}
\arguments{
   \item{data}{a data frame in wide format.}
   \item{grp1}{either the name (given as a character string) or the position (given as a single number) of the first group variable in the data frame.}
   \item{grp2}{either the name (given as a character string) or the position (given as a single number) of the second group variable in the data frame.}
   \item{last}{optional character string to specify which group will be placed in the last column of the matrix (must be one of the groups in the group variables). If not given, the most frequently occurring second group is placed last.}
   \item{shorten}{logical to specify whether the variable names corresponding to the group names should be shortened (the default is \code{FALSE}).}
   \item{minlen}{integer to specify the minimum length of the shortened variable names (the default is 2).}
   \item{check}{logical to specify whether the variables names should be checked to ensure that they are syntactically valid variable names and if not, they are adjusted (by \code{\link{make.names}}) so that they are (the default is \code{TRUE}).}
   \item{append}{logical to specify whether the contrast matrix should be appended to the data frame specified via the \code{data} argument (the default is \code{TRUE}). If \code{append=FALSE}, only the contrast matrix is returned.}
}
\details{
   The function can be used to construct a matrix that indicates which two groups have been contrasted against each other in each row of a data frame (with \code{1} for the first group, \code{-1} for the second group, and \code{0} otherwise).

   The \code{grp1} and \code{grp2} arguments are used to specify the group variables in the dataset (either as character strings or as numbers indicating the column positions of these variables in the dataset). Optional argument \code{last} is used to specify which group will be placed in the last column of the matrix.

   If \code{shorten=TRUE}, the variable names corresponding to the group names are shortened (to at least \code{minlen}; the actual length might be longer to ensure uniqueness of the variable names).

   The examples below illustrate the use of this function.
}
\value{
   A matrix with as many variables as there are groups.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{to.wide}} for a function to create \sQuote{wide} format datasets.

   \code{\link[metadat]{dat.senn2013}}, \code{\link[metadat]{dat.hasselblad1998}}, \code{\link[metadat]{dat.lopez2019}} for illustrative examples.
}
\examples{
### restructure to wide format
dat <- dat.senn2013
dat <- dat[c(1,4,3,2,5,6)]
dat <- to.wide(dat, study="study", grp="treatment", ref="placebo", grpvars=4:6)
dat

### add contrast matrix
dat <- contrmat(dat, grp1="treatment.1", grp2="treatment.2")
dat

### data in long format
dat <- dat.hasselblad1998
dat

### restructure to wide format
dat <- to.wide(dat, study="study", grp="trt", ref="no_contact", grpvars=6:7)
dat

### add contrast matrix
dat <- contrmat(dat, grp1="trt.1", grp2="trt.2", shorten=TRUE, minlen=3)
dat
}
\keyword{manip}
