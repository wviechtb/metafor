\name{rcalc}
\alias{rcalc}
\title{Calculate the Variance-Covariance of Dependent Correlation Coefficients}
\description{
   Function to calculate the variance-covariance matrix of correlation coefficients computed based on the same sample of subjects. \loadmathjax
}
\usage{
rcalc(x, ni, data, rtoz=FALSE, nfun="min", sparse=FALSE, \dots)
}
\arguments{
   \item{x}{a formula of the form \code{ri ~ var1 + var2 | study}. Can also be a correlation matrix or list thereof. See \sQuote{Details}.}
   \item{ni}{vector to specify the sample sizes based on which the correlations were computed.}
   \item{data}{data frame containing the variables specified via the formula (and the sample sizes).}
   \item{rtoz}{logical to specify whether to transform the correlations via Fisher's r-to-z transformation (the default is \code{FALSE}).}
   \item{nfun}{a character string to specify how the \sQuote{common} sample size within each study should be computed. Possible options are \code{"min"} (for the minimum), \code{"harmonic"} (for the harmonic mean), or \code{"mean"} (for the arithmetic mean). Can also be a function. See \sQuote{Details}.}
   \item{sparse}{logical to specify whether the variance-covariance matrix should be returned as a sparse matrix (the default is \code{FALSE}).}
   \item{\dots}{other arguments.}
}
\details{
   A meta-analysis of correlation coefficients may involve multiple correlation coefficients extracted from the same study. When these correlations are computed based on the same sample of subjects, then they are typically not independent. The \code{rcalc} function can be used to create a dataset with the correlation coefficients (possibly transformed with Fisher's r-to-z transformation) and the corresponding variance-covariance matrix. The dataset and variance-covariance matrix can then be further meta-analyzed using the \code{\link{rma.mv}} function.

   When computing the covariance between two correlation coefficients, we can distinguish two cases:
   \enumerate{
   \item In the first case, one of the variables involved in the two correlation coefficients is the same. For example, in \mjseqn{r_{12}} and \mjseqn{r_{13}}, variable 1 is common to both correlation coefficients. This is sometimes called the (partially) \sQuote{overlapping} case. The covariance between the two correlation coefficients, \mjeqn{\text{Cov}[r_{12}, r_{13}]}{Cov[r_{12}, r_{13}]}, then depends on the degree of correlation between variables 2 and 3 (i.e., \mjseqn{r_{23}}).
   \item In the second case, none of the variables are common to both correlation coefficients. For example, this would be the case if we have correlations \mjseqn{r_{12}} and \mjseqn{r_{34}} based on 4 variables. This is sometimes called the \sQuote{non-overlapping} case. The covariance between the two correlation coefficients, \mjeqn{\text{Cov}[r_{12}, r_{34}]}{Cov[r_{12}, r_{34}]}, then depends on \mjseqn{r_{13}}, \mjseqn{r_{14}}, \mjseqn{r_{23}}, and \mjseqn{r_{24}}.
   }

   Equations to compute these covariances can be found, for example, in Steiger (1980) and Olkin and Finn (1990).

   To use the \code{rcalc} function, one needs to construct a data frame that contains a study identifier (say \code{study}), two variable identifiers (say \code{var1} and \code{var2}), the corresponding correlation coefficients (say \code{ri}), and the sample sizes based on which the correlation coefficients were computed (say \code{ni}). Then the first argument should be a formula of the form \code{ri ~ var1 + var2 | study}, argument \code{ni} is set equal to the variable name containing the sample sizes, and the data frame containing these variables is specified via the \code{data} argument. When using the function for a single study, one can leave out the study identifier from the formula.

   When argument \code{rtoz} is set to \code{TRUE}, then the correlations are transformed with Fisher's r-to-z transformation (Fisher, 1921) and the variance-covariance matrix is computed for the transformed values.

   In some cases, the sample size may not be identical within a study (e.g., \mjseqn{r_{12}} may have been computed based on 120 subjects while \mjseqn{r_{13}} was computed based on 118 subjects due to 2 missing values in variable 3). For constructing the variance-covariance matrix, we need to assume a \sQuote{common} sample size for all correlation coefficients within the study. Argument \code{nfun} provides some options for how the common sample size should be computed. Possible options are \code{"min"} (for using the minimum sample size within a study as the common sample size), \code{"harmonic"} (for using the harmonic mean), or \code{"mean"} (for using the arithmetic mean). The default is \code{"min"}, which is a conservative choice (i.e., it will overestimate the sampling variances of coefficients that were computed based on a sample size that was actually larger than the minimum sample size). One can also specify a function via the \code{nfun} argument (which should take a numeric vector as input and return a single value).

   Instead of specifying a formula, one can also pass a correlation matrix to the function via argument \code{x}. Argument \code{ni} then specifies the (common) sample size based on which the elements in the correlation matrix were computed. One can also pass a list of correlation matrices via argument \code{x}, in which case argument \code{ni} should be a vector of sample sizes of the same length as \code{x}.
}
\value{
   A list containing the following components:
   \item{dat}{a data frame with the study identifier, the two variable identifiers, a variable pair identifier, the correlation coefficients (possibly transformed with Fisher's r-to-z transformation), and the (common) sample sizes.}
   \item{V}{corresponding variance-covariance matrix (given as a sparse matrix when \code{sparse=TRUE}).}

   Note that a particular covariance can only be computed when all of the correlation coefficients involved in the covariance equation are included in the dataset. If one or more coefficients needed for the computation are missing, then the resulting covariance will also be missing (i.e., \code{NA}).
}
\note{
   For raw correlation coefficients, the variance-covariance matrix is computed with \mjseqn{n-1} in the denominator (instead of \mjseqn{n} as suggested in Steiger, 1980, and Olkin & Finn, 1990). This is more consistent with the usual equation for computing the sampling variance of a correlation coefficient (which also typically uses \mjseqn{n-1} in the denominator).

   For raw and r-to-z transformed coefficients, the variance-covariance matrix will only be computed when the (common) sample size for a study is at least 5.
}
\author{
   Wolfgang Viechtbauer (\email{wvb@metafor-project.org}, \url{https://www.metafor-project.org}).
}
\references{
   Fisher, R. A. (1921). On the \dQuote{probable error} of a coefficient of correlation deduced from a small sample. \emph{Metron}, \bold{1}, 1--32. \verb{http://hdl.handle.net/2440/15169}

   Olkin, I., & Finn, J. D. (1990). Testing correlated correlations. \emph{Psychological Bulletin}, \bold{108}(2), 330--333. \verb{https://doi.org/10.1037/0033-2909.108.2.330}

   Steiger, J. H. (1980). Tests for comparing elements of a correlation matrix. \emph{Psychological Bulletin}, \bold{87}(2), 245--251. \verb{https://doi.org/10.1037/0033-2909.87.2.245}
}
\seealso{
   \code{\link{rma.mv}} for a model fitting function that can be used to meta-analyze dependent correlation coefficients.

   \code{\link[metadat]{dat.craft2003}} for an illustrative example.
}
\examples{
############################################################################

### copy data into 'dat' and examine the first 12 rows
dat <- dat.craft2003
head(dat, 12)

### construct dataset and var-cov matrix of the correlations
tmp <- rcalc(ri ~ var1 + var2 | study, ni=ni, data=dat)
V <- tmp$V
dat <- tmp$dat

### examine data and var-cov matrix for study 1
dat[dat$study == 1,]
blsplit(V, dat$study, round, 4)$`1`

### examine data and var-cov matrix for study 6
dat[dat$study == 6,]
blsplit(V, dat$study, round, 4)$`6`

### examine data and var-cov matrix for study 17
dat[dat$study == 17,]
blsplit(V, dat$study, round, 4)$`17`

############################################################################

### copy data into 'dat' and examine the first 12 rows
dat <- dat.craft2003
head(dat, 12)

### restructure data from study 1 into a correlation matrix
R1 <- diag(4)
R1[lower.tri(R1)] <- dat$ri[dat$study == 1]
R1[upper.tri(R1)] <- t(R1)[upper.tri(R1)]
rownames(R1) <- colnames(R1) <- c("perf", "acog", "asom", "conf")
R1

### restructure data from study 3 into a correlation matrix
R3 <- diag(4)
R3[lower.tri(R3)] <- dat$ri[dat$study == 3]
R3[upper.tri(R3)] <- t(R3)[upper.tri(R3)]
rownames(R3) <- colnames(R3) <- c("perf", "acog", "asom", "conf")
R3

### an example where a correlation matrix is passed to rcalc()
rcalc(R1, ni=142)

### an example where a list of correlation matrices is passed to rcalc()
tmp <- rcalc(list("1"=R1,"3"=R3), ni=c(142,37))
V <- tmp$V
dat <- tmp$dat

### examine data and var-cov matrix for study 1
dat[dat$id == 1,]
blsplit(V, dat$id, round, 4)$`1`

### examine data and var-cov matrix for study 3
dat[dat$id == 3,]
blsplit(V, dat$id, round, 4)$`3`

############################################################################
}
\keyword{datagen}
