\name{matreg}
\alias{matreg}
\title{Fit Regression Model based on a Correlation Matrix}
\description{
   Function to fit a regression model based on a correlation matrix. \loadmathjax
}
\usage{
matreg(y, x, R, n, V, ztor=FALSE, nearPD=FALSE, level=95, digits)
}
\arguments{
   \item{y}{index of the outcome variable.}
   \item{x}{indices of the predictor variable(s).}
   \item{R}{correlation matrix (or only the lower triangular part including the diagonal).}
   \item{n}{sample size based on which the elements in the correlation matrix were computed.}
   \item{V}{variance-covariance matrix of the lower triangular elements of the correlation matrix. Either \code{V} or \code{n} should be specified, not both. See \sQuote{Details}.}
   \item{ztor}{logical indicating whether \code{R} is a matrix of r-to-z transformed correlations (the default is \code{FALSE}). See \sQuote{Details}.}
   \item{nearPD}{logical indicating whether the \code{\link[Matrix]{nearPD}} function from the \code{Matrix} package should be used when the \mjeqn{R_{x,x}}{R[x,x]} matrix cannot be inverted. See \sQuote{Note}.}
   \item{level}{numerical value between 0 and 100 specifying the confidence interval level (the default is 95).}
   \item{digits}{integer specifying the number of decimal places to which the printed results should be rounded (if unspecified, the default is 4).}
}
\details{
   The function can be used to fit a regression model based on a correlation matrix.

   Let \mjseqn{R} be a \mjeqn{p \times p}{pxp} correlation matrix. Let \mjseqn{y} denote the row/column of the outcome variable and \mjseqn{x} the row(s)/column(s) of the predictor variable(s). Let \mjeqn{R_{x,x}}{R[x,x]} and \mjeqn{R_{x,y}}{R[x,y]} denote the corresponding submatrices of \mjseqn{R}. Then \mjdeqn{b = R_{x,x}^{-1} R_{x,y}}{b = R[x,x]^(-1) R[x,y]} yields the (standardized) regression slopes when regressing the outcome variable on the predictor variable(s).

   The \mjseqn{R} matrix may be computed based on a single sample of \mjseqn{n} subjects. In this case, one should specify the sample size via argument \code{n}. The variance-covariance matrix of the estimated model coefficients is then given by \mjeqn{\mbox{Var}[b] = \mbox{MSE} \times R_{x,x}^{-1}}{Var[b] = MSE * R[x,x]^(-1)}, where \mjeqn{\mbox{MSE} = (1 - b'R_{x,y}) / (n - m)}{MSE = (1 - b'R[x,y]) / (n -m)}, where \mjseqn{m} denotes the number of predictor variables. The standard errors of the estimated model coefficients are then given by the square root of the diagonal elements of \mjeqn{\mbox{Var}[b]}{Var[b]}. Test statistics (in this case, t-statistics) and the corresponding p-values can then be computed as in regular regression.

   Alternatively, \mjseqn{R} may be the result of a meta-analysis of correlation coefficients. In this case, the elements in \mjseqn{R} are really pooled correlation coefficients and there should be a variance-covariance matrix of these pooled coefficients. In this case, one should specify this variance-covariance matrix via argument \code{V}. The order of elements in \code{V} should correspond to the order of elements in the lower triangular part of \mjseqn{R} column-wise. For example, if \mjseqn{R} is a \mjeqn{4 \times 4}{4x4} matrix, then the elements are \mjseqn{r_{2,1}}, \mjseqn{r_{3,1}}, \mjseqn{r_{4,1}}, \mjseqn{r_{3,2}}, \mjseqn{r_{4,2}}, and \mjseqn{r_{4,3}} and hence \code{V} should be a \mjeqn{6 \times 6}{6x6} variance-covariance matrix of these elements. The variance-covariance matrix of the estimated model coefficients (i.e., \mjeqn{\mbox{Var}[b]}{Var[b]}) is then computed as described in Becker (1992) using the multivariate delta method. The standard errors of the estimated model coefficients are then given by the square root of the diagonal elements of \mjeqn{\mbox{Var}[b]}{Var[b]}. Test statistics (in this case, z-statistics) and the corresponding p-values can then be computed in the usual manner.

   In case \mjseqn{R} is the result of a meta-analysis of Fisher r-to-z transformed correlation coefficients (and hence \code{V} is then the corresponding variance-covariance matrix of these pooled transformed coefficients), one should set argument \code{ztor=TRUE}, so that the appropriate back-transformation is then applied to \code{R} (and \code{V}) within the function.
}
\value{
   An object of class \code{"matreg"}. The object is a list containing the following components:
   \item{tab}{a data frame with the estimated model coefficients, standard errors, test statistics, p-values, and lower/upper confidence interval bounds.}
   \item{vb}{the variance-covariance matrix of the estimated model coefficients.}
   \item{\dots}{some additional elements/values.}

   The results are formatted and printed with the \code{print.matreg} function.
}
\note{
   Only the lower triangular part of \code{R} is used in the computations.

   If \mjeqn{R_{x,x}}{R[x,x]} is not invertible, an error will be issued. In this case, one can set argument \code{nearPD=TRUE}, in which case the \code{\link[Matrix]{nearPD}} function from the \code{Matrix} package will be used to find the nearest positive semi-definite matrix, which should be invertible. The results should be treated with caution when this is done.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{http://www.metafor-project.org}
}
\references{
   Becker, B. J. (1992). Using results from replicated studies to estimate linear models. \emph{Journal of Educational Statistics}, \bold{17}(4), 341--362. \verb{https://doi.org/10.3102/10769986017004341}

   Becker, B. J. (1995). Corrections to "Using results from replicated studies to estimate linear models". \emph{Journal of Educational and Behavioral Statistics}, \bold{20}(1), 100--102. \verb{https://doi.org/10.3102/10769986020001100}
}
\seealso{
   \code{\link{rcalc}}, \code{\link{rma.mv}}
}
\examples{
### copy data into 'dat'
dat <- dat.craft2003

### construct dataset and var-cov matrix of the correlations
tmp <- rcalc(ri ~ var1 + var2 | study, ni=ni, data=dat)
V <- tmp$V
dat <- tmp$dat

### turn var1.var2 into a factor with the desired order of levels
dat$var1.var2 <- factor(dat$var1.var2,
   levels=c("acog.perf", "asom.perf", "conf.perf", "acog.asom", "acog.conf", "asom.conf"))

### multivariate random-effects model
res <- rma.mv(yi, V, mods = ~ var1.var2 - 1, random = ~ var1.var2 | study, struct="UN", data=dat)
res

### restructure estimated mean correlations into a 4x4 matrix
R <- matrix(NA, nrow=4, ncol=4)
R[lower.tri(R)] <- coef(res)
rownames(R) <- colnames(R) <- c("perf", "acog", "asom", "conf")
round(R, digits=3)

### check that order in vcov(res) corresponds to order in R
round(vcov(res), digits=4)

### fit regression model with 'perf' as outcome and 'acog', 'asom', and 'conf' as predictors
fit <- matreg(1, 2:4, R=R, V=vcov(res))
fit

\dontrun{
### repeat the above but with r-to-z transformed correlations
dat <- dat.craft2003
tmp <- rcalc(ri ~ var1 + var2 | study, ni=ni, data=dat, rtoz=TRUE)
V <- tmp$V
dat <- tmp$dat
dat$var1.var2 <- factor(dat$var1.var2,
   levels=c("acog.perf", "asom.perf", "conf.perf", "acog.asom", "acog.conf", "asom.conf"))
res <- rma.mv(yi, V, mods = ~ var1.var2 - 1, random = ~ var1.var2 | study, struct="UN", data=dat)
R <- matrix(NA, nrow=4, ncol=4)
R[lower.tri(R)] <- coef(res)
rownames(R) <- colnames(R) <- c("perf", "acog", "asom", "conf")
fit <- matreg(1, 2:4, R=R, V=vcov(res), ztor=TRUE)
fit}
}
\keyword{models}
