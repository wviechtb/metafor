\name{matreg}
\alias{matreg}
\title{Fit Regression Models based on Correlation and Covariance Matrices}
\description{
   Function to fit regression models based on correlation and covariance matrices. \loadmathjax
}
\usage{
matreg(y, x, R, n, V, cov=FALSE, means, ztor=FALSE, nearpd=FALSE, level=95, digits, \dots)
}
\arguments{
   \item{y}{index of the outcome variable.}
   \item{x}{indices of the predictor variables.}
   \item{R}{correlation or covariance matrix (or only the lower triangular part including the diagonal).}
   \item{n}{sample size based on which the elements in the correlation/covariance matrix were computed.}
   \item{V}{variance-covariance matrix of the lower triangular elements of the correlation/covariance matrix. Either \code{V} or \code{n} should be specified, not both. See \sQuote{Details}.}
   \item{cov}{logical to specify whether \code{R} is a covariance matrix (the default is \code{FALSE}).}
   \item{means}{optional vector to specify the means of the variables (only relevant when \code{cov=TRUE}).}
   \item{ztor}{logical to specify whether \code{R} is a matrix of r-to-z transformed correlations and it should be back-transformed to raw correlations (the default is \code{FALSE}). See \sQuote{Details}.}
   \item{nearpd}{logical to specify whether the \code{\link[Matrix]{nearPD}} function from the \href{https://cran.r-project.org/package=Matrix}{Matrix} package should be used when the \mjeqn{R_{x,x}}{R[x,x]} matrix cannot be inverted. See \sQuote{Note}.}
   \item{level}{numeric value between 0 and 100 to specify the confidence interval level (the default is 95).}
   \item{digits}{integer to specify the number of decimal places to which the printed results should be rounded. If unspecified, the default is 4.}
   \item{\dots}{other arguments.}
}
\details{
   Let \mjseqn{R} be a \mjeqn{p \times p}{pxp} correlation or covariance matrix. Let \mjseqn{y} denote the row/column of the outcome variable and \mjseqn{x} the row(s)/column(s) of the predictor variable(s) in this matrix. Let \mjeqn{R_{x,x}}{R[x,x]} and \mjeqn{R_{x,y}}{R[x,y]} denote the corresponding submatrices of \mjseqn{R}. Then \mjdeqn{b = R_{x,x}^{-1} R_{x,y}}{b = R[x,x]^(-1) R[x,y]} yields the standardized or raw regression coefficients (depending on whether \mjseqn{R} is a correlation or covariance matrix, respectively) when regressing the outcome variable on the predictor variable(s).

   The \mjseqn{R} matrix may be computed based on a single sample of \mjseqn{n} subjects. In this case, one should specify the sample size via argument \code{n}. The variance-covariance matrix of the standardized regression coefficients is then given by \mjeqn{\mbox{Var}[b] = \mbox{MSE} \times R_{x,x}^{-1}}{Var[b] = MSE * R[x,x]^(-1)}, where \mjeqn{\mbox{MSE} = (1 - b'R_{x,y}) / (n - m)}{MSE = (1 - b'R[x,y]) / (n -m)} and \mjseqn{m} denotes the number of predictor variables. The standard errors of the regression coefficients are then given by the square root of the diagonal elements of \mjeqn{\mbox{Var}[b]}{Var[b]}. Test statistics (in this case, t-statistics) and the corresponding p-values can then be computed as in a regular regression analysis. When \mjseqn{R} is a covariance matrix, one should set \code{cov=TRUE} and specify the means of the \mjseqn{p} variables via argument \code{means} to obtain raw regression coefficients including the intercept and corresponding standard errors.

   Alternatively, \mjseqn{R} may be the result of a meta-analysis of correlation coefficients. In this case, the elements in \mjseqn{R} are pooled correlation coefficients and the variance-covariance matrix of these pooled coefficients should be specified via argument \code{V}. The order of elements in \code{V} should correspond to the order of elements in the lower triangular part of \mjseqn{R} column-wise. For example, if \mjseqn{R} is a \mjeqn{4 \times 4}{4x4} matrix\ifelse{text}{,}{ of the form: \mjtdeqn{\left[ \begin{array}{cccc} 1 & & & \\\ r_{21} & 1 & & \\\ r_{31} & r_{32} & 1 & \\\ r_{41} & r_{42} & r_{43} & 1 \end{array} \right]}{\begin{bmatrix} 1 & & & \\\\\ r_{21} & 1 & & \\\\\ r_{31} & r_{32} & 1 & \\\\\ r_{41} & r_{42} & r_{43} & 1 \end{bmatrix}}{}} then the elements are \mjseqn{r_{21}}, \mjseqn{r_{31}}, \mjseqn{r_{41}}, \mjseqn{r_{32}}, \mjseqn{r_{42}}, and \mjseqn{r_{43}} and hence \code{V} should be a \mjeqn{6 \times 6}{6x6} variance-covariance matrix of these elements in this order. The variance-covariance matrix of the standardized regression coefficients (i.e., \mjeqn{\mbox{Var}[b]}{Var[b]}) is then computed as a function of \code{V} as described in Becker (1992) using the multivariate delta method. The standard errors of the standardized regression coefficients are then given by the square root of the diagonal elements of \mjeqn{\mbox{Var}[b]}{Var[b]}. Test statistics (in this case, z-statistics) and the corresponding p-values can then be computed in the usual manner.

   In case \mjseqn{R} is the result of a meta-analysis of Fisher r-to-z transformed correlation coefficients (and hence \code{V} is then the corresponding variance-covariance matrix of these pooled transformed coefficients), one should set argument \code{ztor=TRUE}, so that the appropriate back-transformation is then applied to \code{R} (and \code{V}) within the function.

   Finally, \mjseqn{R} may be a covariance matrix based on a meta-analysis (e.g., the estimated variance-covariance matrix of the random effects in a multivariate model). In this case, one should set \code{cov=TRUE} and \code{V} should again be the variance-covariance matrix of the elements in \mjseqn{R}, but now including the diagonal. Hence, if \mjseqn{R} is a \mjeqn{4 \times 4}{4x4} matrix\ifelse{text}{,}{ of the form: \mjtdeqn{\left[ \begin{array}{cccc} \tau_1^2 & & & \\\ \tau_{21} & \tau_2^2 & & \\\ \tau_{31} & \tau_{32} & \tau_3^2 & \\\ \tau_{41} & \tau_{42} & \tau_{43} & \tau_4^2 \end{array} \right]}{\begin{bmatrix} \tau_1^2 & & & \\\\\ \tau_{21} & \tau_2^2 & & \\\\\ \tau_{31} & \tau_{32} & \tau_3^2 & \\\\\ \tau_{41} & \tau_{42} & \tau_{43} & \tau_4^2 \end{bmatrix}}{}} then the elements are \mjseqn{\tau^2_1}, \mjseqn{\tau_{21}}, \mjseqn{\tau_{31}}, \mjseqn{\tau_{41}}, \mjseqn{\tau^2_2}, \mjseqn{\tau_{32}}, \mjseqn{\tau_{42}}, \mjseqn{\tau^2_3}, \mjseqn{\tau_{43}}, and \mjseqn{\tau^2_4}, and hence \code{V} should be a \mjeqn{10 \times 10}{10x10} variance-covariance matrix of these elements in this order. Argument \code{means} can then again be used to specify the means of the variables.
}
\value{
   An object of class \code{"matreg"}. The object is a list containing the following components:
   \item{tab}{a data frame with the estimated model coefficients, standard errors, test statistics, degrees of freedom (only for t-tests), p-values, and lower/upper confidence interval bounds.}
   \item{vb}{the variance-covariance matrix of the estimated model coefficients.}
   \item{\dots}{some additional elements/values.}

   The results are formatted and printed with the \code{print.matreg} function.
}
\note{
   Only the lower triangular part of \code{R} (and \code{V} if it is specified) is used in the computations.

   If \mjeqn{R_{x,x}}{R[x,x]} is not invertible, an error will be issued. In this case, one can set argument \code{nearpd=TRUE}, in which case the \code{\link[Matrix]{nearPD}} function from the \href{https://cran.r-project.org/package=Matrix}{Matrix} package will be used to find the nearest positive semi-definite matrix, which should be invertible. The results should be treated with caution when this is done.

   When \mjseqn{R} is a covariance matrix with \code{V} and \code{means} specified, the means are treated as known constants when estimating the standard error of the intercept.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
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
R <- vec2mat(coef(res))
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
R <- vec2mat(coef(res))
rownames(R) <- colnames(R) <- c("perf", "acog", "asom", "conf")
fit <- matreg(1, 2:4, R=R, V=vcov(res), ztor=TRUE)
fit}
}
\keyword{models}
