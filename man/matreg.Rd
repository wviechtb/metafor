\name{matreg}
\alias{matreg}
\title{Fit Regression Models based on Correlation and Covariance Matrices}
\description{
   Function to fit regression models based on correlation and covariance matrices. \loadmathjax
}
\usage{
matreg(y, x, R, n, V, cov=FALSE, means, ztor=FALSE,
       nearpd=FALSE, level=95, digits, \dots)
}
\arguments{
   \item{y}{index (or name given as a character string) of the outcome variable.}
   \item{x}{indices (or names given as a character vector) of the predictor variables.}
   \item{R}{correlation or covariance matrix (or only the lower triangular part including the diagonal).}
   \item{n}{sample size based on which the elements in the correlation/covariance matrix were computed.}
   \item{V}{variance-covariance matrix of the lower triangular elements of the correlation/covariance matrix. Either \code{V} or \code{n} should be specified, not both. See \sQuote{Details}.}
   \item{cov}{logical to specify whether \code{R} is a covariance matrix (the default is \code{FALSE}).}
   \item{means}{optional vector to specify the means of the variables (only relevant when \code{cov=TRUE}).}
   \item{ztor}{logical to specify whether \code{R} is a matrix of r-to-z transformed correlations and hence should be back-transformed to raw correlations (the default is \code{FALSE}). See \sQuote{Details}.}
   \item{nearpd}{logical to specify whether the \code{\link[Matrix]{nearPD}} function from the \href{https://cran.r-project.org/package=Matrix}{Matrix} package should be used when the \mjeqn{R_{x,x}}{R[x,x]} matrix cannot be inverted. See \sQuote{Note}.}
   \item{level}{numeric value between 0 and 100 to specify the confidence interval level (the default is 95; see \link[=misc-options]{here} for details).}
   \item{digits}{optional integer to specify the number of decimal places to which the printed results should be rounded.}
   \item{\dots}{other arguments.}
}
\details{
   Let \mjseqn{R} be a \mjeqn{p \times p}{pxp} correlation or covariance matrix. Let \mjseqn{y} denote the row/column of the outcome variable and \mjseqn{x} the row(s)/column(s) of the predictor variable(s) in this matrix. Let \mjeqn{R_{x,x}}{R[x,x]} and \mjeqn{R_{x,y}}{R[x,y]} denote the corresponding submatrices of \mjseqn{R}. Then \mjdeqn{b = R_{x,x}^{-1} R_{x,y}}{b = R[x,x]^(-1) R[x,y]} yields the standardized or raw regression coefficients (depending on whether \mjseqn{R} is a correlation or covariance matrix, respectively) when regressing the outcome variable on the predictor variable(s).

   The \mjseqn{R} matrix may be computed based on a single sample of \mjseqn{n} subjects. In this case, one should specify the sample size via argument \code{n}. The variance-covariance matrix of the standardized regression coefficients is then given by \mjeqn{\mbox{Var}[b] = \mbox{MSE} \times R_{x,x}^{-1}}{Var[b] = MSE * R[x,x]^(-1)}, where \mjeqn{\mbox{MSE} = (1 - b'R_{x,y}) / (n - m)}{MSE = (1 - b'R[x,y]) / (n -m)} and \mjseqn{m} denotes the number of predictor variables. The standard errors are then given by the square root of the diagonal elements of \mjeqn{\mbox{Var}[b]}{Var[b]}. Test statistics (in this case, t-statistics) and the corresponding p-values can then be computed as in a regular regression analysis. When \mjseqn{R} is a covariance matrix, one should set \code{cov=TRUE} and specify the means of the \mjseqn{p} variables via argument \code{means} to obtain raw regression coefficients including the intercept and corresponding standard errors.

   Alternatively, \mjseqn{R} may be the result of a meta-analysis of correlation coefficients. In this case, the elements in \mjseqn{R} are pooled correlation coefficients and the variance-covariance matrix of these pooled coefficients should be specified via argument \code{V}. The order of elements in \code{V} should correspond to the order of elements in the lower triangular part of \mjseqn{R} column-wise. For example, if \mjseqn{R} is a \mjeqn{4 \times 4}{4x4} matrix\ifelse{text}{,}{ of the form: \mjtdeqn{\left[ \begin{array}{cccc} 1 & & & \\\ r_{21} & 1 & & \\\ r_{31} & r_{32} & 1 & \\\ r_{41} & r_{42} & r_{43} & 1 \end{array} \right]}{\begin{bmatrix} 1 & & & \\\\\ r_{21} & 1 & & \\\\\ r_{31} & r_{32} & 1 & \\\\\ r_{41} & r_{42} & r_{43} & 1 \end{bmatrix}}{}} then the elements are \mjseqn{r_{21}}, \mjseqn{r_{31}}, \mjseqn{r_{41}}, \mjseqn{r_{32}}, \mjseqn{r_{42}}, and \mjseqn{r_{43}} and hence \code{V} should be a \mjeqn{6 \times 6}{6x6} variance-covariance matrix of these elements in this order. The variance-covariance matrix of the standardized regression coefficients (i.e., \mjeqn{\mbox{Var}[b]}{Var[b]}) is then computed as a function of \code{V} as described in Becker (1992) using the multivariate delta method. The standard errors are then again given by the square root of the diagonal elements of \mjeqn{\mbox{Var}[b]}{Var[b]}. Test statistics (in this case, z-statistics) and the corresponding p-values can then be computed in the usual manner.

   In case \mjseqn{R} is the result of a meta-analysis of Fisher r-to-z transformed correlation coefficients (and hence \code{V} is then the corresponding variance-covariance matrix of these pooled transformed coefficients), one should set argument \code{ztor=TRUE}, so that the appropriate back-transformation is then applied to \code{R} (and \code{V}) within the function.

   Finally, \mjseqn{R} may be a covariance matrix based on a meta-analysis (e.g., the estimated variance-covariance matrix of the random effects in a multivariate model). In this case, one should set \code{cov=TRUE} and \code{V} should again be the variance-covariance matrix of the elements in \mjseqn{R}, but now including the diagonal. Hence, if \mjseqn{R} is a \mjeqn{4 \times 4}{4x4} matrix\ifelse{text}{,}{ of the form: \mjtdeqn{\left[ \begin{array}{cccc} \tau_1^2 & & & \\\ \tau_{21} & \tau_2^2 & & \\\ \tau_{31} & \tau_{32} & \tau_3^2 & \\\ \tau_{41} & \tau_{42} & \tau_{43} & \tau_4^2 \end{array} \right]}{\begin{bmatrix} \tau_1^2 & & & \\\\\ \tau_{21} & \tau_2^2 & & \\\\\ \tau_{31} & \tau_{32} & \tau_3^2 & \\\\\ \tau_{41} & \tau_{42} & \tau_{43} & \tau_4^2 \end{bmatrix}}{}} then the elements are \mjseqn{\tau^2_1}, \mjseqn{\tau_{21}}, \mjseqn{\tau_{31}}, \mjseqn{\tau_{41}}, \mjseqn{\tau^2_2}, \mjseqn{\tau_{32}}, \mjseqn{\tau_{42}}, \mjseqn{\tau^2_3}, \mjseqn{\tau_{43}}, and \mjseqn{\tau^2_4}, and hence \code{V} should be a \mjeqn{10 \times 10}{10x10} variance-covariance matrix of these elements in this order. Argument \code{means} can then again be used to specify the means of the variables.
}
\value{
   An object of class \code{"matreg"}. The object is a list containing the following components:
   \item{tab}{a data frame with the estimated model coefficients, standard errors, test statistics, degrees of freedom (only for t-tests), p-values, and lower/upper confidence interval bounds.}
   \item{vb}{the variance-covariance matrix of the estimated model coefficients.}
   \item{\dots}{some additional elements/values.}

   The results are formatted and printed with the \code{\link[=print.matreg]{print}} function. Extractor functions include \code{\link[=coef.matreg]{coef}} and \code{\link[=vcov.matreg]{vcov}}.
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

   Becker, B. J., & Aloe, A. (2019). Model-based meta-analysis and related approaches. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta-analysis} (3rd ed., pp. 339--363). New York: Russell Sage Foundation.
}
\seealso{
   \code{\link{rma.mv}} for a function to meta-analyze multiple correlation coefficients that can be used to construct an \mjseqn{R} matrix.

   \code{\link{rcalc}} for a function to construct the variance-covariance matrix of dependent correlation coefficients.
}
\examples{
############################################################################

### first an example unrelated to meta-analysis, simply demonstrating that
### one can obtain the same results from lm() and matreg()

### fit a regression model with lm() to the 'mtcars' dataset
res <- lm(mpg ~ hp + wt + am, data=mtcars)
summary(res)

### covariance matrix of the dataset
S <- cov(mtcars)

### fit the same regression model using matreg()
res <- matreg(y="mpg", x=c("hp","wt","am"), R=S, cov=TRUE,
              means=colMeans(mtcars), n=nrow(mtcars))
summary(res)

### copy the 'mtcars' dataset to 'dat' and standardize all variables
dat <- mtcars
dat[] <- scale(dat)

### fit a regression model with lm() to obtain standardized regression coefficients ('betas')
res <- lm(mpg ~ 0 + hp + wt + am, data=dat)
summary(res)

### correlation matrix of the dataset
R <- cor(mtcars)

### fit the same regression model using matreg()
res <- matreg(y="mpg", x=c("hp","wt","am"), R=R, n=nrow(mtcars))
summary(res)

### note: the standard errors of the betas should not be used to construct CIs
### as they assume that the null hypothesis (H0: beta_j = 0) is true

### construct the var-cov matrix of correlations in R
V <- rcalc(R, ni=nrow(mtcars))$V

### fit the same regression model using matreg() but now supply V
res <- matreg(y="mpg", x=c("hp","wt","am"), R=R, V=V)
summary(res)

### the standard errors computed in this way can now be used to construct
### CIs for the betas (here, the difference is relatively small)

############################################################################

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
matreg(1, 2:4, R=R, V=vcov(res))

### can also specify variable names
matreg("perf", c("acog","asom","conf"), R=R, V=vcov(res))

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
matreg(1, 2:4, R=R, V=vcov(res), ztor=TRUE)
}

############################################################################

### a different example based on van Houwelingen et al. (2002)

### create dataset in long format
dat.long <- to.long(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.colditz1994)
dat.long <- escalc(measure="PLO", xi=out1, mi=out2, data=dat.long)
dat.long$tpos <- dat.long$tneg <- dat.long$cpos <- dat.long$cneg <- NULL
levels(dat.long$group) <- c("CON", "EXP")

### fit bivariate model
res <- rma.mv(yi, vi, mods = ~ group - 1, random = ~ group | trial, struct="UN",
              data=dat.long, method="ML")
res

### regression of log(odds)_EXP on log(odds)_CON
matreg(y=2, x=1, R=res$G, cov=TRUE, means=coef(res), n=res$g.levels.comb.k)

### but the SE of the CON coefficient is not computed correctly, since above we treat res$G as if
### it was a var-cov matrix computed from raw data based on res$g.levels.comb.k (= 13) data points

### fit bivariate model and get the var-cov matrix of the estimates in res$G
res <- rma.mv(yi, vi, mods = ~ group - 1, random = ~ group | trial, struct="UN",
              data=dat.long, method="ML", cvvc="varcov", control=list(nearpd=TRUE))

### now use res$vvc as the var-cov matrix of the estimates in res$G
matreg(y=2, x=1, R=res$G, cov=TRUE, means=coef(res), V=res$vvc)
}
\keyword{models}
