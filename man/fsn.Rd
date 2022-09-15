\name{fsn}
\alias{fsn}
\title{Fail-Safe N Analysis (File Drawer Analysis)}
\description{
   Function to compute the fail-safe N (also called a file drawer analysis).
}
\usage{
fsn(yi, vi, sei, data, type="Rosenthal", alpha=.05,
    target, weighted=FALSE, subset, digits, \dots)
}
\arguments{
   \item{yi}{vector with the observed effect sizes or outcomes.}
   \item{vi}{vector with the corresponding sampling variances.}
   \item{sei}{vector with the corresponding standard errors (note: only one of the two, \code{vi} or \code{sei}, needs to be specified).}
   \item{data}{optional data frame containing the variables given to the arguments above.}
   \item{type}{character string to specify the method to use for the calculation of the fail-safe N. Possible options are \code{"Rosenthal"} (the default), \code{"Orwin"}, or \code{"Rosenberg"}. See \sQuote{Details}.}
   \item{alpha}{target alpha level to use for the Rosenthal and Rosenberg methods (the default is .05).}
   \item{target}{target average effect size or outcome to use for the Orwin method. If undefined, then the target average effect size or outcome will be equal to the observed average effect size or outcome divided by 2.}
   \item{weighted}{logical to specify whether Orwin's method should be based on unweighted (the default) or weighted averages (the default is \code{FALSE}).}
   \item{subset}{optional (logical or numeric) vector to specify the subset of studies that should be used for the calculations.}
   \item{digits}{optional integer to specify the number of decimal places to which the printed results should be rounded.}
   \item{\dots}{other arguments.}
}
\details{
   The function can be used in combination with any of the usual effect sizes / outcome measures used in meta-analyses (e.g., log risk ratios, log odds ratios, risk differences, mean differences, standardized mean differences, raw correlation coefficients, correlation coefficients transformed with Fisher's r-to-z transformation), or, more generally, any set of estimates (with corresponding sampling variances) one would like to analyze. Simply specify the observed outcomes via the \code{yi} argument and the corresponding sampling variances via the \code{vi} argument (instead of specifying \code{vi}, one can specify the standard errors via the \code{sei} argument). The \code{\link{escalc}} function can be used to compute a wide variety of effect sizes / outcome measures (and the corresponding sampling variances) based on summary statistics.

   The Rosenthal method (sometimes called a \sQuote{file drawer analysis}) calculates the number of studies averaging null results that would have to be added to the given set of observed outcomes to reduce the combined significance level (p-value) to a particular alpha level (e.g., .05). The calculation is based on Stouffer's method to combine p-values and is described in Rosenthal (1979).

   The Orwin method calculates the number of studies averaging null results that would have to be added to the given set of observed outcomes to reduce the (unweighted or weighted) average outcome to a target value (as specified via the \code{target} argument). The method is described in Orwin (1983). If \code{weighted=FALSE} (the default), the method does not require (or makes use) of \code{vi} (or \code{sei}), so these arguments are then not relevant for this method. If the \code{target} argument is not specified, then the target average outcome will be equal to the observed average outcome divided by 2 (which is quite arbitrary). One should really set \code{target} to a value that reflects an outcome one would consider practically irrelevant. Note that if \code{target} has the opposite sign as the actually observed average outcome, then its sign is automatically flipped.

   The Rosenberg method calculates the number of studies averaging null results that would have to be added to the given set of observed outcomes to reduce the significance level (i.e., p-value) of the weighted average outcome (based on an equal-effects model) to a particular alpha level (e.g., .05). The method is described in Rosenberg (2005).

   If the combined/observed significance level is above the specified alpha level (for \code{type = "Rosenthal"} or \code{type = "Rosenberg"}) or if the observed average outcome is below the target average outcome (for \code{type = "Orwin"}), then the fail-safe N value will be 0.
}
\value{
   An object of class \code{"fsn"}. The object is a list containing the following components:
   \item{type}{the method used.}
   \item{fsnum}{the calculated fail-safe N.}
   \item{alpha}{the specified alpha level.}
   \item{pval}{the p-value of the observed results. \code{NA} for the Orwin method.}
   \item{meanes}{the average outcome of the observed results. \code{NA} for the Rosenthal method.}
   \item{target}{the target value. \code{NA} for the Rosenthal and Rosenberg methods.}

   The results are formatted and printed with the \code{\link[=print.fsn]{print}} function.
}
\note{
   For the Rosenberg method, the p-value is calculated based on a standard normal distribution (instead of a t-distribution, as suggested by Rosenberg, 2005).
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Rosenthal, R. (1979). The "file drawer problem" and tolerance for null results. \emph{Psychological Bulletin}, \bold{86}(3), 638--641. \verb{https://doi.org/10.1037/0033-2909.86.3.638}

   Orwin, R. G. (1983). A fail-safe N for effect size in meta-analysis. \emph{Journal of Educational Statistics}, \bold{8}(2), 157--159. \verb{https://doi.org/10.3102/10769986008002157}

   Rosenberg, M. S. (2005). The file-drawer problem revisited: A general weighted method for calculating fail-safe numbers in meta-analysis. \emph{Evolution}, \bold{59}(2), 464--468. \verb{https://doi.org/10.1111/j.0014-3820.2005.tb01004.x}

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{ranktest}} for the rank correlation test, \code{\link{regtest}} for the regression test, \code{\link{trimfill}} for the trim and fill method, \code{\link{tes}} for the test of excess significance, and \code{\link{selmodel}} for selection models.
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### fit equal-effects model
rma(yi, vi, data=dat, method="EE")

### fail-safe N computations
fsn(yi, vi, data=dat)
fsn(yi, data=dat, type="Orwin", target=log(0.95)) # target corresponds to a 5% risk reduction
fsn(yi, vi, data=dat, type="Orwin", weighted=TRUE, target=log(0.95))
fsn(yi, vi, data=dat, type="Rosenberg")
}
\keyword{htest}
