\name{fsn}
\alias{fsn}
\title{Fail-Safe N Analysis (File Drawer Analysis)}
\description{
   Function to compute the fail-safe N (also called a file drawer analysis). \loadmathjax
}
\usage{
fsn(yi, vi, sei, data, type="Rosenthal", alpha=.05, target,
    method="REML", exact=FALSE, subset, verbose=FALSE, digits, \dots)
}
\arguments{
   \item{yi}{vector with the observed effect sizes or outcomes.}
   \item{vi}{vector with the corresponding sampling variances.}
   \item{sei}{vector with the corresponding standard errors (note: only one of the two, \code{vi} or \code{sei}, needs to be specified).}
   \item{data}{optional data frame containing the variables given to the arguments above.}
   \item{type}{character string to specify the type of method to use for the calculation of the fail-safe N. Possible options are \code{"Rosenthal"} (the default), \code{"Orwin"}, \code{"Rosenberg"}, or \code{"General"}. Can be abbreviated.}
   \item{alpha}{target alpha level for the Rosenthal, Rosenberg, and General methods (the default is .05).}
   \item{target}{target average effect size or outcome for the Orwin and General methods.}
   \item{method}{character string to specify the model fitting method for \code{type="General"} (the default is \code{"REML"}). See \code{\link{rma.uni}} for options.}
   \item{exact}{logical to indicate whether the General method should be based on exact (but slower) calculations or approximate (but faster) ones.}
   \item{subset}{optional (logical or numeric) vector to specify the subset of studies that should be used for the calculations.}
   \item{verbose}{logical to specify whether output should be generated on the progress of the calculations for \code{type="General"} (the default is \code{FALSE}).}
   \item{digits}{optional integer to specify the number of decimal places to which the printed results should be rounded.}
   \item{\dots}{other arguments.}
}
\details{
   The function can be used to calculate the \sQuote{fail-safe N}, that is, the minimum number of studies averaging null results that would have to be added to a given set of \mjseqn{k} studies to change the conclusions of a meta-analysis. If this number is small (in relation to the actual number of studies), then this indicates that the results based on the observed studies are not robust to publication bias (of the form assumed by the method, that is, where a set of studies averaging null results is missing). The method is also called a \sQuote{file drawer analysis} as it assumes that there is a set of studies averaging null results hiding in file drawers, which can overturn the findings from a meta-analysis. Note that \emph{the method does not estimate the number of missing studies}, only how many studies must be hiding in file drawers for the findings to be overturned.

   There are various types of methods that are all based on the same principle, which are described in more detail further below.

   The function can be used in combination with any of the usual effect sizes / outcome measures used in meta-analyses (e.g., log risk ratios, log odds ratios, risk differences, mean differences, standardized mean differences, raw correlation coefficients, correlation coefficients transformed with Fisher's r-to-z transformation), or, more generally, any set of estimates (with corresponding sampling variances) one would like to analyze. Simply specify the observed outcomes of a given set of studies via the \code{yi} argument and the corresponding sampling variances via the \code{vi} argument (instead of specifying \code{vi}, one can specify the standard errors via the \code{sei} argument). The \code{\link{escalc}} function can be used to compute a wide variety of effect sizes / outcome measures (and the corresponding sampling variances) based on summary statistics.

   \subsection{Rosenthal Method}{

      The Rosenthal method (\code{type="Rosenthal"}) calculates the minimum number of studies averaging null results that would have to be added to a given set of studies to reduce the (one-tailed) combined significance level (i.e., p-value) to a particular alpha level, which can be specified via the \code{alpha} argument (.05 by default). The calculation is based on Stouffer's method for combining p-values and is described in Rosenthal (1979). Note that the method is primarily of interest for historical reasons, but the other methods described below are more closely aligned with the way meta-analyses are typically conducted in practice.

   }

   \subsection{Orwin Method}{

      The Orwin method (\code{type="Orwin"}) calculates the minimum number of studies averaging null results that would have to be added to a given set of studies to reduce the (unweighted or weighted) average effect size / outcome to a target value (as specified via the \code{target} argument). The method is described in Orwin (1983). When \code{vi} (or \code{sei}) is not specified, the method is based on the unweighted average of the effect sizes / outcomes; otherwise, the method uses the inverse-variance weighted average. If the \code{target} argument is not specified, then the target value will be equal to the observed average effect size / outcome divided by 2 (which is entirely arbitrary and will always lead to a fail-safe N number that is equal to \mjseqn{k}). One should really set \code{target} to a value that reflects an effect size / outcome one would consider practically irrelevant. Note that if \code{target} has the opposite sign as the actually observed average, then its sign is automatically flipped.

   }

   \subsection{Rosenberg Method}{

      The Rosenberg method (\code{type="Rosenberg"}) calculates the minimum number of studies averaging null results that would have to be added to a given set of studies to reduce the significance level (i.e., p-value) of the average effect size / outcome (as estimated based on an equal-effects model) to a particular alpha level, which can be specified via the \code{alpha} argument (.05 by default). The method is described in Rosenberg (2005). Note that the p-value is calculated based on a standard normal distribution (instead of a t-distribution, as suggested by Rosenberg, 2005), but the difference is typically negligible.

   }

   \subsection{General Method}{

      This method is a generalization of the methods by Orwin and Rosenberg. By default (i.e., when \code{target} is not specified), it calculates the minimum number of studies averaging null results that would have to be added to a given set of studies to reduce the significance level (i.e., p-value) of the average effect size / outcome (as estimated based on a random-effects model) to a particular alpha level, which can be specified via the \code{alpha} argument (.05 by default). The random-effects model is fitted with the method specified via the \code{method} argument, with \code{"REML"} being the default (see \code{\link{rma.uni}} for other options). When setting \code{method="EE"}, then an equal-effects model is used, which yields (essentially) identical results as Rosenberg's method.

      If \code{target} is specified, then the method calculates the minimum number of studies averaging null results that would have to be added to a given set of studies to reduce the average effect size / outcome (as estimated based on a random-effects model) to a target value (as specified via the \code{target} argument). When setting \code{method="EE"}, then an equal-effects model is used, which yields (essentially) identical results as Orwin's method with inverse-variance weights.

      The method uses an iterative algorithm for calculating the fail-safe N, which can be computationally expensive especially when N is large. By default, the method uses approximate (but faster) calculations, but when setting \code{exact=TRUE}, the method uses exact (but slower) calculations. The difference between the two is typically negligible. If N is larger than \mjseqn{10^7}, then the calculated number is given as \code{>1e+07}.

   }

}
\value{
   An object of class \code{"fsn"}. The object is a list containing the following components (some of which may be \code{NA} if they are not applicable to the chosen method):
   \item{type}{the type of method used.}
   \item{fsnum}{the calculated fail-safe N.}
   \item{est}{the average effect size / outcome based on the observed studies.}
   \item{tau2}{the estimated amount of heterogeneity based on the observed studies.}
   \item{pval}{the p-value of the observed results.}
   \item{alpha}{the specified target alpha level.}
   \item{target}{the target average effect size / outcome.}
   \item{est.fsn}{the average effect size / outcome when combining the observed studies with those in the file drawer.}
   \item{tau2}{the estimated amount of heterogeneity when combining the observed studies with those in the file drawer.}
   \item{pval}{the p-value when combining the observed studies with those in the file drawer.}
   \item{\dots}{some additional elements/values.}

   The results are formatted and printed with the \code{\link[=print.fsn]{print}} function.
}
\note{
   If the significance level of the observed studies is already above the specified alpha level or if the average effect size / outcome of the observed studies is already below the target average effect size / outcome, then the fail-safe N value is zero.
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
   \code{\link{regtest}} for the regression test, \code{\link{ranktest}} for the rank correlation test, \code{\link{trimfill}} for the trim and fill method, \code{\link{tes}} for the test of excess significance, and \code{\link{selmodel}} for selection models.
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### fit equal-effects model
rma(yi, vi, data=dat, method="EE")

### fail-safe N computations
fsn(yi, vi, data=dat)
fsn(yi, data=dat, type="Orwin", target=log(0.95)) # target corresponds to a 5% risk reduction
fsn(yi, vi, data=dat, type="Orwin", target=log(0.95)) # Orwin's method with 1/vi weights
fsn(yi, vi, data=dat, type="General", target=log(0.95), method="EE") # like Orwin's method
fsn(yi, vi, data=dat, type="Rosenberg")
fsn(yi, vi, data=dat, type="General", method="EE") # like Rosenberg's method
fsn(yi, vi, data=dat, type="General") # based on a random-effects model
fsn(yi, vi, data=dat, type="General", target=log(0.95)) # based on a random-effects model
}
\keyword{htest}
