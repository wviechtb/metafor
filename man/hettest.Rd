\name{hettest}
\alias{hettest}
\alias{print.hettest}
\title{Test for Heteroscedastic Heterogeneity}
\description{
   Function to test if the heterogeneity in the true effects is heteroscedastic. \loadmathjax
}
\usage{
hettest(x, vi, sei, subset, data, method="REML",
        test="score", boot=TRUE, progbar=TRUE, digits, control, \dots)

\method{print}{hettest}(x, digits=x$digits, \dots)
}
\arguments{
   \item{x}{a vector with the effect size estimates or an object of class \code{"rma"}.}
   \item{vi}{vector with the corresponding sampling variances. Ignored if \code{x} is an object of class \code{"rma"}.}
   \item{sei}{vector with the corresponding standard errors (note: only one of the two, \code{vi} or \code{sei}, needs to be specified).}
   \item{subset}{optional (logical or numeric) vector to specify the subset of studies that should be used for the calculation. Ignored if \code{x} is an object of class \code{"rma"}.}
   \item{data}{optional data frame containing the variables given to the arguments above.}
   \item{method}{character string to specify whether the test should be based on maximum likelihood (\code{"ML"}) or restricted maximum likelihood (\code{"REML"}) estimation. The default is \code{"REML"}. Ignored if \code{x} is an object of class \code{"rma"}.}
   \item{test}{character string (case insensitive) to specify which test to run. Either \code{"lrt"}, \code{"wald"}, or \code{"score"} (the default).}
   \item{boot}{logical to specify whether bootstrapping should be used (the default is \code{TRUE}). Can also be a number to specify the number of bootstrap samples to take.}
   \item{progbar}{logical to specify whether a progress bar should be shown when using bootstrapping (the default is \code{TRUE}).}
   \item{digits}{optional integer to specify the number of decimal places to which the printed results should be rounded.}
   \item{control}{optional list of control values. See \sQuote{Note}.}
   \item{\dots}{other arguments.}
}
\details{

   The standard random-effects model (that can be fitted with the \link{rma.uni} function) assumes that the heterogeneity in the true effects is homoscedastic. This function can be used to test this assumption.

   Via argument \code{x}, one can either pass a vector with the effect size estimates to the function (and then the corresponding sampling variances via \code{vi} or the standard errors via \code{sei}) or an object of class \code{"rma"}. When passing a model object to argument \code{x}, the model must be a standard random-effects model that was fitted with either maximum likelihood (ML) or restricted maximum likelihood (REML) estimation.

   When \code{x} is a vector with the effect size estimates, then argument \code{method} can be used to select whether the test should be based on ML or REML estimation. The latter is the default. On the other hand, if \code{x} is a model object, then the estimation method that was used for fitting the model determines on which estimation method the test is based.

   The test to run is determined by the \code{test} argument and can either be:

   \itemize{
   \item \code{test="lrt"} for a likelihood ratio test,
   \item \code{test="wald"} for a Wald-type test,
   \item \code{test="score"} for a score test.
   }

   The likelihood ratio test compares the fit of the standard random-effects model with a saturated random-effects model where each study has a study-specific \mjseqn{\tau^2_i} value. The Wald-type test directly tests all \mjseqn{k-1} contrasts between the estimated \mjseqn{\tau^2_i} values. Finally, the score test evaluates whether allowing study-specific \mjseqn{\tau^2_i} values would improve model fit based on the score vector and the information matrix evaluated under the standard random-effects model, without requiring estimation of the saturated model. This tends to make this test more stable, yielding better control of the Type I error rate, and therefore is the default option.

   However, the assumed null distribution for all of these tests (a chi-square distribution with \mjseqn{k-1} degrees of freedom) is not correct. Therefore, when \code{boot=TRUE} (the default), the function uses parametric bootstrapping to obtain a p-value for a given test. By default, 1000 bootstrap samples are taken, but one can also specify the number of samples to take via the \code{boot} argument. When \code{progbar=TRUE} (the default), a progress bar is shown to indicate the progress of the bootstrapping.

}
\value{
   An object of class \code{"hettest"}. The object is a list containing the following components:
   \item{statistic}{the test statistic.}
   \item{df}{the degrees of freedom.}
   \item{pval}{the p-value.}
   \item{tau2i}{the estimated \mjseqn{\tau^2_i} values from the saturated random-effects model.}
   \item{se.tau2i}{the corresponding standard errors (only when running a Wald-type test).}
   \item{statistic.boot}{the bootstrap values of the test statistic (only when \code{boot=TRUE}).}
   \item{\dots}{some additional elements/values.}

   The results are formatted and printed with the \code{print} function.
}
\note{
   The saturated random-effects model is fitted with an iterative algorithm, alternating between the computation of the \mjseqn{\tau^2_i} estimates and the estimate of \mjseqn{\mu}. By default, the initial \mjseqn{\tau^2_i} values are set equal to the \mjseqn{\tau^2} estimate from the standard random-effects model. This process continues for a maximum of 1000 iterations until the change in all of the estimated \mjseqn{\tau^2_i} values is smaller than \mjeqn{10^{-8}}{10^(-8)} from one iteration to the next. The starting \mjseqn{\tau^2_i} values, threshold, and maximum number of iterations can be adjusted via the \code{control} argument by setting \code{control=list(tau2i.init=values, threshold=value, maxiter=value)}.

   When bootstrapping, the iterative algorithm must be run for each bootstrap sample. In rare cases, the algorithm may not converge. The output shows for how many bootstrap samples the calculations could be completed (this should usually be \sQuote{1000 / 1000}).
}
\author{
   Wolfgang Viechtbauer (\email{wvb@metafor-project.org}, \url{https://www.metafor-project.org}).
}
\keyword{htest}
\examples{
### copy data from Bangert-Drowns et al. (2004) into 'dat'
dat <- dat.bangertdrowns2004

### fit a standard random-effects model
res <- rma(yi, vi, data=dat)
res

\dontrun{
### test for heteroscedastic heterogeneity
set.seed(1234)
hettest(res)
}
}
