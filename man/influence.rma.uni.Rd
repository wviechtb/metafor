\name{influence.rma.uni}
\alias{influence.rma.uni}
\alias{print.infl.rma.uni}
\alias{cooks.distance.rma.uni}
\alias{dfbetas.rma.uni}
\alias{hatvalues.rma.uni}
\title{Outlier and Influential Case Diagnostics for 'rma.uni' Objects}
\description{
   The functions compute various outlier and influential case diagnostics (some of which indicate the influence of deleting one case/study at a time on the model fit and the fitted/residual values) for objects of class \code{"rma.uni"}.
}
\usage{
\method{influence}{rma.uni}(model, digits, progbar=FALSE, \dots)

\method{print}{infl.rma.uni}(x, digits, \dots)

\method{cooks.distance}{rma.uni}(model, progbar=FALSE, \dots)
\method{dfbetas}{rma.uni}(model, progbar=FALSE, \dots)
\method{hatvalues}{rma.uni}(model, type="diagonal", \dots)
}
\arguments{
   \item{model}{an object of class \code{"rma.uni"}.}
   \item{x}{an object of class \code{"infl.rma.uni"} (for \code{print}).}
   \item{digits}{integer specifying the number of decimal places to which the printed results should be rounded (if unspecified, the default is to take the value from the object).}
   \item{progbar}{logical indicating whether a progress bar should be shown (the default is \code{FALSE}).}
   \item{type}{character string indicating whether to return only the diagonal of the hat matrix (\code{"diagonal"}) or the entire hat matrix (\code{"matrix"}).}
   \item{\dots}{other arguments.}
}
\details{
   The \code{influence} function calculates the following leave-one-out diagnostics for each study:
   \itemize{
   \item externally standardized residual,
   \item DFFITS value,
   \item Cook's distance,
   \item covariance ratio,
   \item the leave-one-out amount of (residual) heterogeneity,
   \item the leave-one-out test statistic for the test of (residual) heterogeneity,
   \item DFBETAS value(s).
   } The diagonal elements of the hat matrix and the weights (in \%) given to the observed effects or outcomes during the model fitting are also provided (except for their scaling, the hat values and weights are the same for models without moderators, but will differ when moderators are included).

   For details on externally standardized residuals, see \code{\link{rstudent.rma.uni}}.

   The DFFITS value essentially indicates how many standard deviations the predicted (average) effect for the \eqn{i}th study changes after excluding the \eqn{i}th study from the model fitting.

   Cook's distance can be interpreted as the Mahalanobis distance between the entire set of predicted values once with the \eqn{i}th study included and once with the \eqn{i}th study excluded from the model fitting.

   The covariance ratio is defined as the determinant of the variance-covariance matrix of the parameter estimates based on the dataset with the \eqn{i}th study removed divided by the determinant of the variance-covariance matrix of the parameter estimates based on the complete dataset. A value below 1 therefore indicates that removal of the \eqn{i}th study yields more precise estimates of the model coefficients.

   The leave-one-out amount of (residual) heterogeneity is the estimated value of \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}} based on the dataset with the \eqn{i}th study removed. Note that this is always equal to \code{0} for fixed-effects models.

   Similarly, the leave-one-out test statistic for the test of (residual) heterogeneity is the value of the test statistic of the test for (residual) heterogeneity calculated based on the dataset with the \eqn{i}th study removed.

   Finally, the DFBETAS value(s) essentially indicate(s) how many standard deviations the estimated coefficient(s) change(s) after excluding the \eqn{i}th study from the model fitting.

   A study may be considered to be \sQuote{influential} if at least one of the following is true:
   \itemize{
   \item The absolute DFFITS value is larger than \eqn{3\sqrt{p/(k-p)}}{3\sqrt(p/(k-p))}, where \eqn{p} is the number of model coefficients and \eqn{k} the number of studies.
   \item The lower tail area of a chi-square distribution with \eqn{p} degrees of freedom cut off by the Cook's distance is larger than 50\%.
   \item The hat value is larger than \eqn{3(p/k)}.
   \item Any DFBETAS value is larger than \eqn{1}.
   } Studies which are considered influential with respect to any of these measures are marked with an asterisk. Note that the chosen cut-offs are (somewhat) arbitrary. Substantively informed judgment should always be used when examining the influence of each study on the results.
}
\value{
   An object of class \code{"infl.rma.uni"}. The object is a list containing the following components:
   \item{inf}{a data frame with columns equal to the externally standardized residuals, DFFITS values, Cook's distances, covariance ratios, leave-one-out \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}} estimates, leave-one-out (residual) heterogeneity test statistics, hat values, and weights.}
   \item{dfbs}{a data frame with columns equal to the DFBETAS values.}
   \item{is.infl}{a logical vector indicating whether a particular study is considered to be \sQuote{influential}.}
   \item{\dots}{some additional elements/values.}

   The results are printed with \code{print.infl.rma.uni} and plotted with \code{\link{plot.infl.rma.uni}}.
}
\note{
   Right now, the leave-one-out diagnostics are calculated by refitting the model \eqn{k} times. Depending on how large \eqn{k} is, it may take a few moments to finish the calculations. There are shortcuts for calculating at least some of these values without refitting the model each time, but these are currently not implemented (and may not exist for all of the leave-one-out diagnostics calculated by the function).

   It may not be possible to fit the model after deletion of the \eqn{i}th study from the dataset. This will result in \code{NA} values for that study.

   Certain relationships between the leave-one-out diagnostics and the (internally or externally) standardized residuals (Belsley, Kuh, & Welsch, 1980; Cook & Weisberg, 1982) no longer hold for the meta-analytic models. Maybe there are other relationships. These remain to be determined.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \cr
   package website: \url{http://www.metafor-project.org/} \cr
   author homepage: \url{http://www.wvbauer.com/}
}
\references{
   Belsley, D. A., Kuh, E., & Welsch, R. E. (1980). \emph{Regression diagnostics}. New York: Wiley.

   Cook, R. D., & Weisberg, S. (1982). \emph{Residuals and influence in regression}. London: Chapman and Hall.

   Hedges, L. V., & Olkin, I. (1985). \emph{Statistical methods for meta-analysis}. San Diego, CA: Academic Press.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.

   Viechtbauer, W., & Cheung, M. W.-L. (2010). Outlier and influence diagnostics for meta-analysis. \emph{Research Synthesis Methods}, \bold{1}, 112--125.
}
\seealso{
   \code{\link{plot.infl.rma.uni}}, \code{\link{rstudent.rma.uni}}, \code{\link{weights.rma.uni}}
}
\examples{
### meta-analysis of the log risk ratios using a mixed-effects model
### with two moderators (absolute latitude and publication year)
res <- rma(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg,
           mods = ~ ablat + year, data=dat.bcg)
influence(res)
plot(influence(res))

cooks.distance(res)
dfbetas(res)
hatvalues(res)
}
\keyword{models}
