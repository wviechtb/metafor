\name{residuals.rma}
\alias{residuals.rma}
\alias{rstandard.rma.uni}
\alias{rstandard.rma.mh}
\alias{rstandard.rma.mv}
\alias{rstandard.rma.peto}
\alias{rstudent.rma.uni}
\alias{rstudent.rma.mh}
\alias{rstudent.rma.peto}
\title{Residual Values based on 'rma' Objects}
\description{
   The \code{residuals}, \code{rstandard}, and \code{rstudent} functions can be used to compute residuals, corresponding standard errors, and standardized residuals for models fitted with the \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}, and \code{\link{rma.mv}} functions.
}
\usage{
\method{residuals}{rma}(object, \dots)

\method{rstandard}{rma.uni}(model, digits, \dots)
\method{rstandard}{rma.mh}(model, digits, \dots)
\method{rstandard}{rma.mv}(model, digits, \dots)
\method{rstandard}{rma.peto}(model, digits, \dots)

\method{rstudent}{rma.uni}(model, digits, progbar=FALSE, \dots)
\method{rstudent}{rma.mh}(model, digits, progbar=FALSE, \dots)
\method{rstudent}{rma.peto}(model, digits, progbar=FALSE, \dots)
}
\arguments{
   \item{object}{an object of class \code{"rma"} (for \code{residuals}).}
   \item{model}{an object of class \code{"rma.uni"}, \code{"rma.mh"}, \code{"rma.peto"}, or \code{"rma.mv"} (for \code{rstandard} and \code{rstudent}).}
   \item{digits}{integer specifying the number of decimal places to which the printed results should be rounded (if unspecified, the default is to take the value from the object).}
   \item{progbar}{logical indicating whether a progress bar should be shown (the default is \code{FALSE}).}
   \item{\dots}{other arguments.}
}
\details{
   The observed residuals (obtained with \code{residuals}) are simply equal to the \sQuote{observed - fitted} values.

   Dividing the observed residuals by their corresponding standard errors yields (internally) standardized residuals. These can be obtained with \code{rstandard}.

   The \code{rstudent} function calculates externally standardized residuals (studentized deleted residuals). The externally standardized residual for the \eqn{i}th case is obtained by deleting the \eqn{i}th case from the dataset, fitting the model based on the remaining cases, calculating the predicted value for the \eqn{i}th case based on the fitted model, taking the difference between the observed and the predicted value for the \eqn{i}th case (the deleted residual), and then standardizing the deleted residual. The standard error of the deleted residual is equal to the square root of the sampling variance of the \eqn{i}th case plus the variance of the predicted value plus the amount of (residual) heterogeneity from the fitted model (for fixed-effects models, this last part is always equal to zero).

   If a particular study fits the model, its standardized residual follows (asymptotically) a standard normal distribution. A large standardized residual for a study therefore may suggest that the study does not fit the assumed model (i.e., it may be an outlier).

   See also \code{\link{influence.rma.uni}} for other leave-one-out diagnostics that are useful for detecting influential cases in models fitted with the \code{\link{rma.uni}} function.
}
\value{
   Either a vector with the observed residuals (for \code{residuals}) or an object of class \code{"list.rma"}, which is a list containing the following components:
   \item{resid}{observed residuals (for \code{rstandard}) or deleted residuals (for \code{rstudent}).}
   \item{se}{corresponding standard errors.}
   \item{z}{standardized residuals (internally standardized for \code{rstandard} or externally standardized for \code{rstudent}).}

   The \code{"list.rma"} object is formated and printed with \code{\link{print.list.rma}}.
}
\note{
   Right now, the externally standardized residuals are calculated by refitting the model \eqn{k} times. Depending on how large \eqn{k} is, it may take a few moments to finish the calculations.

   It may not be possible to fit the model after deletion of the \eqn{i}th case from the dataset. This will result in \code{NA} values for that case when calling \code{rstudent}.

   For objects of class \code{"rma.mh"} and \code{"rma.mh"}, \code{rstandard} computes Pearson (or semi-standardized) residuals.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \cr
   package website: \url{http://www.metafor-project.org/} \cr
   author homepage: \url{http://www.wvbauer.com/}
}
\references{
   Hedges, L. V., & Olkin, I. (1985). \emph{Statistical methods for meta-analysis}. San Diego, CA: Academic Press.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.

   Viechtbauer, W., & Cheung, M. W.-L. (2010). Outlier and influence diagnostics for meta-analysis. \emph{Research Synthesis Methods}, \bold{1}, 112--125.
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}, \code{\link{rma.glmm}}, \code{\link{rma.mv}}, \code{\link{influence.rma.uni}}
}
\examples{
### meta-analysis of the log relative risks using a random-effects model
res <- rma(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
rstudent(res)

### mixed-effects model with absolute latitude as a moderator
res <- rma(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, mods = ~ ablat,
           data=dat.bcg)
rstudent(res)
}
\keyword{models}
