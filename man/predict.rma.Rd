\name{predict.rma}
\alias{predict}
\alias{predict.rma}
\alias{predict.rma.ls}
\title{Predicted Values for 'rma' Objects}
\description{
   The function computes predicted values, corresponding standard errors, confidence intervals, and prediction intervals for objects of class \code{"rma"}. \loadmathjax
}
\usage{
\method{predict}{rma}(object, newmods, intercept, tau2.levels, gamma2.levels, addx=FALSE,
        level, digits, transf, targs, vcov=FALSE, \dots)

\method{predict}{rma.ls}(object, newmods, intercept, addx=FALSE, newscale, addz=FALSE,
        level, digits, transf, targs, vcov=FALSE, \dots)
}
\arguments{
   \item{object}{an object of class \code{"rma"} or \code{"rma.ls"}.}
   \item{newmods}{optional vector or matrix to specify the values of the moderator values for which the predicted values should be calculated. See \sQuote{Details}.}
   \item{intercept}{logical to specify whether the intercept should be included when calculating the predicted values for \code{newmods}. If unspecified, the intercept is automatically added when the original model also included an intercept.}
   \item{tau2.levels}{vector to specify the levels of the inner factor when computing prediction intervals. Only relevant for models of class \code{"rma.mv"} (see \code{\link{rma.mv}}) and when the model includes more than a single \mjseqn{\tau^2} value. See \sQuote{Details}.}
   \item{gamma2.levels}{vector to specify the levels of the inner factor when computing prediction intervals. Only relevant for models of class \code{"rma.mv"} (see \code{\link{rma.mv}}) and when the model includes more than a single \mjseqn{\gamma^2} value. See \sQuote{Details}.}
   \item{addx}{logical to specify whether the values of the moderator variables should be added to the returned object. See \sQuote{Examples}.}
   \item{newscale}{optional vector or matrix to specify the values of the scale variables for which the predicted values should be calculated. Only relevant for location-scale models (see \code{\link{rma}}). See \sQuote{Details}.}
   \item{addz}{logical to specify whether the values of the scale variables should be added to the returned object.}
   \item{level}{numeric value between 0 and 100 to specify the confidence and prediction interval level (if unspecified, the default is to take the value from the object).}
   \item{digits}{integer to specify the number of decimal places to which the printed results should be rounded (if unspecified, the default is to take the value from the object).}
   \item{transf}{optional argument to specify a function that should be used to transform the predicted values and interval bounds (e.g., \code{transf=exp}; see also \link{transf}). If unspecified, no transformation is used.}
   \item{targs}{optional arguments needed by the function specified under \code{transf}.}
   \item{vcov}{logical to specify whether the variance-covariance matrix of the predicted values should also be returned (the default is \code{FALSE}).}
   \item{\dots}{other arguments.}
}
\details{
   For a fixed-effects model, \code{predict(object)} returns the estimated (average) outcome in the set of studies included in the meta-analysis. This is the same as the estimated intercept in the fixed-effects model (i.e., \mjseqn{\hat{\theta}}).

   For a random-effects model, \code{predict(object)} returns the estimated (average) outcome in the hypothetical population of studies from which the set of studies included in the meta-analysis are assumed to be a random selection. This is the same as the estimated intercept in the random-effects model (i.e., \mjseqn{\hat{\mu}}).

   For models including one or more moderators, \code{predict(object)} returns the estimated (average) outcomes for values of the moderator(s) equal to those of the \mjseqn{k} studies included in the meta-analysis (i.e., the \sQuote{fitted values} for the \mjseqn{k} studies).

   For models including \mjseqn{p'} moderator variables, new moderator values (for \mjeqn{k_{new}}{k_new} hypothetical new studies) can be specified by setting \code{newmods} equal to a \mjeqn{k_{new} \times p'}{k_new x p'} matrix with the corresponding new moderator values. If the model object included an intercept, then it should not be explicitly specified under \code{newmods}, as it will be added by default (unless one sets \code{intercept=FALSE}). Also, any factors in the original model get turned into the appropriate contrast variables within the \code{rma} function, so that \code{newmods} should actually include the values for the contrast variables. Examples are shown below.

   For random/mixed-effects models, an approximate prediction interval is also calculated (Riley et al., 2011). The interval estimates where \code{level} \% of the true effect sizes or outcomes fall in the hypothetical population of studies (and hence where the true effect/outcome of a new study from the population of studies should fall in \code{level} \% of the cases).

   For random-effects models that were fitted with the \code{\link{rma.mv}} function, the model may actually include multiple \mjseqn{\tau^2} values (i.e., when the \code{random} argument includes an \sQuote{\code{~ inner | outer}} term and \code{struct="HCS"}, \code{struct="DIAG"}, \code{struct="HAR"}, or \code{struct="UN"}). In that case, the function will provide prediction intervals for each level of the inner factor (since the prediction intervals differ depending on the \mjseqn{\tau^2} value). Alternatively, one can use the \code{tau2.levels} argument to specify for which level(s) the prediction interval should be provided. If the model includes a second \sQuote{\code{~ inner | outer}} term with multiple \mjseqn{\gamma^2} values, prediction intervals for each combination of levels of the inner factors will be provided. Alternatively, one can use the \code{tau2.levels} and \code{gamma2.levels} arguments to specify for which level combination(s) the prediction interval should be provided.

   When using the \code{newmods} argument for mixed-effects models that were fitted with the \code{\link{rma.mv}} function, if the model includes multiple \mjseqn{\tau^2} (and multiple \mjseqn{\gamma^2}) values, then one must use the \code{tau2.levels} (and \code{gamma2.levels}) argument to specify the levels of the inner factor(s) (i.e., a vector of length \mjeqn{k_{new}}{k_new}) to obtain the appropriate prediction interval(s).

   For location-scale models fitted with the \code{\link{rma}} function, one can use \code{newmods} to specify the values of the \mjseqn{p'} moderator variables included in the model and \code{newscale} to specify the values of the \mjseqn{q'} scale variables included in the model. Whenever \code{newmods} is specified, the function computes predicted effects/outcomes for the specified moderators values. To obtain the corresponding prediction intervals, one must also specify the corresponding \code{newscale} values. If only \code{newscale} is specified (and not \code{newmods}), the function computes the predicted log-transformed \mjseqn{\tau^2} values (when using a log link) for to the specified scale values. By setting \code{transf=exp}, one can then obtain the predicted \mjseqn{\tau^2} values.
}
\value{
   An object of class \code{"list.rma"}. The object is a list containing the following components:
   \item{pred}{predicted value(s).}
   \item{se}{corresponding standard error(s).}
   \item{ci.lb}{lower bound of the confidence interval(s).}
   \item{ci.ub}{upper bound of the confidence interval(s).}
   \item{pi.lb}{lower bound of the prediction interval(s) (only for random/mixed-effects models).}
   \item{pi.ub}{upper bound of the prediction interval(s) (only for random/mixed-effects models).}
   \item{tau2.level}{the level(s) of the inner factor (only for models of class \code{"rma.mv"} with multiple \mjseqn{\tau^2} values).}
   \item{gamma2.level}{the level(s) of the inner factor (only for models of class \code{"rma.mv"} with multiple \mjseqn{\gamma^2} values).}
   \item{X}{the moderator value(s) used to calculate the predicted values (only when \code{addx=TRUE}).}
   \item{Z}{the scale value(s) used to calculate the predicted values (only when \code{addz=TRUE} and only for location-scale models).}
   \item{\dots}{some additional elements/values.}

   If \code{vcov=TRUE}, then the returned object is a list with the first element equal to the one as described above and the second element equal to the variance-covariance matrix of the predicted values.

   The object is formatted and printed with \code{\link{print.list.rma}}.
}
\note{
   Confidence and prediction intervals are calculated based on the critical values from a standard normal distribution (i.e., \mjeqn{\pm 1.96}{±1.96} for \code{level=95}). When the model was fitted with \code{test="t"} or \code{test="knha"}, then a t-distribution with \mjseqn{k-p} degrees of freedom is used.

   For a random-effects model (where \mjseqn{p=1}) fitted with the \code{\link{rma.uni}} function, note that this differs slightly from Riley et al. (2001), who suggest to use a t-distribution with \mjseqn{k-2} degrees of freedom for constructing the prediction interval. Neither a normal, nor a t-distribution with \mjseqn{k-1} or \mjseqn{k-2} degrees of freedom is correct; all of these are approximations. The computations are done in the way described above, so that the prediction interval is identical to the confidence interval when \mjeqn{\hat{\tau}^2 = 0}{hat(\tau)^2 = 0}, which could be argued is the logical thing that should happen. If the prediction interval should be computed exactly as described by Riley et al. (2001), then one can use argument \code{pi.type="riley"}.

   The predicted values are based only on the fixed effects of the model. Best linear unbiased predictions (BLUPs) that combine the fitted values based on the fixed effects and the estimated contributions of the random effects can be obtained with \code{\link{blup.rma.uni}} (currently only for objects of class \code{"rma.uni"}).

   When using the \code{transf} option, the transformation is applied to the predicted values and the corresponding interval bounds. The standard errors are omitted from the printed output. Also, \code{vcov=TRUE} is ignored when using the \code{transf} option.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Hedges, L. V., & Olkin, I. (1985). \emph{Statistical methods for meta-analysis}. San Diego, CA: Academic Press.

   Riley, R. D., Higgins, J. P. T., & Deeks, J. J. (2011). Interpretation of random effects meta-analyses. \emph{British Medical Journal}, \bold{342}, d549. \verb{https://doi.org/10.1136/bmj.d549}

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{fitted.rma}}, \code{\link{blup.rma.uni}}
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### fit random-effects model
res <- rma(yi, vi, data=dat)

### average risk ratio with 95\% CI
predict(res, transf=exp)

### fit mixed-effects model with absolute latitude as a moderator
res <- rma(yi, vi, mods = ~ ablat, data=dat)

### predicted average risk ratios for given absolute latitude values
predict(res, transf=exp, addx=TRUE)

### predicted average risk ratios for 10-60 degrees absolute latitude
predict(res, newmods=c(10, 20, 30, 40, 50, 60), transf=exp, addx=TRUE)

### fit mixed-effects model with absolute latitude and publication year as moderators
res <- rma(yi, vi, mods = ~ ablat + year, data=dat)

### predicted average risk ratios for 10 and 60 degrees latitude in 1950 and 1980
predict(res, newmods=cbind(c(10,60,10,60),c(1950,1950,1980,1980)), transf=exp, addx=TRUE)

### fit mixed-effects model with two moderators (one of which is a factor)
res <- rma(yi, vi, mods = ~ ablat + factor(alloc), data=dat)

### examine how the factor was actually coded for the studies in the dataset
predict(res, addx=TRUE)

### predictd average risk ratios at 30 degrees for the three factor levels
### note: the contrast (dummy) variables need to specified explicitly here
predict(res, newmods=c(30, 0, 0), addx=TRUE)   # for alternate  allocation
predict(res, newmods=c(30, 1, 0), addx=TRUE)   # for random     allocation
predict(res, newmods=c(30, 0, 1), addx=TRUE)   # for systematic allocation

### can also use named vector with arbitrary order and abbreviated variable names
predict(res, newmods=c(sys=0, ran=0, abl=30))
predict(res, newmods=c(sys=0, ran=1, abl=30))
predict(res, newmods=c(sys=1, ran=0, abl=30))
}
\keyword{models}
