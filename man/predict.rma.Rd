\name{predict.rma}
\alias{predict.rma}
\title{Predicted Values for 'rma' Objects}
\description{
   The function computes predicted values, corresponding standard errors, confidence intervals, and (approximate) credibility/prediction intervals for objects of class \code{"rma"}.
}
\usage{
\method{predict}{rma}(object, newmods, intercept, tau2.levels, gamma2.levels,
        addx=FALSE, level, digits, transf, targs, \dots)
}
\arguments{
   \item{object}{an object of class \code{"rma"}.}
   \item{newmods}{optional vector or matrix specifying the values of the moderator values for which the predicted values should be calculated. See \sQuote{Details}.}
   \item{intercept}{logical specifying whether the intercept term should be included when calculating the predicted values for \code{newmods}. If unspecified, the intercept is automatically added when the original model also included an intercept.}
   \item{tau2.levels}{vector specifying the levels of the inner factor when computing credibility/prediction intervals. Only relevant for models of class \code{"rma.mv"} (see \code{\link{rma.mv}}) and when the model includes more than a single \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}} value. See \sQuote{Details}.}
   \item{gamma2.levels}{vector specifying the levels of the inner factor when computing credibility/prediction intervals. Only relevant for models of class \code{"rma.mv"} (see \code{\link{rma.mv}}) and when the model includes more than a single \if{latex}{\eqn{\gamma^2}}\if{html}{\eqn{\gamma²}}\if{text}{\eqn{\gamma^2}} value. See \sQuote{Details}.}
   \item{addx}{logical specifying whether the values of the moderator variables should be added to the returned object. See \sQuote{Examples}.}
   \item{level}{numerical value between 0 and 100 specifying the confidence and credibility/prediction interval level (if unspecified, the default is to take the value from the object).}
   \item{digits}{integer specifying the number of decimal places to which the printed results should be rounded (if unspecified, the default is to take the value from the object).}
   \item{transf}{optional argument specifying the name of a function that should be used to transform the predicted values and interval bounds (e.g., \code{transf=exp}; see also \link{transf}). If unspecified, no transformation is used.}
   \item{targs}{optional arguments needed by the function specified under \code{transf}.}
   \item{\dots}{other arguments.}
}
\details{
   For a fixed-effects model, \code{predict(object)} returns the estimated (average) outcome in the set of studies included in the meta-analysis. This is the same as the estimated intercept in the fixed-effects model.

   For a random-effects model, \code{predict(object)} returns the estimated (average) outcome in the hypothetical population of studies from which the set of studies included in the meta-analysis are assumed to be a random selection. This is the same as the estimated intercept in the random-effects model.

   For models including one or more moderators, \code{predict(object)} returns the estimated (average) outcomes for values of the moderator(s) equal to those of the \eqn{k} studies included in the meta-analysis (i.e., the fitted values for the \eqn{k} studies).

   For models including \eqn{q} moderator variables, new moderator values for \eqn{k_{new}}{k_new} new studies can be specified by setting \code{newmods} equal to an \eqn{k_{new} \times q}{k_new x q} matrix with the corresponding new moderator values. If the original model includes an intercept, then it should not be explicitly specified under \code{newmods}, as it will be added by default (unless one sets \code{intercept=FALSE}). Also, any factors in the original model get turned into the appropriate contrast variables within the \code{rma} function, so that \code{newmods} should actually include the values for the contrast variables. Examples are shown below.

   For random/mixed-effects models, an approximate credibility/credible or prediction interval is also calculated (Raudenbush, 2009; Riley et al., 2011). The interval estimates where \code{level} \% of the true outcomes fall in the hypothetical population of studies. Note that this interval is calculated under the assumption that the value of \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}} is known (and not estimated). A method for calculating a credibility/prediction interval that accounts for the uncertainty in the estimate of \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}} may be implemented in the future.

   For random-effects models of class \code{"rma.mv"} (see \code{\link{rma.mv}}) with multiple \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}} values (i.e., when the \code{random} argument includes a \code{~ inner | outer} term and \code{struct="HCS"}, \code{struct="HAR"}, or \code{struct="UN"}), the function will provide credibility/prediction intervals for each level of the inner factor (since the credibility/prediction intervals differ depending on the \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}} value). Alternatively, one can use the \code{tau2.levels} argument to specify for which level(s) the credibility/prediction interval should be provided.

   If the model includes a second \code{~ inner | outer} term with multiple \if{latex}{\eqn{\gamma^2}}\if{html}{\eqn{\gamma²}}\if{text}{\eqn{\gamma^2}} values, credibility/prediction intervals for each combination of levels of the inner factors will be provided. Alternatively, one can use the \code{tau2.levels} and \code{gamma2.levels} arguments to specify for which level combination(s) the credibility/prediction interval should be provided.

   When using the \code{newmods} argument for mixed-effects models of class \code{"rma.mv"} with multiple \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}} (and multiple \if{latex}{\eqn{\gamma^2}}\if{html}{\eqn{\gamma²}}\if{text}{\eqn{\gamma^2}}) values, one must use the \code{tau2.levels} (and \code{gamma2.levels}) argument to specify the levels of the inner factor(s) (i.e., a vector of length \eqn{k_{new}}{k_new}) to obtain the appropriate credibility/prediction interval(s).
}
\value{
   An object of class \code{"list.rma"}. The object is a list containing the following components:
   \item{pred}{predicted value(s).}
   \item{se}{corresponding standard error(s).}
   \item{ci.lb}{lower bound of the confidence interval(s).}
   \item{ci.ub}{upper bound of the confidence interval(s).}
   \item{cr.lb}{lower bound of the credibility/prediction interval(s) (only random/mixed-effects models).}
   \item{cr.ub}{upper bound of the credibility/prediction interval(s) (only random/mixed-effects models).}
   \item{tau2.level}{the level(s) of the inner factor (only for models of class \code{"rma.mv"} with multiple \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}} values).}
   \item{gamma2.level}{the level(s) of the inner factor (only for models of class \code{"rma.mv"} with multiple \if{latex}{\eqn{\gamma^2}}\if{html}{\eqn{\gamma²}}\if{text}{\eqn{\gamma^2}} values).}
   \item{X}{the moderator value(s) used to calculate the predicted values (only when \code{addx=TRUE}).}
   \item{\dots}{some additional elements/values.}

   The \code{"list.rma"} object is formatted and printed with \code{\link{print.list.rma}}.
}
\note{
   The predicted values are based only on the fixed effects of the model. Best linear unbiased predictions (BLUPs) that combine the fitted values based on the fixed effects and the estimated contributions of the random effects can be obtained with \code{\link{blup.rma.uni}} (currently only for objects of class \code{"rma.uni"}).

   When using the \code{transf} option, the transformation is applied to the predicted values and the corresponding interval bounds. The standard errors are then set equal to \code{NA} and are omitted from the printed output.

   Confidence and credibility/prediction intervals are calculated based on the critical values from a standard normal distribution (i.e., \eqn{\pm 1.96}{+- 1.96} for \code{level=95}). When the model was fitted with the Knapp and Hartung (2003) method (i.e., \code{test="knha"} in the \code{\link{rma.uni}} function) or with \code{test="t"} in the \code{\link{rma.glmm}} and \code{\link{rma.mv}} functions, then the t-distribution with \eqn{k-p} degrees of freedom is used.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \cr
   package website: \url{http://www.metafor-project.org/} \cr
   author homepage: \url{http://www.wvbauer.com/}
}
\references{
   Hedges, L. V., & Olkin, I. (1985). \emph{Statistical methods for meta-analysis}. San Diego, CA: Academic Press.

   Raudenbush, S. W. (2009). Analyzing effect sizes: Random effects models. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta-analysis} (2nd ed., pp. 295--315). New York: Russell Sage Foundation.

   Riley, R. D., Higgins, J. P. T., & Deeks, J. J. (2011). Interpretation of random effects meta-analyses. \emph{British Medical Journal}, \bold{342}, d549.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
   \code{\link{fitted.rma}}, \code{\link{blup.rma.uni}}
}
\examples{
### meta-analysis of the log risk ratios using a random-effects model
res <- rma(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### average risk ratio with 95\% CI
predict(res, transf=exp)

### mixed-effects model with absolute latitude as a moderator
res <- rma(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, mods = ~ ablat,
           data=dat.bcg)

### predicted average risk ratios for given absolute latitude values
predict(res, transf=exp, addx=TRUE)

### predicted average risk ratios for 10-60 degrees absolute latitude
predict(res, newmods=c(10, 20, 30, 40, 50, 60), transf=exp)

### mixed-effects model with two moderators (absolute latitude and publication year)
res <- rma(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg,
           mods = ~ ablat + year, data=dat.bcg)

### predicted average risk ratios for 10 and 60 degrees latitude in 1950 and 1980
predict(res, newmods=cbind(c(10,60,10,60),c(1950,1950,1980,1980)), transf=exp, addx=TRUE)

### mixed-effects model with two moderators (one of which is a factor)
res <- rma(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg,
           mods = ~ ablat + factor(alloc), data=dat.bcg)

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
