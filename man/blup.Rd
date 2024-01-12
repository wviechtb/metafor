\name{blup}
\alias{blup}
\alias{blup.rma.uni}
\title{Best Linear Unbiased Predictions for 'rma.uni' Objects}
\description{
   Function to compute best linear unbiased predictions (BLUPs) of the study-specific true effect sizes or outcomes (by combining the fitted values based on the fixed effects and the estimated contributions of the random effects) for objects of class \code{"rma.uni"}. Corresponding standard errors and prediction interval bounds are also provided. \loadmathjax
}
\usage{
blup(x, \dots)

\method{blup}{rma.uni}(x, level, digits, transf, targs, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma.uni"}.}
   \item{level}{numeric value between 0 and 100 to specify the prediction interval level (see \link[=misc-options]{here} for details). If unspecified, the default is to take the value from the object.}
   \item{digits}{optional integer to specify the number of decimal places to which the printed results should be rounded. If unspecified, the default is to take the value from the object.}
   \item{transf}{optional argument to specify a function to transform the predicted values and interval bounds (e.g., \code{transf=exp}; see also \link{transf}). If unspecified, no transformation is used.}
   \item{targs}{optional arguments needed by the function specified under \code{transf}.}
   \item{\dots}{other arguments.}
}
\value{
   An object of class \code{"list.rma"}. The object is a list containing the following components:
   \item{pred}{predicted values.}
   \item{se}{corresponding standard errors.}
   \item{pi.lb}{lower bound of the prediction intervals.}
   \item{pi.ub}{upper bound of the prediction intervals.}
   \item{\dots}{some additional elements/values.}

   The object is formatted and printed with the \code{\link[=print.list.rma]{print}} function. To format the results as a data frame, one can use the \code{\link[=as.data.frame.list.rma]{as.data.frame}} function.
}
\note{
   For best linear unbiased predictions of only the random effects, see \code{\link{ranef}}.

   For predicted/fitted values that are based only on the fixed effects of the model, see \code{\link[=fitted.rma]{fitted}} and \code{\link[=predict.rma]{predict}}.

   For conditional residuals (the deviations of the observed effect sizes or outcomes from the BLUPs), see \code{rstandard.rma.uni} with \code{type="conditional"}.

   Equal-effects models do not contain random study effects. The BLUPs for these models will therefore be equal to the fitted values, that is, those obtained with \code{\link[=fitted.rma]{fitted}} and \code{\link[=predict.rma]{predict}}.

   When using the \code{transf} argument, the transformation is applied to the predicted values and the corresponding interval bounds. The standard errors are then set equal to \code{NA} and are omitted from the printed output.

   By default, a standard normal distribution is used to construct the prediction intervals. When the model was fitted with \code{test="t"}, \code{test="knha"}, \code{test="hksj"}, or \code{test="adhoc"}, then a t-distribution with \mjseqn{k-p} degrees of freedom is used.

   To be precise, it should be noted that the function actually computes empirical BLUPs (eBLUPs), since the predicted values are a function of the estimated value of \mjseqn{\tau^2}.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Kackar, R. N., & Harville, D. A. (1981). Unbiasedness of two-stage estimation and prediction procedures for mixed linear models. Communications in Statistics, Theory and Methods, \bold{10}(13), 1249--1261. \verb{https://doi.org/10.1080/03610928108828108}

   Raudenbush, S. W., & Bryk, A. S. (1985). Empirical Bayes meta-analysis. \emph{Journal of Educational Statistics}, \bold{10}(2), 75--98. \verb{https://doi.org/10.3102/10769986010002075}

   Robinson, G. K. (1991). That BLUP is a good thing: The estimation of random effects. \emph{Statistical Science}, \bold{6}(1), 15--32. \verb{https://doi.org/10.1214/ss/1177011926}

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{rma.uni}} for the function to fit models for which BLUPs can be extracted.

   \code{\link[=predict.rma]{predict}} and \code{\link[=fitted.rma]{fitted}} for functions to compute the predicted/fitted values based only on the fixed effects and \code{\link{ranef}} for a function to compute the BLUPs based only on the random effects.
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### meta-analysis of the log risk ratios using a random-effects model
res <- rma(yi, vi, data=dat)

### BLUPs of the true risk ratios for each study
blup(res, transf=exp)

### illustrate shrinkage of BLUPs towards the (estimated) population average
res <- rma(yi, vi, data=dat)
blups <- blup(res)$pred
plot(NA, NA, xlim=c(.8,2.4), ylim=c(-2,0.5), pch=19,
     xaxt="n", bty="n", xlab="", ylab="Log Risk Ratio")
segments(rep(1,13), dat$yi, rep(2,13), blups, col="darkgray")
points(rep(1,13), dat$yi, pch=19)
points(rep(2,13), blups, pch=19)
axis(side=1, at=c(1,2), labels=c("Observed\nValues", "BLUPs"), lwd=0)
segments(0, res$beta, 2.15, res$beta, lty="dotted")
text(2.3, res$beta, substitute(hat(mu)==muhat, list(muhat=round(res$beta[[1]], 2))), cex=1)
}
\keyword{models}
