\name{emmprep}
\alias{emmprep}
\title{Create a Reference Grid for the 'emmeans' Function}
\description{
   Function to create a reference grid for use with the \code{\link[emmeans]{emmeans}} function from the package of the same name. \loadmathjax
}
\usage{
emmprep(x, verbose=FALSE, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma"}.}
   \item{verbose}{logical to specify whether information on some (extracted) settings should be printed when creating the reference grid (the default is \code{FALSE}).}
   \item{\dots}{other arguments that will be passed on to the \code{\link[emmeans]{qdrg}} function.}
}
\details{
   The \href{https://cran.r-project.org/package=emmeans}{emmeans} package is a popular package that facilitates the computation of 'estimated marginal means'. The function is a wrapper around the \code{\link[emmeans]{qdrg}} function from the \code{emmeans} package to make \code{"rma"} objects compatible with the latter. Unless one needs to pass additional arguments to the \code{\link[emmeans]{qdrg}} function, one simply applies this function to the \code{"rma"} object and then the \code{\link[emmeans]{emmeans}} function (or one of the other functions that can be applied to \code{"emmGrid"} objects) to the resulting object to obtain the desired estimated marginal means.
}
\value{
   An \code{"emmGrid"} object as created by the \code{\link[emmeans]{qdrg}} function from the \code{emmeans} package.

   The resulting object will typically be used in combination with the \code{\link[emmeans]{emmeans}} function.
}
\note{
   When creating the reference grid, the function extracts the degrees of freedom for tests/confidence intervals from the model object (if the model was fitted with \code{test="t"}, \code{test="knha"}, \code{test="hksj"}, or \code{test="adhoc"}; otherwise the degrees of freedom are infinity). In some cases, there is not just a single value for the degrees of freedom, but an entire vector (e.g., for models fitted with \code{\link{rma.mv}}). In this case, the smallest value will be used (as a conservative option). One can set a different/custom value for the degrees of freedom with \code{emmprep(..., df=value)}.

   When the model object contains information about the outcome measure used in the analysis (which should be the case if the observed outcomes were computed with \code{\link{escalc}} or if the \code{measure} argument was set when fitting the model), then information about the appropriate back-transformation (if available) is stored as part of the returned object. If so, the back-transformation is automatically applied when calling \code{\link[emmeans]{emmeans}} with \code{type="response"}.

   The function also tries to extract the estimated value of \mjseqn{\tau^2} (or more precisely, its square root) from the model object (when the model is a random/mixed-effects model). This value is only needed when computing prediction intervals (i.e., when \code{interval="predict"} in \code{\link[emmeans]{predict.emmGrid}}) or when applying the bias adjustment in the back-transformation (i.e., when \code{bias.adjust=TRUE} in \code{\link[emmeans]{summary.emmGrid}}). For some models (e.g., those fitted with \code{\link{rma.mv}}), it is not possible to automatically extract the estimate. In this case, one can manually set the value with \code{emmprep(..., sigma=value)} (note: the argument is called \code{sigma}, following the conventions of \code{\link[emmeans]{summary.emmGrid}} and one must supply the square root of the \mjseqn{\tau^2} estimate).

   By default, the reference grid is created based on the data used for fitting the original model (which is typically the sensible thing to do). One can specify a different dataset with \code{emmprep(..., data=obj)}, where \code{obj} must be a data frame that contains the same variables as used in the original model fitted.

   If the original model fitted involved redundant predictors that were dropped from the model (due to \sQuote{rank deficiencies}), then the function cannot be used. In this case, one should remove any redundancies in the original model fitted before using this function.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### fit meta-regression model with absolute latitude as predictor
res <- rma(yi, vi, mods = ~ ablat, data=dat)
res

### create reference grid
sav <- emmprep(res, verbose=TRUE)

### estimated marginal mean (back-transformed to the risk ratio scale)
if (require(emmeans))
   emmeans(sav, specs="1", type="response")

### same as the predicted effect at the mean absolute latitude
predict(res, newmods = mean(model.matrix(res, asdf=TRUE)$ablat), transf=exp, digits=3)

### fit meta-regression model with allocation factor
res <- rma(yi, vi, mods = ~ alloc, data=dat)
res

### create reference grid
sav <- emmprep(res)

### estimated marginal mean using proportional cell weighting
if (require(emmeans))
   emmeans(sav, specs="1", type="response", weights="proportional")

### estimated marginal mean using equal cell weighting (this is actually the default)
if (require(emmeans))
   emmeans(sav, specs="1", type="response", weights="equal")

### same as the predicted effect using cell proportions as observed in the data
### or using equal proportions for the three groups
predict(res, newmods = colMeans(model.matrix(res))[-1], transf=exp, digits=3)
predict(res, newmods = c(1/3,1/3), transf=exp, digits=3)

### fit meta-regression model with absolute latitude and allocation as predictors
res <- rma(yi, vi, mods = ~ ablat + alloc, data=dat)
res

### create reference grid
sav <- emmprep(res)

### estimated marginal mean using equal cell weighting
if (require(emmeans))
   emmeans(sav, specs="1", type="response")

### same as the predicted effect at the mean absolute latitude and using equal proportions
### for the allocation factor
predict(res, newmods = c(mean(model.matrix(res, asdf=TRUE)$ablat),1/3,1/3), transf=exp, digits=3)

### create reference grid with ablat set equal to 10, 30, and 50 degrees
sav <- emmprep(res, at=list(ablat=c(10,30,50)))

### estimated marginal means at the three ablat values
if (require(emmeans))
   emmeans(sav, specs="1", by="ablat", type="response")

### same as the predicted effect at the chosen absolute latitude values and using equal
### proportions for the allocation factor
predict(res, newmods = cbind(c(10,30,50),1/3,1/3), transf=exp, digits=3)

############################################################################

### copy data into 'dat' and examine data
dat <- dat.mcdaniel1994
head(dat)

### calculate r-to-z transformed correlations and corresponding sampling variances
dat <- escalc(measure="ZCOR", ri=ri, ni=ni, data=dat)

### mixed-effects model with interview type as factor
res <- rma(yi, vi, mods = ~ factor(type), data=dat, test="knha")
res

### create reference grid
sav <- emmprep(res, verbose=TRUE)

### estimated marginal mean (back-transformed to the correlation scale)
if (require(emmeans))
   emmeans(sav, specs="1", type="response")

### same as the predicted correlation using equal cell proportions
predict(res, newmods = c(1/3,1/3), transf=transf.ztor, digits=3)

### estimated marginal means for the three interview types
if (require(emmeans))
   emmeans(sav, specs="type", type="response")

### same as the predicted correlations
predict(res, newmods = rbind(c(0,0), c(1,0), c(0,1)), transf=transf.ztor, digits=3)

### illustrate use of the 'df' and 'sigma' arguments
res <- rma.mv(yi, vi, mods = ~ factor(type), random = ~ 1 | study,
              data=dat, test="t", dfs="contain")
res

### create reference grid
sav <- emmprep(res, verbose=TRUE, df=154, sigma=0.1681)

### estimated marginal mean (back-transformed to the correlation scale)
if (require(emmeans))
   emmeans(sav, specs="1", type="response")
}
\keyword{manip}
