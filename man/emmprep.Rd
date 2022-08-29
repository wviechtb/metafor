\name{emmprep}
\alias{emmprep}
\title{Create a Reference Grid for the 'emmeans' Function}
\description{
   Function to create a reference grid for use with the \code{\link[emmeans]{emmeans}} function from the package of the same name.
}
\usage{
emmprep(x, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma"}.}
   \item{\dots}{other arguments that will be passed on to the \code{\link[emmeans]{qdrg}} function.}
}
\details{
   The function is a wrapper around the \code{\link[emmeans]{qdrg}} function from the \code{emmeans} package to make \code{"rma"} objects compatible with the latter. Unless one needs to pass additional arguments to the \code{\link[emmeans]{qdrg}} function, one simply applies this function to the \code{"rma"} object and then the \code{\link[emmeans]{emmeans}} function to the resulting object to obtain the desired estimated marginal means.
}
\value{
   An \code{"emmGrid"} object as created by the \code{\link[emmeans]{qdrg}} function from the \code{emmeans} package.

   The resulting object will typically be used in combination with the \code{\link[emmeans]{emmeans}} function.
}
\note{
   When calling \code{emmeans} on the resulting object with \code{type="response"}, the appropriate back-transformation should automatically be applied when the model object contains information about the outcome measure used.

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
sav <- emmprep(res)

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

### estimated marginal mean using proportional and equal cell weighting
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
}
\keyword{manip}
