\name{pairwise}
\alias{pairwise}
\title{Construct a Pairwise Contrast Matrix for 'rma' Objects}
\description{
   Functions to construct a matrix of pairwise contrasts for objects of class \code{"rma"}. \loadmathjax
}
\usage{
pairwise(x, btt, btt2, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma"}.}
   \item{btt}{vector of indices to specify for which coefficients pairwise contrasts should be constructed. Can also be a string to \code{\link{grep}} for. See \sQuote{Details}.}
   \item{btt2}{optional argument to specify a second set of coefficients that should also be included in the contrast matrix.}
   \item{\dots}{other arguments.}
}
\value{
   When a meta-regression model includes a categorical moderator variable (i.e., a factor), there is often interest in testing whether the coefficients representing the various levels of the factor differ significantly from each other. The present function constructs the pairwise contrast matrix between all factor levels for a particular factor, which can be used together with the \code{\link[=anova.rma]{anova}} function to carry out such tests and the \code{\link[=predict.rma]{predict}} function to obtain corresponding confidence intervals.

   The \code{x} argument is used to specify a meta-regression model and the \code{btt} argument the indices of the coefficients for which pairwise contrasts should be constructed. For example, with \code{btt=2:4}, contrasts are formed based on the second, third, and fourth coefficient of the model. Instead of specifying the coefficient numbers, one can specify a string for \code{btt}. In that case, \code{\link{grep}} will be used to search for all coefficient names that match the string.

   At times, it may be useful to include a second set of coefficients in the contrast matrix (not as pairwise contrasts, but as \sQuote{main effects}). This can be done via the \code{btt2} argument.

   When using the present function in a call to the \code{\link[=anova.rma]{anova}} or \code{\link[=predict.rma]{predict}} functions, argument \code{x} does not need to specified, as the function will then automatically construct the contrast matrix based on the model object passed to the \code{anova} or \code{predict} function. See below for examples.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.glmm}}, and \code{\link{rma.mv}} for functions to fit meta-regression models for which pairwise contrasts may be useful.

   \code{\link[=anova.rma]{anova}} for a function to carry out tests of the pairwise contrasts and \code{\link[=predict.rma]{predict}} to obtain corresponding confidence/prediction intervals.
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### mixed-effects meta-regression model with the allocation method as a moderator;
### by removing the intercept term, we obtain the estimated average effect for each
### factor level from the model
res <- rma(yi, vi, mods = ~ 0 + alloc, data=dat)
res

### construct the contrast matrix for the 'alloc' factor
pairwise(res, btt=1:3)
pairwise(res, btt="alloc")

### test all pairwise contrasts
anova(res, X=pairwise(btt=1:3))
anova(res, X=pairwise(btt="alloc"))

### obtain the corresponding confidence intervals
predict(res, newmods=pairwise(btt="alloc"))

### test all pairwise contrasts adjusting for multiple testing
anova(res, X=pairwise(btt="alloc"), adjust="bonf")

### fit the same model, but including the intercept term; then 'alternate' is the
### reference level and the coefficients for 'random' and 'systematic' already
### represent pairwise contrasts with this reference level
res <- rma(yi, vi, mods = ~ alloc, data=dat)
res

### in this case, we want to include these coefficients directly in the contrast
### matrix (btt2=2:3) but also include the pairwise contrast between them (btt=2:3)
pairwise(res, btt=2:3, btt2=2:3)
pairwise(res, btt="alloc", btt2="alloc")

### test all pairwise contrasts
anova(res, X=pairwise(btt=2:3, btt2=2:3))
anova(res, X=pairwise(btt="alloc", btt2="alloc"))

### obtain the corresponding confidence intervals
predict(res, newmods=pairwise(btt="alloc", btt2="alloc"))

### meta-regression model with 'ablat' and 'alloc' as moderators
res <- rma(yi, vi, mods = ~ ablat + alloc, data=dat)
res

### test all pairwise contrasts between the 'alloc' levels (while controlling for 'ablat')
anova(res, X=pairwise(btt="alloc", btt2="alloc"))
anova(res, X=pairwise(btt="alloc", btt2="alloc"))

### obtain the corresponding confidence intervals
predict(res, newmods=pairwise(btt="alloc", btt2="alloc"))

### an example of a meta-regression model with more factors levels
dat <- dat.bangertdrowns2004
res <- rma(yi, vi, mods = ~ 0 + factor(grade), data=dat)
res

### test all pairwise contrasts between the 'grade' levels
anova(res, X=pairwise(btt="grade"))

### obtain the corresponding confidence intervals
predict(res, newmods=pairwise(btt="grade"))

### test all pairwise contrasts adjusting for multiple testing
anova(res, X=pairwise(btt="grade"), adjust="bonf")
}
\keyword{models}
