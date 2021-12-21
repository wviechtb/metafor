\name{regplot}
\alias{regplot}
\alias{regplot.rma}
\alias{points.regplot}
\title{Scatter Plots / Bubble Plots}
\description{
   Function to create scatter plots / bubble plots based on meta-regression models. \loadmathjax
}
\usage{
regplot(x, \dots)

\method{regplot}{rma}(x, mod, pred=TRUE, ci=TRUE, pi=FALSE, shade=TRUE,
        xlim, ylim, predlim, olim, xlab, ylab, at, digits=2L,
        transf, atransf, targs, level=x$level,
        pch=21, psize, plim=c(0.5,3), col="black", bg="darkgray",
        grid=FALSE, refline, label=FALSE, offset=c(1,1), labsize=1,
        lcol, lwd, lty, legend=FALSE, xvals, \dots)

\method{points}{regplot}(x, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma.uni"}, \code{"rma.mv"}, or \code{"rma.glmm"} including one or multiple moderators (or an object of class \code{"regplot"} for \code{points}).}
   \item{mod}{either a scalar to specify the position of the moderator variable in the model or a character string to specify the name of the moderator variable.}
   \item{pred}{logical to indicate whether the (marginal) regression line based on the moderator should be added to the plot (the default is \code{TRUE}). Can also be an object from \code{\link[=predict.rma]{predict}}. See \sQuote{Details}.}
   \item{ci}{logical to indicate whether the corresponding confidence interval bounds should be added to the plot (the default is \code{TRUE}).}
   \item{pi}{logical to indicate whether the corresponding prediction interval bounds should be added to the plot (the default is \code{FALSE}).}
   \item{shade}{logical to indicate whether the confidence/prediction interval regions should be shaded (the default is \code{TRUE}). Can also be a two-element character vector to specify the colors for shading the confidence and prediction interval regions (if shading only the former, a single color can also be specified).}
   \item{xlim}{x-axis limits. If unspecified, the function tries to set the x-axis limits to some sensible values.}
   \item{ylim}{y-axis limits. If unspecified, the function tries to set the y-axis limits to some sensible values.}
   \item{predlim}{optional argument to specify the limits of the (marginal) regression line. If unspecified, the limits are based on the range of the moderator variable.}
   \item{olim}{optional argument to specify observation/outcome limits. If unspecified, no limits are used.}
   \item{xlab}{title for the x-axis. If unspecified, the function tries to set an appropriate axis title.}
   \item{ylab}{title for the y-axis. If unspecified, the function tries to set an appropriate axis title.}
   \item{at}{position of the y-axis tick marks and corresponding labels. If unspecified, the function tries to set the tick mark positions/labels to some sensible values.}
   \item{digits}{integer to specify the number of decimal places to which the tick mark labels of the y-axis should be rounded. When specifying an integer (e.g., \code{2L}), trailing zeros after the decimal mark are dropped for the y-axis labels. When specifying a numeric value (e.g., \code{2}), trailing zeros are retained.}
   \item{transf}{optional argument to specify a function that should be used to transform the observed outcomes, predicted values, and confidence/prediction interval bounds (e.g., \code{transf=exp}; see also \link{transf}). If unspecified, no transformation is used.}
   \item{atransf}{optional argument to specify a function that should be used to transform the y-axis labels (e.g., \code{atransf=exp}; see also \link{transf}). If unspecified, no transformation is used.}
   \item{targs}{optional arguments needed by the function specified via \code{transf} or \code{atransf}.}
   \item{level}{numeric value between 0 and 100 to specify the confidence/prediction interval level (the default is to take the value from the object).}
   \item{pch}{plotting symbol to use for the observed outcomes. By default, a filled circle is used. Can also be a vector of values. See \code{\link{points}} for other options.}
   \item{psize}{optional numeric value to specify the point sizes for the observed outcomes. If unspecified, the point sizes are a function of the model weights. Can also be a vector of values. Can also be a character string (either \code{"seinv"} or \code{"vinv"}) to make the point sizes proportional to the inverse standard errors or inverse sampling variances.}
   \item{plim}{numeric vector of length 2 to scale the point sizes (ignored when a numeric value or vector is specified for \code{psize}). See \sQuote{Details}.}
   \item{col}{character string to specify the color to use for plotting the observed outcomes (the default is \code{"black"}). Can also be a vector.}
   \item{bg}{character string to specify the background color for open plot symbols (the default is \code{"darkgray"}). Can also be a vector.}
   \item{grid}{logical to specify whether a grid should be added to the plot. Can also be a color name.}
   \item{refline}{optional numeric value to specify the location of a horizontal reference line that should be added to the plot.}
   \item{label}{argument to control the labeling of the points (the default is \code{FALSE}). See \sQuote{Details}.}
   \item{offset}{argument to control the distance between the points and the corresponding labels. See \sQuote{Details}.}
   \item{labsize}{numeric value to control the size of the labels.}
   \item{lcol}{optional vector of (up to) four elements to specify the color of the regression line, of the confidence interval bounds, of the prediction interval bounds, and of the horizontal reference line.}
   \item{lty}{optional vector of (up to) four elements to specify the line type of the regression line, of the confidence interval bounds, of the prediction interval bounds, and of the horizontal reference line.}
   \item{lwd}{optional vector of (up to) four elements to specify the line width of the regression line, of the confidence interval bounds, of the prediction interval bounds, and of the horizontal reference line.}
   \item{legend}{logical to indicate whether a legend should be added to the plot (the default is \code{FALSE}). Can also be a keyword to indicate the position of the legend (see \code{\link{legend}}).}
   \item{xvals}{optional numeric vector to specify the values of the moderator for which predicted values should be computed. Needs to be specified when passing an object from \code{\link[=predict.rma]{predict}} to the \code{pred} argument. See \sQuote{Details}.}
   \item{\dots}{other arguments.}
}
\details{
   The function draws a scatter plot of the values of a moderator variable in a meta-regression model (on the x-axis) against the observed effect sizes or outcomes (on the y-axis). The regression line from the model (with corresponding confidence interval bounds) is added to the plot by default. These types of plots are also often referred to as \sQuote{bubble plots} as the points are typically drawn in different sizes to reflect their precision or weight in the model.

   By default (i.e., when \code{psize} is not specified), the size of the points is a function of the square root of the model weights. This way, their area is proportional to the the weights. However, the point sizes are rescaled so that the smallest point size is \code{plim[1]} and the largest point size is \code{plim[2]}. As a result, their relative sizes (i.e., areas) no longer exactly correspond to their relative weights. If exactly relative point sizes are desired, one can set \code{plim[2]} to \code{NA}, in which case the points are rescaled so that the smallest point size corresponds to \code{plim[1]} and all other points are scaled accordingly. As a result, the largest point may be very large. Alternatively, one can set \code{plim[1]} to \code{NA}, in which case the points are rescaled so that the largest point size corresponds to \code{plim[2]} and all other points are scaled accordingly. As a result, the smallest point may be very small. To avoid the latter, one can also set \code{plim[3]}, which enforces a minimal point size.

   One can also set \code{psize} to a scalar (e.g., \code{psize=1}) to avoid that the points are drawn in different sizes. One can also specify the point sizes manually by passing a vector of the appropriate length to \code{psize}. Finally, one can also set \code{psize} to either \code{"seinv"} or \code{"vinv"} to make the point sizes proportional to the inverse standard errors or inverse sampling variances.

   For a model with more than one predictor, the regression line reflects the \sQuote{marginal} relationship between the chosen moderator and the effect sizes or outcomes (i.e., all other moderators except the one being plotted are held constant at their means).

   With the \code{label} argument, one can control whether points in the plot will be labeled. If \code{label="all"} (or \code{label=TRUE}), all points in the plot will be labeled. If \code{label="ciout"} or \code{label="piout"}, points falling outside of the confidence/prediction interval will be labeled. Alternatively, one can set this argument to a logical or numeric vector to specify which points should be labeled. The labels are placed above the points when they fall above the regression line and otherwise below. With the \code{offset} argument, one can adjust the distance between the labels and the corresponding points. This can either be a single numeric value, which is used as a multiplicative factor for the point sizes (so that the distance between labels and points is larger for larger points) or a numeric vector with two values, where the first is used as an additive factor independent of the point sizes and the second again as a multiplicative factor for the point sizes. The values are given as percentages of the y-axis range. It may take some trial and error to find two values for the \code{offset} argument so that the labels are placed right next to the boundary of the points. With \code{labsize}, one can control the size of the labels.

   One can also pass an object from \code{\link[=predict.rma]{predict}} to the \code{pred} argument. This can be useful when the meta-regression model reflects a more complex relationship between the moderator variable and the effect sizes or outcomes (e.g., when using polynomials or splines) or when the model involves interactions. In this case, one also needs to specify the \code{xvals} argument. See \sQuote{Examples}.
}
\note{
   For certain types of models, it may not be possible to draw the prediction interval bounds (if this is the case, a warning will be issued).

   When specifying vectors for \code{pch}, \code{psize}, \code{col}, \code{bg}, and/or \code{label}, the variables specified are assumed to be of the same length as the data passed to the model fitting function. Any subsetting and removal of studies with missing values is automatically applied to the variables specified via these arguments.

   If the outcome measure used for creating the plot is bounded (e.g., correlations are bounded between -1 and +1, proportions are bounded between 0 and 1), one can use the \code{olim} argument to enforce those limits (the observed outcomes and confidence/prediction intervals cannot exceed those bounds then).
}
\value{
   An object of class \code{"regplot"} with components:
   \item{slab}{the study labels}
   \item{ids}{the study ids}
   \item{xi}{the x-axis coordinates of the points that were plotted.}
   \item{yi}{the y-axis coordinates of the points that were plotted.}
   \item{pch}{the plotting symbols of the points that were plotted.}
   \item{psize}{the point sizes of the points that were plotted.}
   \item{col}{the colors of the points that were plotted.}
   \item{bg}{the background colors of the points that were plotted.}
   \item{label}{logical vector indicating whether a point was labeled or not.}

   Note that the object is returned invisibly. Using \code{points.regplot}, one can redraw the points (and labels) in case one wants to superimpose the points on top of any elements that were added manually to the plot (see \sQuote{Examples}).
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Thompson, S. G., & Higgins, J. P. T. (2002). How should meta-regression analyses be undertaken and interpreted? \emph{Statistics in Medicine}, \bold{21}(11), 1559--1573. \verb{https://doi.org/10.1002/sim.1187}

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.glmm}}, and \code{\link{rma.mv}} for functions to fit models for which scatter plots / bubble plots can be drawn.
}
\examples{
### copy BCG vaccine data into 'dat'
dat <- dat.bcg

### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat)

############################################################################

### fit mixed-effects model with absolute latitude as a moderator
res <- rma(yi, vi, mods = ~ ablat, data=dat)
res

### draw plot
regplot(res, mod="ablat", xlab="Absolute Latitude")

### adjust x-axis limits and back-transform to risk ratios
regplot(res, mod="ablat", xlab="Absolute Latitude", xlim=c(0,60), transf=exp)

### also extend the prediction limits for the regression line
regplot(res, mod="ablat", xlab="Absolute Latitude", xlim=c(0,60), predlim=c(0,60), transf=exp)

### add the prediction interval to the plot, add a reference line at 1, and add a legend
regplot(res, mod="ablat", pi=TRUE, xlab="Absolute Latitude",
        xlim=c(0,60), predlim=c(0,60), transf=exp, refline=1, legend=TRUE)

### label points outside of the prediction interval
regplot(res, mod="ablat", pi=TRUE, xlab="Absolute Latitude",
        xlim=c(0,60), predlim=c(0,60), transf=exp, refline=1, legend=TRUE,
        label="piout", labsize=0.8)

############################################################################

### fit mixed-effects model with absolute latitude and publication year as moderators
res <- rma(yi, vi, mods = ~ ablat + year, data=dat)
res

### plot the marginal relationships
regplot(res, mod="ablat", xlab="Absolute Latitude")
regplot(res, mod="year",  xlab="Publication Year")

############################################################################

### fit a quadratic polynomial meta-regression model
res <- rma(yi, vi, mods = ~ ablat + I(ablat^2), data=dat)
res

### compute predicted values using predict()
xs <- seq(0,60,length=601)
tmp <- predict(res, newmods=cbind(xs, xs^2))

### can now pass these results to the 'pred' argument (and have to specify xvals accordingly)
regplot(res, mod="ablat", pred=tmp, xlab="Absolute Latitude", xlim=c(0,60), xvals=xs)

### back-transform to risk ratios and add reference line
regplot(res, mod="ablat", pred=tmp, xlab="Absolute Latitude", xlim=c(0,60), xvals=xs,
        transf=exp, refline=1)

############################################################################

### fit a model with an interaction between a quantitative and a categorical predictor
### (note: just for illustration purposes; this model is too complex for this dataset)
res <- rma(yi, vi, mods = ~ ablat * alloc, data=dat)
res

### draw bubble plot but do not add regression line or CI
tmp <- regplot(res, mod="ablat", xlab="Absolute Latitude", xlim=c(0,60), pred=FALSE, ci=FALSE)

### add regression lines for the three alloc levels
xs <- seq(0, 60, length=100)
preds <- predict(res, newmods=cbind(xs, 0, 0, 0, 0))
lines(xs, preds$pred, lwd=3)
preds <- predict(res, newmods=cbind(xs, 1, 0, xs, 0))
lines(xs, preds$pred, lwd=3)
preds <- predict(res, newmods=cbind(xs, 0, 1, 0, xs))
lines(xs, preds$pred, lwd=3)

### add points back to the plot (so they are on top of the lines)
points(tmp)
}
\keyword{hplot}
