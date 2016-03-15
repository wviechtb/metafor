\name{forest.rma}
\alias{forest.rma}
\title{Forest Plots (Method for 'rma' Objects)}
\description{Function to create forest plots for objects of class \code{"rma"}.}
\usage{
\method{forest}{rma}(x, annotate=TRUE, addfit=TRUE, addcred=FALSE, showweights=FALSE,
       xlim, alim, clim, ylim, at, steps=5,
       level=x$level, refline=0, digits=2, width,
       xlab, slab, mlab, ilab, ilab.xpos, ilab.pos,
       order, transf, atransf, targs, rows,
       efac=1, pch=15, psize, col, border, lty,
       cex, cex.lab, cex.axis, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma"}.}
   \item{annotate}{logical specifying whether annotations should be added to the plot (the default is \code{TRUE}).}
   \item{addfit}{logical specifying whether the summary estimate (for models without moderators) or fitted values (for models with moderators) should be added to the plot. See \sQuote{Details}.}
   \item{addcred}{logical specifying whether the bounds of the credibility/prediction interval should be added to the plot. See \sQuote{Details}.}
   \item{showweights}{logical specifying whether the annotations should also include the weights given to the observed effects or outcomes during the model fitting (the default is \code{FALSE}).}
   \item{xlim}{horizontal limits of the plot region. If unspecified, the function tries to set the horizontal plot limits to some sensible values.}
   \item{alim}{the actual x-axis limits. If unspecified, the function tries to set the x-axis limits to some sensible values.}
   \item{clim}{limits for the confidence/credibility/prediction intervals. If unspecified, no limits are used.}
   \item{ylim}{the y limits of the plot. If unspecified, the function tries to set the y-axis limits to some sensible values.}
   \item{at}{position of the x-axis tick marks and corresponding labels. If unspecified, the function tries to set the tick mark positions/labels to some sensible values.}
   \item{steps}{the number of tick marks for the x-axis (the default is 5). Ignored when the user specifies the positions via the \code{at} argument.}
   \item{level}{numerical value between 0 and 100 specifying the confidence interval level (the default is to take the value from the object).}
   \item{refline}{value at which a vertical \sQuote{reference} line should be drawn (the default is 0). The line can be suppressed by setting this argument to \code{NA}.}
   \item{digits}{integer specifying the number of decimal places to which the tick mark labels of the x-axis and the annotations should be rounded (the default is 2). Can also be a vector of two integers, the first specifying the number of decimal places for the annotations, the second for the x-axis labels.}
   \item{width}{optional integer to manually adjust the width of the columns for the annotations.}
   \item{xlab}{title for the x-axis. If unspecified, the function tries to set an appropriate axis title.}
   \item{slab}{optional vector with labels for the \eqn{k} studies. If unspecified, the labels are either taken from the object (if study labels were specified) or simple labels are created within the function. To suppress labels, set this argument to \code{NA}.}
   \item{mlab}{optional character string giving a label to the summary estimate from a fixed- or random-effects model. If unspecified, the label is created within the function.}
   \item{ilab}{optional vector or matrix of character strings providing additional information about the studies.}
   \item{ilab.xpos}{vector of numerical value(s) specifying the x-axis position(s) of the character vector(s) given via \code{ilab} (must be specified if \code{ilab} is specified).}
   \item{ilab.pos}{integer(s) (either 1, 2, 3, or 4) specifying the alignment of the character vector(s) given via \code{ilab} (2 means right, 4 mean left aligned). If unspecified, the default is to center the labels.}
   \item{order}{optional character string specifying how the studies should be ordered. See \sQuote{Details}.}
   \item{transf}{optional argument specifying the name of a function that should be used to transform the observed effect sizes, summary estimates, fitted values, and confidence interval bounds (e.g., \code{transf=exp}). If unspecified, no transformation is used.}
   \item{atransf}{optional argument specifying the name of a function that should be used to transform the x-axis labels and annotations (e.g., \code{transf=exp}). If unspecified, no transformation is used.}
   \item{targs}{optional arguments needed by the function specified via \code{transf} or \code{atransf}.}
   \item{rows}{optional vector specifying the rows (or more generally, the horizontal positions) for plotting the outcomes. If unspecified, the function sets this value automatically. Can also be a single value specifying the row (horizontal position) of the first outcome (the remaining outcomes are then plotted below this starting row).}
   \item{efac}{vertical expansion factor for confidence interval limits, arrows, and the symbol used to denote summary estimates. The default value of 1 should usually work okay. Can also be a vector of two numbers, the first for CI limits and arrows, the second for summary estimates. Can also be a vector of three numbers, the first for CI limits, the second for arrows, the third for summary estimates.}
   \item{pch}{plotting symbol to use for the observed effect sizes or outcomes. By default, a filled square is used. See \code{\link{points}} for other options. Can also be a vector of values.}
   \item{psize}{optional vector with point sizes for the observed effect sizes or outcomes. If unspecified, the point sizes are drawn proportional to the precision of the estimates.}
   \item{cex}{optional character and symbol expansion factor. If unspecified, the function tries to set this to a sensible value.}
   \item{cex.lab}{optional expansion factor for the x-axis titel. If unspecified, the function tries to set this to a sensible value.}
   \item{cex.axis}{optional expansion factor for the x-axis labels. If unspecified, the function tries to set this to a sensible value.}
   \item{col}{optional character string specifying the name of a color to use for the summary polygon or fitted values. If unspecified, the function sets a default color.}
   \item{border}{optional character string specifying the name of a color to use for the border of the summary polygon or fitted values. If unspecified, the function sets a default color.}
   \item{lty}{optional character string specifying the line type for the confidence intervals (if unspecified, the function sets this to \code{"solid"} by default).}
   \item{\dots}{other arguments.}
}
\details{
   The plot shows the individual observed effect sizes or outcomes with corresponding confidence intervals.

   For fixed- and random-effects models (i.e., for models without moderators), a polygon is added to the bottom of the forest plot, showing the summary estimate based on the model (with the outer edges of the polygon indicating the confidence interval limits). The \code{col} and \code{border} arguments can be used to adjust the (border) color of the polygon.

   For random-effects models and if \code{addcred=TRUE}, a dotted line indicates the (approximate) bounds of the credibility/credible or prediction interval (the interval indicates where \code{level} \% of the true effects are expected to fall) (Riley et al., 2011). For random-effects models of class \code{"rma.mv"} (see \code{\link{rma.mv}}) with multiple \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}} values, the \code{addcred} argument can be used to specify for which level of the inner factor the credibility/prediction intervals should be provided (since the credibility/prediction intervals differ depending on the \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}} value). If the model should also contain multiple multiple \if{latex}{\eqn{\gamma^2}}\if{html}{\eqn{\gamma²}}\if{text}{\eqn{\gamma^2}} values, the \code{addcred} argument should then be of length 2 to specify the levels of both inner factors.

   For models involving moderators, the fitted value for each study is added as a polygon to the plot. By default, the width of the polygons corresponds to the confidence interval limits for the fitted values. By setting \code{addcred=TRUE}, the width reflects the credibility/prediction interval limits. Again, the \code{col} and \code{border} arguments can be used to adjust the (border) color of the polygons. These polygons can be suppressed by setting \code{addfit=FALSE}.

   With the \code{transf} argument, the observed effect sizes or outcomes, summary estimate, fitted values, confidence interval bounds, and credibility/prediction interval bounds can be transformed with some suitable function. For example, when plotting log odds ratios, one could use \code{transf=exp} to obtain a forest plot showing the odds ratios. Alternatively, one can use the \code{atransf} argument to transform the x-axis labels and annotations (e.g., \code{atransf=exp}). The examples below illustrate the use of these arguments.

   By default, the studies are ordered from top to bottom (i.e., the first study in the dataset will be placed in row \eqn{k}, the second study in row \eqn{k-1}, and so on, until the last study, which is placed in the first row). The studies can be reordered with the \code{order} argument:
   \itemize{
   \item \code{order="obs"}: the studies are ordered by the observed effect sizes,
   \item \code{order="fit"}: the studies are ordered by the fitted values,
   \item \code{order="prec"}: the studies are ordered by their sampling variances,
   \item \code{order="resid"}: the studies are ordered by the size of their residuals,
   \item \code{order="rstandard"}: the studies are ordered by the size of their standardized residuals,
   \item \code{order="abs.resid"}: the studies are ordered by the size of their absolute residuals,
   \item \code{order="abs.rstandard"}: the studies are ordered by the size of their absolute standardized residuals.
   } Alternatively, it is also possible to set \code{order} equal to a vector with indices specifying the desired order (see examples below).

   Additional summary estimates can also be added to the plot with the \code{\link{addpoly}} function. See the documentation for that function for examples.
}
\note{
   The function tries to set some sensible values for the optional arguments, but it may be necessary to tweak these in certain circumstances.

   The function actually returns some information about the chosen defaults invisibly. Printing this information is useful as a starting point to make adjustments to the plot (see examples below).

   If the number of studies is quite large, the labels, annotations, and symbols may become quite small and impossible to read. Stretching the plot window vertically may then provide a more readable figure (one should call the function again after adjusting the window size, so that the label/symbol sizes can be properly adjusted). Also, the \code{cex}, \code{cex.lab}, and \code{cex.axis} arguments are then useful to adjust the symbol and text sizes.

   If the horizontal plot and/or x-axis limits are set by the user, then the horizontal plot limits (\code{xlim}) must be at least as wide as the x-axis limits (\code{alim}). This restriction is enforced inside the function.

   If the outcome measure used for creating the plot has absolute limits (e.g., correlations are bounded between \code{-1} and \code{1}, proportions are bounded between \code{0} and \code{1}), one can use the \code{clim} argument to enforce those limits (confidence/credibility/prediction intervals cannot exceed those bounds then).

   The models without moderators, the \code{col} argument can also be a vector of two elements, the first for specifying the color of the summary polygon, the second for specifying the color of the line for the credibility/prediction interval.

   The \code{lty} argument can also be a vector of up to three elements, the first for specifying the line type of the individual CIs (\code{"solid"} by default), the second for the line type of the credibility/prediction interval (\code{"dotted"} by default), the third for the line type of the horizontal lines that are automatically added to the plot (\code{"solid"} by default; set to \code{"blank"} to remove them).
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \cr
   package website: \url{http://www.metafor-project.org/} \cr
   author homepage: \url{http://www.wvbauer.com/}
}
\references{
   Lewis, S., & Clarke, M. (2001). Forest plots: Trying to see the wood and the trees. \emph{British Medical Journal}, \bold{322}, 1479--1480.

   Riley, R. D., Higgins, J. P. T., & Deeks, J. J. (2011). Interpretation of random effects meta-analyses. \emph{British Medical Journal}, \bold{342}, d549.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
   \code{\link{forest}}, \code{\link{forest.default}}, \code{\link{addpoly}}
}
\examples{
### meta-analysis of the log relative risks using a random-effects model
res <- rma(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg,
           slab=paste(author, year, sep=", "))

### default forest plot of the log relative risks and summary estimate
forest(res)

### summary estimate in row -1; studies in rows k=13 through 1; horizontal
### lines in rows 0 and k+1; and two extra lines of space at the top for
### headings or other annotations
op <- par(xpd=TRUE)
text(x=-8.5, y=-1:16, -1:16, pos=4, cex=.5)
par(op)

### can also inspect defaults chosen
defaults <- forest(res)
defaults

### several forest plots illustrating the use of various arguments
forest(res, cex=.8)
forest(res, cex=.8, addcred=TRUE)
forest(res, cex=.8, alim=c(-3,3))
forest(res, cex=.8, order="prec", alim=c(-3,3))
forest(res, cex=.8, order=order(dat.bcg$ablat), addcred=TRUE)

### adjust xlim values to see how that changes the plot
forest(res)
par("usr")[1:2] ### this shows what xlim values were chosen by default
forest(res, xlim=c(-16,14))
forest(res, xlim=c(-18,10))
forest(res, xlim=c(-10,10))

### illustrate transf argument
forest(res, transf=exp, at=c(0,1,2,4,6), xlim=c(-8,12), cex=.8, refline=1)

### illustrate atransf argument
forest(res, atransf=exp, at=log(c(.05,.25,1,4,20)), xlim=c(-8,7), cex=.8)

### showweights argument
forest(res, atransf=exp, at=log(c(.05,.25,1,4,20)), xlim=c(-8,8),
       order="prec", showweights=TRUE, cex=.8)

### forest plot with extra annotations (may need to stretch plotting device)
forest(res, atransf=exp, at=log(c(.05, .25, 1, 4)), xlim=c(-16,6),
       ilab=cbind(dat.bcg$tpos, dat.bcg$tneg, dat.bcg$cpos, dat.bcg$cneg),
       ilab.xpos=c(-9.5,-8,-6,-4.5), cex=.75)
op <- par(cex=.75, font=2)
text(c(-9.5,-8,-6,-4.5), 15, c("TB+", "TB-", "TB+", "TB-"))
text(c(-8.75,-5.25),     16, c("Vaccinated", "Control"))
text(-16,                15, "Author(s) and Year",     pos=4)
text(6,                  15, "Relative Risk [95\% CI]", pos=2)
par(op)

### mixed-effects model with absolute latitude in the model
res <- rma(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, mods = ~ ablat,
           data=dat.bcg, slab=paste(author, year, sep=", "))

### forest plot with observed and fitted values
forest(res, xlim=c(-9,5), order="fit", cex=.8, ilab=dat.bcg$ablat,
       ilab.xpos=-3.5, atransf=exp, at=log(c(.05,.25,1,4)))
op <- par(font=2)
text(-9,   15, "Author(s) and Year",    pos=4, cex=.8)
text( 5,   15, "Observed RR [95\% CI]",  pos=2, cex=.8)
text(-3.5, 15, "Latitude",                     cex=.8)
par(op)

### meta-analysis of the log relative risks using a random-effects model
res <- rma(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg,
           slab=paste(author, year, sep=", "))

### for more complicated plots, the ylim and rows arguments may be useful
forest(res)
forest(res, ylim=c(-1.5, 16)) ### the default
forest(res, ylim=c(-1.5, 20)) ### extra space in plot
forest(res, ylim=c(-1.5, 20), rows=c(17:15, 12:6, 3:1)) ### set positions

### forest plot with subgrouping of studies (may need to stretch plotting device)
forest(res, xlim=c(-16, 6), at=log(c(.05, .25, 1, 4)), atransf=exp,
       ilab=cbind(dat.bcg$tpos, dat.bcg$tneg, dat.bcg$cpos, dat.bcg$cneg),
       ilab.xpos=c(-9.5,-8,-6,-4.5), cex=.75, ylim=c(-1, 21),
       order=order(dat.bcg$alloc), rows=c(1:2,5:11,14:17))
op <- par(cex=.75, font=4)
text(-16, c(18,12,3), c("Systematic Allocation", "Random Allocation",
                        "Alternate Allocation"), pos=4)
par(cex=.75, font=2)
text(c(-9.5,-8,-6,-4.5), 20, c("TB+", "TB-", "TB+", "TB-"))
text(c(-8.75,-5.25),     21, c("Vaccinated", "Control"))
text(-16,                20, "Author(s) and Year",     pos=4)
text(6,                  20, "Relative Risk [95\% CI]", pos=2)
par(op)

### see also the addpoly.rma function for an example where summaries
### for the three subgroups are added to such a forest plot

### illustrate use of clim argument with a meta-analysis of raw correlation
### coefficients (data from Pritz, 1997); without clim=c(0,1), some of the
### CIs would have upper bounds larger than 1
dat <- escalc(measure="PR", xi=xi, ni=ni, data=dat.pritz1997)
res <- rma(yi, vi, data=dat, slab=paste0(study, ") ", authors))
forest(res, xlim=c(-1,2), alim=c(0,1), psize=1, refline=coef(res), clim=c(0,1))
text(-1, 16, "Author(s)", pos=4)
text( 2, 16, "Proportion [95\% CI]", pos=2)
}
\keyword{hplot}
