\name{forest.default}
\alias{forest.default}
\title{Forest Plots (Default Method)}
\description{
   Function to create forest plots for a given set of data. \loadmathjax
}
\usage{
\method{forest}{default}(x, vi, sei, ci.lb, ci.ub,
       annotate=TRUE, showweights=FALSE, header=FALSE,
       xlim, alim, olim, ylim, at, steps=5,
       level=95, refline=0, digits=2L, width,
       xlab, slab, ilab, ilab.lab, ilab.xpos, ilab.pos,
       order, subset, transf, atransf, targs, rows,
       efac=1, pch, psize, plim=c(0.5,1.5), col,
       shade, colshade, lty, fonts, cex, cex.lab, cex.axis, \dots)
}
\arguments{
   \item{x}{vector of length \mjseqn{k} with the observed effect sizes or outcomes.}
   \item{vi}{vector of length \mjseqn{k} with the corresponding sampling variances.}
   \item{sei}{vector of length \mjseqn{k} with the corresponding standard errors (note: only one of the two, \code{vi} or \code{sei}, needs to be specified).}
   \item{ci.lb}{vector of length \mjseqn{k} with the corresponding lower confidence interval bounds. Not needed if \code{vi} or \code{sei} is specified. See \sQuote{Details}.}
   \item{ci.ub}{vector of length \mjseqn{k} with the corresponding upper confidence interval bounds. Not needed if \code{vi} or \code{sei} is specified. See \sQuote{Details}.}
   \item{annotate}{logical to specify whether annotations should be added to the plot (the default is \code{TRUE}).}
   \item{showweights}{logical to specify whether the annotations should also include the inverse variance weights (the default is \code{FALSE}).}
   \item{header}{logical to specify whether column headings should be added to the plot (the default is \code{FALSE}). Can also be a character vector to specify the left and right headings  (or only the left one).}
   \item{xlim}{horizontal limits of the plot region. If unspecified, the function sets the horizontal plot limits to some sensible values.}
   \item{alim}{the x-axis limits. If unspecified, the function sets the x-axis limits to some sensible values.}
   \item{olim}{optional argument to specify observation/outcome limits. If unspecified, no limits are used.}
   \item{ylim}{the y-axis limits of the plot. If unspecified, the function sets the y-axis limits to some sensible values. Can also be a single value to set the lower bound (while the upper bound is still set automatically).}
   \item{at}{position of the x-axis tick marks and corresponding labels. If unspecified, the function sets the tick mark positions/labels to some sensible values.}
   \item{steps}{the number of tick marks for the x-axis (the default is 5). Ignored when the positions are specified via the \code{at} argument.}
   \item{level}{numeric value between 0 and 100 to specify the confidence interval level (the default is 95; see \link[=misc-options]{here} for details).}
   \item{refline}{numeric value to specify the location of the vertical \sQuote{reference} line (the default is 0). The line can be suppressed by setting this argument to \code{NA}. Can also be a vector to add multiple lines.}
   \item{digits}{integer to specify the number of decimal places to which the tick mark labels of the x-axis and the annotations should be rounded (the default is \code{2L}). Can also be a vector of two integers, the first to specify the number of decimal places for the annotations, the second for the x-axis labels (when \code{showweights=TRUE}, can also specify a third value for the weights). When specifying an integer (e.g., \code{2L}), trailing zeros after the decimal mark are dropped for the x-axis labels. When specifying a numeric value (e.g., \code{2}), trailing zeros are retained.}
   \item{width}{optional integer to manually adjust the width of the columns for the annotations (either a single integer or a vector of the same length as the number of annotation columns).}
   \item{xlab}{title for the x-axis. If unspecified, the function sets an appropriate axis title. Can also be a vector of three/two values (to also/only add labels at the end points of the x-axis limits).}
   \item{slab}{optional vector with labels for the \mjseqn{k} studies. If unspecified, the function tries to extract study labels from \code{x} and otherwise simple labels are created within the function. To suppress labels, set this argument to \code{NA}.}
   \item{ilab}{optional vector, matrix, or data frame providing additional information about the studies that should be added to the plot.}
   \item{ilab.lab}{optional character vector with (column) labels for the variable(s) given via \code{ilab}.}
   \item{ilab.xpos}{optional numeric vector to specify the horizontal position(s) of the variable(s) given via \code{ilab}.}
   \item{ilab.pos}{integer(s) (either 1, 2, 3, or 4) to specify the alignment of the variable(s) given via \code{ilab} (2 means right, 4 means left aligned). If unspecified, the default is to center the values.}
   \item{order}{optional character string to specify how the studies should be ordered. Can also be a variable based on which the studies will be ordered. See \sQuote{Details}.}
   \item{subset}{optional (logical or numeric) vector to specify the subset of studies that should be included in the plot.}
   \item{transf}{optional argument to specify a function to transform the observed outcomes and corresponding confidence interval bounds (e.g., \code{transf=exp}; see also \link{transf}). If unspecified, no transformation is used.}
   \item{atransf}{optional argument to specify a function to transform the x-axis labels and annotations (e.g., \code{atransf=exp}; see also \link{transf}). If unspecified, no transformation is used.}
   \item{targs}{optional arguments needed by the function specified via \code{transf} or \code{atransf}.}
   \item{rows}{optional vector to specify the rows (or more generally, the horizontal positions) for plotting the outcomes. Can also be a single value to specify the row (horizontal position) of the first outcome (the remaining outcomes are then plotted below this starting row).}
   \item{efac}{vertical expansion factor for confidence interval limits and arrows. The default value of 1 should usually work okay. Can also be a vector of two numbers, the first for CI limits, the second for arrows.}
   \item{pch}{plotting symbol to use for the observed outcomes. By default, a filled square is used. See \code{\link{points}} for other options. Can also be a vector of values.}
   \item{psize}{optional numeric value to specify the point sizes for the observed outcomes. If unspecified, the point sizes are a function of the precision of the estimates. Can also be a vector of values.}
   \item{plim}{numeric vector of length 2 to scale the point sizes (ignored when \code{psize} is specified). See \sQuote{Details}.}
   \item{col}{optional character string to specify the color of the observed outcomes. Can also be a vector.}
   \item{shade}{optional character string or a (logical or numeric) vector for shading rows of the plot. See \sQuote{Details}.}
   \item{colshade}{optional argument to specify the color for the shading.}
   \item{lty}{optional character string to specify the line type for the confidence intervals. If unspecified, the function sets this to \code{"solid"} by default.}
   \item{fonts}{optional character string to specify the font for the study labels, annotations, and the extra information (if specified via \code{ilab}). If unspecified, the default font is used.}
   \item{cex}{optional character and symbol expansion factor. If unspecified, the function sets this to a sensible value.}
   \item{cex.lab}{optional expansion factor for the x-axis title. If unspecified, the function sets this to a sensible value.}
   \item{cex.axis}{optional expansion factor for the x-axis labels. If unspecified, the function sets this to a sensible value.}
   \item{\dots}{other arguments.}
}
\details{
   The plot shows the observed effect sizes or outcomes with corresponding confidence intervals. To use the function, one should specify the observed outcomes (via the \code{x} argument) together with the corresponding sampling variances (via the \code{vi} argument) or with the corresponding standard errors (via the \code{sei} argument). Alternatively, one can specify the observed outcomes together with the corresponding confidence interval bounds (via the \code{ci.lb} and \code{ci.ub} arguments).

   With the \code{transf} argument, the observed outcomes and corresponding confidence interval bounds can be transformed with some suitable function. For example, when plotting log odds ratios, then one could use \code{transf=exp} to obtain a forest plot showing the odds ratios. Alternatively, one can use the \code{atransf} argument to transform the x-axis labels and annotations (e.g., \code{atransf=exp}). See also \link{transf} for some other useful transformation functions in the context of a meta-analysis. The examples below illustrate the use of these arguments.

   By default, the studies are ordered from top to bottom (i.e., the first study in the dataset will be placed in row \mjseqn{k}, the second study in row \mjseqn{k-1}, and so on, until the last study, which is placed in the first row). The studies can be reordered with the \code{order} argument:

   \itemize{
   \item \code{order="obs"}: the studies are ordered by the observed outcomes,
   \item \code{order="prec"}: the studies are ordered by their sampling variances.
   }

   Alternatively, it is also possible to set \code{order} equal to a variable based on which the studies will be ordered (see \sQuote{Examples}).

   Additional columns with information about the studies can be added to the plot via the \code{ilab} argument. This can either be a single variable or an entire matrix / data frame (with as many rows as there are studies in the forest plot). The \code{ilab.xpos} argument can be used to specify the horizontal position of the variables specified via \code{ilab}. The \code{ilab.pos} argument can be used to specify how the variables should be aligned.

   Summary estimates can be added to the plot with the \code{\link{addpoly}} function. See the documentation for that function for examples.

   By default (i.e., when \code{psize} is not specified), the point sizes are a function of the precision (i.e., inverse standard errors) of the outcomes. This way, more precise estimates are visually more prominent in the plot. By making the point sizes a function of the inverse standard errors of the estimates, their areas are proportional to the inverse sampling variances, which corresponds to the weights they would receive in an equal-effects model. However, the point sizes are rescaled so that the smallest point size is \code{plim[1]} and the largest point size is \code{plim[2]}. As a result, their relative sizes (i.e., areas) no longer exactly correspond to their relative weights in such a model. If exactly relative point sizes are desired, one can set \code{plim[2]} to \code{NA}, in which case the points are rescaled so that the smallest point size corresponds to \code{plim[1]} and all other points are scaled accordingly. As a result, the largest point may be very large. Alternatively, one can set \code{plim[1]} to \code{NA}, in which case the points are rescaled so that the largest point size corresponds to \code{plim[2]} and all other points are scaled accordingly. As a result, the smallest point may be very small and essentially indistinguishable from the confidence interval line. To avoid the latter, one can also set \code{plim[3]}, which enforces a minimal point size.

   With the \code{shade} argument, one can shade rows of the plot. The argument can be set to one of the following character strings: \code{"zebra"} (same as \code{shade=TRUE}) or \code{"zebra2"} to use zebra-style shading (starting either at the first or second study) or to \code{"all"} in which case all rows are shaded. Alternatively, the argument can be set to a logical or numeric vector to specify which rows should be shaded. The \code{colshade} argument can be used to set the color of shaded rows.
}
\section{Note}{
   The function sets some sensible values for the optional arguments, but it may be necessary to adjust these in certain circumstances.

   The function actually returns some information about the chosen values invisibly. Printing this information is useful as a starting point to make adjustments to the plot.

   If the number of studies is quite large, the labels, annotations, and symbols may become quite small and impossible to read. Stretching the plot window vertically may then provide a more readable figure (one should call the function again after adjusting the window size, so that the label/symbol sizes can be properly adjusted). Also, the \code{cex}, \code{cex.lab}, and \code{cex.axis} arguments are then useful to adjust the symbol and text sizes.

   If the outcome measure used for creating the plot is bounded (e.g., correlations are bounded between -1 and +1, proportions are bounded between 0 and 1), one can use the \code{olim} argument to enforce those limits (the observed outcomes and confidence intervals cannot exceed those bounds then).

   The \code{lty} argument can also be a vector of two elements, the first for specifying the line type of the individual CIs (\code{"solid"} by default), the second for the line type of the horizontal line that is automatically added to the plot (\code{"solid"} by default; set to \code{"blank"} to remove it).
}
\section{Additional Optional Arguments}{
   There are some additional optional arguments that can be passed to the function via \code{...} (hence, they cannot be abbreviated):
   \describe{
   \item{top}{single numeric value to specify the amount of space (in terms of number of rows) to leave empty at the top of the plot (e.g., for adding headers). The default is 3.}
   \item{annosym}{vector of length 3 to select the left bracket, separation, and right bracket symbols for the annotations. The default is \code{c(" [", ", ", "]")}. Can also include a 4th element to adjust the look of the minus symbol, for example to use a proper minus sign (\ifelse{latex}{\mjseqn{-}}{\enc{−}{-}}) instead of a hyphen-minus (-). Can also include a 5th element that should be a space-like symbol (e.g., an \sQuote{en space}) that is used in place of numbers (only relevant when trying to line up numbers exactly). For example, \code{annosym=c(" [", ", ", "]", "\u2212", "\u2002")} would use a proper minus sign and an \sQuote{en space} for the annotations. The decimal point character can be adjusted via the \code{OutDec} argument of the \code{\link{options}} function before creating the plot (e.g., \code{options(OutDec=",")}).}
   \item{tabfig}{single numeric value (either a 1, 2, or 3) to set \code{annosym} automatically to a vector that will exactly align the numbers in the annotations when using a font that provides \sQuote{tabular figures}. Value 1 corresponds to using \code{"\u2212"} (a minus) and \code{"\u2002"} (an \sQuote{en space}) in \code{annoyym} as shown above. Value 2 corresponds to \code{"\u2013"} (an \sQuote{en dash}) and \code{"\u2002"} (an \sQuote{en space}). Value 3 corresponds to \code{"\u2212"} (a minus) and \code{"\u2007"} (a \sQuote{figure space}). The appropriate value for this argument depends on the font used. For example, for fonts Calibri and Carlito, 1 or 2 should work; for fonts Source Sans 3 and Palatino Linotype, 1, 2, and 3 should all work; for Computer/Latin Modern and Segoe UI, 2 should work; for Lato, Roboto, and Open Sans (and maybe Arial), 3 should work. Other fonts may work as well, but this is untested.}
   \item{textpos}{numeric vector of length 2 to specify the placement of the study labels and the annotations. The default is to use the horizontal limits of the plot region, i.e., the study labels to the right of \code{xlim[1]} and the annotations to the left of \code{xlim[2]}.}
   \item{rowadj}{numeric vector of length 3 to vertically adjust the position of the study labels, the annotations, and the extra information (if specified via \code{ilab}). This is useful for fine-tuning the position of text added with different positional alignments (i.e., argument \code{pos} in the \code{\link{text}} function).}
   }
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Lewis, S., & Clarke, M. (2001). Forest plots: Trying to see the wood and the trees. \emph{British Medical Journal}, \bold{322}(7300), 1479--1480. \verb{https://doi.org/10.1136/bmj.322.7300.1479}

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{forest}} for an overview of the various \code{forest} functions and especially \code{\link{forest.rma}} for the function draw forest plots including a summary polygon.

   \code{\link{addpoly}} for a function to add polygons to forest plots.
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg,
              data=dat.bcg, slab=paste(author, year, sep=", "))

### default forest plot of the observed log risk ratios
forest(dat$yi, dat$vi, header=TRUE)

### the with() function can be used to avoid having to retype dat$... over and over
with(dat, forest(yi, vi, header=TRUE))

### forest plot of the observed risk ratios (transform outcomes)
with(dat, forest(yi, vi, transf=exp, alim=c(0,2), steps=5,
                 xlim=c(-2.5,4), refline=1, cex=0.9, header=TRUE, top=2))

### forest plot of the observed risk ratios (transformed x-axis)
with(dat, forest(yi, vi, atransf=exp, at=log(c(0.05,0.25,1,4,20)),
                 xlim=c(-10,8), cex=0.9, header=TRUE, top=2))

### forest plot of the observed risk ratios with studies ordered by the RRs
with(dat, forest(yi, vi, atransf=exp, at=log(c(0.05,0.25,1,4,20)),
                 xlim=c(-10,8), cex=0.9, header=TRUE, top=2, order="obs"))

### forest plot of the observed risk ratios with studies ordered by absolute latitude
with(dat, forest(yi, vi, atransf=exp, at=log(c(0.05,0.25,1,4,20)),
                 xlim=c(-10,8), cex=0.9, header=TRUE, top=2, order=ablat))

### see also examples for the forest.rma function
}
\keyword{hplot}
