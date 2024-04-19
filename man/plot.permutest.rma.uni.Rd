\name{plot.permutest.rma.uni}
\alias{plot.permutest.rma.uni}
\title{Plot Method for 'permutest.rma.uni' Objects}
\description{
   Function to plot objects of class \code{"permutest.rma.uni"}.
}
\usage{
\method{plot}{permutest.rma.uni}(x, beta, alpha, QM=FALSE, QS=FALSE,
     breaks="Scott", freq=FALSE, col, border, col.out, col.ref, col.density,
     trim=0, adjust=1, lwd=c(2,0,0,4), layout, legend=FALSE, \dots)
}
\arguments{
   \item{x}{an object of class \code{"permutest.rma.uni"} obtained with \code{\link{permutest}}.}
   \item{beta}{optional vector of indices to specify which (location) coefficients should be plotted.}
   \item{alpha}{optional vector of indices to specify which scale coefficients should be plotted. Only relevant for location-scale models (see \code{\link{rma.uni}}).}
   \item{QM}{logical to specify whether the permutation distribution of the omnibus test of the (location) coefficients should be plotted (the default is \code{FALSE}).}
   \item{QS}{logical to specify whether the permutation distribution of the omnibus test of the scale coefficients should be plotted (the default is \code{FALSE}). Only relevant for location-scale models (see \code{\link{rma.uni}}).}
   \item{breaks}{argument to be passed on to the corresponding argument of \code{\link{hist}} to set (the method for determining) the (number of) breakpoints.}
   \item{freq}{logical to specify whether frequencies or probability densities should be plotted (the default is \code{FALSE} to plot densities).}
   \item{col}{optional character string to specify the color of the histogram bars.}
   \item{border}{optional character string to specify the color of the borders around the bars.}
   \item{col.out}{optional character string to specify the color of the bars that are more extreme than the observed test statistic (the default is a semi-transparent shade of red).}
   \item{col.ref}{optional character string to specify the color of the theoretical reference/null distribution that is superimposed on top of the histogram (the default is a dark shade of gray).}
   \item{col.density}{optional character string to specify the color of the kernel density estimate of the permutation distribution that is superimposed on top of the histogram (the default is blue).}
   \item{trim}{the fraction (up to 0.5) of observations to be trimmed from the tails of each permutation distribution before its histogram is plotted.}
   \item{adjust}{numeric value to be passed on to the corresponding argument of \code{\link{density}} (for adjusting the bandwidth of the kernel density estimate).}
   \item{lwd}{numeric vector to specify the width of the vertical lines corresponding to the value of the observed test statistic, of the theoretical reference/null distribution, of the density estimate, and of the vertical line at 0 (note: by default, the theoretical reference/null distribution and the density estimate both have a line width of 0 and are therefore not plotted).}
   \item{layout}{optional vector of two numbers to specify the number of rows and columns for the layout of the figure.}
   \item{legend}{logical to specify whether a legend should be added to the plot (the default is \code{FALSE}). Can also be a keyword to specify the position of the legend (see \code{\link{legend}}).}
   \item{\dots}{other arguments.}
}
\details{
   The function plots the permutation distribution of each model coefficient as a histogram.

   For models with moderators, one can choose via argument \code{beta} which coefficients to plot (by default, all permutation distributions except that of the intercept are plotted). One can also choose to plot the permutation distribution of the omnibus test of the model coefficients (by setting \code{QM=TRUE}).

   Arguments \code{breaks}, \code{freq}, \code{col}, and \code{border} are passed on to the \code{\link{hist}} function for the plotting.

   Argument \code{trim} can be used to trim away a certain fraction of observations from the tails of each permutation distribution before its histogram is plotted. By setting this to a value above 0, one can quickly remove some of the extreme values that might lead to the bulk of the distribution getting squished together at the center (typically, a small value such as \code{trim=0.01} is sufficient for this purpose).

   The observed test statistic is indicated as a vertical dashed line (in both tails for a two-sided test).

   Argument \code{col.out} is used to specify the color for the bars in the histogram that are more extreme than the observed test statistic. The p-value of a permutation test corresponds to the area of these bars.

   One can superimpose the theoretical reference/null distribution on top of the histogram (i.e., the distribution as assumed by the model). The p-value for the standard (i.e., non-permutation) test is the area that is more extreme than the observed test statistic under this reference/null distribution.

   A kernel density estimate of the permutation distribution can also be superimposed on top of the histogram (as a smoothed representation of the permutation distribution).

   Note that the theoretical reference/null distribution and the kernel density estimate of the permutation distribution are only shown when setting the line width for these elements greater than 0 via the \code{lwd} argument (e.g., \code{lwd=c(2,2,2,4)}).

   For location-scale models (see \code{\link{rma.uni}} for details), one can also use arguments \code{alpha} and \code{QS} to specify which scale coefficients to plot and whether to also plot the permutation distribution of the omnibus test of the scale coefficients (by setting \code{QS=TRUE}).
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link[=permutest.rma.uni]{permutest}} for the function to create \code{permutest.rma.uni} objects.
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### random-effects model
res <- rma(yi, vi, data=dat)
res

\dontrun{
### permutation test (exact)
permres <- permutest(res, exact=TRUE)
permres

### plot of the permutation distribution
### dashed horizontal line: the observed value of the test statistic (in both tails)
### black curve: standard normal density (theoretical reference/null distribution)
### blue curve: kernel density estimate of the permutation distribution
plot(permres, lwd=c(2,3,3,4))

### mixed-effects model with two moderators (absolute latitude and publication year)
res <- rma(yi, vi, mods = ~ ablat + year, data=dat)
res

### permutation test (approximate)
set.seed(1234) # for reproducibility
permres <- permutest(res, iter=10000)
permres

### plot of the permutation distribution for absolute latitude
### note: the tail area under the permutation distribution is larger
### than under a standard normal density (hence, the larger p-value)
plot(permres, beta=2, lwd=c(2,3,3,4), xlim=c(-5,5))
}
}
\keyword{hplot}
