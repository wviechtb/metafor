\name{labbe}
\alias{labbe}
\alias{labbe.rma}
\title{L'Abbe Plots for 'rma' Objects}
\description{Function to create L'Abbé plots for objects of class \code{"rma"}.}
\usage{
labbe(x, \dots)

\method{labbe}{rma}(x, xlim, ylim, xlab, ylab,
      add=x$add, to=x$to, transf, targs,
      pch=21, psize, bg="gray", \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma"}. See \sQuote{Details}.}
   \item{xlim}{x-axis limits. If unspecified, the function tries to set the x-axis limits to some sensible values.}
   \item{ylim}{y-axis limits. If unspecified, the function tries to set the y-axis limits to some sensible values.}
   \item{xlab}{title for the x-axis. If unspecified, the function tries to set an appropriate axis title.}
   \item{ylab}{title for the y-axis. If unspecified, the function tries to set an appropriate axis title.}
   \item{add}{See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{to}{See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{transf}{optional argument specifying the name of a function that should be used to transform the outcomes (e.g., \code{transf=exp}). If unspecified, no transformation is used.}
   \item{targs}{optional arguments needed by the function specified under \code{transf}.}
   \item{pch}{plotting symbol to use for the outcomes. By default, a filled circle is used. Can also be a vector of values. See \code{\link{points}} for other options.}
   \item{psize}{optional vector with point sizes for the outcomes. If unspecified, the point sizes are drawn proportional to the precision of the estimates.}
   \item{bg}{color to use for filling the plotting symbol (default is \code{"gray"}). Can also be a vector of values. Set to \code{NA} to make the plotting symbols transparent.}
   \item{\dots}{other arguments.}
}
\details{
   The model specified via \code{x} must be a model without moderators (i.e., either a fixed- or a random-effects model) fitted with either the \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}, or \code{\link{rma.glmm}} function. Moreover, the model must be fitted with \code{measure} set equal to \code{"RD"} (for risk differences), \code{"RR"} (for relative risks), \code{"OR"} (for odds ratios), \code{"AS"} (for arcsine square-root transformed risk differences), \code{"IRR"} (for incidence rate ratios), \code{"IRD"} (for incidence rate differences), or \code{"IRSD"} (for square-root transformed incidence rate differences).

   The function calculates the arm-level outcomes for the two experimental groups (e.g., treatment and control groups) and plots them against each other. In particular, the function plots the raw proportions of the two groups again each other when analyzing risk differences, the log of the proportions when analyzing (log) relative risks, the log odds when analyzing (log) odds ratios, the arcsine square-root transformed proportions when analyzing arcsine square-root transformed risk differences, the raw incidence rates when analyzing incidence rate differences, the log of the incidence rates when analyzing (log) incidence rate ratios, and the square-root transformed incidence rates when analyzing square-root transformed incidence rate differences. The \code{transf} argument can be used to transform these values (for example, \code{transf=exp} to transform the log of the proportions back to raw proportions).

   As described under the documentation for the \code{\link{escalc}} function, zero cells can lead to problems when calculating particular outcomes. Adding a small constant to the cells of the \eqn{2 \times 2}{2x2} tables is a common solution to this problem. By default, the functions adopts the same method for handling zero cells as was done when fitting the model.

   The size of the points is drawn proportional to the precision (inverse standard error) of the outcomes. The solid line corresponds to identical outcomes in the two groups (i.e., the absence of a difference between the two groups). The dashed line indicates the estimated effect based on the fitted model.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \cr
   package website: \url{http://www.metafor-project.org/} \cr
   author homepage: \url{http://www.wvbauer.com/}
}
\references{
   L'Abbé, K. A., Detsky, A. S., & O'Rourke, K. (1987). Meta-analysis in clinical research. \emph{Annals of Internal Medicine}, \bold{107}, 224--233.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}, \code{\link{rma.glmm}}
}
\examples{
### load BCG vaccine data
data(dat.bcg)

### meta-analysis of the log relative risks using a random-effects model
res <- rma(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### default plot
labbe(res)

### funnel plot with risk values on the x- and y-axis
labbe(res, transf=exp)
}
\keyword{hplot}
