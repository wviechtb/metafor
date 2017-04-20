\name{radial}
\alias{radial}
\alias{galbraith}
\alias{radial.rma}
\title{Radial (Galbraith) Plots for 'rma' Objects}
\description{Function to create radial (also called Galbraith) plots for objects of class \code{"rma"}.}
\usage{
radial(x, \dots)
galbraith(x, \dots)

\method{radial}{rma}(x, center=FALSE, xlim, zlim, xlab, zlab,
       atz, aty, steps=7, level=x$level, digits=2, back="lightgray",
       transf, targs, pch=19, arc.res=100, cex, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma"}.}
   \item{center}{logical to indicate whether the plot should be centered horizontally at the model estimate (the default is \code{FALSE}).}
   \item{xlim}{x-axis limits. If unspecified, the function tries to set the x-axis limits to some sensible values.}
   \item{zlim}{z-axis limits. If unspecified, the function tries to set the z-axis limits to some sensible values (note that the z-axis limits are the actual vertical limit of the plotting region).}
   \item{xlab}{title for the x-axis. If unspecified, the function tries to set an appropriate axis title.}
   \item{zlab}{title for the z-axis. If unspecified, the function tries to set an appropriate axis title.}
   \item{atz}{position for the z-axis tick marks and labels. If unspecified, these values are set by the function.}
   \item{aty}{position for the y-axis tick marks and labels. If unspecified, these values are set by the function.}
   \item{steps}{the number of tick marks for the y-axis (the default is 7). Ignored when argument \code{aty} is used.}
   \item{level}{numerical value between 0 and 100 specifying the level of the z-axis error region (the default is to take the value from the object).}
   \item{digits}{integer specifying the number of decimal places to which the tick mark labels of the y-axis should be rounded (the default is 2).}
   \item{back}{color of the z-axis error region. Set to \code{NA} to suppress shading of the region.}
   \item{transf}{optional argument specifying the name of a function that should be used to transform the y-axis labels (e.g., \code{transf=exp}; see also \link{transf}). If unspecified, no transformation is used.}
   \item{targs}{optional arguments needed by the function specified via \code{transf}.}
   \item{pch}{plotting symbol. By default, a solid circle is used. See \code{\link{points}} for other options.}
   \item{arc.res}{integer specifying the number of line segments to use when drawing the y-axis and confidence interval arcs (the default is 100).}
   \item{cex}{optional character and symbol expansion factor. If unspecified, the function tries to set this to a sensible value.}
   \item{\dots}{other arguments.}
}
\details{
   For a fixed-effects model, the plot shows the inverse of the standard errors on the horizontal axis against the individual observed effect sizes or outcomes standardized by their corresponding standard errors on the vertical axis. Since the vertical axis corresponds to standardized values, it is referred to as the z-axis within this function. On the right hand side of the plot, an arc is drawn (referred to as the y-axis within this function) corresponding to the individual observed effect sizes or outcomes. A line projected from (0,0) through a particular point within the plot onto this arc indicates the value of the individual observed effect size or outcome for that point.

   For a random-effects model, the function uses \if{latex}{\eqn{1/\sqrt{v_i + \tau^2}}}\if{html}{\eqn{1/\sqrt(vᵢ + \tau²)}}\if{text}{\eqn{1/\sqrt(v_i + \tau^2)}} for the horizontal axis, where \if{latex}{\eqn{v_i}}\if{html}{\eqn{vᵢ}}\if{text}{\eqn{v_i}} is the sampling variance of the observed effect size or outcome and \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}} is the amount of heterogeneity as estimated based on the model. For the z-axis, \if{latex}{\eqn{\sqrt{v_i + \tau^2}}}\if{html}{\eqn{\sqrt(vᵢ + \tau²)}}\if{text}{\eqn{\sqrt(v_i + \tau^2)}} is used to standardize the individual observed effect sizes or outcomes.

   If the model contains moderators, the function returns an error.
}
\value{
   A data frame with components:
   \item{x}{the x coordinates of the points that were plotted.}
   \item{y}{the y coordinates of the points that were plotted.}
   \item{slab}{study labels of the points that were plotted.}
   Note that the data frame is returned invisibly.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \cr
   package website: \url{http://www.metafor-project.org/} \cr
   author homepage: \url{http://www.wvbauer.com/}
}
\references{
   Galbraith, R. F. (1988). Graphical display of estimates having differing standard errors. \emph{Technometrics}, \bold{30}, 271--281.

   Galbraith, R. F. (1988). A note on graphical presentation of estimated odds ratios from several clinical trials. \emph{Statistics in Medicine}, \bold{7}, 889--894.

   Galbraith, R. F (1994). Some applications of radial plots. \emph{Journal of the American Statistical Association}, \bold{89}, 1232--1242.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}, \code{\link{rma.glmm}}, \code{\link{rma.mv}}
}
\examples{
### meta-analysis of the log risk ratios using a fixed-effects model
res <- rma(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, method="FE")
radial(res)

### line from (0,0) with slope equal to the log risk ratio from the 4th study
abline(a=0, b=-1.44155119, lty="dotted")

### meta-analysis of the log risk ratios using a random-effects model
res <- rma(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
radial(res)
}
\keyword{hplot}
