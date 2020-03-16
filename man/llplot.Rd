\name{llplot}
\alias{llplot}
\title{Likelihood Plot of a Parameter Corresponding to an Effect Size or Outcome Measure}
\description{
   Function to plot the likelihood of a certain parameter corresponding to an effect size or outcome measure given the study data.
}
\usage{
llplot(measure, yi, vi, sei, ai, bi, ci, di, n1i, n2i, data, subset, drop00=TRUE,
       xvals=1000, xlim, ylim, xlab, ylab, scale=TRUE,
       lty, lwd, col, level=99.99, refline=0, \dots)
}
\arguments{
   \item{measure}{a character string indicating for which effect size or outcome measure the likelihoods should be calculated. See \sQuote{Details} for possible options and how the data should be specified.}
   \item{yi}{vector with the observed effect sizes or outcomes.}
   \item{vi}{vector with the corresponding sampling variances.}
   \item{sei}{vector to specify the corresponding standard.}
   \item{ai}{vector to specify the \eqn{2 \times 2}{2x2} table frequencies (upper left cell).}
   \item{bi}{vector to specify the \eqn{2 \times 2}{2x2} table frequencies (upper right cell).}
   \item{ci}{vector to specify the \eqn{2 \times 2}{2x2} table frequencies (lower left cell).}
   \item{di}{vector to specify the \eqn{2 \times 2}{2x2} table frequencies (lower right cell).}
   \item{n1i}{vector to specify the group sizes or row totals (first group/row).}
   \item{n2i}{vector to specify the group sizes or row totals (second group/row).}
   \item{data}{optional data frame containing the variables given to the arguments above.}
   \item{subset}{optional vector indicating the subset of studies that should be used. This can be a logical vector or a numeric vector indicating the indices of the studies to include.}
   \item{drop00}{logical indicating whether studies with no cases (or only cases) in both groups should be dropped. See \sQuote{Details}.}
   \item{xvals}{integer specifying for how many distinct values the likelihood should be evaluated.}
   \item{xlim}{x-axis limits. If unspecified, the function tries to set the x-axis limits to some sensible values.}
   \item{ylim}{y-axis limits. If unspecified, the function tries to set the y-axis limits to some sensible values.}
   \item{xlab}{title for the x-axis. If unspecified, the function tries to set an appropriate axis title.}
   \item{ylab}{title for the y-axis. If unspecified, the function tries to set an appropriate axis title.}
   \item{scale}{logical indicating whether the likelihood values should be scaled, so that the total area under each curve is (approximately) equal to 1.}
   \item{lty}{the line types (either a single value or a vector of length \eqn{k}). If unspecified, the function sets the line types according to some characteristics of the likelihood function. See \sQuote{Details}.}
   \item{lwd}{the line widths (either a single value or a vector of length \eqn{k}). If unspecified, the function sets the widths according to the sampling variances (so that the line is thicker for more precise studies and vice-versa).}
   \item{col}{the line colors (either a single value or a vector of length \eqn{k}). If unspecified, the function uses various shades of gray according to the sampling variances (so that darker shades are used for more precise studies and vice-versa).}
   \item{level}{numerical value between 0 and 100 specifying the plotting limits for each likelihood line in terms of the confidence interval (the default is 99.99).}
   \item{refline}{value at which a vertical \sQuote{reference} line should be drawn (the default is 0). The line can be suppressed by setting this argument to \code{NA}.}
   \item{\dots}{other arguments.}
}
\details{
   At the moment, the function only accepts \code{measure="GEN"} or \code{measure="OR"}.

   For \code{measure="GEN"}, one must specify arguments \code{yi} for the observed effect sizes or outcomes and \code{vi} for the corresponding sampling variances (instead of specifying \code{vi}, one can specify the standard errors via the \code{sei} argument). The function then plots the likelihood of the true outcome based on a normal sampling distribution with observed outcome as given by \code{yi} and variance as given by \code{vi} for each study.

   For \code{measure="OR"}, one must specify arguments \code{ai}, \code{bi}, \code{ci}, and \code{di}, which denote the cell frequencies of the \eqn{2 \times 2}{2x2} tables. Alternatively, one can specify \code{ai}, \code{ci}, \code{n1i}, and \code{n2i}. See \code{\link{escalc}} function for more details. The function then plots the likelihood of the true log odds ratio based on the non-central hypergeometric distribution for each \eqn{2 \times 2}{2x2} table. Since studies with no cases (or only cases) in both groups have a flat likelihood and are not informative about the odds ratio, they are dropped by default (i.e., \code{drop00=TRUE}) and are hence not drawn (if \code{drop00=FALSE}, these likelihood are indicated by dotted lines). For studies that have a single zero count, the MLE of the odds ratio is infinite and these likelihoods are indicated by dashed lines.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{http://www.metafor-project.org/}
}
\references{
   van Houwelingen, H. C., Zwinderman, K. H., & Stijnen, T. (1993). A bivariate approach to meta-analysis. \emph{Statistics in Medicine}, \bold{12}, 2273--2284.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{https://www.jstatsoft.org/v036/i03}.
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.glmm}}
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### draw likelihoods
llplot(measure="GEN", yi=yi, vi=vi, data=dat, lwd=1, refline=NA, xlim=c(-3,2))

### create plot (Figure 2 in van Houwelingen, Zwinderman, & Stijnen, 1993)
llplot(measure="OR", ai=b.xci, n1i=nci, ci=b.xti, n2i=nti, data=dat.collins1985a,
       lwd=1, refline=NA, xlim=c(-4,4), drop00=FALSE)
}
\keyword{hplot}
