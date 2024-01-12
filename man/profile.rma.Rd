\name{profile.rma}
\alias{profile}
\alias{profile.rma}
\alias{profile.rma.uni}
\alias{profile.rma.mv}
\alias{profile.rma.uni.selmodel}
\alias{profile.rma.ls}
\alias{print.profile.rma}
\alias{plot.profile.rma}
\title{Profile Likelihood Plots for 'rma' Objects}
\description{
   Functions to profile the (restricted) log-likelihood for objects of class \code{"rma.uni"}, \code{"rma.mv"}, \code{"rma.uni.selmodel"}, and \code{"rma.ls"}. \loadmathjax
}
\usage{
\method{profile}{rma.uni}(fitted, xlim, ylim, steps=20, lltol=1e-03,
        progbar=TRUE, parallel="no", ncpus=1, cl, plot=TRUE, \dots)

\method{profile}{rma.mv}(fitted, sigma2, tau2, rho, gamma2, phi, xlim, ylim, steps=20, lltol=1e-03,
        progbar=TRUE, parallel="no", ncpus=1, cl, plot=TRUE, \dots)

\method{profile}{rma.uni.selmodel}(fitted, tau2, delta, xlim, ylim, steps=20, lltol=1e-03,
        progbar=TRUE, parallel="no", ncpus=1, cl, plot=TRUE, \dots)

\method{profile}{rma.ls}(fitted, alpha, xlim, ylim, steps=20, lltol=1e-03,
        progbar=TRUE, parallel="no", ncpus=1, cl, plot=TRUE, \dots)

\method{print}{profile.rma}(x, \dots)
\method{plot}{profile.rma}(x, xlim, ylim, pch=19, xlab, ylab, main, refline=TRUE, cline=FALSE, \dots)
}
\arguments{
   \item{fitted}{an object of class \code{"rma.uni"}, \code{"rma.mv"}, \code{"rma.uni.selmodel"}, or \code{"rma.ls"}.}
   \item{x}{an object of class \code{"profile.rma"} (for \code{plot} and \code{print}).}
   \item{sigma2}{optional integer to specify for which \mjseqn{\sigma^2} parameter the likelihood should be profiled.}
   \item{tau2}{optional integer to specify for which \mjseqn{\tau^2} parameter the likelihood should be profiled.}
   \item{rho}{optional integer to specify for which \mjseqn{\rho} parameter the likelihood should be profiled.}
   \item{gamma2}{optional integer to specify for which \mjseqn{\gamma^2} parameter the likelihood should be profiled.}
   \item{phi}{optional integer to specify for which \mjseqn{\phi} parameter the likelihood should be profiled.}
   \item{delta}{optional integer to specify for which \mjseqn{\delta} parameter the likelihood should be profiled.}
   \item{alpha}{optional integer to specify for which \mjseqn{\alpha} parameter the likelihood should be profiled.}
   \item{xlim}{optional vector to specify the lower and upper limit of the parameter over which the profiling should be done. If unspecified, the function sets these limits automatically.}
   \item{ylim}{optional vector to specify the y-axis limits when plotting the profiled likelihood. If unspecified, the function sets these limits automatically.}
   \item{steps}{number of points between \code{xlim[1]} and \code{xlim[2]} (inclusive) for which the likelihood should be evaluated (the default is 20). Can also be a numeric vector of length 2 or longer to specify for which parameter values the likelihood should be evaluated (in this case, \code{xlim} is automatically set to \code{range(steps)} if unspecified).}
   \item{lltol}{numerical tolerance used when comparing values of the profiled log-likelihood with the log-likelihood of the fitted model (the default is 1e-03).}
   \item{progbar}{logical to specify whether a progress bar should be shown (the default is \code{TRUE}).}
   \item{parallel}{character string to specify whether parallel processing should be used (the default is \code{"no"}). For parallel processing, set to either \code{"snow"} or \code{"multicore"}. See \sQuote{Details}.}
   \item{ncpus}{integer to specify the number of processes to use in the parallel processing.}
   \item{cl}{optional cluster to use if \code{parallel="snow"}. If unspecified, a cluster on the local machine is created for the duration of the call.}
   \item{plot}{logical to specify whether the profile plot should be drawn after profiling is finished (the default is \code{TRUE}).}
   \item{pch}{plotting symbol to use. By default, a filled circle is used. See \code{\link{points}} for other options.}
   \item{refline}{logical to specify whether the value of the parameter estimate should be indicated by a dotted vertical line and its log-likelihood value by a dotted horizontal line (the default is \code{TRUE}).}
   \item{cline}{logical to specify whether a horizontal reference line should be added to the plot that indicates the log-likelihood value corresponding to the 95\% profile confidence interval (the default is \code{FALSE}). Can also be a numeric value between 0 and 100 to specify the confidence interval level.}
   \item{xlab}{title for the x-axis. If unspecified, the function sets an appropriate axis title.}
   \item{ylab}{title for the y-axis. If unspecified, the function sets an appropriate axis title.}
   \item{main}{title for the plot. If unspecified, the function sets an appropriate title.}
   \item{\dots}{other arguments.}
}
\details{
   The function fixes a particular parameter of the model and then computes the maximized (restricted) log-likelihood over the remaining parameters of the model. By doing this for a range of values for the parameter that was fixed, a profile of the (restricted) log-likelihood is constructed.

   \subsection{Selecting the Parameter(s) to Profile}{

      The parameters that can be profiled depend on the model object:

      \itemize{
      \item For objects of class \code{"rma.uni"} obtained with the \code{\link{rma.uni}} function, the function profiles over \mjseqn{\tau^2} (not for equal-effects models).
      \item For objects of class \code{"rma.mv"} obtained with the \code{\link{rma.mv}} function, profiling is done by default over all (non-fixed) variance and correlation components of the model. Alternatively, one can use the \code{sigma2}, \code{tau2}, \code{rho}, \code{gamma2}, or \code{phi} arguments to specify over which parameter the profiling should be done. Only one of these arguments can be used at a time. A single integer is used to specify the number of the parameter.
      \item For selection model objects of class \code{"rma.uni.selmodel"} obtained with the \code{\link{selmodel}} function, profiling is done by default over \mjseqn{\tau^2} (for models where this is an estimated parameter) and all (non-fixed) selection model parameters. Alternatively, one can choose to profile only \mjseqn{\tau^2} by setting \code{tau2=TRUE} or one can select one of the selection model parameters to profile by specifying its number via the \code{delta} argument.
      \item For location-scale model objects of class \code{"rma.ls"} obtained with the \code{\link{rma.uni}} function, profiling is done by default over all (non-fixed) \mjseqn{\alpha} parameters that are part of the scale model.
      }

   }

   \subsection{Interpreting a Likelihood Profile}{

      A profile plot should show a single peak at the corresponding ML/REML estimate (assuming that the model was fitted with ML/REML estimation). If \code{refline=TRUE} (the default), the value of the parameter estimate is indicated by a dotted vertical line and its log-likelihood value by a dotted horizontal line. Hence, the intersection of these two lines should correspond to the peak.

      When profiling a variance component (or some other parameter that cannot be negative), the peak may be at zero (if this corresponds to the ML/REML estimate of the parameter). In this case, the profiled log-likelihood should be a monotonically decreasing function of the parameter. Similarly, when profiling a correlation component, the peak may be at -1 or +1.

      If the profiled log-likelihood has multiple peaks, this indicates that the likelihood surface is not unimodal. In such cases, the ML/REML estimate may correspond to a local optimum (when the intersection of the two dotted lines is not at the highest peak).

      If the profile is flat (over the entire parameter space or large portions of it), then this suggests that at least some of the parameters of the model are not identifiable (and the parameter estimates obtained are to some extent arbitrary). See Raue et al. (2009) for some further discussion of parameter identifiability (structurally and practically) and the use of profile likelihoods to check for this.

      The function checks whether any profiled log-likelihood value is actually larger than the log-likelihood of the fitted model (using a numerical tolerance of \code{lltol}). If so, a warning is issued as this might indicate that the optimizer did not identify the actual ML/REML estimate of the parameter profiled.

   }

   \subsection{Parallel Processing}{

      Profiling requires repeatedly refitting the model, which can be slow when \mjseqn{k} is large and/or the model is complex (the latter especially applies to \code{"rma.mv"} objects and also to certain \code{"rma.uni.selmodel"} or \code{"rma.ls"} objects). On machines with multiple cores, one can try to speed things up by delegating the model fitting to separate worker processes, that is, by setting \code{parallel="snow"} or \code{parallel="multicore"} and \code{ncpus} to some value larger than 1. Parallel processing makes use of the \code{\link[parallel]{parallel}} package, using the \code{\link[parallel]{makePSOCKcluster}} and \code{\link[parallel]{parLapply}} functions when \code{parallel="snow"} or using \code{\link[parallel]{mclapply}} when \code{parallel="multicore"} (the latter only works on Unix/Linux-alikes). With \code{parallel::detectCores()}, one can check on the number of available cores on the local machine.

   }

}
\value{
   An object of class \code{"profile.rma"}. The object is a list (or list of such lists) containing the following components:

   One of the following (depending on the parameter that was actually profiled):
   \item{sigma2}{values of \mjseqn{\sigma^2} over which the likelihood was profiled.}
   \item{tau2}{values of \mjseqn{\tau^2} over which the likelihood was profiled.}
   \item{rho}{values of \mjseqn{\rho} over which the likelihood was profiled.}
   \item{gamma2}{values of \mjseqn{\gamma^2} over which the likelihood was profiled.}
   \item{phi}{values of \mjseqn{\phi} over which the likelihood was profiled.}
   \item{delta}{values of \mjseqn{\delta} over which the likelihood was profiled.}
   \item{alpha}{values of \mjseqn{\alpha} over which the likelihood was profiled.}

   In addition, the following components are included:
   \item{ll}{(restricted) log-likelihood values at the corresponding parameter values.}
   \item{beta}{a matrix with the estimated model coefficients at the corresponding parameter values.}
   \item{ci.lb}{a matrix with the lower confidence interval bounds of the model coefficients at the corresponding parameter values.}
   \item{ci.ub}{a matrix with the upper confidence interval bounds of the model coefficients at the corresponding parameter values.}
   \item{\dots}{some additional elements/values.}
   Note that the list is returned invisibly.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Raue, A., Kreutz, C., Maiwald, T., Bachmann, J., Schilling, M., Klingmuller, U., & Timmer, J. (2009). Structural and practical identifiability analysis of partially observed dynamical models by exploiting the profile likelihood. \emph{Bioinformatics}, \bold{25}(15), 1923--1929. \verb{https://doi.org/10.1093/bioinformatics/btp358}

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}

   Viechtbauer, W., & \enc{López-López}{Lopez-Lopez}, J. A. (2022). Location-scale models for meta-analysis. \emph{Research Synthesis Methods}. \bold{13}(6), 697--715. \verb{https://doi.org/10.1002/jrsm.1562}
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mv}}, and \code{\link[=selmodel.rma.uni]{selmodel}} for functions to fit models for which profile likelihood plots can be drawn.

   \code{\link[=confint.rma]{confint}} for functions to compute corresponding profile likelihood confidence intervals.
}
\examples{
### calculate log odds ratios and corresponding sampling variances
dat <- escalc(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### fit random-effects model using rma.uni()
res <- rma(yi, vi, data=dat)

### profile over tau^2
profile(res, progbar=FALSE)

### adjust xlim
profile(res, progbar=FALSE, xlim=c(0,1))

### specify tau^2 values at which to profile the likelihood
profile(res, progbar=FALSE, steps=c(seq(0,0.2,length=20),seq(0.3,1,by=0.1)))

### change data into long format
dat.long <- to.long(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, append=FALSE)

### set levels of group variable ("exp" = experimental/vaccinated; "con" = control/non-vaccinated)
levels(dat.long$group) <- c("exp", "con")

### set "con" to reference level
dat.long$group <- relevel(dat.long$group, ref="con")

### calculate log odds and corresponding sampling variances
dat.long <- escalc(measure="PLO", xi=out1, mi=out2, data=dat.long)
dat.long

### fit bivariate random-effects model using rma.mv()
res <- rma.mv(yi, vi, mods = ~ group, random = ~ group | study, struct="UN", data=dat.long)
res

### profile over tau^2_1, tau^2_2, and rho
### note: for rho, adjust region over which profiling is done ('zoom in' on area around estimate)
\dontrun{
par(mfrow=c(2,2))
profile(res, tau2=1)
profile(res, tau2=2)
profile(res, rho=1, xlim=c(0.90, 0.98))
par(mfrow=c(1,1))
}

### an example where the peak of the likelihood profile is at 0
dat <- escalc(measure="RD", n1i=n1i, n2i=n2i, ai=ai, ci=ci, data=dat.hine1989)
res <- rma(yi, vi, data=dat)
profile(res, progbar=FALSE)
}
\keyword{hplot}
