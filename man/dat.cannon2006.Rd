\name{dat.cannon2006}
\docType{data}
\alias{dat.cannon2006}
\title{Studies on the Effectiveness of Intensive Versus Moderate Statin Therapy for Preventing Coronary Death or Myocardial Infarction}
\description{Results from 4 trials examining the effectiveness of intensive (high dose) versus moderate (standard dose) statin therapy for preventing coronary death or myocardial infarction.}
\usage{dat.cannon2006}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{trial} \tab \code{character} \tab trial name \cr
\bold{pop}   \tab \code{character} \tab study population (post-ACS: post acute coronary syndrome; stable CAD: stable coronary artery disease) \cr
\bold{nt}    \tab \code{numeric}   \tab number of patients in the high dose group \cr
\bold{nc}    \tab \code{numeric}   \tab number of patients in the standard dose group \cr
\bold{ep1t}  \tab \code{numeric}   \tab number of events in the high dose group for end point 1: coronary death or non-fatal myocardial infarction \cr
\bold{ep1c}  \tab \code{numeric}   \tab number of events in the standard dose group for end point 1: coronary death or non-fatal myocardial infarction \cr
\bold{ep2t}  \tab \code{numeric}   \tab number of events in the high dose group for end point 2: coronary death or any cardiovascular event (MI, stroke, hospitalization for unstable angina, or revascularization) \cr
\bold{ep2c}  \tab \code{numeric}   \tab number of events in the standard dose group for end point 2: coronary death or any cardiovascular event (MI, stroke, hospitalization for unstable angina, or revascularization) \cr
\bold{ep3t}  \tab \code{numeric}   \tab number of events in the high dose group for end point 3: cardiovascular death \cr
\bold{ep3c}  \tab \code{numeric}   \tab number of events in the standard dose group for end point 3: cardiovascular death \cr
\bold{ep4t}  \tab \code{numeric}   \tab number of events in the high dose group for end point 4: non-cardiovascular death \cr
\bold{ep4c}  \tab \code{numeric}   \tab number of events in the standard dose group for end point 4: non-cardiovascular death \cr
\bold{ep5t}  \tab \code{numeric}   \tab number of events in the high dose group for end point 5: deaths (all-cause mortality) \cr
\bold{ep5c}  \tab \code{numeric}   \tab number of events in the standard dose group for end point 5: deaths (all-cause mortality) \cr
\bold{ep6t}  \tab \code{numeric}   \tab number of events in the high dose group for end point 6: stroke \cr
\bold{ep6c}  \tab \code{numeric}   \tab number of events in the standard dose group for end point 6: stroke
}
}
\details{
   The data were obtained from Figures 2, 3, 4, and 5 in Cannon et al. (2006). The authors used the Mantel-Haenszel method for combining the results from the 4 trials. This approach is implemented in the \code{\link{rma.mh}} function.
}
\source{
   Cannon, C. P., Steinberg, B. A., Murphy, S. A., Mega, J. L., & Braunwald, E. (2006). Meta-analysis of cardiovascular outcomes trials comparing intensive versus moderate statin therapy. \emph{Journal of the American College of Cardiology}, \bold{48}(3), 438--445. \verb{https://doi.org/10.1016/j.jacc.2006.04.070}
}
\examples{
### copy data into 'dat' and examine data
dat <- dat.cannon2006
dat

### meta-analysis of log odds ratios using the MH method for endpoint 1
res <- rma.mh(measure="OR", ai=ep1t, n1i=nt, ci=ep1c, n2i=nc, data=dat, slab=trial)
print(res, digits=2)

### forest plot
forest(res, xlim=c(-.8,.8), atransf=exp, at=log(c(2/3, 1, 3/2)),
       header=TRUE, top=2, cex=1.2, xlab="Odds Ratio")
mtext("(high dose better)", side=1, line=par("mgp")[1]-0.5, at=log(2/3), cex=1.2, font=3)
mtext("(standard dose better)", side=1, line=par("mgp")[1]-0.5, at=log(3/2), cex=1.2, font=3)
}
\keyword{datasets}
