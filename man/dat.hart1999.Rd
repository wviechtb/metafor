\name{dat.hart1999}
\docType{data}
\alias{dat.hart1999}
\title{Studies on the Effectiveness of Warfarin for Preventing Strokes}
\description{Results from 6 clinical trials examining the effectiveness of adjusted-dose warfarin for preventing strokes in patients with atrial fibrillation.}
\usage{dat.hart1999}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{trial}    \tab \code{numeric}   \tab trial number \cr
\bold{study}    \tab \code{character} \tab study name (abbreviated) \cr
\bold{year}     \tab \code{numeric}   \tab publication year \cr
\bold{x1i}      \tab \code{numeric}   \tab number of strokes in the warfarin group \cr
\bold{n1i}      \tab \code{numeric}   \tab number of patients in the warfarin group \cr
\bold{t1i}      \tab \code{numeric}   \tab total person-time (in years) in the warfarin group \cr
\bold{x2i}      \tab \code{numeric}   \tab number of strokes in the placebo/control group \cr
\bold{n2i}      \tab \code{numeric}   \tab number of patients in the placebo/control group \cr
\bold{t2i}      \tab \code{numeric}   \tab total person-time (in years) in the placebo/control group \cr
\bold{compgrp}  \tab \code{character} \tab type of comparison group (placebo or control) \cr
\bold{prevtype} \tab \code{character} \tab type of prevention (primary or secondary) \cr
\bold{trinr}    \tab \code{character} \tab target range for the international normalized ratio (INR)
}
}
\details{
   The 6 studies provide data with respect to the number of strokes in the warfarin and the comparison (placebo or control) group. In addition, the number of patients and the total person-time (in years) is provided for the two groups. The goal of the meta-analysis was to examine the effectiveness of adjusted-dose warfarin for preventing strokes in patients with atrial fibrillation.
}
\source{
   Hart, R. G., Benavente, O., McBride, R., & Pearce, L. A. (1999). Antithrombotic therapy to prevent stroke in patients with atrial fibrillation: A meta-analysis. \emph{Annals of Internal Medicine}, \bold{131}(7), 492--501. \url{https://doi.org/10.7326/0003-4819-131-7-199910050-00003}
}
\examples{
### copy data into 'dat'
dat <- dat.hart1999

### calculate log incidence rate ratios and corresponding sampling variances
dat <- escalc(measure="IRR", x1i=x1i, x2i=x2i, t1i=t1i, t2i=t2i, data=dat)
dat

### meta-analysis of log incidence rate ratios using a random-effects model
res <- rma(yi, vi, data=dat)
res

### average incidence rate ratio with 95\% CI
predict(res, transf=exp)

### forest plot with extra annotations
par(mar=c(5,4,1,2))
forest(res, xlim=c(-11, 5), at=log(c(.05, .25, 1, 4)), atransf=exp,
       slab=paste0(dat$study, " (", dat$year, ")"),
       ilab=cbind(paste(dat$x1i, "/", dat$t1i, sep=" "),
       paste(dat$x2i, "/", dat$t2i, sep=" ")),
       ilab.xpos=c(-6.5,-4), cex=.85, header="Study (Year)")
op <- par(cex=.85, font=2)
text(c(-6.5,-4), 8.5, c("Warfarin", "Control"))
text(c(-6.5,-4), 7.5, c("Strokes / PT", "Strokes / PT"))
segments(x0=-8, y0=8, x1=-2.75, y1=8)
par(op)

### meta-analysis of incidence rate differences using a random-effects model
res <- rma(measure="IRD", x1i=x1i, x2i=x2i, t1i=t1i, t2i=t2i, data=dat)
res
}
\keyword{datasets}
