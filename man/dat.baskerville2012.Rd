\name{dat.baskerville2012}
\docType{data}
\alias{dat.baskerville2012}
\title{Studies on the Effectiveness of Practice Facilitation Interventions}
\description{Results from 23 studies on the effectiveness of practice facilitation interventions within the primary care practice setting.}
\usage{dat.baskerville2012}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{author}    \tab \code{character} \tab study author(s) \cr
\bold{year}      \tab \code{numeric}   \tab publication year \cr
\bold{score}     \tab \code{numeric}   \tab quality score (0 to 12 scale) \cr
\bold{design}    \tab \code{character} \tab study design (cct = controlled clinical trial, rct = randomized clinical trial, crct = cluster randomized clinical trial) \cr
\bold{alloconc}  \tab \code{numeric}   \tab allocation concealed (0 = no, 1 = yes) \cr
\bold{blind}     \tab \code{numeric}   \tab single- or double-blind study (0 = no, 1 = yes) \cr
\bold{itt}       \tab \code{numeric}   \tab intention to treat analysis (0 = no, 1 = yes) \cr
\bold{fumonths}  \tab \code{numeric}   \tab follow-up months \cr
\bold{retention} \tab \code{numeric}   \tab retention (in percent) \cr
\bold{country}   \tab \code{character} \tab country where study was conducted \cr
\bold{outcomes}  \tab \code{numeric}   \tab number of outcomes assessed \cr
\bold{duration}  \tab \code{numeric}   \tab duration of intervention \cr
\bold{pperf}     \tab \code{numeric}   \tab practices per facilitator \cr
\bold{meetings}  \tab \code{numeric}   \tab (average) number of meetings \cr
\bold{hours}     \tab \code{numeric}   \tab (average) hours per meeting \cr
\bold{tailor}    \tab \code{numeric}   \tab intervention tailored to the context and needs of the practice (0 = no, 1 = yes) \cr
\bold{smd}       \tab \code{numeric}   \tab standardized mean difference \cr
\bold{se}        \tab \code{numeric}   \tab corresponding standard error
}
}
\details{
   "Outreach or practice facilitation is a multifaceted approach that involves skilled individuals who enable others, through a range of intervention components and approaches, to address the challenges in implementing evidence-based care guidelines within the primary care setting." (Baskerville et al., 2012). The studies included in this dataset examined the effectiveness of such practice facilitation interventions for improving some relevant evidence-based practice behavior. The effect was quantified in terms of a standardized mean difference, comparing the change (from pre- to post-intervention) in the intervention versus the comparison group (or the difference from baseline in prospective cohort studies).
}
\source{
   Baskerville, N. B., Liddy, C., & Hogg, W. (2012). Systematic review and meta-analysis of practice facilitation within primary care settings. \emph{Annals of Family Medicine}, \bold{10}(1), 63--74. \verb{https://doi.org/10.1370/afm.1312}
}
\examples{
### copy data into 'dat' and examine data
dat <- dat.baskerville2012
dat

### random-effects model
res <- rma(smd, sei=se, data=dat, method="DL")
print(res, digits=2)

### funnel plot
funnel(res, xlab="Standardized Mean Difference", ylim=c(0,0.6))

### rank and regression tests for funnel plot asymmetry
ranktest(res)
regtest(res)

### meta-regression analyses examining various potential moderators
rma(smd, sei=se, mods = ~ score, data=dat, method="DL")
rma(smd, sei=se, mods = ~ alloconc, data=dat, method="DL")
rma(smd, sei=se, mods = ~ blind,    data=dat, method="DL")
rma(smd, sei=se, mods = ~ itt,      data=dat, method="DL")
rma(smd, sei=se, mods = ~ duration, data=dat, method="DL")
rma(smd, sei=se, mods = ~ tailor,   data=dat, method="DL")
rma(smd, sei=se, mods = ~ pperf,    data=dat, method="DL")
rma(smd, sei=se, mods = ~ I(meetings * hours), data=dat, method="DL")
}
\keyword{datasets}
