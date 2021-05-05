\name{dat.landenberger2005}
\docType{data}
\alias{dat.landenberger2005}
\title{Studies on the Effectiveness of CBT for Reducing Recidivism}
\description{Results from 58 studies on the effectiveness of cognitive-behavioral therapy (CBT) for reducing recidivism in juvenile and adult offenders.}
\usage{dat.landenberger2005}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{study}            \tab \code{character} \tab (first) author and year \cr
\bold{pubtype}          \tab \code{character} \tab publication type (book chapter, journal article, report, or thesis) \cr
\bold{country}          \tab \code{character} \tab country where study was carried out (Canada, New Zealand, UK, or USA) \cr
\bold{design}           \tab \code{character} \tab study design (matched groups, nonequivalent groups, or randomized trial) \cr
\bold{program}          \tab \code{character} \tab purpose of setting up the CBT program (for demonstration, practice, or research purposes) \cr
\bold{setting}          \tab \code{character} \tab treatment setting (community or prison) \cr
\bold{designprob}       \tab \code{character} \tab indication of study design problems (no, favors the control group, or favors the treatment group) \cr
\bold{n.ctrl.rec}       \tab \code{numeric}   \tab number of recidivists in the control group \cr
\bold{n.ctrl.non}       \tab \code{numeric}   \tab number of non-recidivists in the control group \cr
\bold{n.cbt.rec}        \tab \code{numeric}   \tab number of recidivists in the CBT group \cr
\bold{n.cbt.non}        \tab \code{numeric}   \tab number of non-recidivists in the CBT group \cr
\bold{interval}         \tab \code{numeric}   \tab recidivism interval (in months) \cr
\bold{group}            \tab \code{numeric}   \tab study group (adults or juveniles) \cr
\bold{age}              \tab \code{numeric}   \tab mean age of the study group \cr
\bold{male}             \tab \code{numeric}   \tab percentage of males in the study group \cr
\bold{minority}         \tab \code{numeric}   \tab percentage of minorities in the study group \cr
\bold{length}           \tab \code{numeric}   \tab treatment length (in weeks) \cr
\bold{sessions}         \tab \code{numeric}   \tab number of CBT sessions per week \cr
\bold{hrs_week}         \tab \code{numeric}   \tab treatment hours per week \cr
\bold{hrs_total}        \tab \code{numeric}   \tab total hours of treatment \cr
\bold{cbt.cogskills}    \tab \code{character} \tab CBT component: cognitive skills (yes, no) \cr
\bold{cbt.cogrestruct}  \tab \code{character} \tab CBT component: cognitive restructuring (yes, no) \cr
\bold{cbt.intpprbsolv}  \tab \code{character} \tab CBT component: interpersonal problem solving (yes, no) \cr
\bold{cbt.socskills}    \tab \code{character} \tab CBT component: social skills (yes, no) \cr
\bold{cbt.angerctrl}    \tab \code{character} \tab CBT component: anger control (yes, no) \cr
\bold{cbt.victimimpact} \tab \code{character} \tab CBT component: victim impact (yes, no) \cr
\bold{cbt.subabuse}     \tab \code{character} \tab CBT component: substance abuse (yes, no) \cr
\bold{cbt.behavmod}     \tab \code{character} \tab CBT component: behavior modification (yes, no) \cr
\bold{cbt.relapseprev}  \tab \code{character} \tab CBT component: relapse prevention (yes, no) \cr
\bold{cbt.moralrsng}    \tab \code{character} \tab CBT component: moral reasoning (yes, no) \cr
\bold{cbt.roletaking}   \tab \code{character} \tab CBT component: role taking (yes, no) \cr
\bold{cbt.other}        \tab \code{character} \tab CBT component: other (yes, no)
}
}
\details{
   Landenberger and Lipsey (2005) conducted a meta-analysis of 58 experimental and quasi-experimental studies of the effects of cognitive-behavioral therapy (CBT) on the recidivism rates of adult and juvenile offenders (see also Lipsey et al., 2007). The present dataset includes the results of these studies and a range of potential moderator variables to identify factors associated with variation in treatment effects.
}
\source{
   Personal communication.
}
\references{
   Landenberger, N. A., & Lipsey, M. W. (2005). The positive effects of cognitive-behavioral programs for offenders: A meta-analysis of factors associated with effective treatment. \emph{Journal of Experimental Criminology}, \bold{1}, 451--476. \verb{https://doi.org/10.1007/s11292-005-3541-7}

   Lipsey, M. W., Landenberger, N. A., & Wilson, S. J. (2007). Effects of cognitive-behavioral programs for criminal offenders. \emph{Campbell Systematic Reviews}, \bold{3}(1), 1--27. \verb{https://doi.org/10.4073/csr.2007.6}
}
\examples{
### copy data into 'dat' and examine data
dat <- dat.landenberger2005
head(dat)

### calculate log odds ratios (for non-recidivism in CBT vs. control groups) and sampling variances
dat <- escalc(measure="OR", ai=n.cbt.non, bi=n.cbt.rec, ci=n.ctrl.non, di=n.ctrl.rec, data=dat)

### fit random-effects model
res <- rma(yi, vi, data=dat)
res

### estimated average OR and corresponding 95% CI/PI
predict(res, transf=exp, digits=2)

### examine if number of treatment sessions per week is a potential moderator
res <- rma(yi, vi, mods = ~ sessions, data=dat)
res

### predicted ORs for 1, 2, 5, or 10 sessions per week
predict(res, newmods=c(1,2,5,10), transf=exp, digits=2)
}
\keyword{datasets}
