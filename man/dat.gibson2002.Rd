\name{dat.gibson2002}
\docType{data}
\alias{dat.gibson2002}
\title{Studies on the Effectiveness of Self-Management Education and Regular Medical Review for Adults with Asthma}
\description{Results from 15 trials examining the effectiveness of self-management education and regular medical review for adults with asthma.}
\usage{dat.gibson2002}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{author}   \tab \code{character} \tab first author of study \cr
\bold{year}     \tab \code{numeric}   \tab publication year \cr
\bold{n1i}      \tab \code{numeric}   \tab number of participants in the intervention group \cr
\bold{m1i}      \tab \code{numeric}   \tab mean number of days off work/school in the intervention group \cr
\bold{sd1i}     \tab \code{numeric}   \tab standard deviation of the number of days off work/school in the intervention group \cr
\bold{n2i}      \tab \code{numeric}   \tab number of participants in the control/comparison group \cr
\bold{m2i}      \tab \code{numeric}   \tab mean number of days off work/school in the control/comparison group \cr
\bold{sd2i}     \tab \code{numeric}   \tab standard deviation of the number of days off work/school in the control/comparison group \cr
\bold{ai}       \tab \code{numeric}   \tab number of participants who had one or more days off work/school in the intervention group \cr
\bold{bi}       \tab \code{numeric}   \tab number of participants who no days off work/school in the intervention group \cr
\bold{ci}       \tab \code{numeric}   \tab number of participants who had one or more days off work/school in the control/comparison group \cr
\bold{di}       \tab \code{numeric}   \tab number of participants who no days off work/school in the control/comparison group \cr
\bold{type}     \tab \code{numeric}   \tab numeric code for the intervention type (see below)
}
}
\details{
   Asthma management guidelines typically recommend for patients to receive education and regular medical review. While self-management programs have been shown to increase patient knowledge, it is less clear to what extent they actually impact health outcomes. The systematic review by Gibson et al. (2002) examined the effectiveness of self-management education and regular medical review for adults with asthma. In each study, participants receiving a certain management intervention were compared against those in a control/comparison group with respect to a variety of health outcomes. One of the outcomes examined in a number of studies was the number of days off work/school.

   The majority of studies reporting this outcome provided means and standard deviations allowing a meta-analysis of standardized mean differences. Seven studies also reported the number of participants who had one or more days off work/school in each group. These studies could be meta-analyzed using, for example, (log) risk ratios. Finally, one could also consider a combined analysis based on standardized mean differences computed from the means and standard deviations where available and using probit transformed risk differences (which also provide estimates of the standardized mean difference) for the remaining studies.

   Some degree of patient education was provided in all studies. In addition, the \code{type} variable indicates what additional intervention components were included in each study:
   \enumerate{
   \item optimal self-management (writing action plan, self-monitoring, regular medical review),
   \item self-monitoring and regular medical review,
   \item self-monitoring only,
   \item regular medical review only,
   \item written action plan only.
   }
}
\source{
   Gibson, P. G., Powell, H., Wilson, A., Abramson, M. J., Haywood, P., Bauman, A., Hensley, M. J., Walters, E. H., & Roberts, J. J. L. (2002). Self-management education and regular practitioner review for adults with asthma. \emph{Cochrane Database of Systematic Reviews}, \bold{3}, CD001117. \url{https://doi.org/10.1002/14651858.CD001117}
}
\examples{
### copy data into 'dat' and examine data
dat <- dat.gibson2002
dat

### fixed-effects model analysis of the standardized mean differences
dat <- escalc(measure="SMD", m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i, data=dat)
res <- rma(yi, vi, data=dat, method="FE")
print(res, digits=2)

### fixed-effects model analysis of the (log) risk ratios
dat <- escalc(measure="RR", ai=ai, bi=bi, ci=ci, di=di, data=dat)
res <- rma(yi, vi, data=dat, method="FE")
print(res, digits=2)
predict(res, transf=exp, digits=2)

### fixed-effects model analysis of the standardized mean differences and the probit transformed
### risk differences (which also provide estimates of the standardized mean difference)
dat <- escalc(measure="SMD", m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i, data=dat)
dat <- escalc(measure="PBIT", ai=ai, bi=bi, ci=ci, di=di, data=dat, replace=FALSE)
dat
res <- rma(yi, vi, data=dat, method="FE")
print(res, digits=2)
}
\keyword{datasets}
