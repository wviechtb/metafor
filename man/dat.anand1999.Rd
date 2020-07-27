\name{dat.anand1999}
\docType{data}
\alias{dat.anand1999}
\title{Studies on the Effectiveness of Oral Anticoagulants in Patients with Coronary Artery Disease}
\description{Results from 34 trials examining the effectiveness of oral anticoagulants in patients with coronary artery disease.}
\usage{
dat.anand1999
}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{study}     \tab \code{character} \tab author(s) or trial name \cr
\bold{year}      \tab \code{numeric}   \tab publication year \cr
\bold{intensity} \tab \code{character} \tab intensity of anticoagulation (low, medium, or high) \cr
\bold{asp.t}     \tab \code{numeric}   \tab concomitant use of aspirin in the treatment group (0 = no, 1 = yes) \cr
\bold{asp.c}     \tab \code{numeric}   \tab concomitant use of aspirin in the control group (0 = no, 1 = yes) \cr
\bold{ai}        \tab \code{numeric}   \tab number of deaths in the treatment group \cr
\bold{n1i}       \tab \code{numeric}   \tab number of patients in the treatment group \cr
\bold{ci}        \tab \code{numeric}   \tab number of deaths in the control group \cr
\bold{n2i}       \tab \code{numeric}   \tab number of patients in the control group
}
}
\details{
   The dataset includes the results from 34 randomized clinical trials that examined the effectiveness of oral anticoagulants in patients with coronary artery disease. The results given here are focused on the total mortality in the treatment versus control groups.
}
\note{
   Strictly speaking, there are only 31 trials, since Breddin et al. (1980) and ATACS (1990) are multiarm trials.

   According to a correction, \code{dat.anand1999$ci[29]} should be 1. But then \code{dat.anand1999$ci[21]} would also have to be 1 (if these data indeed refer to the same control group). This appears contradictory, so this correction was not made.
}
\source{
   Anand, S. S., & Yusuf, S. (1999). Oral anticoagulant therapy in patients with coronary artery disease: A meta-analysis. \emph{Journal of the American Medical Association}, \bold{282}(21), 2058--2067. \verb{https://doi.org/10.1001/jama.282.21.2058}
}
\examples{
### copy data into 'dat' and examine data
dat <- dat.anand1999
dat

### High-Intensity OA vs Control
rma.mh(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat,
       subset=(intensity=="high" & asp.t==0 & asp.c==0), digits=2)

### High- or Moderate-Intensity OA vs Aspirin
rma.mh(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat,
       subset=(intensity \%in\% c("high","moderate") & asp.t==0 & asp.c==1), digits=2)

### Moderate-Intensity OA vs Control
rma.mh(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat,
       subset=(intensity=="moderate" & asp.t==0 & asp.c==0), digits=2)

### High- or Moderate-Intensity OA and Aspirin vs Aspirin
rma.mh(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat,
       subset=(intensity \%in\% c("high","moderate") & asp.t==1 & asp.c==1), digits=2)

### Low-Intensity OA and Aspirin vs Aspirin
rma.mh(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat,
       subset=(intensity=="low" & asp.t==1 & asp.c==1), digits=2)
}
\keyword{datasets}
