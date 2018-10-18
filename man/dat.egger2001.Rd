\name{dat.egger2001}
\docType{data}
\alias{dat.egger2001}
\title{Studies on the Effectiveness of Intravenous Magnesium in Acute Myocardial Infarction}
\description{Results from 16 trials examining the effectiveness of intravenous magnesium in the prevention of death following acute myocardial infarction.}
\usage{
dat.egger2001
}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{id}    \tab \code{numeric}   \tab trial id number \cr
\bold{study} \tab \code{character} \tab first author or trial name \cr
\bold{year}  \tab \code{numeric}   \tab publication year \cr
\bold{ai}    \tab \code{numeric}   \tab number of deaths in the magnesium group \cr
\bold{n1i}   \tab \code{numeric}   \tab number of patients in the magnesium group \cr
\bold{ci}    \tab \code{numeric}   \tab number of deaths in the control group \cr
\bold{n2i}   \tab \code{numeric}   \tab number of patients in the control group
}
}
\details{
   The dataset includes the results from 16 randomized clinical trials that examined the effectiveness of intravenous magnesium in the prevention of death following acute myocardial infarction. Studies 1-7 were included in the meta-analyses by Teo et al. (1991) and Horner (1992) and were combined with the results from the LIMIT-2 trial (Woods et al., 1992) in Yusuf et al. (1993), suggesting that magnesium is an effective treatment for reducing mortality. However, the results from the ISIS-4 mega trial (ISIS-4 Collaborative Group, 1995) indicated no reduction in mortality with magnesium treatment. Publication bias has been suggested as one possible explanation for the conflicting findings (Egger & Davey Smith, 1995).

   The present dataset includes some additional trials and are based on Table 18.2 from Egger, Davey Smith, and Altman (2001).
}
\source{
   Egger, M., Davey Smith, G., & Altman, D. G. (Eds.) (2001). \emph{Systematic reviews in health care: Meta-analysis in context} (2nd ed.). London: BMJ Books.
}
\references{
   Egger, M., & Davey Smith, G. (1995). Misleading meta-analysis: Lessons from \dQuote{an effective, safe, simple} intervention that wasn't. \emph{British Medical Journal}, \bold{310}, 752--754.

   Horner, S. M. (1992). Efficacy of intravenous magnesium in acute myocardial infarction in reducing arrhythmias and mortality: Meta-analysis of magnesium in acute myocardial infarction. \emph{Circulation}, \bold{86}, 774--779.

   ISIS-4 Collaborative Group (1995). ISIS-4: A randomised factorial trial assessing early oral captopril, oral mononitrate, and intravenous magnesium sulphate in 58,050 patients with suspected acute myocardial infarction. \emph{Lancet}, \bold{345}, 669--685.

   Teo, K. K., Yusuf, S., Collins, R., Held, P. H., & Peto, R. (1991). Effects of intravenous magnesium in suspected acute myocardial infarction: Overview of randomised trials. \emph{British Medical Journal}, 303, 1499--1503.

   Woods, K. L., Fletcher, S., Roffe, C., & Haider, Y. (1992). Intravenous magnesium sulphate in suspected acute myocardial infarction: Results of the second Leicester Intravenous Magnesium Intervention Trial (LIMIT-2). \emph{Lancet}, \bold{339}, 1553--1558.

   Yusuf, S., Teo, K., & Woods, K. (1993). Intravenous magnesium in acute myocardial infarction: An effective, safe, simple, and inexpensive treatment. \emph{Circulation}, \bold{87}, 2043--2046.
}
\seealso{
   \code{\link{dat.li2007}}
}
\examples{
### load data
dat <- get(data(dat.egger2001))

### meta-analysis of trials 1-7 using Peto's method (as in Teo et al., 1991)
res <- rma.peto(ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, subset=1:7)
print(res, digits=2)

### meta-analysis of trials 1-7 and LIMIT-2 (as in Yusuf et al., 1993)
res <- rma.peto(ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, subset=c(1:7,14))
print(res, digits=2)

### meta-analysis of all trials except ISIS-4
res <- rma.peto(ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, subset=-16)
print(res, digits=2)
predict(res, transf=exp, digits=2)

### meta-analysis of all trials including ISIS-4
res <- rma.peto(ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat)
print(res, digits=2)
predict(res, transf=exp, digits=2)

### contour-enhanced funnel plot centered at 0
funnel(res, refline=0, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"))
}
\keyword{datasets}
