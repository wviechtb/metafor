\name{dat.hine1989}
\docType{data}
\alias{dat.hine1989}
\title{Studies on Prophylactic Use of Lidocaine After a Heart Attack}
\description{Results from 6 studies evaluating mortality from prophylactic use of lidocaine in acute myocardial infarction.}
\usage{dat.hine1989}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{study}  \tab \code{numeric}   \tab study number \cr
\bold{source} \tab \code{character} \tab source of data \cr
\bold{n1i}    \tab \code{numeric}   \tab number of patients in lidocaine group \cr
\bold{n2i}    \tab \code{numeric}   \tab number of patients in control group \cr
\bold{ai}     \tab \code{numeric}   \tab number of deaths in lidocaine group \cr
\bold{ci}     \tab \code{numeric}   \tab number of deaths in control group \cr
}
}
\details{
   Hine et al. (1989) conducted a meta-analysis of death rates in randomized controlled trials in which prophylactic lidocaine was administered to patients with confirmed or suspected acute myocardial infarction. The dataset describes the mortality at the end of the assigned treatment period for control and intravenous lidocaine treatment groups for six studies. The question of interest is whether there is a detrimental effect of lidocaine. Because the studies were conducted to compare rates of arrhythmias following a heart attack, the studies, taken individually, are too small to detect important differences in mortality rates.

   The data in this dataset were obtained from Table I in Normand (1999, p. 322).
}
\source{
   Normand, S. T. (1999). Meta-analysis: Formulating, evaluating, combining, and reporting. \emph{Statistics in Medicine}, \bold{18}, 321--359.
}
\references{
   Hine, L. K., Laird, N., Hewitt, P., & Chalmers, T. C. (1989). Meta-analytic evidence against prophylactic use of lidocaine in acute myocardial infarction. \emph{Archives of Internal Medicine}, \bold{149}, 2694--2698.
}
\examples{
### load data
dat <- get(data(dat.hine1989))

### calculate risk differences and corresponding sampling variances
dat <- escalc(measure="RD", n1i=n1i, n2i=n2i, ai=ai, ci=ci, data=dat)
dat

### meta-analysis of risk differences using a random-effects model
res <- rma(yi, vi, data=dat)
res
}
\keyword{datasets}
