\name{dat.normand1999}
\docType{data}
\alias{dat.normand1999}
\title{Studies on the Length of Hospital Stay of Stroke Patients}
\description{Results from 9 studies on the length of the hospital stay of stroke patients under specialized care and under conventional/routine (non-specialist) care.}
\usage{dat.normand1999}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{study}    \tab \code{numeric}   \tab study number \cr
\bold{source}   \tab \code{character} \tab source of data \cr
\bold{n1i}      \tab \code{numeric}   \tab number of patients under specialized care \cr
\bold{m1i}      \tab \code{numeric}   \tab mean length of stay (in days) under specialized care \cr
\bold{sd1i}     \tab \code{numeric}   \tab standard deviation of the length of stay under specialized care \cr
\bold{n2i}      \tab \code{numeric}   \tab number of patients under routine care \cr
\bold{m2i}      \tab \code{numeric}   \tab mean length of stay (in days) under routine care \cr
\bold{sd2i}     \tab \code{numeric}   \tab standard deviation of the length of stay under routine care
}
}
\details{
   The 9 studies provide data in terms of the mean length of the hospital stay (in days) of stroke patients under specialized care and under conventional/routine (non-specialist) care. The goal of the meta-analysis was to examine the hypothesis whether specialist stroke unit care will result in a shorter length of hospitalization compared to routine management.
}
\source{
   Normand, S. T. (1999). Meta-analysis: Formulating, evaluating, combining, and reporting. \emph{Statistics in Medicine}, \bold{18}(3), 321--359. \url{https://doi.org/10.1002/(sici)1097-0258(19990215)18:3<321::aid-sim28>3.0.co;2-p}
}
\examples{
### copy data into 'dat'
dat <- dat.normand1999

### calculate mean differences and corresponding sampling variances
dat <- escalc(measure="MD", m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i, data=dat)
dat

### meta-analysis of mean differences using a random-effects model
res <- rma(yi, vi, data=dat)
res

### meta-analysis of standardized mean differences using a random-effects model
res <- rma(measure="SMD", m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i,
           data=dat, slab=source)
res

### draw forest plot
forest(res, xlim=c(-7,5), alim=c(-3,1), cex=.8, header="Study/Source")
}
\keyword{datasets}
