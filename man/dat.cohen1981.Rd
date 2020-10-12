\name{dat.cohen1981}
\docType{data}
\alias{dat.cohen1981}
\title{Studies on the Relationship between Course Ratings and Student Achievement}
\description{Results from 20 studies on the correlation between overall course ratings and student achievement.}
\usage{dat.cohen1981}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{study}     \tab \code{character} \tab study author(s) and year \cr
\bold{sample}    \tab \code{character} \tab course type \cr
\bold{control}   \tab \code{character} \tab ability control \cr
\bold{ni}        \tab \code{numeric}   \tab sample size of the study (number of sections) \cr
\bold{ri}        \tab \code{numeric}   \tab observed correlation
}
}
\details{
   The studies included in this dataset examined to what extent students' ratings of a course correlated with their achievement in the course. Instead of correlating individual ratings and achievement scores, the studies were carried out in multisection courses, in which the sections had different instructors but all sections used a common achievement measure (e.g., a final exam). The correlation coefficients reflects the correlation between the mean overall course rating of each section and the mean achievement score.

   Note that this dataset (extracted from Table A.3 in Cooper & Hedges, 1994) only contains studies with at least 10 sections.
}
\source{
   Cooper, H., & Hedges, L. V. (1994). Appendix A: Data Sets. In H. Cooper & L. V. Hedges (Eds.), \emph{The handbook of research synthesis} (pp. 543-547). New York: Russell Sage Foundation.
}
\references{
   Cohen, P. A. (1981). Student ratings of instruction and student achievement: A meta-analysis of multisection validity studies. \emph{Review of Educational Research}, \bold{51}(3), 281--309. \verb{https://doi.org/10.3102/00346543051003281}
}
\examples{
### copy data into 'dat'
dat <- dat.cohen1981

### calculate r-to-z transformed correlations and corresponding sampling variances
dat <- escalc(measure="ZCOR", ri=ri, ni=ni, data=dat[c(1,4,5)])
dat

### meta-analysis of the transformed correlations using a random-effects model
res <- rma(yi, vi, data=dat, digits=2)
res

### predicted average correlation with 95\% CI
predict(res, transf=transf.ztor)
}
\keyword{datasets}
