\name{dat.bornmann2007}
\docType{data}
\alias{dat.bornmann2007}
\title{Studies on Gender Differences in Grant and Fellowship Awards}
\description{Results from 21 studies on gender differences in grant and fellowship awards.}
\usage{dat.bornmann2007}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{study}      \tab \code{character} \tab study reference \cr
\bold{obs}        \tab \code{numeric}   \tab observation within study \cr
\bold{doctype}    \tab \code{character} \tab document type \cr
\bold{gender}     \tab \code{character} \tab gender of the study authors \cr
\bold{year}       \tab \code{numeric}   \tab (average) cohort year \cr
\bold{org}        \tab \code{character} \tab funding organization / program \cr
\bold{country}    \tab \code{character} \tab country of the funding organization / program \cr
\bold{type}       \tab \code{character} \tab fellowship or grant application \cr
\bold{discipline} \tab \code{character} \tab discipline / field \cr
\bold{waward}     \tab \code{numeric}   \tab number of women who received a grant/fellowship award \cr
\bold{wtotal}     \tab \code{numeric}   \tab number of women who applied for an award \cr
\bold{maward}     \tab \code{numeric}   \tab number of men who received a grant/fellowship award \cr
\bold{mtotal}     \tab \code{numeric}   \tab number of men who applied for an award
}
}
\details{
   The studies in this dataset examine whether the chances of receiving a grant or fellowship award differs for men and women. Note that many studies provide multiple comparisons (e.g., for different years / cohorts / disciplines). A multilevel meta-analytic model can be used to account for the multilevel structure in these data.
}
\source{
   Bornmann, L., Mutz, R., & Daniel, H. (2007). Gender differences in grant peer review: A meta-analysis. \emph{Journal of Informetrics}, \bold{1}(3), 226--238. \verb{https://doi.org/10.1016/j.joi.2007.03.001}
}
\references{
   Marsh, H. W., Bornmann, L., Mutz, R., Daniel, H.-D., & O'Mara, A. (2009). Gender effects in the peer reviews of grant proposals: A comprehensive meta-analysis comparing traditional and multilevel approaches. \emph{Review of Educational Research}, \bold{79}(3), 1290--1326. \verb{https://doi.org/10.3102/0034654309334143}
}
\examples{
### copy data into 'dat' and examine data
dat <- dat.bornmann2007
dat

### calculate log odds ratios and corresponding sampling variances
dat <- escalc(measure="OR", ai=waward, n1i=wtotal, ci=maward, n2i=mtotal, data=dat)

### fit multilevel meta-analytic model
res <- rma.mv(yi, vi, random = ~ 1 | study/obs, data=dat)
res

### estimated average odds ratio (with 95% CI/PI)
predict(res, transf=exp, digits=2)

### test for a difference between fellowship and grant applications
res <- rma.mv(yi, vi, mods = ~ type, random = ~ 1 | study/obs, data=dat)
res
predict(res, newmods=0:1, transf=exp, digits=2)
}
\keyword{datasets}
