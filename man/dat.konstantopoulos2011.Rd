\name{dat.konstantopoulos2011}
\docType{data}
\alias{dat.konstantopoulos2011}
\title{Studies on the Effects of Modified School Calendars on Student Achievement}
\description{Results from 56 studies on the effects of modified school calendars on student achievement.}
\usage{dat.konstantopoulos2011}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{district} \tab \code{numeric}   \tab district id number \cr
\bold{school}   \tab \code{numeric}   \tab school id number (within district) \cr
\bold{study}    \tab \code{numeric}   \tab study id number \cr
\bold{yi}       \tab \code{numeric}   \tab standardized mean difference \cr
\bold{vi}       \tab \code{numeric}   \tab corresponding sampling variance \cr
\bold{year}     \tab \code{numeric}   \tab year of the study
}
}
\details{
   Instead of following the more traditional school calendar with a long summer break (in addition to a short winter and spring break), some schools have switched to a modified school calendar comprising more frequent but shorter intermittent breaks (e.g., 9 weeks of school followed by 3 weeks off), while keeping the total number of days at school approximately the same. The effects of using such a modified calendar on student achievement have been examined in a number of studies and were meta-analyzed by Cooper et al. (2003).

   The dataset (taken from Konstantopoulos, 2011) contains the results from 56 studies, each comparing the level of academic achievement in a group of students following a modified school calendar with that of a group of students following a more traditional school calendar. The difference between the two groups was quantified in terms of a standardized mean difference (with positive values indicating a higher mean level of achievement in the group following the modified school calendar).
   
   The studies were conducted at various schools that were clustered within districts. The data therefore have a multilevel structure, with schools nested within districts. A multilevel meta-analysis of these data can be used to estimate and account for the amount of heterogeneity between districts and between schools within districts.
}
\source{
   Konstantopoulos, S. (2011). Fixed effects and variance components estimation in three-level meta-analysis. \emph{Research Synthesis Methods}, \bold{2}, 61--76.
}
\references{
   Cooper, H., Valentine, J. C., Charlton, K., & Melson, A. (2003). The effects of modified school calendars on student achievement and on school and community attitudes. \emph{Review of Educational Research}, \bold{73}, 1--52.
}
\examples{
### load data
dat <- get(data(dat.konstantopoulos2011))

### regular random-effects model
res <- rma(yi, vi, data=dat)
print(res, digits=3)

### regular random-effects model using rma.mv()
res <- rma.mv(yi, vi, random = ~ 1 | study, data=dat)
print(res, digits=3)

### multilevel random-effects model
res.ml <- rma.mv(yi, vi, random = ~ 1 | district/school, data=dat)
print(res.ml, digits=3)

### profile variance components
par(mfrow=c(2,1))
profile(res.ml, sigma2=1)
profile(res.ml, sigma2=2)

### multivariate parameterization of the model
res.mv <- rma.mv(yi, vi, random = ~ factor(school) | district, data=dat)
print(res.mv, digits=3)

### tau^2 from multivariate model = sum of the two variance components from the multilevel model
round(sum(res.ml$sigma2), 3)

### rho from multivariate model = intraclass correlation coefficient based on the multilevel model
round(res.ml$sigma2[1] / sum(res.ml$sigma2), 3)
}
\keyword{datasets}
