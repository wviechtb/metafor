\name{dat.pritz1997}
\docType{data}
\alias{dat.pritz1997}
\title{Studies on the Effectiveness of Hyperdynamic Therapy for Treating Cerebral Vasospasm}
\description{Results from 14 studies on the effectiveness of hyperdynamic therapy for treating cerebral vasospasm.}
\usage{dat.pritz1997}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{study}    \tab \code{numeric}   \tab study number \cr
\bold{authors}  \tab \code{character} \tab study authors \cr
\bold{xi}       \tab \code{numeric}   \tab number of patients that improved with hyperdynamic therapy \cr
\bold{ni}       \tab \code{numeric}   \tab total number of patients treated \cr
}
}
\details{
   As described in Zhou et al. (1999), "hyperdynamic therapy refers to induced hypertension and hypervolaemia (volume expansion) to treat ischaemic symptoms due to vasospasm, and the success of this therapy is defined as clinical improvement in terms of neurologic deficits." For each study that was included in the meta-analysis, the dataset includes information on the number of patients that improved under this form of therapy and the total number of patients that were treated. The goal of the meta-analysis is to estimate the true (average) success rate of hyperdynamic therapy.
}
\source{
   Zhou, X.-H., Brizendine, E. J., & Pritz, M. B. (1999). Methods for combining rates from several studies. \emph{Statistics in Medicine}, \bold{18}, 557--566.
}
\references{
   Pritz, M. B. (1997). Treatment of cerebral vasospasm due to aneurysmal subarachnoid hemorrhage: Past, present, and future of hyperdynamic therapy. \emph{Neurosurgery Quarterly}, \bold{7}, 273--285.
}
\examples{
### load data
dat <- get(data(dat.pritz1997))

### computation of "weighted average" in Zhou et al. (1999), Table IV
dat <- escalc(measure="PR", xi=xi, ni=ni, data=dat, add=0)
theta.hat    <- sum(dat$ni * dat$yi) / sum(dat$ni)
se.theta.hat <- sqrt(sum(dat$ni^2 * dat$vi) / sum(dat$ni)^2)
ci.lb        <- theta.hat - 1.96*se.theta.hat
ci.ub        <- theta.hat + 1.96*se.theta.hat
round(c(estimate = theta.hat, se = se.theta.hat, ci.lb = ci.lb, ci.ub = ci.ub), 4)

### this is identical to a FE model with sample size weights
rma(yi, vi, weights=ni, method="FE", data=dat)

### random-effects model with raw proportions
dat <- escalc(measure="PR", xi=xi, ni=ni, data=dat)
res <- rma(yi, vi, data=dat)
predict(res)

### random-effects model with logit transformed proportions
dat <- escalc(measure="PLO", xi=xi, ni=ni, data=dat)
res <- rma(yi, vi, data=dat)
predict(res, transf=transf.ilogit)

### mixed-effects logistic regression model
res <- rma.glmm(measure="PLO", xi=xi, ni=ni, data=dat)
predict(res, transf=transf.ilogit)
}
\keyword{datasets}
