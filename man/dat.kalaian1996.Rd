\name{dat.kalaian1996}
\docType{data}
\alias{dat.kalaian1996}
\title{Studies on the Effectiveness of Coaching for the SAT}
\description{Results from studies examining the effectiveness of coaching on the performance on the Scholastic Aptitude Test (SAT).}
\usage{dat.kalaian1996}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{id}          \tab \code{numeric}   \tab row (effect) id \cr
\bold{study}       \tab \code{character} \tab study identifier \cr
\bold{year}        \tab \code{numeric}   \tab publication year \cr
\bold{n1i}         \tab \code{numeric}   \tab number of participants in the coached group \cr
\bold{n2i}         \tab \code{numeric}   \tab number of participants in the uncoached group \cr
\bold{outcome}     \tab \code{character} \tab subtest (verbal or math) \cr
\bold{yi}          \tab \code{numeric}   \tab standardized mean difference \cr
\bold{vi}          \tab \code{numeric}   \tab corresponding sampling variance \cr
\bold{hrs}         \tab \code{numeric}   \tab hours of coaching \cr
\bold{ets}         \tab \code{numeric}   \tab study conducted by the Educational Testing Service (ETS) (0 = no, 1 = yes) \cr
\bold{homework}    \tab \code{numeric}   \tab assignment of homework outside of the coaching course (0 = no, 1 = yes) \cr
\bold{type}        \tab \code{numeric}   \tab study type (1 = randomized study, 2 = matched study, 3 = nonequivalent comparison study)
}
}
\details{
   The effectiveness of coaching for the Scholastic Aptitude Test (SAT) has been examined in numerous studies. This dataset contains standardized mean differences comparing the performance of a coached versus uncoached group on the verbal and/or math subtest of the SAT. Studies may report a standardized mean difference for the verbal subtest, the math subtest, or both. In the latter case, the two standardized mean differences are not independent (since they were measured in the same group of subjects). The number of hours of coaching (variable \code{hrs}), whether the study was conducted by the Educational Testing Service (variable \code{ets}), whether homework was assigned outside of the coaching course (variable \code{homework}), and the study type (variable \code{type}) may be potential moderators of the treatment effect.
}
\note{
   The dataset was obtained from Table 1 in Kalaian and Raudenbush (1996). However, there appear to be some inconsistencies between the data in the table and those that were actually used for the analyses (see \sQuote{Examples}).
}
\source{
   Kalaian, H. A., & Raudenbush, S. W. (1996). A multivariate mixed linear model for meta-analysis. \emph{Psychological Methods}, \bold{1}(3), 227--235. \verb{https://doi.org/10.1037/1082-989X.1.3.227}
}
\examples{
### copy data into 'dat'
dat <- dat.kalaian1996

### check ranges
range(dat$yi[dat$outcome == "verbal"]) # -0.35 to 0.74 according to page 230
range(dat$yi[dat$outcome == "math"])   # -0.53 to 0.60 according to page 231

### comparing this with Figure 1 in the paper reveals some discrepancies
par(mfrow=c(1,2), mar=c(5,4,1,1))
plot(log(dat$hrs[dat$outcome == "verbal"]), dat$yi[dat$outcome == "verbal"],
     pch=19, xlab="Log(Coaching Hours)", ylab="Effect Size (verbal)",
     xlim=c(1,6), ylim=c(-0.5,1), xaxs="i", yaxs="i")
abline(h=c(-0.5,0,0.5), lty="dotted")
abline(v=log(c(5,18)),  lty="dotted")
plot(log(dat$hrs[dat$outcome == "math"]), dat$yi[dat$outcome == "math"],
     pch=19, xlab="Log(Coaching Hours)", ylab="Effect Size (math)",
     xlim=c(1,6), ylim=c(-1.0,1), xaxs="i", yaxs="i")
abline(h=c(-0.5,0,0.5), lty="dotted")
abline(v=log(c(5,18)),  lty="dotted")

### construct variance-covariance matrices assuming rho = 0.66
vcalc <- function(v, rho) {
   S <- diag(sqrt(v), nrow=length(v), ncol=length(v))
   R <- matrix(rho, nrow=length(v), ncol=length(v))
   diag(R) <- 1
   S \%*\% R \%*\% S
}
V <- lapply(split(dat$v, dat$study), vcalc, rho=0.66)
V <- bldiag(V, order=dat$study)

### fit multivariate random-effects model
res <- rma.mv(yi, V, mods = ~ outcome - 1,
              random = ~ outcome | study, struct="UN",
              data=dat, digits=3)
res

### test whether the effect differs for the math and verbal subtest
anova(res, X=c(1,-1))

### log-transform and mean center the hours of coaching variable
dat$loghrs <- log(dat$hrs) - mean(log(dat$hrs), na.rm=TRUE)

### fit multivariate model with log(hrs) as moderator
res <- rma.mv(yi, V, mods = ~ outcome + outcome:loghrs - 1,
              random = ~ outcome | study, struct="UN",
              data=dat, digits=3)
res

### fit model with tau2 = 0 for outcome verbal (which also constrains rho = 0)
res <- rma.mv(yi, V, mods = ~ outcome + outcome:loghrs - 1,
              random = ~ outcome | study, struct="UN", tau2=c(NA,0),
              data=dat, digits=3)
res
}
\keyword{datasets}
