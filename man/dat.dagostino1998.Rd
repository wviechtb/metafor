\name{dat.dagostino1998}
\docType{data}
\alias{dat.dagostino1998}
\title{Studies on the Effectiveness of Antihistamines in Reducing Symptoms of the Common Cold}
\description{Results from 9 studies on the effectiveness of antihistamines in reducing the severity of runny nose and sneezing in the common cold.}
\usage{dat.dagostino1998}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{study}    \tab \code{numeric}   \tab study id \cr
\bold{cold}     \tab \code{character} \tab natural or induced cold study \cr
\bold{scale.rn} \tab \code{character} \tab scale for measuring runny nose severity \cr
\bold{scale.sn} \tab \code{character} \tab scale for measuring sneezing severity \cr
\bold{drug}     \tab \code{character} \tab type of antihistamine studied \cr
\bold{tnt}      \tab \code{numeric}   \tab total sample size of the treatment group \cr
\bold{tnc}      \tab \code{numeric}   \tab total sample size of the control (placebo) group \cr
\bold{outcome}  \tab \code{character} \tab outcome variable (see \sQuote{Details})\cr
\bold{mt}       \tab \code{numeric}   \tab mean in the treatment group \cr
\bold{sdt}      \tab \code{numeric}   \tab SD in the treatment group \cr
\bold{mc}       \tab \code{numeric}   \tab mean in the control group \cr
\bold{sdc}      \tab \code{numeric}   \tab SD in the control group \cr
\bold{xt}       \tab \code{numeric}   \tab number of patients reaching the therapy goal in the treatment group \cr
\bold{xc}       \tab \code{numeric}   \tab number of patients reaching the therapy goal in the control (placebo) group \cr
\bold{nt}       \tab \code{numeric}   \tab sample size of the treatment group for measuring the outcome \cr
\bold{nc}       \tab \code{numeric}   \tab sample size of the control group for measuring the outcome
}
}
\details{
   The studies for this meta-analysis were assembled to examine the effectiveness of antihistamines in reducing the severity of runny nose and sneezing in the common cold. Effectiveness was measured after one and two days of treatment in terms of 4 different outcome variables:
   \enumerate{
   \item \code{rnic1} and \code{rnic2} (continuous): incremental change (improvement) in runny nose severity at day 1 and day 2,
   \item \code{rngoal1} and \code{rngoal2} (dichotomous): reaching the goal of therapy (of at least a 50\% reduction in runny nose severity) at day 1 and day 2,
   \item \code{snic1} and \code{snic2} (continuous): incremental change (improvement) in sneezing severity at day 1 and day 2, and
   \item \code{rngoal1} and \code{rngoal2} (dichotomous): reaching the goal of therapy (of at least a 50\% reduction in sneezing severity) at day 1 and day 2.
   }
   For the continuous outcomes, standardized mean differences can be computed to quantify the difference between the treatment and control groups. For the dichotomous outcomes, one can compute (log) odds ratios to quantify the difference between the treatment and control groups.
}
\source{
   D'Agostino, R. B., Sr., Weintraub, M., Russell, H. K., Stepanians, M., D'Agostino, R. B., Jr., Cantilena, L. R., Jr., Graumlich, J. F., Maldonado, S., Honig, P., & Anello, C. (1998). The effectiveness of antihistamines in reducing the severity of runny nose and sneezing: A meta-analysis. \emph{Clinical Pharmacology & Therapeutics}, \bold{64}(6), 579--596. \verb{https://doi.org/10.1016/S0009-9236(98)90049-2}
}
\examples{
### copy data into 'dat' and examine data
dat <- dat.dagostino1998
dat

### compute standardized mean differences and corresponding sampling variances
dat <- escalc(measure="SMD", m1i=mt, m2i=mc, sd1i=sdt, sd2i=sdc, n1i=nt, n2i=nc, data=dat,
              add.measure=TRUE)

### compute log odds ratios and corresponding sampling variances
dat <- escalc(measure="OR",  ai=xt, ci=xc, n1i=nt, n2i=nc, data=dat,
              replace=FALSE, add.measure=TRUE, add=1/2, to="all")

### inspect data for the first study
head(dat, 8)

### fit a random-effects model for incremental change in runny nose severity at day 1
res <- rma(yi, vi, data=dat, subset=outcome=="rnic1")
res

### fit a random-effects model for reaching the goal of therapy for runny nose severity at day 1
res <- rma(yi, vi, data=dat, subset=outcome=="rngoal1")
res
predict(res, transf=exp)

### construct approximate V matrix assuming a correlation of 0.7 for sampling errors within studies
V <- lapply(split(dat$vi, dat$study), function(v) {
   S <- diag(sqrt(v), nrow=length(v), ncol=length(v))
   R <- matrix(0.7, nrow=length(v), ncol=length(v))
   diag(R) <- 1
   S \%*\% R \%*\% S
})
V <- bldiag(V, order=dat$study)

### fit a model for incremental change in runny nose severity at day 1 and at day 2, allowing for
### correlated sampling errors (no random effects added, since there does not appear to be any
### noteworthy heterogeneity in these data)
res <- rma.mv(yi, V, mods = ~ outcome - 1, data=dat, subset=outcome \%in\% c("rnic1","rnic2"))
res

### test if there is a difference in effects at day 1 and day 2
anova(res, X=c(1,-1))
}
\keyword{datasets}
