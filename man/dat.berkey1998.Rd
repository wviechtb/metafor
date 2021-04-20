\name{dat.berkey1998}
\docType{data}
\alias{dat.berkey1998}
\title{Studies on Treatments for Periodontal Disease}
\description{Results from 5 trials comparing surgical and non-surgical treatments for medium-severity periodontal disease one year after treatment. \loadmathjax}
\usage{dat.berkey1998}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{trial}      \tab \code{numeric}    \tab trial number \cr
\bold{author}     \tab \code{character}  \tab study author(s) \cr
\bold{year}       \tab \code{numeric}    \tab publication year \cr
\bold{ni}         \tab \code{numeric}    \tab number of patients \cr
\bold{outcome}    \tab \code{character}  \tab outcome (PD = probing depth; AL = attachment level) \cr
\bold{yi}         \tab \code{numeric}    \tab observed mean difference in outcome (surgical versus non-surgical) \cr
\bold{vi}         \tab \code{numeric}    \tab corresponding sampling variance \cr
\bold{v1i}        \tab \code{numeric}    \tab variances and covariances of the observed effects \cr
\bold{v2i}        \tab \code{numeric}    \tab variances and covariances of the observed effects
}
}
\details{
   The dataset includes the results from 5 trials that compared surgical and non-surgical methods for the treatment of medium-severity periodontal disease. Reported outcomes include the change in probing depth (PD) and attachment level (AL) one year after the treatment. The outcome measure used for this meta-analysis was the (raw) mean difference, calculated in such a way that positive values indicate that surgery was more effective than non-surgical treatment in decreasing the probing depth and increasing the attachment level (so, the results from the various trials indicate that surgery is preferable for reducing the probing depth, while non-surgical treatment is preferable for increasing the attachment level). Since each trial provides effect size estimates for both outcomes, the estimates are correlated. A multivariate model can be used to meta-analyze the two outcomes simultaneously.

   The \code{v1i} and \code{v2i} values are the variances and covariances of the observed effects. In particular, for each study, variables \code{v1i} and \code{v2i} form a \mjeqn{2 \times 2}{2x2} variance-covariance matrix of the observed effects, with the diagonal elements corresponding to the sampling variances of the mean differences (the first for probing depth, the second for attachment level) and the off-diagonal value corresponding to the covariance of the two mean differences. Below, the full (block diagonal) variance-covariance for all studies is constructed from these two variables.
}
\source{
   Berkey, C. S., Antczak-Bouckoms, A., Hoaglin, D. C., Mosteller, F., & Pihlstrom, B. L. (1995). Multiple-outcomes meta-analysis of treatments for periodontal disease. \emph{Journal of Dental Research}, \bold{74}(4), 1030--1039. \verb{https://doi.org/10.1177/00220345950740040201}

   Berkey, C. S., Hoaglin, D. C., Antczak-Bouckoms, A., Mosteller, F., & Colditz, G. A. (1998). Meta-analysis of multiple outcomes by regression with random effects. \emph{Statistics in Medicine}, \bold{17}(22), 2537--2550. \verb{https://doi.org/10.1002/(sici)1097-0258(19981130)17:22<2537::aid-sim953>3.0.co;2-c}
}
\examples{
### copy data into 'dat' and examine data
dat <- dat.berkey1998
dat

### construct list with the variance-covariance matrices of the observed outcomes for the studies
V <- lapply(split(dat[c("v1i", "v2i")], dat$trial), as.matrix)

### construct block diagonal matrix
V <- bldiag(V)

### fit multiple outcomes (meta-regression) model (with REML estimation)
res <- rma.mv(yi, V, mods = ~ outcome - 1, random = ~ outcome | trial, struct="UN", data=dat)
print(res, digits=3)

### test/estimate difference between the two outcomes
anova(res, X=c(1,-1))

### fit model including publication year as moderator for both outcomes (with ML estimation)
res <- rma.mv(yi, V, mods = ~ outcome + outcome:I(year - 1983) - 1,
              random = ~ outcome | trial, struct="UN", data=dat, method="ML")
print(res, digits=3)
}
\keyword{datasets}
