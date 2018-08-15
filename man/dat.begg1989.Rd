\name{dat.begg1989}
\docType{data}
\alias{dat.begg1989}
\title{Studies on Bone-Marrow Transplantation versus Chemotherapy for the Treatment of Leukemia}
\description{Results from controlled and uncontrolled studies on the effectiveness of allogeneic bone-marrow transplantation (BMT) and conventional chemotherapy (CMO) in the treatment of acute nonlymphocytic leukemia.}
\usage{dat.berkey1989}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{study}      \tab \code{numeric}    \tab study number \cr
\bold{trt}        \tab \code{character}  \tab treatment (BMT or CMO) \cr
\bold{arms}       \tab \code{numeric}    \tab number of arms in the study (1 = uncontrolled studies; 2 = controlled studies) \cr
\bold{yi}         \tab \code{numeric}    \tab 2-year disease-free survival \cr
\bold{sei}        \tab \code{numeric}    \tab corresponding standard error \cr
\bold{vi}         \tab \code{numeric}    \tab corresponding sampling variance
}
}
\details{
   The dataset includes the results from controlled and uncontrolled studies on the 2-year disease-free survival in patients with acute nonlymphocytic leukemia receiving either allogeneic bone-marrow transplantation (BMT) or conventional chemotherapy (CMO). In the controlled (two-arm) studies (studies 1-4), a cohort of patients in complete remission and potentially eligible for BMT was assembled, and those who consented and for whom a donor could be found received BMT, with the remaining patients used as controls (receiving CMO). In the uncontrolled (one-arm) studies (studies 5-16), only a single group was studied, receiving either BMT or CMO.

   The data in this dataset were obtained from Table 1 in Begg & Pilote (1991, p. 902).
}
\source{
   Begg, C. B., & Pilote, L. (1991). A model for incorporating historical controls into a meta-analysis. \emph{Biometrics}, \bold{47}, 899--906.
}
\references{
   Begg, C. B., Pilote, L., & McGlave, P. B. (1989). Bone marrow transplantation versus chemotherapy in acute non-lymphocytic leukemia: A meta-analytic review. \emph{European Journal of Cancer and Clinical Oncology}, \bold{25}, 1519--1523.
}
\examples{
### load data
dat <- get(data(dat.begg1989))

### turn trt and arms into factors and set reference levels
dat$trt  <- relevel(factor(dat$trt), ref="CMO")
dat$arms <- relevel(factor(dat$arms), ref="2")

### create data frame with the treatment differences for the controlled studies
dat2 <- data.frame(yi = dat$yi[c(1,3,5,7)] - dat$yi[c(2,4,6,8)],
                   vi = dat$vi[c(1,3,5,7)] + dat$vi[c(2,4,6,8)])

### DerSimonian and Laird method using the treatment differences
res <- rma(yi, vi, data=dat2, method="DL", digits=2)
res

### Begg & Pilote (1991) model incorporating the uncontrolled studies
res <- rma.mv(yi, vi, mods = ~ trt, random = ~ 1 | study,
              data=dat, method="ML", digits=2)
res

### model involving bias terms for the uncontrolled studies
res <- rma.mv(yi, vi, mods = ~ trt + trt:arms, random = ~ 1 | study,
              data=dat, method="ML", digits=2)
res

### model with random treatment effect
res <- rma.mv(yi, vi, mods = ~ trt, random = list(~ 1 | study, ~ trt | study),
              struct="UN", tau2=c(0,NA), rho=0, data=dat, method="ML", digits=2)
res

### model with random treatment effect, but with equal variances in both arms
res <- rma.mv(yi, vi, mods = ~ trt, random = list(~ 1 | study, ~ trt | study),
              struct="CS", rho=0, data=dat, method="ML", digits=2)
res
}
\keyword{datasets}
