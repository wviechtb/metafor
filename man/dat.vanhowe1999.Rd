\name{dat.vanhowe1999}
\docType{data}
\alias{dat.vanhowe1999}
\title{Studies on the Association between Circumcision and HIV Infection}
\description{Results from 33 studies examining the association between male circumcision and HIV infection. \loadmathjax}
\usage{
dat.vanhowe1999
}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{study}    \tab \code{character} \tab study author \cr
\bold{category} \tab \code{character} \tab study type (high-risk group, partner study, or population survey) \cr
\bold{non.pos}  \tab \code{numeric}   \tab number of non-circumcised HIV positive cases \cr
\bold{non.neg}  \tab \code{numeric}   \tab number of non-circumcised HIV negative cases \cr
\bold{cir.pos}  \tab \code{numeric}   \tab number of circumcised HIV positive cases \cr
\bold{cir.neg}  \tab \code{numeric}   \tab number of circumcised HIV negative cases
}
}
\details{
   The 33 studies provide data in terms of \mjeqn{2 \times 2}{2x2} tables in the form:
   \tabular{lcc}{
                   \tab HIV positive   \tab HIV negative \cr
   non-circumcised \tab \code{non.pos} \tab \code{non.neg} \cr
   circumcised     \tab \code{cir.pos} \tab \code{cir.neg}
   } The goal of the meta-analysis was to examine if the risk of an HIV infection differs between non-circumcised versus circumcised men.

   The dataset is interesting because it can be used to illustrate the difference between naively pooling results by summing up the counts across studies and then computing the odds ratio based on the aggregated table (as was done by Van Howe, 1999) and conducting a proper meta-analysis (as illustrated by O'Farrell & Egger, 2000). In fact, a proper meta-analysis shows that the HIV infection risk is on average higher in non-circumcised men, which is the opposite of what the naive pooling approach yields (which makes this an illustration of Simpson's paradox).
}
\source{
   Van Howe, R. S. (1999). Circumcision and HIV infection: Review of the literature and meta-analysis. \emph{International Journal of STD & AIDS}, \bold{10}(1), 8--16. \verb{https://doi.org/10.1258/0956462991913015}
}
\references{
   O'Farrell, N., & Egger, M. (2000). Circumcision in men and the prevention of HIV infection: A 'meta-analysis' revisited. \emph{International Journal of STD & AIDS}, \bold{11}(3), 137--142. \verb{https://doi.org/10.1258/0956462001915480}
}
\examples{
### copy data into 'dat'
dat <- dat.vanhowe1999

### naive pooling by summing up the counts within categories and then
### computing the odds ratios and corresponding confidence intervals
cat1 <- with(dat[dat$category=="high-risk group",],
   escalc(measure="OR", ai=sum(non.pos), bi=sum(non.neg), ci=sum(cir.pos), di=sum(cir.neg)))
cat2 <- with(dat[dat$category=="partner study",],
   escalc(measure="OR", ai=sum(non.pos), bi=sum(non.neg), ci=sum(cir.pos), di=sum(cir.neg)))
cat3 <- with(dat[dat$category=="population survey",],
   escalc(measure="OR", ai=sum(non.pos), bi=sum(non.neg), ci=sum(cir.pos), di=sum(cir.neg)))
summary(cat1, transf=exp, digits=2)
summary(cat2, transf=exp, digits=2)
summary(cat3, transf=exp, digits=2)

### naive pooling across all studies
all <- escalc(measure="OR", ai=sum(dat$non.pos), bi=sum(dat$non.neg),
                            ci=sum(dat$cir.pos), di=sum(dat$cir.neg))
summary(all, transf=exp, digits=2)

### calculate log odds ratios and corresponding sampling variances
dat <- escalc(measure="OR", ai=non.pos, bi=non.neg, ci=cir.pos, di=cir.neg, data=dat)
dat

### random-effects model
res <- rma(yi, vi, data=dat, method="DL")
res
predict(res, transf=exp, digits=2)

### random-effects model within subgroups
res <- rma(yi, vi, data=dat, method="DL", subset=category=="high-risk group")
predict(res, transf=exp, digits=2)
res <- rma(yi, vi, data=dat, method="DL", subset=category=="partner study")
predict(res, transf=exp, digits=2)
res <- rma(yi, vi, data=dat, method="DL", subset=category=="population survey")
predict(res, transf=exp, digits=2)
}
\keyword{datasets}
