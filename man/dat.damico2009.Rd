\name{dat.damico2009}
\docType{data}
\alias{dat.damico2009}
\title{Studies on Topical plus Systemic Antibiotics to Prevent Respiratory Tract Infections}
\description{Results from 16 studies examining the effectiveness of topical plus systemic antibiotics to prevent respiratory tract infections (RTIs).}
\usage{
dat.damico2009
}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{study}   \tab \code{character} \tab first author \cr
\bold{year}    \tab \code{numeric}   \tab publication year \cr
\bold{xt}      \tab \code{numeric}   \tab number of RTIs in the treatment group \cr
\bold{nt}      \tab \code{numeric}   \tab number of patients in the treatment group \cr
\bold{xc}      \tab \code{numeric}   \tab number of RTIs in the control group \cr
\bold{nc}      \tab \code{numeric}   \tab number of patients in the control group \cr
\bold{conceal} \tab \code{numeric}   \tab allocation concealment (0 = not adequate, 1 = adequate) \cr
\bold{blind}   \tab \code{numeric}   \tab blinding (0 = open, 1 = double-blind)
}
}
\details{
   The dataset includes the results from 16 studies that examined the effectiveness of topical plus systemic antibiotics versus no prophylaxis to prevent respiratory tract infections (RTIs).
}
\source{
   D'Amico, R., Pifferi, S., Torri, V., Brazzi, L., Parmelli, E., & Liberati, A. (2009). Antibiotic prophylaxis to reduce respiratory tract infections and mortality in adults receiving intensive care. \emph{Cochrane Database of Systematic Reviews}, \bold{4}, CD000022. \verb{https://doi.org/10.1002/14651858.CD000022.pub3}
}
\examples{
### copy data into 'dat' and examine data
dat <- dat.damico2009
dat

### meta-analysis of the (log) odds ratios using the Mantel-Haenszel method
rma.mh(measure="OR", ai=xt, n1i=nt, ci=xc, n2i=nc, data=dat, digits=2)

### compute log odds ratios and corresponding sampling variances
dat <- escalc(measure="OR", ai=xt, n1i=nt, ci=xc, n2i=nc, data=dat)

### meta-analysis using a random-effects model
res <- rma(yi, vi, data=dat, method="DL")
res
predict(res, transf=exp, digits=2)
}
\keyword{datasets}
