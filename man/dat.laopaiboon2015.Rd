\name{dat.laopaiboon2015}
\docType{data}
\alias{dat.laopaiboon2015}
\title{Studies on the Effectiveness of Azithromycin for Treating Lower Respiratory Tract Infections}
\description{Results from 15 studies the effectiveness of azithromycin versus amoxycillin or amoxycillin/clavulanic acid (amoxyclav) in the treatment of acute lower respiratory tract infections.}
\usage{
dat.laopaiboon2015
}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{author}  \tab \code{character} \tab author(s) \cr
\bold{year}    \tab \code{numeric}   \tab publication year \cr
\bold{ai}      \tab \code{numeric}   \tab number of clinical failures in the group treated with azithromycin \cr
\bold{n1i}     \tab \code{numeric}   \tab number of patients in the group treated with azithromycin \cr
\bold{ci}      \tab \code{numeric}   \tab number of clinical failures in the group treated with amoxycillin or amoxyclav \cr
\bold{n2i}     \tab \code{numeric}   \tab number of patients in the group treated with amoxycillin or amoxyclav \cr
\bold{age}     \tab \code{character} \tab whether the trial included adults or children \cr
\bold{diag.ab} \tab \code{numeric}   \tab trial included patients with a diagnosis of acute bacterial bronchitis \cr
\bold{diag.cb} \tab \code{numeric}   \tab trial included patients with a diagnosis of chronic bronchitis with acute exacerbation \cr
\bold{diag.pn} \tab \code{numeric}   \tab trial included patients with a diagnosis of pneumonia \cr
\bold{ctrl}    \tab \code{character} \tab antibiotic in control group (amoxycillin or amoxyclav) \cr
}
}
\details{
   Azithromycin is an antibiotic useful for the treatment of a number of bacterial infections. Laopaiboon et al. (2015) conducted a meta-analysis of trials comparing the effectiveness of azithromycin versus amoxycillin or amoxycillin/clavulanic acid (amoxyclav) in the treatment of acute lower respiratory tract infections, including acute bacterial bronchitis, acute exacerbations of chronic bronchitis, and pneumonia. The results from 15 trials are included in this dataset.
}
\source{
   Laopaiboon, M., Panpanich, R., & Swa Mya, K. (2015). Azithromycin for acute lower respiratory tract infections. \emph{Cochrane Database of Systematic Reviews}, \bold{3}, CD001954.
}
\examples{
### load data
dat <- get(data(dat.laopaiboon2015))

### analysis using the Mantel-Haenszel method
rma.mh(measure="RR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, digits=3)

### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat)

### random-effects model
res <- rma(yi, vi, data=dat)
res

### average risk ratio with 95\% CI
predict(res, transf=exp)
}
\keyword{datasets}
