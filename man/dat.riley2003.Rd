\name{dat.riley2003}
\docType{data}
\alias{dat.riley2003}
\title{Studies on MYC-N as a Prognostic Marker for Neuroblastoma}
\description{Results from 81 studies examining overall and disease-free survival in neuroblastoma patients with amplified versus normal MYC-N protein levels.}
\usage{dat.riley2003}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{study}      \tab \code{numeric}    \tab study number \cr
\bold{yi}         \tab \code{numeric}    \tab log hazard ratio of the outcome in those with amplified versus normal MYC-N protein levels \cr
\bold{sei}        \tab \code{numeric}    \tab standard error of the log hazard ratio \cr
\bold{outcome}    \tab \code{character}  \tab outcome (OS = overall survival; DFS = disease-free survival)
}
}
\details{
   The meta-analysis by Riley et al. (2003) examined a variety of prognostic markers for overall and disease-free survival in patients with neuroblastoma. One of the markers examined was amplified levels of the MYC-N protein, with is associated with poorer outcomes.

   The dataset given here was extracted from Riley (2011) and has been used in several other publications (e.g., Riley et al., 2004, 2007). The dataset provides the (log) hazard ratios (and corresponding standard errors) with respect to these two outcomes in 81 studies, with positive values indicating a greater risk of death (for OS) or disease recurrence/death (for DFS) for patients with high MYC-N levels compared to those with normal/low levels. Note that information on both outcomes could only be extracted from 17 studies (39 studies only provided sufficient information to extract the OS estimate, while 25 studies only allowed for extraction of the DFS estimate).
}
\source{
   Riley, R. D., Sutton, A. J., Abrams, K. R., & Lambert, P. C. (2004). Sensitivity analyses allowed more appropriate and reliable meta-analysis conclusions for multiple outcomes when missing data was present. \emph{Journal of Clinical Epidemiology}, \bold{57}, 911--924.

   Riley, R. D., Abrams, K. R., Lambert, P. C., Sutton, A. J., & Thompson, J. R. (2007). An evaluation of bivariate random-effects meta-analysis for the joint synthesis of two correlated outcomes. \emph{Statistics in Medicine}, \bold{26}, 78--97.

   Riley, R. D. (2011). Erratum. \emph{Statistics in Medicine}, \bold{30}, 400.
}
\references{
   Riley, R. D., Burchill, S. A., Abrams, K. R., Heney, D., Lambert, P. C., Jones, D. R., et al. (2003). A systematic review and evaluation of the use of tumour markers in paediatric oncology: Ewing's sarcoma and neuroblastoma. \emph{Health Technology Assessment}, \bold{7}, 1--162.
}
\examples{
### copy data into 'dat' and examine data
dat <- dat.riley2003
dat

### random-effects model analysis for outcome DFS
res <- rma(yi, sei=sei, data=dat, subset=(outcome == "DFS"), method="DL")
res
predict(res, transf=exp, digits=2)

### random-effects model analysis for outcome OS
res <- rma(yi, sei=sei, data=dat, subset=(outcome == "OS"), method="DL")
res
predict(res, transf=exp, digits=2)
}
\keyword{datasets}
