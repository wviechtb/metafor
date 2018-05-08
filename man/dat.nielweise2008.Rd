\name{dat.nielweise2008}
\docType{data}
\alias{dat.nielweise2008}
\title{Studies on Anti-Infective-Treated Central Venous Catheters for Prevention of Catheter-Related Bloodstream Infections}
\description{Results from 18 studies comparing the risk of catheter-related bloodstream infection when using anti-infective-treated versus standard catheters for total parenteral nutrition or chemotherapy.}
\usage{dat.nielweise2008}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{study}   \tab \code{numeric}   \tab study number \cr
\bold{authors} \tab \code{character} \tab study authors \cr
\bold{year}    \tab \code{numeric}   \tab publication year \cr
\bold{x1i}     \tab \code{numeric}   \tab number of CRBSIs in patients receiving an anti-infective catheter \cr
\bold{t1i}     \tab \code{numeric}   \tab total number of catheter days for patients receiving an anti-infective catheter \cr
\bold{x2i}     \tab \code{numeric}   \tab number of CRBSIs in patients receiving a standard catheter \cr
\bold{t2i}     \tab \code{numeric}   \tab total number of catheter days for patients receiving a standard catheter \cr
}
}
\details{
   The use of a central venous catheter may lead to a catheter-related bloodstream infection (CRBSI), which in turn increases the risk of morbidity and mortality. Anti-infective-treated catheters have been developed that are meant to reduce the risk of CRBSIs. Niel-Weise et al. (2008) conducted a meta-analysis of studies comparing infection risk when using anti-infective-treated versus standard catheters for total parenteral nutrition or chemotherapy. The results from 9 such studies are included in this dataset.

   The dataset was used in the article by Stijnen et al. (2010) to illustrate various generalized linear mixed-effects models for the meta-analysis of incidence rates and incidence rate ratios (see \sQuote{References}).
}
\source{
   Niel-Weise, B. S., Stijnen, T., & van den Broek, P. J. (2008). Anti-infective-treated central venous catheters for total parenteral nutrition or chemotherapy: A systematic review. \emph{Journal of Hospital Infection}, \bold{69}, 114--123.
}
\references{
   Stijnen, T., Hamza, T. H., & Ozdemir, P. (2010). Random effects meta-analysis of event outcome in the framework of the generalized linear mixed model with applications in sparse data. \emph{Statistics in Medicine}, \bold{29}, 3046--3067.
}
\examples{
### load data
dat <- get(data(dat.nielweise2008))

### standard (inverse-variance) random-effects model
res <- rma(measure="IRR", x1i=x1i, t1i=t1i, x2i=x2i, t2i=t2i, data=dat)
print(res, digits=3)
predict(res, transf=exp, digits=2)

### random-effects conditional Poisson model
\dontrun{
res <- rma.glmm(measure="IRR", x1i=x1i, t1i=t1i, x2i=x2i, t2i=t2i, data=dat, model="CM.EL")
print(res, digits=3)
predict(res, transf=exp, digits=2)}
}
\keyword{datasets}
