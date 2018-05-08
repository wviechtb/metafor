\name{dat.nielweise2007}
\docType{data}
\alias{dat.nielweise2007}
\title{Studies on Anti-Infective-Treated Central Venous Catheters for Prevention of Catheter-Related Bloodstream Infections}
\description{Results from 18 studies comparing the risk of catheter-related bloodstream infection when using anti-infective-treated versus standard catheters in the acute care setting.}
\usage{dat.nielweise2007}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{study}  \tab \code{numeric}   \tab study number \cr
\bold{author} \tab \code{character} \tab (first) author \cr
\bold{year}   \tab \code{numeric}   \tab publication year \cr
\bold{ai}     \tab \code{numeric}   \tab number of CRBSIs in patients receiving an anti-infective catheter \cr
\bold{n1i}    \tab \code{numeric}   \tab number of patients receiving an anti-infective catheter \cr
\bold{ci}     \tab \code{numeric}   \tab number of CRBSIs in patients receiving a standard catheter \cr
\bold{n2i}    \tab \code{numeric}   \tab number of patients receiving a standard catheter \cr
}
}
\details{
   The use of a central venous catheter may lead to a catheter-related bloodstream infection (CRBSI), which in turn increases the risk of morbidity and mortality. Anti-infective-treated catheters have been developed that are meant to reduce the risk of CRBSIs. Niel-Weise et al. (2007) conducted a meta-analysis of studies comparing infection risk when using anti-infective-treated versus standard catheters in the acute care setting. The results from 18 such studies are included in this dataset.

   The dataset was used in the article by Stijnen et al. (2010) to illustrate various generalized linear mixed-effects models for the meta-analysis of proportions and odds ratios (see \sQuote{References}).
}
\source{
   Niel-Weise, B. S., Stijnen, T., & van den Broek, P. J. (2007). Anti-infective-treated central venous catheters: A systematic review of randomized controlled trials. \emph{Intensive Care Medicine}, \bold{33}, 2058--2068.
}
\references{
   Stijnen, T., Hamza, T. H., & Ozdemir, P. (2010). Random effects meta-analysis of event outcome in the framework of the generalized linear mixed model with applications in sparse data. \emph{Statistics in Medicine}, \bold{29}, 3046--3067.
}
\examples{
### load data
dat <- get(data(dat.nielweise2007))

### standard (inverse-variance) random-effects model
res <- rma(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, drop00=TRUE)
print(res, digits=3)
predict(res, transf=exp, digits=2)

### random-effects conditional logistic model
\dontrun{
res <- rma.glmm(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, model="CM.EL")
print(res, digits=3)
predict(res, transf=exp, digits=2)}
}
\keyword{datasets}
