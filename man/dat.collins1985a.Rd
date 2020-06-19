\name{dat.collins1985a}
\docType{data}
\alias{dat.collins1985a}
\title{Studies on the Treatment of Upper Gastrointestinal Bleeding by a Histamine H2 Antagonist}
\description{Results from studies examining the effectiveness of histamine H2 antagonists (cimetidine or ranitidine) in treating patients with acute upper gastrointestinal hemorrhage.}
\usage{dat.collins1985a}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{id}    \tab \code{numeric}   \tab study number \cr
\bold{trial} \tab \code{character} \tab first author of trial \cr
\bold{year}  \tab \code{numeric}   \tab year of publication \cr
\bold{ref}   \tab \code{numeric}   \tab reference number \cr
\bold{trt}   \tab \code{character} \tab C = cimetidine, R = ranitidine \cr
\bold{ctrl}  \tab \code{character} \tab P = placebo, AA = antacids, UT = usual treatment \cr
\bold{nti}   \tab \code{numeric}   \tab number of patients in treatment group \cr
\bold{b.xti} \tab \code{numeric}   \tab number of patients in treatment group with persistent or recurrent bleedings \cr
\bold{o.xti} \tab \code{numeric}   \tab number of patients in treatment group in need of operation \cr
\bold{d.xti} \tab \code{numeric}   \tab number of patients in treatment group that died \cr
\bold{nci}   \tab \code{numeric}   \tab number of patients in control group \cr
\bold{b.xci} \tab \code{numeric}   \tab number of patients in control group with persistent or recurrent bleedings \cr
\bold{o.xci} \tab \code{numeric}   \tab number of patients in control group in need of operation \cr
\bold{d.xci} \tab \code{numeric}   \tab number of patients in control group that died
}
}
\details{
   The data were obtained from Tables 1 and 2 in Collins and Langman (1985). The authors used Peto's (one-step) method for meta-analyzing the 27 trials. This approach is implemented in the \code{\link{rma.peto}} function. Using the same dataset, van Houwelingen, Zwinderman, and Stijnen (1993) describe some alternative approaches for analyzing these data, including fixed and random-effects conditional logistic models. Those are implemented in the \code{\link{rma.glmm}} function.
}
\source{
   Collins, R., & Langman, M. (1985). Treatment with histamine H2 antagonists in acute upper gastrointestinal hemorrhage. \emph{New England Journal of Medicine}, \bold{313}(11), 660--666. \verb{https://doi.org/10.1056/NEJM198509123131104}
}
\references{
   van Houwelingen, H. C., Zwinderman, K. H., & Stijnen, T. (1993). A bivariate approach to meta-analysis. \emph{Statistics in Medicine}, \bold{12}(24), 2273--2284. \verb{https://doi.org/10.1002/sim.4780122405}
}
\examples{
### copy data into 'dat' and examine data
dat <- dat.collins1985a
dat

### meta-analysis of log ORs using Peto's method (outcome: persistent or recurrent bleedings)
res <- rma.peto(ai=b.xti, n1i=nti, ci=b.xci, n2i=nci, data=dat)
print(res, digits=2)

\dontrun{
### meta-analysis of log ORs using a conditional logistic regression model (FE model)
res <- rma.glmm(measure="OR", ai=b.xti, n1i=nti, ci=b.xci, n2i=nci, data=dat,
                model="CM.EL", method="FE")
summary(res)
predict(res, transf=exp, digits=2)

### plot the likelihoods of the odds ratios
llplot(measure="OR", ai=b.xti, n1i=nti, ci=b.xci, n2i=nci, data=dat,
       lwd=1, refline=NA, xlim=c(-4,4), drop00=FALSE)

### meta-analysis of log odds ratios using a conditional logistic regression model (RE model)
res <- rma.glmm(measure="OR", ai=b.xti, n1i=nti, ci=b.xci, n2i=nci, data=dat,
                model="CM.EL", method="ML")
summary(res)
predict(res, transf=exp, digits=2)
}

### meta-analysis of log ORs using Peto's method (outcome: need for surgery)
res <- rma.peto(ai=o.xti, n1i=nti, ci=o.xci, n2i=nci, data=dat)
print(res, digits=2)

### meta-analysis of log ORs using Peto's method (outcome: death)
res <- rma.peto(ai=d.xti, n1i=nti, ci=d.xci, n2i=nci, data=dat)
print(res, digits=2)
}
\keyword{datasets}
