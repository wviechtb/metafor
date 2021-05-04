\name{dat.knapp2017}
\docType{data}
\alias{dat.knapp2017}
\title{Studies on Differences in Planning Performance in Schizophrenia Patients versus Healthy Controls}
\description{Results from 31 studies examining differences in planning performance in schizophrenia patients versus healthy controls.}
\usage{dat.knapp2017}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{author}     \tab \code{character} \tab study author(s) \cr
\bold{year}       \tab \code{numeric}   \tab publication year \cr
\bold{study}      \tab \code{numeric}   \tab study id number \cr
\bold{yi}         \tab \code{numeric}   \tab standardized mean difference for planning performance \cr
\bold{vi}         \tab \code{numeric}   \tab corresponding sampling variance \cr
\bold{difficulty} \tab \code{numeric}   \tab task difficulty \cr
\bold{n_sz}       \tab \code{numeric}   \tab number of schizophrenic patients \cr
\bold{n_hc}       \tab \code{numeric}   \tab number of healthy controls \cr
\bold{comp}       \tab \code{numeric}   \tab id for comparisons within studies \cr
\bold{yi}         \tab \code{numeric}   \tab standardized mean difference for IQ \cr
\bold{vi}         \tab \code{numeric}   \tab corresponding sampling variance
}
}
\details{
   The studies included in this dataset examined differences between schizophrenia patients and healthy controls with respect to their performance on the tower of London test (\url{https://en.wikipedia.org/wiki/Tower_of_London_test}) or a similar cognitive tasks measuring planning ability. The outcome measure for this meta-analysis was the standardized mean difference (with positive values indicating better performance in the healthy controls compared to the schizophrenia patients).

   The dataset has a more complex structure for two reasons:
   \enumerate{
   \item Studies 2, 3, 9, and 20 included more than schizophrenia patient group and the standardized mean differences were computed by comparing these groups against a single healthy control group.
   \item Studies 6, 12, 14, 15, 18, 19, 22, and 26 had the patients and controls complete different tasks of varying complexity (essentially the average number of moves required to complete a task).
   }
   Both of these issues lead to correlated sampling errors, which should be taken into consideration in the analysis.
}
\source{
   Knapp, F., Viechtbauer, W., Leonhart, R., Nitschke, K., & Kaller, C. P. (2017). Planning performance in schizophrenia patients: A meta-analysis of the influence of task difficulty and clinical and sociodemographic variables. \emph{Psychological Medicine}, \bold{47}(11), 2002--2016. \verb{https://doi.org/10.1017/S0033291717000459}
}
\examples{
### copy data into 'dat' and examine data
dat <- dat.knapp2017
dat

### fit a standard random-effects model ignoring correlated sampling errors
res <- rma(yi, vi, data=dat)
res

### fit a multilevel model with random effects for studies and comparisons within studies
res <- rma.mv(yi, vi, random = ~ 1 | study/comp, data=dat)
res

### construct an approximate V matrix assuming a correlation of 0.4 for the sampling errors
### of different comparisons within the same study
V <- lapply(split(dat$vi, dat$study), function(v) {
   S <- diag(sqrt(v), nrow=length(v), ncol=length(v))
   R <- matrix(0.4, nrow=length(v), ncol=length(v))
   diag(R) <- 1
   S \%*\% R \%*\% S
})
V <- bldiag(V, order=dat$study)

### fit the same multilevel model, but now use this V matrix in the model
res <- rma.mv(yi, V, random = ~ 1 | study/comp, data=dat)
res

### use cluster-robust inference methods based on this model
robust(res, cluster=dat$study)

### examine if task difficulty is a potential moderator of the effect
res <- rma.mv(yi, V, mods = ~ difficulty, random = ~ 1 | study/comp, data=dat)
res
sav <- robust(res, cluster=dat$study)
sav

### draw bubble plot
regplot(sav, xlab="Task Difficulty", ylab="Standardized Mean Difference", las=1, digits=1, bty="l")
}
\keyword{datasets}
