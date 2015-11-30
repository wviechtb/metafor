\name{dat.yusuf1985}
\docType{data}
\alias{dat.yusuf1985}
\title{Studies of Beta Blockers During and After Myocardial Infarction}
\description{Results from studies examining the effectiveness of beta blockers for reducing mortality and reinfarction.}
\usage{dat.yusuf1985}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{table}  \tab \code{character} \tab table number \cr
\bold{id}     \tab \code{character} \tab trial id number \cr
\bold{trial}  \tab \code{character} \tab trial name or first author \cr
\bold{ai}     \tab \code{numeric}   \tab number of deaths/reinfarctions in treatment group \cr
\bold{n1i}    \tab \code{numeric}   \tab number of patients in treatment group \cr
\bold{ci}     \tab \code{numeric}   \tab number of deaths/reinfarctions in control group \cr
\bold{n2i}    \tab \code{numeric}   \tab number of patients in control group
}
}
\details{
   The dataset contains table 6 (total mortality from short-term trials of oral beta blockers), 9 (total mortality at one week from trials with an initial IV dose of a beta blocker), 10 (total mortality from long-term trials with treatment starting late and mortality from day 8 onwards in long-term trials that began early and continued after discharge), 11 (nonfatal reinfarction from long-term trials of beta blockers), 12a (sudden death in long-term beta blocker trials), and 12b (nonsudden death in long-term beta blocker trials) from the meta-analysis by Yusuf et al. (1985) on the effectiveness of of beta blockers for reducing mortality and reinfarction.

   The article also describes what is sometimes called Peto's one-step method for meta-analyzing \eqn{2 \times 2}{2x2} table data. This method is implemented in the \code{\link{rma.peto}} function.
}
\source{
   Yusuf, S., Peto, R., Lewis, J., Collins, R., & Sleight, P. (1985). Beta blockade during and after myocardial infarction: An overview of the randomized trials. \emph{Progress in Cardiovascular Disease}, \bold{27}, 335--371.
}
\examples{
### load data
dat <- get(data(dat.yusuf1985))

### to select a table for the analysis
tab <- "6" ### either: 6, 9, 10, 11, 12a, 12b

### to double-check total counts as reported in article
apply(dat[dat$table==tab,4:7], 2, sum, na.rm=TRUE)

### meta-analysis using Peto's one-step method
res <- rma.peto(ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, subset=(table==tab))
res
predict(res, transf=exp, digits=2)
}
\keyword{datasets}
