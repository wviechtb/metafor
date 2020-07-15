\name{dat.lau1992}
\docType{data}
\alias{dat.lau1992}
\title{Studies on Intravenous Streptokinase for Acute Myocardial Infarction}
\description{Results from 33 trials comparing intravenous streptokinase versus placebo or no therapy in patients who had been hospitalized for acute myocardial infarction.}
\usage{dat.lau1992}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{trial} \tab \code{character} \tab trial name \cr
\bold{year}  \tab \code{numeric}   \tab publication year \cr
\bold{ai}    \tab \code{numeric}   \tab number of deaths in the streptokinase group \cr
\bold{n1i}   \tab \code{numeric}   \tab number of patients in the streptokinase group \cr
\bold{ci}    \tab \code{numeric}   \tab number of deaths in the control group \cr
\bold{n2i}   \tab \code{numeric}   \tab number of patients in the control group
}
}
\details{
   In the paper by Lau et al. (1992), the data are used to illustrate the idea of a cumulative meta-analysis, where the results are updated as each trial is added to the dataset. See \sQuote{Examples} for code that replicates the results and shows corresponding forest plots.
}
\source{
   Lau, J., Antman, E. M., Jimenez-Silva, J., Kupelnick, B., Mosteller, F., & Chalmers, T. C. (1992). Cumulative meta-analysis of therapeutic trials for myocardial infarction. \emph{New England Journal of Medicine}, \bold{327}(4), 248--254. \verb{https://doi.org/10.1056/NEJM199207233270406}
}
\examples{
### copy data into 'dat' and examine data
dat <- dat.lau1992
dat

### meta-analysis of log odds ratios using the MH method
res <- rma.mh(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, slab=trial)
print(res, digits=2)

### forest plot
forest(res, xlim=c(-10,9), atransf=exp, at=log(c(.01, 0.1, 1, 10, 100)),
       header=TRUE, top=2, ilab=dat$year, ilab.xpos=-6)
text(-6, 35, "Year", font=2)

### cumulative meta-analysis
sav <- cumul(res)

### forest plot of the cumulative results
forest(sav, xlim=c(-5,4), atransf=exp, at=log(c(0.1, 0.5, 1, 2, 10)),
       header=TRUE, top=2, ilab=dat$year, ilab.xpos=-3)
text(-3, 35, "Year", font=2)
id <- c(4, 8, 15, 33) # rows for which the z/p-values should be shown (as in Lau et al., 1992)
text(1.1, (res$k:1)[id], paste0("z = ", formatC(sav$zval[id], format="f", digits=2),
                                ", p = ", formatC(sav$pval[id], format="f", digits=4)))
}
\keyword{datasets}
