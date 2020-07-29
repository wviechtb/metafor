\name{dat.dorn2007}
\docType{data}
\alias{dat.dorn2007}
\title{Studies on Complementary and Alternative Medicine for Irritable Bowel Syndrome}
\description{Results from 19 trials examining complementary and alternative medicine (CAM) for irritable bowel syndrome (IBS).}
\usage{
dat.dorn2007
}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{id}       \tab \code{numeric}   \tab trial id number \cr
\bold{study}    \tab \code{character} \tab (first) author \cr
\bold{year}     \tab \code{numeric}   \tab publication year \cr
\bold{country}  \tab \code{character} \tab country where trial was conducted \cr
\bold{ibs.crit} \tab \code{character} \tab IBS diagnostic criteria (Manning, Rome I, Rome II, or Other) \cr
\bold{days}     \tab \code{numeric}   \tab number of treatment days \cr
\bold{visits}   \tab \code{numeric}   \tab number of practitioner visits \cr
\bold{jada}     \tab \code{numeric}   \tab Jadad score \cr
\bold{x.a}      \tab \code{numeric}   \tab number of responders in the active treatment group \cr
\bold{n.a}      \tab \code{numeric}   \tab number of participants in the active treatment group \cr
\bold{x.p}      \tab \code{numeric}   \tab number of responders in the placebo group \cr
\bold{n.p}      \tab \code{numeric}   \tab number of participants in the placebo group
}
}
\details{
   The dataset includes the results from 19 randomized clinical trials that examined the effectiveness of complementary and alternative medicine (CAM) for irritable bowel syndrome (IBS).
}
\note{
   The data were extracted from Table I in Dorn et al. (2009). Comparing the funnel plot in Figure 1 with the one obtained below indicates that the data for study 5 (Davis et al., 2006) in the table were not the ones that were used in the actual analyses.
}
\source{
   Dorn, S. D., Kaptchuk, T. J., Park, J. B., Nguyen, L. T., Canenguez, K., Nam, B. H., Woods, K. B., Conboy, L. A., Stason, W. B., & Lembo, A. J. (2007). A meta-analysis of the placebo response in complementary and alternative medicine trials of irritable bowel syndrome. \emph{Neurogastroenterology & Motility}, \bold{19}(8), 630--637. \verb{https://doi.org/10.1111/j.1365-2982.2007.00937.x}
}
\examples{
### copy data into 'dat' and examine data
dat <- dat.dorn2007
dat

### compute log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=x.a, n1i=n.a, ci=x.p, n2i=n.p, data=dat)

### random-effects model
res <- rma(yi, vi, data=dat, digits=2, method="DL")
res

### estimated average risk ratio
predict(res, transf=exp)

### funnel plot with study 5 highlighted in red
funnel(res, atransf=exp, at=log(c(.1, .2, .5, 1, 2, 5, 10)),
       ylim=c(0,1), steps=6, las=1, col=ifelse(dat$id == 5, "red", "black"))

### change log risk ratio for study 5
dat$yi[5] <- -0.44

### results are now more in line with what is reported in the paper
### (although the CI in the paper is not wide enough)
res <- rma(yi, vi, data=dat, digits=2, method="DL")
predict(res, transf=exp)

### funnel plot with study 5 highlighted in red
funnel(res, atransf=exp, at=log(c(.1, .2, .5, 1, 2, 5, 10)),
       ylim=c(0,1), steps=6, las=1, col=ifelse(dat$id == 5, "red", "black"))
}
\keyword{datasets}
