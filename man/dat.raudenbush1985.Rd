\name{dat.raudenbush1985}
\docType{data}
\alias{dat.raudenbush1985}
\title{Studies on Assessing the Effects of Teacher Expectations on Pupil IQ}
\description{Results from 19 studies examining how teachers' expectations about their pupils can influence actual IQ levels.}
\usage{dat.raudenbush1985}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{study}      \tab \code{numeric}    \tab study number \cr
\bold{author}     \tab \code{character}  \tab study author(s) \cr
\bold{year}       \tab \code{numeric}    \tab publication year \cr
\bold{weeks}      \tab \code{numeric}    \tab weeks of contact prior to expectancy induction \cr
\bold{setting}    \tab \code{character}  \tab whether tests were group or individually administered \cr
\bold{tester}     \tab \code{character}  \tab whether test administrator was aware or blind \cr
\bold{n1i}        \tab \code{numeric}    \tab sample size of experimental group \cr
\bold{n2i}        \tab \code{numeric}    \tab sample size of control group \cr
\bold{yi}         \tab \code{numeric}    \tab standardized mean difference \cr
\bold{vi}         \tab \code{numeric}    \tab corresponding sampling variance
}
}
\details{
   In the so-called \sQuote{Pygmalion study} (Rosenthal & Jacobson, 1968), \dQuote{all of the predominantly poor children in the so-called Oak elementary school were administered a test pretentiously labeled the \sQuote{Harvard Test of Inflected Acquisition.} After explaining that this newly designed instrument had identified those children most likely to show dramatic intellectual growth during the coming year, the experimenters gave the names of these \sQuote{bloomers} to the teachers. In truth, the test was a traditional IQ test and the \sQuote{bloomers} were a randomly selected 20\% of the student population. After retesting the children 8 months later, the experimenters reported that those predicted to bloom had in fact gained significantly more in total IQ (nearly 4 points) and reasoning IQ (7 points) than the control group children. Further, at the end of the study, the teachers rated the experimental children as intellectually more curious, happier, better adjusted, and less in need of approval than their control group peers} (Raudenbush, 1984).

   In the following years, a series of studies were conducted attempting to replicate this rather controversial finding. However, the great majority of those studies were unable to demonstrate a statistically significant difference between the two experimental groups in terms of IQ scores. Raudenbush (1984) conducted a meta-analysis based on 19 such studies to further examine the evidence for the existence of the \sQuote{Pygmalion effect}. The dataset includes the results from these studies.

   The outcome measure used for the meta-analysis was the standardized mean difference (\code{yi}), with positive values indicating that the supposed \sQuote{bloomers} had, on average, higher IQ scores than those in the control group. The \code{weeks} variable indicates the number of weeks of prior contact between teachers and students before the expectancy induction. Testing was done either in a group setting or individually, which is indicated by the \code{setting} variable. Finally, the \code{tester} variable indicates whether the test administrators were either aware or blind to the researcher-provided designations of the children's intellectual potential.

   The data in this dataset were obtained from Raudenbush and Bryk (1985) with information on the \code{setting} and \code{tester} variables extracted from Raudenbush (1984).
}
\source{
   Raudenbush, S. W. (1984). Magnitude of teacher expectancy effects on pupil IQ as a function of the credibility of expectancy induction: A synthesis of findings from 18 experiments. \emph{Journal of Educational Psychology}, \bold{76}(1), 85--97. \verb{https://doi.org/10.1037/0022-0663.76.1.85}

   Raudenbush, S. W., & Bryk, A. S. (1985). Empirical Bayes meta-analysis. \emph{Journal of Educational Statistics}, \bold{10}(2), 75--98. \verb{https://doi.org/10.3102/10769986010002075}
}
\examples{
### copy data into 'dat' and examine data
dat <- dat.raudenbush1985
dat

### random-effects model
res <- rma(yi, vi, data=dat)
res

### create weeks variable where values larger than 3 are set to 3
dat$weeks.c <- ifelse(dat$week > 3, 3, dat$week)

### mixed-effects model with weeks.c variable as moderator
res <- rma(yi, vi, mods=~weeks.c, data=dat, digits=3)
res
}
\keyword{datasets}
