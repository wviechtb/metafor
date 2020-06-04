\name{dat.colditz1994}
\docType{data}
\alias{dat.colditz1994}
\alias{dat.bcg}
\title{Studies on the Effectiveness of the BCG Vaccine Against Tuberculosis}
\description{Results from 13 studies examining the effectiveness of the Bacillus Calmette-Guerin (BCG) vaccine against tuberculosis. \loadmathjax}
\usage{
dat.colditz1994
dat.bcg
}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{trial}  \tab \code{numeric}   \tab trial number \cr
\bold{author} \tab \code{character} \tab author(s) \cr
\bold{year}   \tab \code{numeric}   \tab publication year \cr
\bold{tpos}   \tab \code{numeric}   \tab number of TB positive cases in the treated (vaccinated) group \cr
\bold{tneg}   \tab \code{numeric}   \tab number of TB negative cases in the treated (vaccinated) group \cr
\bold{cpos}   \tab \code{numeric}   \tab number of TB positive cases in the control (non-vaccinated) group \cr
\bold{cneg}   \tab \code{numeric}   \tab number of TB negative cases in the control (non-vaccinated) group \cr
\bold{ablat}  \tab \code{numeric}   \tab absolute latitude of the study location (in degrees) \cr
\bold{alloc}  \tab \code{character} \tab method of treatment allocation (random, alternate, or systematic assignment)
}
}
\details{
   The 13 studies provide data in terms of \mjeqn{2 \times 2}{2x2} tables in the form:
   \tabular{lcc}{
                    \tab TB positive \tab TB negative \cr
   vaccinated group \tab \code{tpos} \tab \code{tneg} \cr
   control group    \tab \code{cpos} \tab \code{cneg}
   } The goal of the meta-analysis was to examine the overall effectiveness of the BCG vaccine for preventing tuberculosis and to examine moderators that may potentially influence the size of the effect.

   The dataset has been used in several publications to illustrate meta-analytic methods (see \sQuote{References}).
}
\source{
   Colditz, G. A., Brewer, T. F., Berkey, C. S., Wilson, M. E., Burdick, E., Fineberg, H. V., & Mosteller, F. (1994). Efficacy of BCG vaccine in the prevention of tuberculosis: Meta-analysis of the published literature. \emph{Journal of the American Medical Association}, \bold{271}, 698--702.
}
\references{
   Berkey, C. S., Hoaglin, D. C., Mosteller, F., & Colditz, G. A. (1995). A random-effects regression model for meta-analysis. \emph{Statistics in Medicine}, \bold{14}, 395--411.

   van Houwelingen, H. C., Arends, L. R., & Stijnen, T. (2002). Advanced methods in meta-analysis: Multivariate approach and meta-regression. \emph{Statistics in Medicine}, \bold{21}, 589--624.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{https://www.jstatsoft.org/v036/i03}.
}
\examples{
### copy data into 'dat'
dat <- dat.bcg

### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat)
dat

### random-effects model
res <- rma(yi, vi, data=dat)
res

### average risk ratio with 95\% CI
predict(res, transf=exp)

### mixed-effects model with absolute latitude and publication year as moderators
res <- rma(yi, vi, mods = ~ ablat + year, data=dat)
res

### predicted average risk ratios for 10-60 degrees absolute latitude
### holding the publication year constant at 1970
predict(res, newmods=cbind(seq(from=10, to=60, by=10), 1970), transf=exp)

### note: the interpretation of the results is difficult because absolute
### latitude and publication year are strongly correlated (the more recent
### studies were conducted closer to the equator)
plot(dat$ablat, dat$year)
cor(dat$ablat, dat$year)
}
\keyword{datasets}
