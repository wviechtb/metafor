\name{dat.craft2003}
\docType{data}
\alias{dat.craft2003}
\title{Studies on the Relationship between the Competitive State Anxiety Inventory-2 and Sport Performance}
\description{Results from 10 studies on the relationship between the Competitive State Anxiety Inventory-2 (CSAI-2) and sport performance.}
\usage{dat.craft2003}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{study} \tab \code{numeric}   \tab study number \cr
\bold{ni}    \tab \code{numeric}   \tab sample size \cr
\bold{sport} \tab \code{character} \tab type of sport (T = team sport, I = individual sport) \cr
\bold{ri}    \tab \code{numeric}   \tab correlation coefficient \cr
\bold{var1}  \tab \code{character} \tab variable 1 of the correlation coefficient (see \sQuote{Details}) \cr
\bold{var2}  \tab \code{character} \tab variable 2 of the correlation coefficient (see \sQuote{Details})
}
}
\details{
   The 10 studies included in this dataset are a subset of the studies included in the meta-analysis by Craft et al. (2003) on the relationship between the Competitive State Anxiety Inventory-2 (CSAI-2) and sport performance.

   The CSAI-2 has three subscales: cognitive anxiety (\code{acog}), somatic anxiety (\code{asom}), and self-confidence (\code{conf}). The studies included in this dataset administered the CSAI-2 prior to some sport competition and then measured sport performance based on the competition. Most studies provided all 6 correlations (3 for the correlations among the 3 subscales and 3 for the correlations between the subscales and sport performance), but 2 studies (with study numbers 6 and 17) only provided a subset.
}
\source{
   Becker, B. J., & Aloe, A. M. (2019). Model-based meta-analysis and related approaches. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta-analysis} (3nd ed., pp. 339--363). New York: Russell Sage Foundation.
}
\references{
   Craft, L. L., Magyar, T. M., Becker, B. J., & Feltz, D. L. (2003). The relationship between the Competitive State Anxiety Inventory-2 and sport performance: A meta-analysis. \emph{Journal of Sport and Exercise Psychology}, \bold{25}(1), 44--65. \verb{https://doi.org/10.1123/jsep.25.1.44}
}

\examples{
### copy data into 'dat'
dat <- dat.craft2003
dat

### construct dataset and var-cov matrix of the correlations
tmp <- rcalc(ri ~ var1 + var2 | study, ni=ni, data=dat)
V <- tmp$V
dat <- tmp$dat

### examine data for study 1
dat[dat$study == 1,]
V[dat$study == 1, dat$study == 1]

### examine data for study 6
dat[dat$study == 6,]
V[dat$study == 6, dat$study == 6]

### examine data for study 17
dat[dat$study == 17,]
V[dat$study == 17, dat$study == 17]

### multivariate random-effects model
res <- rma.mv(yi, V, mods = ~ var1.var2 - 1, random = ~ var1.var2 | study, struct="UN", data=dat)
res
}
\keyword{datasets}
