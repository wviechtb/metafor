\name{dat.maire2019}
\docType{data}
\alias{dat.maire2019}
\title{Studies on Temporal Trends in Fish Community Structures in French Rivers}
\description{Results from studies examining changes in the abundance of fish species in French rivers.}
\usage{
dat.maire2019
}
\format{The object is a list containing a data frame called \code{dat} that contains the following columns and distance matrix called \code{dmat}:
\tabular{lll}{
\bold{site}         \tab \code{character} \tab study site \cr
\bold{station}      \tab \code{character} \tab sampling station at site \cr
\bold{site_station} \tab \code{character} \tab site and station combined \cr
\bold{s1}           \tab \code{numeric}   \tab Mann-Kendal trend statistic for relative abundance of non-local species \cr
\bold{vars1}        \tab \code{numeric}   \tab corresponding sampling variance (corrected for temporal autocorrelation) \cr
\bold{s2}           \tab \code{numeric}   \tab Mann-Kendal trend statistic for relative abundance of northern species \cr
\bold{vars2}        \tab \code{numeric}   \tab corresponding sampling variance (corrected for temporal autocorrelation) \cr
\bold{s3}           \tab \code{numeric}   \tab Mann-Kendal trend statistic for relative abundance of non-native species \cr
\bold{vars3}        \tab \code{numeric}   \tab corresponding sampling variance (corrected for temporal autocorrelation) \cr
\bold{const}        \tab \code{numeric}   \tab constant value of 1
}
}
\details{
   The dataset includes the results from 35 sampling stations (at 11 sites along various French rivers) examining the abundance of various fish species over time (i.e., over 19-37 years, all until 2015). The temporal trend in these abundance data was quantified in terms of Mann-Kendal trend statistics, with positive values indicating monotonically increasing trends. The corresponding sampling variances were corrected for the temporal autocorrelation in the data (Hamed & Rao, 1998).

   The distance matrix \code{dmat} indicates the distance of the sampling stations (1-423 river-km). For stations not connected through the river network, a high distance value of 10,000 river-km was set (effectively forcing the spatial correlation to be 0 for such stations).

   The dataset can be used to illustrate a meta-analysis allowing for spatial correlation in the outcomes.
}
\source{
   Maire, A., Thierry, E., Viechtbauer, W., & Daufresne, M. (2019). Poleward shift in large-river fish communities detected with a novel meta-analysis framework. \emph{Freshwater Biology}, \bold{64}(6), 1143--1156. \verb{https://doi.org/10.1111/fwb.13291}
}
\references{
   Hamed, K. H., & Rao, A. R. (1998). A modified Mann-Kendall trend test for autocorrelated data. \emph{Journal of Hydrology}, \bold{204}(1-4), 182--196. \verb{https://doi.org/10.1016/S0022-1694(97)00125-X}
}
\examples{
### copy data into 'dat' and examine data
dat <- dat.maire2019$dat
dat

### copy distance matrix into 'dmat' and examine first 5 rows/columns
dmat <- dat.maire2019$dmat
dmat[1:5,1:5]

### fit a standard random-effects model ignoring spatial correlation
res1 <- rma.mv(s1, vars1, random = ~ 1 | site_station, data=dat)
res1

### fit model allowing for spatial correlation
res2 <- rma.mv(s1, vars1, random = ~ site_station | const, struct="SPGAU",
               data=dat, dist=list(dmat), control=list(rho.init=10))
res2

### add random effects for sites and stations within sites
res3 <- rma.mv(s1, vars1, random = list(~ 1 | site/station, ~ site_station | const), struct="SPGAU",
               data=dat, dist=list(dmat), control=list(rho.init=10))
res3

### likelihood ratio tests comparing the models
anova(res1, res2)
anova(res2, res3)

### profile likelihood plots for model res2
profile(res2, cline=TRUE)

### effective range (river-km for which the spatial correlation is >= .05)
sqrt(3) * res2$rho
}
\keyword{datasets}
