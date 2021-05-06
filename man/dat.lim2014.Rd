\name{dat.lim2014}
\docType{data}
\alias{dat.lim2014}
\title{Studies on the association between maternal size, offspring size, and number of offsprings}
\description{Results from studies examining the association between maternal size, offspring size, and number of offsprings.}
\usage{
dat.lim2014
}
\format{The object is a list containing data frames \code{m_o_size}, \code{m_o_fecundity}, \code{o_o_unadj}, and \code{o_o_adj} that contain the following columns and the corresponding phylogenetic trees called \code{m_o_size_tree}, \code{m_o_fecundity_tree}, \code{o_o_unadj_tree}, and \code{o_o_adj_tree}:
\tabular{lll}{
\bold{article}     \tab \code{numeric}   \tab article id \cr
\bold{author}      \tab \code{character} \tab study author(s) \cr
\bold{year}        \tab \code{numeric}   \tab publication year \cr
\bold{species}     \tab \code{character} \tab species \cr
\bold{amniotes}    \tab \code{character} \tab whether the species was amniotic \cr
\bold{environment} \tab \code{character} \tab whether the species were wild or captive \cr
\bold{reprounit}   \tab \code{character} \tab whether the data were based on lifetime reproductive output or a single reproductive event (only in \code{m_o_size} and \code{m_o_fecundity}) \cr
\bold{ri}          \tab \code{numeric}   \tab correlation coefficient \cr
\bold{ni}          \tab \code{numeric}   \tab sample size
}
}
\details{
   The object \code{dat.lim2014} includes 4 datasets:
   \tabular{ll}{
   \code{m_o_size}      \tab on the correlation between maternal size and offspring size \cr
   \code{m_o_fecundity} \tab on the correlation between maternal size and number of offsprings \cr
   \code{o_o_unadj}     \tab on the correlation between offspring size and number of offsprings \cr
   \code{o_o_adj}       \tab on the correlation between offspring size and number of offsprings adjusted for maternal size
   }
   Objects \code{m_o_size_tree}, \code{m_o_fecundity_tree}, \code{o_o_unadj_tree}, and \code{o_o_adj_tree} are the corresponding phylogenetic trees for the species included in each of these datasets.
}
\source{
   Lim, J. N., Senior, A. M., & Nakagawa, S. (2014). Heterogeneity in individual quality and reproductive trade-offs within species. \emph{Evolution}, \bold{68}(8), 2306--2318. \verb{https://doi.org/10.1111/evo.12446}
}
\references{
   Hadfield, J. D., & Nakagawa, S. (2010). General quantitative genetic methods for comparative biology: Phylogenies, taxonomies and multi-trait models for continuous and categorical characters. \emph{Journal of Evolutionary Biology}, \bold{23}(3), 494--508. \verb{https://doi.org/10.1111/j.1420-9101.2009.01915.x}

   Nakagawa, S., & Santos, E. S. A. (2012). Methodological issues and advances in biological meta-analysis. \emph{Evolutionary Ecology}, \bold{26}(5), 1253--1274. \verb{https://doi.org/10.1007/s10682-012-9555-5}
}
\examples{
### copy data into 'dat' and examine data
dat <- dat.lim2014$o_o_unadj
head(dat)

### calculate r-to-z transformed correlations and corresponding sampling variances
dat <- escalc(measure="ZCOR", ri=ri, ni=ni, data=dat)

\dontrun{
### load 'ape' package
require(ape)

### copy tree to 'tree'
tree <- dat.lim2014$o_o_unadj_tree

### compute branch lengths
tree <- compute.brlen(tree)

### compute phylogenetic correlation matrix
A <- vcv(tree, corr=TRUE)

### make copy of the species variable
dat$species.phy <- dat$species

### create effect size id variable
dat$esid <- 1:nrow(dat)

### fit multilevel phylogenetic meta-analytic model
res <- rma.mv(yi, vi,
   random = list(~ 1 | article, ~ 1 | esid, ~ 1 | species, ~ 1 | species.phy),
   R=list(species.phy=A), data=dat)
res
}
}
\keyword{datasets}
