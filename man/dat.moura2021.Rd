\name{dat.moura2021}
\docType{data}
\alias{dat.moura2021}
\title{Studies on Assortative Mating}
\description{Results from 457 studies on assortative mating in various species.}
\usage{
dat.moura2021
}
\format{The object is a list containing a data frame called \code{dat} that contains the following columns and a phylogenetic tree called \code{tree}:
\tabular{lll}{
\bold{study.id}               \tab \code{character} \tab study ID \cr
\bold{effect.size.id}         \tab \code{numeric}   \tab effect size ID \cr
\bold{species}                \tab \code{character} \tab species \cr
\bold{species.id}             \tab \code{character} \tab species ID (as in the Open Tree of Life reference taxonomy) \cr
\bold{subphylum}              \tab \code{character} \tab the subphyla of the species \cr
\bold{phylum}                 \tab \code{character} \tab the phyla of the species \cr
\bold{assortment.trait}       \tab \code{character} \tab the measure of body size \cr
\bold{trait.dimensions}       \tab \code{character} \tab dimensionality of the measure \cr
\bold{field.collection}       \tab \code{character} \tab whether data were collected in the field \cr
\bold{publication.year}       \tab \code{numeric}   \tab publication year of the study \cr
\bold{pooled.data}            \tab \code{character} \tab whether data were pooled either spatially and/or temporally \cr
\bold{spatially.pooled.data}  \tab \code{character} \tab whether data were pooled spatially \cr
\bold{temporally.pooled.data} \tab \code{character} \tab whether data were pooled temporally \cr
\bold{ri}                     \tab \code{numeric}   \tab correlation coefficient \cr
\bold{ni}                     \tab \code{numeric}   \tab sample size
}
Blah.
}
\details{
   The 457 studies included in this dataset provide 1828 correlation coefficients describing the similarity in some measure of body size in mating couples in 341 different species.
}
\source{
   Rios Moura, R., Oliveira Gonzaga, M., Silva Pinto, N., Vasconcellos-Neto, J., & Requena, G. S. (2021). Assortative mating in space and time: Patterns and biases. \emph{Ecology Letters}, \bold{24}(5), 1089--1102. \verb{https://doi.org/10.1111/ele.13690}
}
\references{
   Hadfield, J. D., & Nakagawa, S. (2010). General quantitative genetic methods for comparative biology: Phylogenies, taxonomies and multi-trait models for continuous and categorical characters. \emph{Journal of Evolutionary Biology}, \bold{23}(3), 494--508. \verb{https://doi.org/10.1111/j.1420-9101.2009.01915.x}

   Nakagawa, S., & Santos, E. S. A. (2012). Methodological issues and advances in biological meta-analysis. \emph{Evolutionary Ecology}, \bold{26}(5), 1253--1274. \verb{https://doi.org/10.1007/s10682-012-9555-5}
}
\examples{
### copy data into 'dat' and examine data
dat <- dat.moura2021$dat
head(dat)

### calculate r-to-z transformed correlations and corresponding sampling variances
dat <- escalc(measure="ZCOR", ri=ri, ni=ni, data=dat)

\dontrun{
### load 'ape' package
require(ape)

### copy tree to 'tree'
tree <- dat.moura2021$tree

### turn tree into an ultrametric one
tree <- compute.brlen(tree)

### compute phylogenetic correlation matrix
A <- vcv(tree, corr=TRUE)

### make copy of the species.id variable
dat$species.id.phy <- dat$species.id

### fit multilevel phylogenetic meta-analytic model
res <- rma.mv(yi, vi,
   random = list(~ 1 | study.id, ~ 1 | effect.size.id, ~ 1 | species.id, ~ 1 | species.id.phy),
   R=list(species.id.phy=A), data=dat)
res

### examine if spatial and/or temporal pooling of data tends to yield larger correlations
res <- rma.mv(yi, vi,
   mods = ~ spatially.pooled.data * temporally.pooled.data,
   random = list(~ 1 | study.id, ~ 1 | effect.size.id, ~ 1 | species.id, ~ 1 | species.id.phy),
   R=list(species.id.phy=A), data=dat)
res

### estimated average correlation without pooling, when pooling spatially,
### when pooling temporally, and when pooling spatially and temporally
predict(res, newmods = rbind(c(0,0,0),c(1,0,0),c(0,1,0),c(1,1,1)), transf=transf.ztor, digits=2)
}
}
\keyword{datasets}
