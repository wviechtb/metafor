\name{dat.hasselblad1998}
\docType{data}
\alias{dat.hasselblad1998}
\title{Studies on the Effectiveness of Counseling for Smoking Cessation}
\description{Results from 24 studies on the effectiveness of various counseling types for smoking cessation.}
\usage{dat.hasselblad1998}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{id}         \tab \code{numeric}    \tab id number for each treatment arm \cr
\bold{study}      \tab \code{numeric}    \tab study id number \cr
\bold{authors}    \tab \code{character}  \tab study author(s) \cr
\bold{year}       \tab \code{numeric}    \tab publication year \cr
\bold{trt}        \tab \code{character}  \tab intervention group \cr
\bold{xi}         \tab \code{numeric}    \tab number of individuals abstinent \cr
\bold{ni}         \tab \code{numeric}    \tab number of individuals in group
}
}
\details{
   The dataset includes the results from 24 studies on the effectiveness of various counseling types for smoking cessation (i.e., self-help, individual counseling, group counseling, and no contact). The dataset indicates the total number of individuals within each study arm and the number that were abstinent from 6 to 12 months. The majority of the studies compared two interventions types against each other, while 2 studies compared three types against each other simultaneously.

   The data can be used for a \sQuote{network meta-analysis} (also called \sQuote{mixed treatment comparison meta-analysis}). The code below shows how such an analysis can be conducted using an arm-based and a contrast-based model (see Salanti et al., 2008, for more details).
}
\source{
   Hasselblad, V. (1998). Meta-analysis of multitreatment studies. \emph{Medical Decision Making}, \bold{18}(1), 37--43. \verb{https://doi.org/10.1177/0272989X9801800110}
}
\references{
   Gleser, L. J., & Olkin, I. (2009). Stochastically dependent effect sizes. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta-analysis} (2nd ed., pp. 357--376). New York: Russell Sage Foundation.

   Law, M., Jackson, D., Turner, R., Rhodes, K., & Viechtbauer, W. (2016). Two new methods to fit models for network meta-analysis with random inconsistency effects. \emph{BMC Medical Research Methodology}, \bold{16}, 87. \verb{https://doi.org/10.1186/s12874-016-0184-5}

   Salanti, G., Higgins, J. P. T., Ades, A. E., & Ioannidis, J. P. A. (2008). Evaluation of networks of randomized trials. \emph{Statistical Methods in Medical Research}, \bold{17}(3), 279--301. \verb{https://doi.org/10.1177/0962280207080643}
}
\examples{
### copy data into 'dat'
dat <- dat.hasselblad1998
dat

### create network graph ('igraph' package must be installed)
\dontrun{
require(igraph)
pairs <- data.frame(do.call(rbind,
   sapply(split(dat$trt, dat$study), function(x) t(combn(x,2)))), stringsAsFactors=FALSE)
lvls <- c("no_contact", "self_help", "ind_counseling", "grp_counseling")
pairs$X1 <- factor(pairs$X1, levels=lvls)
pairs$X2 <- factor(pairs$X2, levels=lvls)
tab <- table(pairs[,1], pairs[,2])
tab # adjacency matrix
g <- graph_from_adjacency_matrix(tab, mode = "plus", weighted=TRUE, diag=FALSE)
vertex_attr(g, "name") <- c("No Contact", "Self-Help",
                            "Individual\nCounseling", "Group\nCounseling")
plot(g, edge.curved=FALSE, edge.width=E(g)$weight, layout=layout_on_grid,
     vertex.size=45, vertex.color="lightgray", vertex.label.color="black", vertex.label.font=2)}

### calculate log odds for each study arm
dat <- escalc(measure="PLO", xi=xi, ni=ni, add=1/2, to="all", data=dat)
dat

### convert trt variable to factor with desired ordering of levels
dat$trt <- factor(dat$trt, levels=c("no_contact", "self_help", "ind_counseling", "grp_counseling"))

### add a space before each level (this makes the output a bit more legible)
levels(dat$trt) <- paste0(" ", levels(dat$trt))

### network meta-analysis using an arm-based model with fixed study effects
### by setting rho=1/2, tau^2 reflects the amount of heterogeneity for all treatment comparisons
res <- rma.mv(yi, vi, mods = ~ factor(study) + trt - 1,
              random = ~ trt | study, rho=1/2, data=dat, btt="trt")
res

### all pairwise odds ratios of interventions versus no contact
predict(res, newmods=cbind(matrix(0, nrow=3, ncol=24), diag(3)),
        intercept=FALSE, transf=exp, digits=2)

### all pairwise odds ratios comparing interventions (ic vs sh, gc vs sh, and gc vs ic)
predict(res, newmods=cbind(matrix(0, nrow=3, ncol=24), rbind(c(-1,1,0), c(-1,0,1), c(0,-1,1))),
        intercept=FALSE, transf=exp, digits=2)

### forest plot of ORs of interventions versus no contact
dev.new(width=7, height=4)
par(mar=c(5,4,1,2))
forest(c(0,res$beta[25:27]), sei=c(0,res$se[25:27]), psize=1, xlim=c(-3,4), digits=c(2,1), efac=2,
       slab=c("No Contact", "Self-Help", "Individual Counseling", "Group Counseling"),
       atransf=exp, at=log(c(.5, 1, 2, 4, 8)), xlab="Odds Ratio for Intervention vs. No Contact",
       header=c("Intervention", "Odds Ratio [95\% CI]"))

############################################################################

### restructure dataset to a contrast-based format
dat <- to.wide(dat.hasselblad1998, study="study", grp="trt", ref="no_contact", grpvars=6:7)

### calculate log odds ratios for each treatment comparison
dat <- escalc(measure="OR", ai=xi.1, n1i=ni.1,
                            ci=xi.2, n2i=ni.2, add=1/2, to="all", data=dat)
dat

### calculate the variance-covariance matrix of the log odds ratios for multitreatment studies
### see Gleser & Olkin (2009), equation (19.11), for the covariance equation
calc.v <- function(x) {
   v <- matrix(1/(x$xi.2[1]+1/2) + 1/(x$ni.2[1] - x$xi.2[1] + 1/2), nrow=nrow(x), ncol=nrow(x))
   diag(v) <- x$vi
   v
}
V <- bldiag(lapply(split(dat, dat$study), calc.v))

### add contrast matrix to dataset
dat <- contrmat(dat, grp1="trt.1", grp2="trt.2")
dat

### network meta-analysis using a contrast-based random-effects model
### by setting rho=1/2, tau^2 reflects the amount of heterogeneity for all treatment comparisons
res <- rma.mv(yi, V, mods = ~ self_help + ind_counseling + grp_counseling - 1,
              random = ~ comp | study, rho=1/2, data=dat)
res

### predicted odds ratios of interventions versus no contact
predict(res, newmods=diag(3), transf=exp, digits=2)

### fit random inconsistency effects model (see Law et al., 2016)
res <- rma.mv(yi, V, mods = ~ self_help + ind_counseling + grp_counseling - 1,
              random = list(~ comp | study, ~ comp | design), rho=1/2, phi=1/2, data=dat)
res
}
\keyword{datasets}
