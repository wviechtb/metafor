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
\bold{ni}         \tab \code{numeric}    \tab number of individuals in group \cr
}
}
\details{
   The dataset includes the results from 24 studies on the effectiveness of various counseling types for smoking cessation (i.e., self-help, individual counseling, group counseling, and no contact). The dataset indicates the total number of individuals within each study arm and the number that were abstinent from 6 to 12 months. The majority of the studies compared two interventions types against each other, while 2 studies compared three types against each other simultaneously.

   The data can be used for a \sQuote{network meta-analysis} (also called \sQuote{mixed treatment comparison meta-analysis}). The code below shows how such an analysis can be conducted using an arm-based and a contrast-based model (see Salanti et al., 2008, for more details).
}
\source{
   Hasselblad, V. (1998). Meta-analysis of multitreatment studies. \emph{Medical Decision Making}, \bold{18}, 37--43.
}
\references{
   Gleser, L. J., & Olkin, I. (2009). Stochastically dependent effect sizes. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta-analysis} (2nd ed., pp. 357--376). New York: Russell Sage Foundation.

   Salanti, G., Higgins, J. P. T., Ades, A. E., & Ioannidis, J. P. A. (2008). Evaluation of networks of randomized trials. \emph{Statistical Methods in Medical Research}, \bold{17}, 279--301.
}
\examples{
### load data
dat <- get(data(dat.hasselblad1998))

### calculate log odds for each study arm
dat <- escalc(measure="PLO", xi=xi, ni=ni, add=1/2, to="all", data=dat)

### create network graph ('plyr' and 'igraph' packages must be installed)
\dontrun{
require(plyr)
require(igraph)
pairs <- do.call(rbind, sapply(split(dat$trt, dat$study), function(x) t(combn(x,2))))
pairs <- ddply(data.frame(pairs), .(X1, X2), count)
g <- graph.edgelist(as.matrix(pairs[,1:2]), directed=FALSE)
plot(g, edge.curved=FALSE, edge.width=pairs$freq, vertex.label.dist=.7,
vertex.label=c("Individual\nCounseling", "Group\nCounseling", "No Contact", "Self-Help"))}

### convert trt variable to factor with desired ordering of levels
dat$trt <- factor(dat$trt, levels=c("no_contact", "self_help", "ind_counseling", "grp_counseling"))

### add a space before each level (this makes the output a bit more legible)
levels(dat$trt) <- paste0(" ", levels(dat$trt))

### network meta-analysis using an arm-based model with fixed study effects
### by setting rho=1/2, tau^2 reflects the amount of heterogeneity for all treatment comparisons
res <- rma.mv(yi, vi, mods = ~ factor(study) + trt - 1,
              random = ~ trt | study, rho=1/2, data=dat, btt=25:27)
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
       atransf=exp, at=log(c(.5, 1, 2, 4, 8)), xlab="Odds Ratio for Intervention vs. No Contact")
text(-3, 6, "Intervention", pos=4)
text( 4, 6, "Odds Ratio [95\% CI]",  pos=2)

### restructure dataset to a contrast-based format
dat.c <- lapply(split(dat.hasselblad1998, dat.hasselblad1998$study),
                function(x) cbind(x[rep(1,nrow(x)-1),], x[-1,5:7]))
dat.c <- do.call(rbind, dat.c)
dat.c <- dat.c[,c(1:5,8,6:7,9:10)]
names(dat.c)[5:10] <- c("trt1", "trt2", "ai", "n1i", "ci", "n2i")
rownames(dat.c) <- 1:nrow(dat.c)
dat.c$id <- 1:nrow(dat.c)
dat.c

### calculate log odds ratios for each treatment comparison
dat.c <- escalc(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, add=1/2, to="all", data=dat.c)
dat.c

### calculate the variance-covariance matrix of the log odds ratios for multitreatment studies
### see Gleser & Olkin (2009), equation (19.11), for the covariance equation
calc.v <- function(x) {
   v <- matrix(1/(x$ai[1]+1/2) + 1/(x$n1i[1] - x$ai[1] + 1/2), nrow=nrow(x), ncol=nrow(x))
   diag(v) <- x$vi
   v
}
V <- bldiag(lapply(split(dat.c, dat.c$study), calc.v))

### convert trt1 and trt2 variables to factors with desired ordering of levels
lvls <- c("no_contact", "self_help", "ind_counseling", "grp_counseling")
dat.c$trt1 <- factor(dat.c$trt1, levels=lvls)
dat.c$trt2 <- factor(dat.c$trt2, levels=lvls)

### create variables to indicate the contrast examined
dat.c <- cbind(dat.c, model.matrix(~ dat.c$trt1 - 1) - model.matrix(~ dat.c$trt2 - 1))
names(dat.c)[(ncol(dat.c)-3):ncol(dat.c)] <- lvls

### network meta-analysis using a contrast-based random-effects model
### by setting rho=1/2, tau^2 reflects the amount of heterogeneity for all treatment comparisons
res <- rma.mv(yi, V, mods = ~ self_help + ind_counseling + grp_counseling - 1,
              random = ~ factor(id) | study, rho=1/2, data=dat.c)
res
}
\keyword{datasets}
