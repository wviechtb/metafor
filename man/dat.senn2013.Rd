\name{dat.senn2013}
\docType{data}
\alias{dat.senn2013}
\title{Studies on the Effectiveness of Glucose-Lowering Agents}
\description{Results from 26 trials examining the effectiveness of glucose-lowering agents in patients with type 2 diabetes}
\usage{dat.senn2013}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{study}     \tab \code{character} \tab (first) author and year of study \cr
\bold{ni}        \tab \code{numeric}   \tab sample size of the study arm \cr
\bold{treatment} \tab \code{character} \tab treatment given \cr
\bold{comment}   \tab \code{character} \tab whether figures given are based on raw values at outcome or on change from baseline \cr
\bold{mi}        \tab \code{numeric}   \tab mean score \cr
\bold{sdi}       \tab \code{numeric}   \tab standard deviation
}
}
\details{
   The dataset includes the results from 26 randomized controlled trials examining the effectiveness of adding various oral glucose-lowering agents to a baseline sulfonylurea therapy in patients with type 2 diabetes. The outcome measured in the studies was either the mean HbA1c level at follow-up or the mean change in HbA1c level from baseline to follow-up. A total of 10 different treatment types were examined in these studies: acarbose, benfluorex, metformin, miglitol, pioglitazone, placebo, rosiglitazone, sitagliptin, sulfonylurea alone, and vildagliptin. One study included three treatment arms (Willms, 1999), while the rest of the studies included two treatment arms (hence, the dataset includes the results from 53 treatment arms).

   The data can be used for a network meta-analysis, either using an arm-based or a contrast-based model. See \sQuote{Examples} below.
}
\source{
   Law, M., Jackson, D., Turner, R., Rhodes, K., & Viechtbauer, W. (2016). Two new methods to fit models for network meta-analysis with random inconsistency effects. \emph{BMC Medical Research Methodology}, \bold{16}, 87. \verb{https://doi.org/10.1186/s12874-016-0184-5}

   Senn, S., Gavini, F., Magrez, D., & Scheen, A. (2013). Issues in performing a network meta-analysis. \emph{Statistical Methods in Medical Research}, \bold{22}(2), 169--189. \verb{https://doi.org/10.1177/0962280211432220}
}
\examples{
### copy data into 'dat' and examine data
dat <- dat.senn2013
dat

### create network graph ('igraph' package must be installed)
\dontrun{
require(igraph)
pairs <- data.frame(do.call(rbind,
   sapply(split(dat$treatment, dat$study), function(x) t(combn(x,2)))), stringsAsFactors=FALSE)
pairs$X1 <- factor(pairs$X1, levels=sort(unique(dat$treatment)))
pairs$X2 <- factor(pairs$X2, levels=sort(unique(dat$treatment)))
tab <- table(pairs[,1], pairs[,2])
tab # adjacency matrix
g <- graph_from_adjacency_matrix(tab, mode = "plus", weighted=TRUE, diag=FALSE)
plot(g, edge.curved=FALSE, edge.width=E(g)$weight, layout=layout_as_star(g, center="placebo"),
     vertex.size=45, vertex.color="lightgray", vertex.label.color="black", vertex.label.font=2)}

### table of studies versus treatments examined
print(addmargins(table(dat$study, dat$treatment)), zero.print="")

### table of frequencies with which treatment pairs were studied
print(as.table(crossprod(table(dat$study, dat$treatment))), zero.print="")

### add means and sampling variances of the means to the dataset
dat <- escalc(measure="MN", mi=mi, sdi=sdi, ni=ni, data=dat)

### turn treatment variable into factor and set reference level
dat$treatment <- relevel(factor(dat$treatment), ref="placebo")

### add a space before each level (this makes the output a bit more legible)
levels(dat$treatment) <- paste0(" ", levels(dat$treatment))

### network meta-analysis using an arm-based fixed-effects model with fixed study effects
res.fe <- rma.mv(yi, vi, mods = ~ study + treatment - 1, data=dat, slab=paste0(study, treatment))
res.fe

### test if treatment factor as a whole is significant
anova(res.fe, btt="treatment")

### forest plot of the contrast estimates (treatments versus placebos)
forest(tail(coef(res.fe), 9), tail(diag(vcov(res.fe)), 9), slab=levels(dat$treatment)[-1],
       xlim=c(-2.5, 2.0), alim=c(-1.5, 0.5), psize=1, xlab="Estimate", header="Treatment")

### weight matrix for the estimation of the fixed effects (leaving out study effects)
w <- t(tail(vcov(res.fe) \%*\% t(model.matrix(res.fe)) \%*\% weights(res.fe, type="matrix"), 9))
rownames(w) <- res.fe$slab

### create shade plot for the diabetes network with placebo as the reference treatment
### negative values in blue shades, positive values in red shades
cols <- colorRampPalette(c("blue", "gray95", "red"))(9)
heatmap(w, Rowv=NA, Colv=NA, scale="none", margins=c(6,11), col=cols,
        cexRow=.7, cexCol=1, labCol=levels(dat$treatment)[-1])

### network meta-analysis using an arm-based random-effects model with fixed study effects
### by setting rho=1/2, tau^2 reflects the amount of heterogeneity for all treatment comparisons
res.re <- rma.mv(yi, vi, mods = ~ study + treatment - 1, random = ~ treatment | study, rho=1/2,
                 data=dat, slab=paste0(study, treatment))
res.re

### test if treatment factor as a whole is significant
anova(res.re, btt="treatment")

### forest plot of the contrast estimates (treatments versus placebos)
forest(tail(coef(res.re), 9), tail(diag(vcov(res.re)), 9), slab=levels(dat$treatment)[-1],
       xlim=c(-3.0, 2.5), alim=c(-1.5, 0.5), psize=1, xlab="Estimate", header="Treatment")

### compute contribution of each study to the overall Q-test value
qi <- sort(by((resid(res.fe) / sqrt(dat$vi))^2, dat$study, sum))

### check that the values add up
sum(qi)
res.fe$QE

### plot the values
s <- length(qi)
par(mar=c(5,10,2,1))
plot(qi, 1:s, pch=19, xaxt="n", yaxt="n", xlim=c(0,40), xlab="Chi-Square Contribution", ylab="")
axis(side=1)
axis(side=2, at=1:s, labels=names(qi), las=1, tcl=0)
segments(rep(0,s), 1:s, qi, 1:s)

############################################################################

### restructure dataset to a contrast-based format
dat <- dat.senn2013[c(1,4:2,5:6)] # reorder variables first
dat <- to.wide(dat, study="study", grp="treatment", ref="placebo", grpvars=4:6)
dat

### calculate mean difference and corresponding sampling variance for each treatment comparison
dat <- escalc(measure="MD", m1i=mi.1, sd1i=sdi.1, n1i=ni.1,
                            m2i=mi.2, sd2i=sdi.2, n2i=ni.2, data=dat)
dat

### calculate the variance-covariance matrix of the mean differences for the multitreatment studies
calc.v <- function(x) {
   v <- matrix(x$sdi.2[1]^2 / x$ni.2[1], nrow=nrow(x), ncol=nrow(x))
   diag(v) <- x$vi
   v
}
V <- bldiag(lapply(split(dat, dat$study), calc.v))

### add contrast matrix to dataset
dat <- contrmat(dat, grp1="treatment.1", grp2="treatment.2")
dat

### network meta-analysis using a contrast-based random-effects model
### by setting rho=1/2, tau^2 reflects the amount of heterogeneity for all treatment comparisons
### the treatment left out (placebo) becomes the reference level for the treatment comparisons
res <- rma.mv(yi, V, mods = ~ acarbose + benfluorex + metformin + miglitol + pioglitazone +
                              rosiglitazone + sitagliptin + sulfonylurea + vildagliptin - 1,
              random = ~ comp | study, rho=1/2, data=dat)
res

### estimate all pairwise differences between treatments (using the 'multcomp' package)
if (require(multcomp)) {
   contr <- contrMat(setNames(rep(1,res$p), colnames(res$X)), type="Tukey")
   sav <- predict(res, newmods=contr)
   sav[["slab"]] <- rownames(contr)
   sav
}

### fit random inconsistency effects model (see Law et al., 2016)
res <- rma.mv(yi, V, mods = ~ acarbose + benfluorex + metformin + miglitol + pioglitazone +
                              rosiglitazone + sitagliptin + sulfonylurea + vildagliptin - 1,
              random = list(~ comp | study, ~ comp | design), rho=1/2, phi=1/2, data=dat)
res
}
\keyword{datasets}
