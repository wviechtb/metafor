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
   Senn, S., Gavini, F., Magrez, D., & Scheen, A. (2013). Issues in performing a network meta-analysis. \emph{Statistical Methods in Medical Research}, \bold{22}, 169--189.
}
\examples{
### copy data into 'dat' and examine data
dat <- dat.senn2013
dat

### create network graph ('plyr' and 'igraph' packages must be installed)
\dontrun{
require(plyr)
require(igraph)
pairs <- do.call(rbind, sapply(split(dat$treatment, dat$study), function(x) t(combn(x,2))))
pairs <- ddply(data.frame(pairs), .(X1, X2), count)
g <- graph.edgelist(as.matrix(pairs[,1:2]), directed=FALSE)
plot(g, edge.curved=FALSE, edge.width=pairs$freq, vertex.label.dist=.8)}

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
dev.new(width=7, height=5)
forest(tail(coef(res.fe), 9), tail(diag(vcov(res.fe)), 9), slab=levels(dat$treatment)[-1],
       xlim=c(-2.5, 2.0), alim=c(-1.5, 0.5), psize=1, xlab="Estimate")
text(-2.5, 11, "Treatment", pos=4, font=2)
text( 2.0, 11, "Estimate [95\% CI]", pos=2, font=2)

### weight matrix for the estimation of the fixed effects (leaving out study effects)
w <- t(tail(vcov(res.fe) \%*\% t(model.matrix(res.fe)) \%*\% weights(res.fe, type="matrix"), 9))
rownames(w) <- res.fe$slab

### create shade plot for the diabetes network with placebo as the reference treatment
### negative values in blue shades, positive values in red shades
dev.new()
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
       xlim=c(-3.0, 2.5), alim=c(-1.5, 0.5), psize=1, xlab="Estimate")
text(-3.0, 11, "Treatment", pos=4, font=2)
text( 2.5, 11, "Estimate [95\% CI]", pos=2, font=2)

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
res.cb <- rma.mv(yi, V, mods = ~ acarbose + benfluorex + metformin + miglitol + pioglitazone +
                                 rosiglitazone + sitagliptin + sulfonylurea + vildagliptin - 1,
                 random = ~ factor(id) | study, rho=1/2, data=dat)
res.cb

### estimate all pairwise differences between treatments (using the 'multcomp' package)
if (require(multcomp)) {
   contr <- contrMat(setNames(rep(1,res.cb$p), colnames(res.cb$X)), type="Tukey")
   sav <- predict(res.cb, newmods=contr)
   sav[["slab"]] <- rownames(contr)
   sav
}
}
\keyword{datasets}
