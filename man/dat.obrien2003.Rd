\name{dat.obrien2003}
\docType{data}
\alias{dat.obrien2003}
\title{Studies on the Relationship Between BMI and Risk of Preeclampsia}
\description{Results from 13 studies on the relationship between maternal body mass index (BMI) and the risk of preeclampsia.}
\usage{dat.obrien2003}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{study}      \tab \code{numeric}    \tab study id \cr
\bold{author}     \tab \code{character}  \tab (first) author of the study \cr
\bold{year}       \tab \code{numeric}    \tab publication year \cr
\bold{ref}        \tab \code{numeric}    \tab reference number \cr
\bold{ch}         \tab \code{character}  \tab exclusion due to chronic hypertension (yes/no) \cr
\bold{dm}         \tab \code{character}  \tab exclusion due to diabetes mellitus (yes/no) \cr
\bold{mg}         \tab \code{character}  \tab exclusion due to multiple gestation (yes/no) \cr
\bold{bmi.lb}     \tab \code{numeric}    \tab lower bound of the BMI interval \cr
\bold{bmi.ub}     \tab \code{numeric}    \tab upper bound of the BMI interval \cr
\bold{bmi}        \tab \code{numeric}    \tab midpoint of the BMI interval \cr
\bold{cases}      \tab \code{numeric}    \tab number of preeclampsia cases in the BMI group \cr
\bold{total}      \tab \code{numeric}    \tab number of individuals in the BMI group
}
}
\details{
   The dataset includes the results from 13 studies examining the relationship between maternal body mass index (BMI) and the risk of preeclampsia. For each study, results are given in terms of the number of preeclampsia cases within two or more groups defined by the lower and upper BMI bounds as shown in the dataset (\code{NA} means that the interval is either open to the left or right). The \code{bmi} variable is the interval midpoint as defined by O'Brien et al. (2003).
}
\source{
   O'Brien, T. E., Ray, J. G., & Chan, W.-S. (2003). Maternal body mass index and the risk of preeclampsia: A systematic overview. \emph{Epidemiology}, \bold{14}(3), 368--374. \verb{https://doi.org/10.1097/00001648-200305000-00020}
}
\examples{
### copy data into 'dat'
dat <- dat.obrien2003
dat

### restructure the data into a wide format
dat2 <- to.wide(dat, study="study", grp="grp", ref=1, grpvars=c("bmi","cases","total"),
                addid=FALSE, adddesign=FALSE, postfix=c(1,2))
dat2

### calculate log risk ratios and corresponding sampling variances
dat2 <- escalc(measure="RR", ai=cases1, n1i=total1, ci=cases2, n2i=total2, data=dat2)
dat2

### forest plot of the risk ratios
dd <- c(0,diff(dat2$study))
dd[dd > 0] <- 1
rows <- (1:nrow(dat2)) + cumsum(dd)
rows <- 1 + max(rows) - rows
slabs <- mapply(function(x,y,z) as.expression(bquote(.(x)^.(y)~.(z))),
                dat2$author, dat2$ref, dat2$year)
with(dat2, forest(yi, vi, header=TRUE, slab=slabs, xlim=c(-7,5.5), fonts="mono", cex=0.8,
   psize=1, pch=19, efac=0, rows=rows, ylim=c(0,max(rows)+3), yaxs="i",
   atransf=exp, at=log(c(.05,0.1,0.2,0.5,1,2,5,10,20)), ilab=comp, ilab.xpos=-4, ilab.pos=4))
text(-4.4, max(rows)+2, "Comparison", font=2, cex=0.8, pos=4)

### within-study mean center the BMI variable
dat$bmicent <- dat$bmi - ave(dat$bmi, dat$study)

### compute the proportion of preeclampsia cases and corresponding sampling variances
dat <- escalc(measure="PR", xi=cases, ni=total, data=dat)

### convert the proportions to percentages (and convert the variances accordingly)
dat$yi <- dat$yi*100
dat$vi <- dat$vi*100^2

### fit multilevel meta-regression model to examine the relationship between the
### (centered) BMI variable and the risk of preeclampsia
res <- rma.mv(yi, vi, mods = ~ bmicent, random = ~ 1 | study/grp, data=dat)
res

### draw scatterplot with regression line
res$slab <- dat$ref
regplot(res, xlab=expression("Within-Study Mean Centered BMI"~(kg/m^2)),
        ylab="Preeclampsia Prevalence (\%)", las=1, bty="l",
        at=seq(0,18,by=2), olim=c(0,100), psize=2, bg="gray90",
        label=TRUE, offset=0, labsize=0.6)
}
\keyword{datasets}
