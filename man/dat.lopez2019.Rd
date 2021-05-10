\name{dat.lopez2019}
\docType{data}
\alias{dat.lopez2019}
\title{Studies on the Effectiveness of CBT for Depression}
\description{Results from 76 studies examining the effectiveness of cognitive behavioral therapy (CBT) for depression in adults.}
\usage{dat.lopez2019}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{study}     \tab \code{character} \tab (first) author and year of study \cr
\bold{treatment} \tab \code{character} \tab treatment provided (see \sQuote{Details}) \cr
\bold{scale}     \tab \code{character} \tab scale used to measure depression symptoms \cr
\bold{n}         \tab \code{numeric}   \tab group size \cr
\bold{diff}      \tab \code{numeric}   \tab standardized mean change \cr
\bold{se}        \tab \code{numeric}   \tab corresponding standard error \cr
\bold{group}     \tab \code{numeric}   \tab type of therapy (0 = individual, 1 = group therapy) \cr
\bold{tailored}  \tab \code{numeric}   \tab whether the intervention was tailored to each patient (0 = no, 1 = yes) \cr
\bold{sessions}  \tab \code{numeric}   \tab number of sessions \cr
\bold{length}    \tab \code{numeric}   \tab average session length (in minutes) \cr
\bold{intensity} \tab \code{numeric}   \tab product of sessions and length \cr
\bold{multi}     \tab \code{numeric}   \tab intervention included multimedia elements (0 = no, 1 = yes) \cr
\bold{cog}       \tab \code{numeric}   \tab intervention included cognitive techniques (0 = no, 1 = yes) \cr
\bold{ba}        \tab \code{numeric}   \tab intervention included behavioral activation (0 = no, 1 = yes) \cr
\bold{psed}      \tab \code{numeric}   \tab intervention included psychoeducation (0 = no, 1 = yes) \cr
\bold{home}      \tab \code{numeric}   \tab intervention included homework (0 = no, 1 = yes) \cr
\bold{prob}      \tab \code{numeric}   \tab intervention included problem solving (0 = no, 1 = yes) \cr
\bold{soc}       \tab \code{numeric}   \tab intervention included social skills training (0 = no, 1 = yes) \cr
\bold{relax}     \tab \code{numeric}   \tab intervention included relaxation (0 = no, 1 = yes) \cr
\bold{goal}      \tab \code{numeric}   \tab intervention included goal setting (0 = no, 1 = yes) \cr
\bold{final}     \tab \code{numeric}   \tab intervention included a final session (0 = no, 1 = yes) \cr
\bold{mind}      \tab \code{numeric}   \tab intervention included mindfulness (0 = no, 1 = yes) \cr
\bold{act}       \tab \code{numeric}   \tab intervention included acceptance and commitment therapy (0 = no, 1 = yes)
}
}
\details{
   The dataset includes the results from 76 studies examining the effectiveness of cognitive behavioral therapy (CBT) for treating depression in adults. Studies included two or more of the following treatments/conditions:
   \enumerate{
   \item treatment as usual (TAU),
   \item no treatment,
   \item wait list,
   \item psychological or attention placebo,
   \item face-to-face CBT,
   \item multimedia CBT,
   \item hybrid CBT (i.e., multimedia CBT with one or more face-to-face sessions).
   }
   Multimedia CBT was defined as CBT delivered via self-help books, audio/video recordings, telephone, computer programs, apps, e-mail, or text messages.

   Variable \code{diff} is the standardized mean change within each group, with negative values indicating a decrease in depression symptoms.
}
\source{
   Personal communication.
}
\references{
   López-López, J. A., Davies, S. R., Caldwell, D. M., Churchill, R., Peters, T. J., Tallon, D., Dawson, S., Wu, Q., Li, J., Taylor, A., Lewis, G., Kessler, D. S., Wiles, N., & Welton, N. J. (2019). The process and delivery of CBT for depression in adults: A systematic review and network meta-analysis. \emph{Psychological Medicine}, \bold{49}(12), 1937--1947. \verb{https://doi.org/10.1017/S003329171900120X}
}
\examples{
### copy data into 'dat' and examine data
dat <- dat.lopez2019
dat[1:10,1:6]

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
plot(g, edge.curved=FALSE, edge.width=E(g)$weight/2,
     layout=layout_in_circle(g, order=c("Wait list", "No treatment", "TAU", "Multimedia CBT",
                                        "Hybrid CBT", "F2F CBT", "Placebo")),
     vertex.size=45, vertex.color="lightgray", vertex.label.color="black", vertex.label.font=2)}

### restructure data into wide format
dat <- to.wide(dat, study="study", grp="treatment", ref="TAU",
               grpvars=c("diff","se","n"), postfix=c("1","2"))

### compute contrasts between treatment pairs and corresponding sampling variances
dat$yi <- with(dat, diff1 - diff2)
dat$vi <- with(dat, se1^2 + se2^2)

### calculate the variance-covariance matrix for multitreatment studies
calc.v <- function(x) {
   v <- matrix(x$se2[1]^2, nrow=nrow(x), ncol=nrow(x))
   diag(v) <- x$vi
   v
}
V <- bldiag(lapply(split(dat, dat$study), calc.v))

### add contrast matrix to the dataset
dat <- contrmat(dat, grp1="treatment1", grp2="treatment2")

### network meta-analysis using a contrast-based random-effects model
### by setting rho=1/2, tau^2 reflects the amount of heterogeneity for all treatment comparisons
### the treatment left out (TAU) becomes the reference level for the treatment comparisons
res <- rma.mv(yi, V, data=dat,
         mods = ~ No.treatment + Wait.list + Placebo + F2F.CBT + Hybrid.CBT + Multimedia.CBT - 1,
         random = ~ comp | study, rho=1/2)
res

### forest plot of the contrast estimates (treatments versus TAU)
forest(coef(res), diag(vcov(res)), slab=sub(".", " ", names(coef(res)), fixed=TRUE),
       xlim=c(-5,5), alim=c(-3,3), psize=1, header="Treatment",
       xlab="Difference in Standardized Mean Change (compared to TAU)")

### fit random inconsistency effects model
res <- rma.mv(yi, V, data=dat,
         mods = ~ No.treatment + Wait.list + Placebo + F2F.CBT + Hybrid.CBT + Multimedia.CBT - 1,
         random = list(~ comp | study, ~ comp | design), rho=1/2, phi=1/2)
res
}
\keyword{datasets}
