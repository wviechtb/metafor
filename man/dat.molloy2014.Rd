\name{dat.molloy2014}
\docType{data}
\alias{dat.molloy2014}
\title{Studies on the Relationship between Conscientiousness and Medication Adherence}
\description{Results from 16 studies on the correlation between conscientiousness and medication adherence.}
\usage{dat.molloy2014}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{authors}   \tab \code{character} \tab study authors \cr
\bold{year}      \tab \code{numeric}   \tab publication year \cr
\bold{ni}        \tab \code{numeric}   \tab sample size of the study \cr
\bold{ri}        \tab \code{numeric}   \tab observed correlation \cr
\bold{controls}  \tab \code{character} \tab number of variables controlled for \cr
\bold{design}    \tab \code{character} \tab whether a cross-sectional or prospective design was used \cr
\bold{a_measure} \tab \code{character} \tab type of adherence measure (self-report or other) \cr
\bold{c_measure} \tab \code{character} \tab type of conscientiousness measure (NEO or other) \cr
\bold{meanage}   \tab \code{numeric}   \tab mean age of the sample \cr
\bold{quality}   \tab \code{numeric}   \tab methodological quality
}
}
\details{
   Conscientiousness, one of the big-5 personality traits, can be defined as \dQuote{socially prescribed impulse control that facilitates task- and goal-directed behaviour, such as thinking before acting, delaying gratification, following norms and rules and planning, organising and prioritising tasks} (John & Srivastava, 1999). Conscientiousness has been shown to be related to a number of health-related behaviors (e.g., tobacco/alcohol/drug use, diet and activity patterns, risky behaviors). A recent meta-analysis by Molloy et al. (2014) examined to what extent conscientiousness is related to medication adherence, that is, the extent to which (typically chronically ill) patients follow a prescribed medication regimen (e.g., taking a daily dose of a cholesterol lowering drug in patients with high LDL serum cholesterol levels). The results from the 16 studies included in this meta-analysis are provided in this dataset.

   Variable \code{a_measure} indicates whether adherence was measured based on self-reports or a more \sQuote{objective} measure (e.g., electronic monitoring of pill bottle openings, pill counts). Variable \code{c_measure} indicates whether conscientiousness was measured with some version of the NEO personality inventory or some other scale. Methodological quality was scored by the authors on a 1 to 4 scale with higher scores indicating higher quality (see article for details on how this score was derived).
}
\source{
   Molloy, G. J., O'Carroll, R. E., & Ferguson, E. (2014). Conscientiousness and medication adherence: A meta-analysis. \emph{Annals of Behavioral Medicine, 47,} 92--101.
}
\references{
   John, O. P., & Srivastava, S. (1999). The Big Five Trait taxonomy: History, measurement, and theoretical perspectives. In L. A. Pervin & O. P. John (Eds.), \emph{Handbook of personality: Theory and research} (2nd ed., pp. 102-138). New York: Guilford Press.
}
\examples{
### load data
dat <- get(data(dat.molloy2014))

### calculate r-to-z transformed correlations and corresponding sampling variances
dat <- escalc(measure="ZCOR", ri=ri, ni=ni, data=dat, slab=paste(authors, year, sep=", "))
dat

### meta-analysis of the transformed correlations using a random-effects model
res <- rma(yi, vi, data=dat)
res

### average correlation with 95\% CI
predict(res, digits=3, transf=transf.ztor)

### forest plot
forest(res, addcred=TRUE, xlim=c(-1.6,1.6), atransf=transf.ztor,
       at=transf.rtoz(c(-.4,-.2,0,.2,.4,.6)), digits=c(2,1), cex=.8)
text(-1.6, 18, "Author(s), Year", pos=4, cex=.8)
text( 1.6, 18, "Correlation [95\% CI]", pos=2, cex=.8)

### funnel plot
funnel(res)
}
\keyword{datasets}
