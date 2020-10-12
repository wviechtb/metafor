\name{dat.hahn2001}
\docType{data}
\alias{dat.hahn2001}
\title{Studies on the Effectiveness of Different Rehydration Solutions for the Prevention of Unscheduled Intravenous Infusion in Children with Diarrhoea}
\description{Results from 12 trials examining the effectiveness of a reduced versus standard rehydration solution for the prevention of unscheduled intravenous infusion in children with diarrhoea.}
\usage{
dat.hahn2001
}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{study} \tab \code{character} \tab trial name and year \cr
\bold{ai}    \tab \code{numeric}   \tab number of children requiring unscheduled intravenous infusion in the reduced rehydration solution group \cr
\bold{n1i}   \tab \code{numeric}   \tab number of children in the reduced rehydration solution group \cr
\bold{ci}    \tab \code{numeric}   \tab number of children requiring unscheduled intravenous infusion in the standard rehydration solution group \cr
\bold{n2i}   \tab \code{numeric}   \tab number of children in the standard rehydration solution group
}
}
\details{
   The dataset includes the results from 12 randomized clinical trials that examined the effectiveness of a reduced osmolarity oral rehydration solution (total osmolarity <250 mmol/l with reduced sodium) with a standard WHO oral rehydration solution (sodium 90 mmol/l, glucose 111mmol/l, total osmolarity 311 mmol/l) for the prevention of unscheduled intravenous infusion in children with diarrhoea.
}
\source{
   Hahn, S., Kim, Y., & Garner, P. (2001). Reduced osmolarity oral rehydration solution for treating dehydration due to diarrhoea in children: Systematic review. \emph{British Medical Journal}, \bold{323}(7304), 81--85. \verb{https://doi.org/10.1136/bmj.323.7304.81}
}
\examples{
### copy data into 'dat' and examine data
dat <- dat.hahn2001
dat

### meta-analysis of (log) odds rations using the Mantel-Haenszel method
res <- rma.mh(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, digits=2, slab=study)
res

### forest plot (also show studies that were excluded from the analysis)
options(na.action="na.pass")
forest(res, atransf=exp, at=log(c(.01, .1, 1, 10, 100)), header=TRUE)
options(na.action="na.omit")
}
\keyword{datasets}
