\name{dat.li2007}
\docType{data}
\alias{dat.li2007}
\title{Studies on the Effectiveness of Intravenous Magnesium in Acute Myocardial Infarction.}
\description{Results from 22 trials examining the effectiveness of intravenous magnesium in the prevention of death following acute myocardial infarction.}
\usage{
dat.li2007
}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{id}    \tab \code{numeric}   \tab trial id number \cr
\bold{study} \tab \code{character} \tab first author or trial name \cr
\bold{year}  \tab \code{numeric}   \tab publication year \cr
\bold{ai}    \tab \code{numeric}   \tab number of deaths in the magnesium group \cr
\bold{n1i}   \tab \code{numeric}   \tab number of patients in the magnesium group \cr
\bold{ci}    \tab \code{numeric}   \tab number of deaths in the control group \cr
\bold{n2i}   \tab \code{numeric}   \tab number of patients in the control group \cr
}
}
\details{
   The dataset includes the results from 22 randomized clinical trials that examined the effectiveness of intravenous magnesium in the prevention of death following acute myocardial infarction. It is similar to the dataset \code{\link{dat.egger2001}}, with some slight differences in the included trials and data used.
}
\source{
   Li, J., Zhang, Q., Zhang, M., & Egger, M. (2007). Intravenous magnesium for acute myocardial infarction. \emph{Cochrane Database of Systematic Reviews}, \bold{2}, CD002755.
}
\seealso{
   \code{\link{dat.egger2001}}
}
\examples{
### load data
dat <- get(data(dat.li2007))

### meta-analysis of all trials except ISIS-4
res <- rma(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, method="FE", subset=-14)
print(res, digits=2)
predict(res, transf=exp, digits=2)

### meta-analysis of all trials including ISIS-4
res <- rma(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, method="FE")
print(res, digits=2)
predict(res, transf=exp, digits=2)

### contour-enhanced funnel plot centered at 0
funnel(res, refline=0, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"))
}
\keyword{datasets}
