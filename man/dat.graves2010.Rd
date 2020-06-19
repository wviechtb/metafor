\name{dat.graves2010}
\docType{data}
\alias{dat.graves2010}
\title{Studies on the Effectiveness of Injected Cholera Vaccines}
\description{Results from 17 studies on the effectiveness of injected vaccines against cholera.}
\usage{
dat.graves2010
}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{study}   \tab \code{character} \tab author/study name and publication year \cr
\bold{ai}      \tab \code{numeric}   \tab number of cholera cases in the vaccinated group \cr
\bold{n1i}     \tab \code{numeric}   \tab number of individuals in the vaccinated group \cr
\bold{ci}      \tab \code{numeric}   \tab number of cholera cases in the placebo group \cr
\bold{n2i}     \tab \code{numeric}   \tab number of individuals in the placebo group
}
}
\details{
   Cholera is an infection caused by certain strains of the bacterium \emph{Vibrio cholerae}. When untreated, mortality rates can be as high as 50-60\%. Proper sanitation practices are usually effective in preventing outbreaks, but a number of oral and injectable vaccines have also been developed. The Cochrane review by Graves et al. (2010) examined the effectiveness of injectable vaccines for preventing cholera cases and death. The present dataset includes results from 17 studies that reported the number of cholera cases in vaccinated and placebo/comparison groups up to 7 months after the treatment.
}
\source{
   Graves, P. M., Deeks, J. J., Demicheli, V., & Jefferson, T. (2010). Vaccines for preventing cholera: Killed whole cell or other subunit vaccines (injected). \emph{Cochrane Database of Systematic Reviews}, \bold{8}, CD000974. \verb{https://doi.org/10.1002/14651858.CD000974.pub2}
}
\examples{
### copy data into 'dat' and examine data
dat <- dat.graves2010
dat

### analysis using the Mantel-Haenszel method
rma.mh(measure="RR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, digits=2)
}
\keyword{datasets}
