\name{cmicalc}
\alias{cmicalc}
\title{Calculate the Bias Correction Factor for SMDs and Related Measures}
\description{
   Function to calculate the bias correction factor for standardized mean differences and related measures. \loadmathjax
}
\usage{
cmicalc(mi, method="exact")
}
\arguments{
   \item{mi}{numeric vector.}
   \item{method}{character string to specify if the exact (\code{"exact"}) or approximate (\code{"approx"}) formula should be used.}
}
\details{
   Standardized mean differences (Cohen's d values) have a slight positive bias (i.e., they tend to be too large on average). Hedges (1981) derived the bias correction factor to make standardized mean differences unbiased (yielding what is commonly referred to as Hedges' g values).

   The function computes the bias correction factor, using \mjseqn{m_i} as input. For standardized mean differences, this is equal to \mjseqn{m_i = n_{1i} + n_{2i} - 2}, where \mjseqn{n_{1i}} and \mjseqn{n_{2i}} are the number of subjects in the two groups. The same bias correction factor can be used for various other related measures, such as standardized mean changes, except that \mjseqn{m_i} may need to be calculated differently.

   By default (i.e., when \code{method="exact"}), the method uses the exact formula to calculate the bias correction factor (equation 6e in Hedges, 1981), that is, \mjsdeqn{c(m_i) = \frac{\Gamma(\frac{m_i}{2})}{\sqrt{\frac{m_i}{2}}\Gamma(\frac{m_i-1}{2})}.} When \code{method="approx"}, the (very accurate) approximation \mjsdeqn{c(m_i) = 1 - \frac{3}{4m_i -1}} is used instead (first equation on page 114 in Hedges, 1981).
}
\value{
   A numeric vector of length equal to the length of \code{mi}.
}
\note{
   The notation \mjseqn{c(m_i)} is based on that used in Hedges (1981). Later publications and references often denote the bias correction factor \mjseqn{J} or \mjseqn{J(m_i)}.
}
\author{
   Wolfgang Viechtbauer (\email{wvb@metafor-project.org}, \url{https://www.metafor-project.org}).
}
\references{
   Hedges, L. V. (1981). Distribution theory for Glass's estimator of effect size and related estimators. \emph{Journal of Educational Statistics}, \bold{6}(2), 107--128. \verb{https://doi.org/10.3102/10769986006002107}
}
\seealso{
   \code{\link{escalc}} for a function that makes use of this bias correction.
}
\examples{
### copy dataset to dat
dat <- dat.furukawa2003

### calculate the bias correction factor for each study
dat$mi <- dat$Ne + dat$Nc - 2
dat$cmi <- cmicalc(dat$mi)
dat
}
\keyword{manip}
