\name{print.anova.rma}
\alias{print.anova.rma}
\title{Print Method for 'anova.rma' Objects}
\description{
   Print method for objects of class \code{"anova.rma"}.
}
\usage{
\method{print}{anova.rma}(x, digits, \dots)
}
\arguments{
   \item{x}{an object of class \code{"anova.rma"}.}
   \item{digits}{integer specifying the number of decimal places to which the printed results should be rounded (if unspecified, the default is to take the value from the object).}
   \item{\dots}{other arguments.}
}
\details{
   The output includes:
   \itemize{
   \item the number of parameters in the full and the reduced model.
   \item the AIC, BIC, AICc, and log-likelihood of the full and the reduced model.
   \item the value of the likelihood ratio test statistic.
   \item the p-value for the likelihood ratio test.
   \item the test statistic for the test of (residual) heterogeneity for the full and the reduced model.
   \item the estimate of \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}} from the full and the reduced model. Suppressed for fixed-effects models.
   \item{R2}{amount of (residual) heterogeneity in the reduced model that is accounted for in the full model (in percent). \code{NA} for fixed-effects models, if the amount of heterogeneity in the reduced model is equal to zero, or for \code{"rma.mv"} objects. This can be regarded as a pseudo \if{latex}{\eqn{R^2}}\if{html}{\eqn{R²}}\if{text}{\eqn{R^2}} statistic (Raudenbush, 2009). Note that the value may not be very accurate unless \eqn{k} is large (Lopez-Lopez et al., 2014).}
   }
   The last two items are not provided when comparing \code{"rma.mv"} models.
}
\value{
   The function does not return an object.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{http://www.metafor-project.org/}
}
\references{
   \enc{López-López}{Lopez-Lopez}, J. A., \enc{Marín-Martínez}{Marin-Martinez}, F., \enc{Sánchez-Meca}{Sanchez-Meca}, J., Van den Noortgate, W., & Viechtbauer, W. (2014). Estimation of the predictive power of the model in mixed-effects meta-regression: A simulation study. \emph{British Journal of Mathematical and Statistical Psychology}, \bold{67}, 30--48.

   Raudenbush, S. W. (2009). Analyzing effect sizes: Random effects models. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta-analysis} (2nd ed., pp. 295--315). New York: Russell Sage Foundation.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
   \code{\link{anova.rma}}
}
\keyword{print}
