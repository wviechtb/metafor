\name{print.anova.rma}
\alias{print.anova.rma}
\alias{print.list.anova.rma}
\title{Print Methods for 'anova.rma' and 'list.anova.rma' Objects}
\description{
   Functions to print objects of class \code{"anova.rma"} and \code{"list.anova.rma"}. \loadmathjax
}
\usage{
\method{print}{anova.rma}(x, digits=x$digits, \dots)
\method{print}{list.anova.rma}(x, digits=x[[1]]$digits, \dots)
}
\arguments{
   \item{x}{an object of class \code{"anova.rma"} or \code{"list.anova.rma"} obtained with \code{\link[=anova.rma]{anova}}.}
   \item{digits}{integer to specify the number of decimal places to which the printed results should be rounded (the default is to take the value from the object).}
   \item{\dots}{other arguments.}
}
\details{
   For a Wald-type test of one or multiple model coefficients, the output includes the test statistic (either a chi-square or F-value) and the corresponding p-value.

   When testing one or multiple contrasts, the output includes the estimated value of the contrast, its standard error, test statistic (either a z- or a t-value), and the corresponding p-value.

   When comparing two model objects, the output includes:
   \itemize{
   \item the number of parameters in the full and the reduced model.
   \item the AIC, BIC, AICc, and log-likelihood of the full and the reduced model.
   \item the value of the likelihood ratio test statistic.
   \item the corresponding p-value.
   \item the test statistic of the test for (residual) heterogeneity for the full and the reduced model.
   \item the estimate of \mjseqn{\tau^2} from the full and the reduced model. Suppressed for equal-effects models.
   \item amount (in percent) of heterogeneity in the reduced model that is accounted for in the full model (\code{NA} for \code{"rma.mv"} objects). This can be regarded as a pseudo \mjseqn{R^2} statistic (Raudenbush, 2009). Note that the value may not be very accurate unless \mjseqn{k} is large (Lopez-Lopez et al., 2014).
   }
   The last two items are not provided when comparing \code{"rma.mv"} models.
}
\value{
   The function does not return an object.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   \enc{López-López}{Lopez-Lopez}, J. A., \enc{Marín-Martínez}{Marin-Martinez}, F., \enc{Sánchez-Meca}{Sanchez-Meca}, J., Van den Noortgate, W., & Viechtbauer, W. (2014). Estimation of the predictive power of the model in mixed-effects meta-regression: A simulation study. \emph{British Journal of Mathematical and Statistical Psychology}, \bold{67}(1), 30--48. \verb{https://doi.org/10.1111/bmsp.12002}

   Raudenbush, S. W. (2009). Analyzing effect sizes: Random effects models. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta-analysis} (2nd ed., pp. 295--315). New York: Russell Sage Foundation.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link[=anova.rma]{anova}} for the function to create \code{anova.rma} objects.
}
\keyword{print}
