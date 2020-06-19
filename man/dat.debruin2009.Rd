\name{dat.debruin2009}
\docType{data}
\alias{dat.debruin2009}
\title{Studies on Standard Care Quality and HAART-Adherence}
\description{Results from 13 trials providing information about standard care quality and HAART-adherence in control groups.}
\usage{dat.debruin2009}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{author}      \tab \code{character} \tab (first) author of study \cr
\bold{year}        \tab \code{numeric}   \tab publication year \cr
\bold{scq}         \tab \code{numeric}   \tab standard care quality \cr
\bold{ni}          \tab \code{numeric}   \tab number of patients in the standard care group \cr
\bold{xi}          \tab \code{numeric}   \tab number of patients with an undetectable viral load in standard care group \cr
\bold{mi}          \tab \code{numeric}   \tab number of patients with a detectable viral load in standard care group \cr
\bold{ethnicity}   \tab \code{character} \tab dominant ethnicity of the patients in the standard care group \cr
\bold{patients}    \tab \code{character} \tab inclusion of patients continuing or starting (a new) treatment \cr
\bold{select}      \tab \code{character} \tab baseline selection of patients with adherence problems or no selection \cr
\bold{sens}        \tab \code{character} \tab sensitivity of viral load assessments (<400 vs. >=400 copies/ml)
}
}
\details{
   Highly active antiretroviral therapy (HAART) refers to a combination of multiple antiretroviral drugs that can effectively suppress the HIV virus. However, achieving viral suppression (to the point that the virus becomes essentially undetectable in a blood sample) requires high levels of adherence to an often complicated medication regimen. A number of trials have examined various interventions that aim to increase adherence levels. In each trial, patients receiving the intervention are compared to patients in a control group receiving standard care (often referred to as \sQuote{care as usual}). However, the quality of standard care can vary substantially between these studies. de Bruin et al. (2009) assessed the quality of standard care provided (based on a quantification of the number of behavior change techniques applied) and examined to what extent the quality of standard care was related to the proportion of patients achieving effective viral suppression in the control groups.
}
\source{
   de Bruin, M., Viechtbauer, W., Hospers, H. J., Schaalma, H. P., & Kok, G. (2009). Standard care quality determines treatment outcomes in control groups of HAART-adherence intervention studies: Implications for the interpretation and comparison of intervention effects. \emph{Health Psychology}, \bold{28}(6), 668--674. \verb{https://doi.org/10.1037/a0015989}
}
\examples{
### copy data into 'dat'
dat <- dat.debruin2009

### calculate proportions and corresponding sampling variances
dat <- escalc(measure="PR", xi=xi, ni=ni, data=dat)
dat

### random-effects model
res <- rma(yi, vi, data=dat)
print(res, digits=2)

### mixed-effects meta-regression model with all predictors/covariates
res <- rma(yi, vi, mods = ~ scq + ethnicity + patients + select + sens, data=dat)
print(res, digits=3)

### mixed-effects meta-regression model with scq and ethnicity as predictors/covariates
res <- rma(yi, vi, mods = ~ scq + ethnicity, data=dat)
print(res, digits=3)
}
\keyword{datasets}
