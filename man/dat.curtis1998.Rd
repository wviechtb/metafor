\name{dat.curtis1998}
\docType{data}
\alias{dat.curtis1998}
\title{Studies on the Effects of Elevated CO2 Levels on Woody Plant Mass}
\description{Results from studies examining the effects of elevated CO2 levels on woody plant mass.}
\usage{dat.curtis1998}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{id}       \tab \code{numeric}   \tab observation number \cr
\bold{paper}    \tab \code{numeric}   \tab paper number \cr
\bold{genus}    \tab \code{character} \tab genus name \cr
\bold{species}  \tab \code{character} \tab species name \cr
\bold{fungrp}   \tab \code{character} \tab plant functional group \cr
\bold{co2.ambi} \tab \code{numeric}   \tab ambient CO2 level (control group) \cr
\bold{co2.elev} \tab \code{numeric}   \tab elevated CO2 level (treatment group) \cr
\bold{units}    \tab \code{character} \tab units for CO2 exposure levels \cr
\bold{time}     \tab \code{numeric}   \tab maximum length of time (days) of CO2 exposure \cr
\bold{pot}      \tab \code{character} \tab growing method (see below) \cr
\bold{method}   \tab \code{character} \tab CO2 exposure facility (see below) \cr
\bold{stock}    \tab \code{character} \tab planting stock code \cr
\bold{xtrt}     \tab \code{character} \tab interacting treatment code (see below) \cr
\bold{level}    \tab \code{character} \tab interacting treatment level codes (see below) \cr
\bold{m1i}      \tab \code{numeric}   \tab mean plant mass under elevated CO2 level (treatment group) \cr
\bold{sd1i}     \tab \code{numeric}   \tab standard deviation of plant mass underelevated CO2 level (treatment group) \cr
\bold{n1i}      \tab \code{numeric}   \tab number of observations under elevated CO2 level (treatment group) \cr
\bold{m2i}      \tab \code{numeric}   \tab mean plant mass under ambient CO2 level (control group) \cr
\bold{sd2i}     \tab \code{numeric}   \tab standard deviation of plant mass under ambient CO2 level (control group) \cr
\bold{n2i}      \tab \code{numeric}   \tab number of observations under ambient CO2 level (control group)
}
}
\details{
   The studies included in this dataset compared the total above- plus below-ground biomass (in grams) for plants that were either exposed to ambient (around 35 Pa) and elevated CO2 levels (around twice the ambient level). The \code{co2.ambi} and \code{co2.elev} variables indicate the CO2 levels in the control and treatment groups, respectively (with the \code{units} variable specifying the units for the CO2 exposure levels). Many of the studies also varied one or more additional environmental variables (defined by the \code{xtrt} and \code{level} variables):
   \itemize{
   \item NONE = no additional treatment factor
   \item FERT = soil fertility (either a \code{CONTROL}, \code{HIGH}, or \code{LOW} level)
   \item LIGHT = light treatment (always a \code{LOW} light level)
   \item FERT+L = soil fertility and light (a \code{LOW} light and soil fertility level)
   \item H2O = well watered vs drought (either a \code{WW} or \code{DRT} level)
   \item TEMP = temperature treatment (either a \code{HIGH} or \code{LOW} level)
   \item OZONE = ozone exposure (either a \code{HIGH} or \code{LOW} level)
   \item UVB = ultraviolet-B radiation exposure (either a \code{HIGH} or \code{LOW} level)
   } In addition, the studies differed with respect to various design variables, including CO2 exposure duration (\code{time}), growing method (\code{pot}: number = pot size in liters; \code{GRND} = plants rooted in ground; \code{HYDRO} = solution or aeroponic culture), CO2 exposure facility (\code{method}: \code{GC} = growth chamber; \code{GH} = greenhouse; \code{OTC} = field-based open-top chamber), and planting stock (\code{stock}: \code{SEED} = plants started from seeds; \code{SAP} = plants started from cuttings). The goal of the meta-analysis was to examine the effects of elevated CO2 levels on plant physiology and growth and the interacting effects of the environmental (and design) variables.
}
\source{
   Hedges, L. V., Gurevitch, J., & Curtis, P. S. (1999). The meta-analysis of response ratios in experimental ecology. \emph{Ecology}, \bold{80}, 1150--1156. (data obtained from \emph{Ecological Archives}, E080-008-S1, at: \url{http://www.esapubs.org/archive/ecol/E080/008/})
}
\references{
   Curtis, P. S., & Wang, X. (1998). A meta-analysis of elevated CO2 effects on woody plant mass, form, and physiology. \emph{Oecologia}, \bold{113}, 299--313.
}
\examples{
### copy data into 'dat'
dat <- dat.curtis1998

### calculate log ratio of means and corresponding sampling variances
dat <- escalc(measure="ROM", m1i=m1i, sd1i=sd1i, n1i=n1i,
                             m2i=m2i, sd2i=sd2i, n2i=n2i, data=dat)
dat

### meta-analysis of log ratio of means using a random-effects model
res <- rma(yi, vi, method="DL", data=dat)
res

### average ratio of means with 95\% CI
predict(res, transf=exp, digits=2)

### meta-analysis for plants grown under nutrient stress
res <- rma(yi, vi, method="DL", data=dat, subset=(xtrt=="FERT" & level=="LOW"))
predict(res, transf=exp, digits=2)

### meta-analysis for plants grown under low light conditions
res <- rma(yi, vi, method="DL", data=dat, subset=(xtrt=="LIGHT" & level=="LOW"))
predict(res, transf=exp, digits=2)
}
\keyword{datasets}
