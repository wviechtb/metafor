\name{vif}
\alias{vif}
\alias{vif.rma}
\alias{print.vif.rma}
\title{Variance Inflation Factors for 'rma' Objects}
\description{
   Function to compute (generalized) variance inflation factors (VIFs) for objects of class \code{"rma"}. \loadmathjax
}
\usage{
vif(x, \dots)

\method{vif}{rma}(x, btt, att, table=FALSE, reestimate=FALSE, sim=FALSE, progbar=TRUE,
    seed=NULL, parallel="no", ncpus=1, cl, digits, \dots)

\method{print}{vif.rma}(x, digits=x$digits, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma"} (for \code{vif}) or \code{"vif.rma"} (for \code{print}).}
   \item{btt}{optional vector of indices (or list thereof) to specify a set of coefficients for which a generalized variance inflation factor (GVIF) should be computed. Can also be a string to \code{\link{grep}} for.}
   \item{att}{optional vector of indices (or list thereof) to specify a set of scale coefficients for which a generalized variance inflation factor (GVIF) should be computed. Can also be a string to \code{\link{grep}} for. Only relevant for location-scale models (see \code{\link{rma.uni}}).}
   \item{table}{logical to specify whether the VIFs should be added to the model coefficient table (the default is \code{FALSE}). Only relevant when \code{btt} (or \code{att}) is not specified.}
   \item{reestimate}{logical to specify whether the model should be reestimated when removing moderator variables from the model for computing a (G)VIF (the default is \code{FALSE}).}
   \item{sim}{logical to specify whether the distribution of each (G)VIF under independence should be simulated (the default is \code{FALSE}). Can also be an integer to specify how many values to simulate (when \code{sim=TRUE}, the default is \code{1000}).}
   \item{progbar}{logical to specify whether a progress bar should be shown when \code{sim=TRUE} (the default is \code{TRUE}).}
   \item{seed}{optional value to specify the seed of the random number generator when \code{sim=TRUE} (for reproducibility).}
   \item{parallel}{character string to specify whether parallel processing should be used (the default is \code{"no"}). For parallel processing, set to either \code{"snow"} or \code{"multicore"}. See \sQuote{Note}.}
   \item{ncpus}{integer to specify the number of processes to use in the parallel processing.}
   \item{cl}{optional cluster to use if \code{parallel="snow"}. If unspecified, a cluster on the local machine is created for the duration of the call.}
   \item{digits}{optional integer to specify the number of decimal places to which the printed results should be rounded. If unspecified, the default is to take the value from the object.}
   \item{\dots}{other arguments.}
}
\details{
   The function computes (generalized) variance inflation factors (VIFs) for meta-regression models. Hence, the model specified via argument \code{x} must include moderator variables (and more than one for this to be useful, as the VIF for a model with a single moderator variable will always be equal to 1).

   \subsection{VIFs for Individual Coefficients}{

      By default (i.e., if \code{btt} is not specified), VIFs are computed for the individual model coefficients.

      Let \mjseqn{b_j} denote the estimate of the \mjeqn{j\textrm{th}}{jth} model coefficient of a particular meta-regression model and \mjeqn{\text{Var}[b_j]}{Var[b_j]} its variance (i.e., the corresponding diagonal element from the matrix obtained with the \code{\link[=vcov.rma]{vcov}} function). Moreover, let \mjseqn{b'_j} denote the estimate of the same model coefficient if the other moderator variables in the model had \emph{not} been included in the model and \mjeqn{\text{Var}[b'_j]}{Var[b'_j]} the corresponding variance. Then the VIF for the model coefficient is given by \mjdeqn{\text{VIF}[b_j] = \frac{\text{Var}[b_j]}{\text{Var}[b'_j]},}{VIF[b_j] = Var[b_j] / Var[b'_j],} which indicates the inflation in the variance of the estimated model coefficient due to potential collinearity of the \mjeqn{j\textrm{th}}{jth} moderator variable with the other moderator variables in the model. Taking the square root of a VIF gives the corresponding standard error inflation factor (SIF).

   }

   \subsection{GVIFs for Sets of Coefficients}{

      If the model includes factors (coded in terms of multiple dummy variables) or other sets of moderator variables that belong together (e.g., for polynomials or cubic splines), one may want to examine how much the variance in all of the coefficients in the set is jointly impacted by collinearity with the other moderator variables in the model. For this, we can compute a generalized variance inflation factor (GVIF) (Fox & Monette, 1992) by setting the \code{btt} argument equal to the indices of those coefficients for which the GVIF should be computed. The square root of a GVIF indicates the inflation in the confidence ellipse/(hyper)ellipsoid for the set of coefficients corresponding to the set due to collinearity. However, to make this value more directly comparable to SIFs (based on single coefficients), the function computes the generalized standard error inflation factor (GSIF) by raising the GVIF to the power of \mjseqn{1/(2m)} (where \mjseqn{m} denotes the number of coefficients in the set). One can also specify a list of indices/strings, in which case GVIFs/GSIFs of all list elements will be computed. See \sQuote{Examples}.

      For location-scale models fitted with the \code{\link{rma.uni}} function, one can use the \code{att} argument in an analogous manner to specify the indices of the scale coefficients for which a GVIF should be computed.

   }

   \subsection{Re-Estimating the Model}{

      The way the VIF is typically computed for a particular model coefficient (or a set thereof for a GVIF) makes use of some clever linear algebra to avoid having to re-estimate the model when removing the other moderator variables from the model. This speeds up the computations considerably. However, this assumes that the other moderator variables do not impact other aspects of the model, in particular the amount of residual heterogeneity (or, more generally, any variance/correlation components in a more complex model, such as those that can be fitted with the \code{\link{rma.mv}} function).

      For a more accurate (but slower) computation of each (G)VIF, one can set \code{reestimate=TRUE}, in which case the model is refitted to account for the impact that removal of the other moderator variables has on all aspects of the model. Note that refitting may fail, in which case the corresponding (G)VIF will be missing.

   }

   \subsection{Interpreting the Size of (G)VIFs}{

      A large VIF value suggests that the precision with which we can estimate a particular model coefficient (or a set thereof for a GVIF) is negatively impacted by multicollinearity among the moderator variables. However, there is no specific cutoff for determining whether a particular (G)VIF is \sQuote{large}. Sometimes, values such as 5 or 10 have been suggested as rules of thumb, but such cutoffs are essentially arbitrary.

   }

   \subsection{Simulating the Distribution of (G)VIFs Under Independence}{

      As a more principled approach, we can simulate the distribution of a particular (G)VIF under independence and then examine how extreme the actually observed (G)VIF value is under this distribution. The distribution is simulated by randomly reshuffling the columns of the model matrix (to break any dependence between the moderators) and recomputing the (G)VIF. When setting \code{sim=TRUE}, this is done 1000 times (but one can also set \code{sim} to an integer to specify how many (G)VIF values should be simulated).

      The way the model matrix is reshuffled depends on how the model was fitted. When the model was specified as a formula via the \code{mods} argument and the data was supplied via the \code{data} argument, then each column of the data frame specified via \code{data} is reshuffled and the formula is evaluated within the reshuffled data (creating the corresponding reshuffled model matrix). This way, factor/character variables are properly reshuffled and derived terms (e.g., interactions, polynomials, splines) are correct constructed. This is the recommended approach.

      On the other hand, if the model matrix was directly supplied to the \code{mods} argument, then each column of the matrix is directly reshuffled. This is not recommended, since this approach cannot account for any inherent relationships between variables in the model matrix (e.g., an interaction term is the product of two variables and should not be reshuffled by itself).

      Once the distribution of a (G)VIF under independence has been simulated, the proportion of simulated values that are smaller than the actually observed (G)VIF value is computed. If the proportion is close to 1, then this indicates that the actually observed (G)VIF value is extreme.

      The general principle underlying the simulation approach is the same as that underlying Horn's parallel analysis (1965) for determining the number of components / factors to keep in a principal component / factor analysis.

   }

}
\value{
   An object of class \code{"vif.rma"}. The object is a list containing the following components:
   \item{vif}{a list of data frames with the (G)VIFs and (G)SIFs and some additional information.}
   \item{vifs}{a vector with just the (G)VIFs.}
   \item{table}{the model coefficient table (only when \code{table=TRUE}).}
   \item{sim}{a matrix with the simulated (G)VIF values (only when \code{sim=TRUE}).}
   \item{prop}{a vector with the proportions of simulated values that are smaller than the actually observed (G)VIF values (only when \code{sim=TRUE}).}
   \item{\dots}{some additional elements/values.}

   When \code{x} was a location-scale model object and (G)VIFs can be computed for both the location and the scale coefficients, then the object is a list with elements \code{beta} and \code{alpha}, where each element is an \code{"vif.rma"} object as described above.

   The results are formatted and printed with the \code{print} function. To format the results as a data frame, one can use the \code{\link[=as.data.frame.vif.rma]{as.data.frame}} function. When \code{sim=TRUE}, the distribution of each (G)VIF can be plotted with the \code{\link[=plot.vif.rma]{plot}} function.
}
\note{
   If the original model fitted involved redundant predictors that were dropped from the model, then \code{sim=TRUE} cannot be used. In this case, one should remove any redundancies in the original model fitted before using this method.

   When using \code{sim=TRUE}, the model needs to be refitted (by default) 1000 times. When \code{sim=TRUE} is combined with \code{reestimate=TRUE}, then this value needs to be multiplied by the total number of (G)VIF values that are computed by the function. Hence, the combination of \code{sim=TRUE} with \code{reestimate=TRUE} is computationally expensive, especially for more complex models where model fitting can be slow.

   When refitting the model fails, the simulated (G)VIF value(s) will be missing. It can also happen that one or multiple model coefficients become inestimable due to redundancies in the model matrix after the reshuffling. In this case, the corresponding simulated (G)VIF value(s) will be set to \code{Inf} (as that is the value of (G)VIFs in the limit as we approach perfect multicollinearity).

   On machines with multiple cores, one can try to speed things up by delegating the model fitting to separate worker processes, that is, by setting \code{parallel="snow"} or \code{parallel="multicore"} and \code{ncpus} to some value larger than 1. Parallel processing makes use of the \code{\link[parallel]{parallel}} package, using the \code{\link[parallel]{makePSOCKcluster}} and \code{\link[parallel]{parLapply}} functions when \code{parallel="snow"} or using \code{\link[parallel]{mclapply}} when \code{parallel="multicore"} (the latter only works on Unix/Linux-alikes). With \code{parallel::detectCores()}, one can check on the number of available cores on the local machine.
}
\author{
   Wolfgang Viechtbauer (\email{wvb@metafor-project.org}, \url{https://www.metafor-project.org}).
}
\references{
   Belsley, D. A., Kuh, E., & Welsch, R. E. (1980). \emph{Regression diagnostics}. New York: Wiley.

   Fox, J., & Monette, G. (1992). Generalized collinearity diagnostics. \emph{Journal of the American Statistical Association}, \bold{87}(417), 178--183. \verb{https://doi.org/10.2307/2290467}

   Horn, J. L. (1965). A rationale and test for the number of factors in factor analysis. \emph{Psychometrika}, \bold{30}(2), 179--185. \verb{https://doi.org/10.1007/BF02289447}

   Stewart, G. W. (1987). Collinearity and least squares regression. \emph{Statistical Science}, \bold{2}(1), 68--84. \verb{https://doi.org/10.1214/ss/1177013439}

   Wax, Y. (1992). Collinearity diagnosis for a relative risk regression-analysis: An application to assessment of diet cancer relationship in epidemiologic studies. \emph{Statistics in Medicine}, \bold{11}(10), 1273--1287. \verb{https://doi.org/10.1002/sim.4780111003}

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}

   Viechtbauer, W., & \enc{López-López}{Lopez-Lopez}, J. A. (2022). Location-scale models for meta-analysis. \emph{Research Synthesis Methods}. \bold{13}(6), 697--715. \verb{https://doi.org/10.1002/jrsm.1562}
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.glmm}}, and \code{\link{rma.mv}} for functions to fit models for which variance inflation factors can be computed.

   \code{\link[=plot.vif.rma]{plot}} for the plot method and \code{\link[=as.data.frame.vif.rma]{as.data.frame}} for the method to format the results as a data frame.
}
\examples{
### copy data from Bangert-Drowns et al. (2004) into 'dat'
dat <- dat.bangertdrowns2004

### fit mixed-effects meta-regression model
res <- rma(yi, vi, mods = ~ length + wic + feedback + info + pers + imag + meta, data=dat)

### get variance inflation factors
vif(res)

### use the simulation approach to analyze the size of the VIFs
\dontrun{
vif(res, sim=TRUE, seed=1234)
}

### get variance inflation factors using the re-estimation approach
vif(res, reestimate=TRUE)

### show that VIFs are not influenced by scaling of the predictors
u <- scale # to standardize the predictors
res <- rma(yi, vi, mods = ~ u(length) + u(wic) + u(feedback) + u(info) +
                            u(pers) + u(imag) + u(meta), data=dat)
vif(res, reestimate=TRUE)

### get full table
vif(res, reestimate=TRUE, table=TRUE)

### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### fit meta-regression model where one predictor (alloc) is a three-level factor
res <- rma(yi, vi, mods = ~ ablat + alloc + year, data=dat)

### get variance inflation factors for all individual coefficients
vif(res, table=TRUE)

### generalized variance inflation factor for the 'alloc' factor
vif(res, btt=3:4)

### can also specify a string to grep for
vif(res, btt="alloc")

### can also specify a list for the 'btt' argument (and use the simulation approach)
\dontrun{
vif(res, btt=list(2,3:4,5), sim=TRUE, seed=1234)
}
}
\keyword{models}
