\name{predict.matreg}
\alias{predict.matreg}
\title{Predicted Values for 'matreg' Objects}
\description{
   The function computes predicted values, corresponding standard errors, and confidence intervals for objects of class \code{"matreg"}. \loadmathjax
}
\usage{
\method{predict}{matreg}(object, newmods, intercept, addx=FALSE,
        level, adjust=FALSE, digits, transf, targs, vcov=FALSE, \dots)
}
\arguments{
   \item{object}{an object of class \code{"matreg"}.}
   \item{newmods}{vector or matrix to specify the values of the moderator/predictor values for which the predicted values should be calculated. See \sQuote{Details}.}
   \item{intercept}{logical to specify whether the intercept should be included when calculating the predicted values for \code{newmods}. If unspecified, the intercept is automatically added when the model also included an intercept.}
   \item{addx}{logical to specify whether the values of the moderator/predictor variables should be added to the returned object. See \sQuote{Examples}.}
   \item{level}{numeric value between 0 and 100 to specify the confidence interval level (see \link[=misc-options]{here} for details). If unspecified, the default is to take the value from the object.}
   \item{adjust}{logical to specify whether the width of confidence intervals should be adjusted using a Bonferroni correction (the default is \code{FALSE}).}
   \item{digits}{optional integer to specify the number of decimal places to which the printed results should be rounded.}
   \item{transf}{optional argument to specify a function to transform the predicted values and interval bounds (e.g., \code{transf=exp}; see also \link{transf}). If unspecified, no transformation is used.}
   \item{targs}{optional arguments needed by the function specified under \code{transf}.}
   \item{vcov}{logical to specify whether the variance-covariance matrix of the predicted values should also be returned (the default is \code{FALSE}).}
   \item{\dots}{other arguments.}
}
\details{
   For models including \mjseqn{p'} moderator/predictor variables, new moderator/predictor values (for \mjeqn{k_{new}}{k_new} hypothetical new studies/cases) must be specified by setting \code{newmods} equal to a \mjeqn{k_{new} \times p'}{k_new x p'} matrix with the corresponding new moderator/predictor values (if \code{newmods} is a vector, then only a single predicted value is computed unless the model only includes a single moderator/predictor, in which case predicted values corresponding to all the vector values are computed). If the model object includes an intercept (so that the model has \mjseqn{p' + 1} coefficients), then it will be automatically added to \code{newmods} unless one sets \code{intercept=FALSE}; alternatively, if \code{newmods} is a \mjeqn{k_{new} \times (p'+1)}{k_new x (p'+1)} matrix, then the \code{intercept} argument is ignored and the first column of the matrix determines whether the intercept is included when computing the predicted values or not. If the matrix specified via \code{newmods} has row names, then these are used to label the predicted values in the output.

   When computing multiple predicted values, one can set \code{adjust=TRUE} to obtain confidence intervals whose width is adjusted based on a Bonferroni correction (e.g., instead of 95\% CIs, the function provides (100-5/\mjeqn{k_{new}}{k_new})\% CIs, where \mjeqn{k_{new}}{k_new} denotes the number of predicted values computed).
}
\value{
   An object of class \code{c("predict.matreg","list.rma")}. The object is a list containing the following components:
   \item{pred}{predicted value(s).}
   \item{se}{corresponding standard error(s).}
   \item{ci.lb}{lower bound of the confidence interval(s).}
   \item{ci.ub}{upper bound of the confidence interval(s).}
   \item{X}{the moderator/predictor value(s) used to calculate the predicted values (only when \code{addx=TRUE}).}
   \item{\dots}{some additional elements/values.}

   If \code{vcov=TRUE}, then the returned object is a list with the first element equal to the one as described above and the second element equal to the variance-covariance matrix of the predicted values.

   The object is formatted and printed with the \code{\link[=print.list.rma]{print}} function. To format the results as a data frame, one can use the \code{\link[=as.data.frame.list.rma]{as.data.frame}} function.
}
\note{
   Under the \sQuote{Regular \mjseqn{R} Matrix} case (see \code{\link{matreg}}), confidence intervals are constructed based on the critical values from a t-distribution with \mjseqn{k-p} degrees of freedom, where \mjseqn{p} denotes the total number of coefficients (i.e., including the intercept term if the model includes one). Otherwise, critical values from a standard normal distribution (i.e., \mjeqn{\pm 1.96}{±1.96} for \code{level=95}) are used.

   When using the \code{transf} option, the transformation is applied to the predicted values and the corresponding interval bounds. The standard errors are omitted from the printed output. Also, \code{vcov=TRUE} is ignored when using the \code{transf} option.
}
\author{
   Wolfgang Viechtbauer (\email{wvb@metafor-project.org}, \url{https://www.metafor-project.org}).
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\examples{
### fit a regression model with lm() to the 'mtcars' dataset
res <- lm(mpg ~ hp + wt + am, data=mtcars)

### obtain a predicted value
predict(res, newdata=data.frame(hp=120, wt=4.2, am=1), interval="confidence")

### covariance matrix of the dataset
S <- cov(mtcars)

### fit the same regression model using matreg()
res <- matreg(mpg ~ hp + wt + am, R=S, cov=TRUE,
              means=colMeans(mtcars), n=nrow(mtcars))

### obtain the same predicted value
predict(res, newmods=c(hp=120, wt=4.2, am=1))
}
\keyword{models}
