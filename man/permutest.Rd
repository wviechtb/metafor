\name{permutest}
\alias{permutest}
\alias{permutest.rma.uni}
\alias{permutest.rma.ls}
\title{Permutation Tests for 'rma.uni' Objects}
\description{
   The function carries out permutation tests for objects of class \code{"rma.uni"}. \loadmathjax
}
\usage{
permutest(x, \dots)

\method{permutest}{rma.uni}(x, exact=FALSE, iter=1000, permci=FALSE,
          progbar=TRUE, digits, control, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma.uni"}.}
   \item{exact}{logical to specify whether an exact permutation test should be carried out (the default is \code{FALSE}). See \sQuote{Details}.}
   \item{iter}{integer to specify the number of iterations for the permutation test when not doing an exact test (the default is \code{1000}).}
   \item{permci}{logical to specify whether permutation-based confidence intervals (CIs) should also be constructed (the default is \code{FALSE}). Can also be a vector of indices to specify for which coefficients a permutation-based CI should be obtained.}
   \item{progbar}{logical to specify whether a progress bar should be shown (the default is \code{TRUE}).}
   \item{digits}{optional integer to specify the number of decimal places to which the printed results should be rounded. If unspecified, the default is to take the value from the object.}
   \item{control}{list of control values for numerical comparisons (\code{comptol}) and for \code{\link{uniroot}} (i.e., \code{tol} and \code{maxiter}). The latter is only relevant when \code{permci=TRUE}. See \sQuote{Note}.}
   \item{\dots}{other arguments.}
}
\details{
   For models without moderators, the permutation test is carried out by permuting the signs of the observed effect sizes or outcomes. The (two-sided) p-value of the permutation test is then equal to the proportion of times that the absolute value of the test statistic under the permuted data is as extreme or more extreme than under the actually observed data. See Follmann and Proschan (1999) for more details.

   For models with moderators, the permutation test is carried out by permuting the rows of the model matrix (i.e., \mjseqn{X}). The (two-sided) p-value for a particular model coefficient is then equal to the proportion of times that the absolute value of the test statistic for the coefficient under the permuted data is as extreme or more extreme than under the actually observed data. Similarly, for the omnibus test, the p-value is the proportion of times that the test statistic for the omnibus test is as extreme or more extreme than the actually observed one. See Higgins and Thompson (2004) and Viechtbauer et al. (2015) for more details.

   \subsection{Exact versus Approximate Permutation Tests}{

      If \code{exact=TRUE}, the function will try to carry out an exact permutation test. An exact permutation test requires fitting the model to each possible permutation. However, the number of possible permutations increases rapidly with the number of outcomes/studies (i.e., \mjseqn{k}). For models without moderators, there are \mjseqn{2^k} possible permutations of the signs. Therefore, for \mjseqn{k=5}, there are 32 possible permutations, for \mjseqn{k=10}, there are already 1024, and for \mjseqn{k=20}, there are over one million such permutations.

      For models with moderators, the increase in the number of possible permutations may be even more severe. The total number of possible permutations of the model matrix is \mjseqn{k!}. Therefore, for \mjseqn{k=5}, there are 120 possible permutations, for \mjseqn{k=10}, there are 3,628,800, and for \mjseqn{k=20}, there are over \mjeqn{10^{18}}{10^18} permutations of the model matrix.

      Therefore, going through all possible permutations may become infeasible. Instead of using an exact permutation test, one can set \code{exact=FALSE} (which is also the default). In that case, the function approximates the exact permutation-based p-value(s) by going through a smaller number (as specified by the \code{iter} argument) of \emph{random} permutations. Therefore, running the function twice on the same data can yield (slightly) different p-values. Setting \code{iter} sufficiently large ensures that the results become stable. For full reproducibility, one can also set the seed of the random number generator before running the function (see \sQuote{Examples}). Note that if \code{exact=FALSE} and \code{iter} is actually larger than the number of iterations required for an exact permutation test, then an exact test will automatically be carried out.

      For models with moderators, the exact permutation test actually only requires fitting the model to each \emph{unique} permutation of the model matrix. The number of unique permutations will be smaller than \mjseqn{k!} when the model matrix contains recurring rows. This may be the case when only including categorical moderators (i.e., factors) in the model or when any quantitative moderators included in the model can only take on a small number of unique values. When \code{exact=TRUE}, the function therefore uses an algorithm to restrict the test to only the unique permutations of the model matrix, which may make the use of the exact test feasible even when \mjseqn{k} is large.

      When using random permutations, the function ensures that the very first permutation will always correspond to the original data. This avoids p-values equal to 0.

   }

   \subsection{Permutation-Based Confidence Intervals}{

      When \code{permci=TRUE}, the function also tries to obtain permutation-based confidence intervals (CIs) of the model coefficient(s). This is done by shifting the observed effect sizes or outcomes by some amount and finding the most extreme values for this amount for which the permutation-based test would just lead to non-rejection. This is computationally expensive and may take a long time to complete. For models with moderators, one can also set \code{permci} to a vector of indices to specify for which coefficient(s) a permutation-based CI should be obtained. When the algorithm fails to determine a particular CI bound, it will be shown as \code{NA} in the output.

   }

   \subsection{Permutation Tests for Location-Scale Models}{

      The function also works with location-scale models (see \code{\link{rma.uni}} for details on such models). Permutation tests will then be carried out for both the location and scale parts of the model. However, note that permutation-based CIs are not available for location-scale models.

   }

}
\value{
   An object of class \code{"permutest.rma.uni"}. The object is a list containing the following components:
   \item{pval}{p-value(s) based on the permutation test.}
   \item{QMp}{p-value for the omnibus test of moderators based on the permutation test.}
   \item{zval.perm}{values of the test statistics of the coefficients under the various permutations.}
   \item{b.perm}{the model coefficients under the various permutations.}
   \item{QM.perm}{the test statistic of the omnibus test of moderators under the various permutations.}
   \item{ci.lb}{lower bound of the confidence intervals for the coefficients (permutation-based when \code{permci=TRUE}).}
   \item{ci.ub}{upper bound of the confidence intervals for the coefficients (permutation-based when \code{permci=TRUE}).}
   \item{\dots}{some additional elements/values are passed on.}

   The results are formatted and printed with the \code{\link[=print.permutest.rma.uni]{print}} function. One can also use \code{\link[=coef.permutest.rma.uni]{coef}} to obtain the table with the model coefficients, corresponding standard errors, test statistics, p-values, and confidence interval bounds. The permutation distribution(s) can be plotted with the \code{\link[=plot.permutest.rma.uni]{plot}} function.
}
\note{
   The p-values obtained with permutation tests cannot reach conventional levels of statistical significance (i.e., \mjseqn{p \le .05}) when \mjseqn{k} is very small. In particular, for models without moderators, the smallest possible (two-sided) p-value is .0625 when \mjseqn{k=5} and .03125 when \mjseqn{k=6}. Therefore, the permutation test is only able to reject the null hypothesis at \mjseqn{\alpha=.05} when \mjseqn{k} is at least equal to 6. For models with moderators, the smallest possible (two-sided) p-value for a particular model coefficient is .0833 when \mjseqn{k=4} and .0167 when \mjseqn{k=5} (assuming that each row in the model matrix is unique). Therefore, the permutation test is only able to reject the null hypothesis at \mjseqn{\alpha=.05} when \mjseqn{k} is at least equal to 5. Consequently, permutation-based CIs can also only be obtained when \mjseqn{k} is sufficiently large.

   When the number of permutations required for the exact test is so large as to be essentially indistinguishable from infinity (e.g., \code{factorial(200)}), the function will terminate with an error.

   Determining whether a test statistic under the permuted data is as extreme or more extreme than under the actually observed data requires making \code{>=} or \code{<=} comparisons. To avoid problems due to the finite precision with which computers generally represent numbers (see \href{https://cran.r-project.org/doc/FAQ/R-FAQ.html#Why-doesn_0027t-R-think-these-numbers-are-equal_003f}{this} FAQ for details), the function uses a numerical tolerance (\code{control} argument \code{comptol}, which is set equal to \code{.Machine$double.eps^0.5} by default) when making such comparisons (e.g., instead of \code{sqrt(3)^2 >= 3}, which may evaluate to \code{FALSE}, we use \code{sqrt(3)^2 >= 3 - .Machine$double.eps^0.5}, which should evaluate to \code{TRUE}).

   When obtaining permutation-based CIs, the function makes use of \code{\link{uniroot}}. By default, the desired accuracy is set equal to \code{.Machine$double.eps^0.25} and the maximum number of iterations to \code{100}. The desired accuracy and the maximum number of iterations can be adjusted with the \code{control} argument (i.e., \code{control=list(tol=value, maxiter=value)}). Also, the interval searched for the CI bounds may be too narrow, leading to \code{NA} for a bound. In this case, one can try setting \code{control=list(distfac=value)} with a value larger than 1 to extend the interval (the value indicating a multiplicative factor by which to extend the width of the interval searched) or \code{control=list(extendInt="yes")} to allow \code{\link{uniroot}} to extend the interval dynamically.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Follmann, D. A., & Proschan, M. A. (1999). Valid inference in random effects meta-analysis. \emph{Biometrics}, \bold{55}(3), 732--737. \verb{https://doi.org/10.1111/j.0006-341x.1999.00732.x}

   Good, P. I. (2009). \emph{Permutation, parametric, and bootstrap tests of hypotheses} (3rd ed.). New York: Springer.

   Higgins, J. P. T., & Thompson, S. G. (2004). Controlling the risk of spurious findings from meta-regression. \emph{Statistics in Medicine}, \bold{23}(11), 1663--1682. \verb{https://doi.org/10.1002/sim.1752}

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}

   Viechtbauer, W., \enc{López-López}{Lopez-Lopez}, J. A., \enc{Sánchez-Meca}{Sanchez-Meca}, J., & \enc{Marín-Martínez}{Marin-Martinez}, F. (2015). A comparison of procedures to test for moderators in mixed-effects meta-regression models. \emph{Psychological Methods}, \bold{20}(3), 360--374. \verb{https://doi.org/10.1037/met0000023}
}
\seealso{
   \code{\link{rma.uni}} for the function to fit models for which permutation tests can be conducted.

   \code{\link{print.permutest.rma.uni}} and \code{\link{plot.permutest.rma.uni}} for the print and plot methods and \code{\link{coef.permutest.rma.uni}} for a method to extract the model results table.
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### random-effects model
res <- rma(yi, vi, data=dat)
res

\dontrun{
### permutation test (approximate and exact)
set.seed(1234) # for reproducibility
permutest(res)
permutest(res, exact=TRUE)
}

### mixed-effects model with two moderators (absolute latitude and publication year)
res <- rma(yi, vi, mods = ~ ablat + year, data=dat)
res

### permutation test (approximate only; exact not feasible)
\dontrun{

### permutation test (approximate)
set.seed(1234) # for reproducibility
permres <- permutest(res, iter=10000)
permres

### plot of the permutation distribution for absolute latitude
### dashed horizontal line: the observed value of the test statistic (in both tails)
### black curve: standard normal density (theoretical reference/null distribution)
### blue curve: kernel density estimate of the permutation distribution
### note: the tail area under the permutation distribution is larger
### than under a standard normal density (hence, the larger p-value)
plot(permres, beta=2, lwd=c(2,3,3,4), xlim=c(-5,5))
}
}
\keyword{models}
