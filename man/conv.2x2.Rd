\name{conv.2x2}
\alias{conv.2x2}
\title{Reconstruct Cell Frequencies of \mjeqn{2 \times 2}{2x2} Tables}
\description{
   Function to reconstruct the cell frequencies of \mjeqn{2 \times 2}{2x2} tables based on other summary statistics. \loadmathjax
}
\usage{
conv.2x2(ori, ri, x2i, ni, n1i, n2i, correct=TRUE, data, include,
         var.names=c("ai","bi","ci","di"), append=TRUE, replace="ifna")
}
\arguments{
   \item{ori}{optional vector with the odds ratios corresponding to the tables.}
   \item{ri}{optional vector with the phi coefficients corresponding to the tables.}
   \item{x2i}{optional vector with the (signed) chi-square statistics corresponding to the tables.}
   \item{ni}{vector with the total sample sizes.}
   \item{n1i}{vector with the marginal counts for the outcome of interest on the first variable.}
   \item{n2i}{vector with the marginal counts for the outcome of interest on the second variable.}
   \item{correct}{optional logical (or vector thereof) to indicate whether chi-square statistics were computed using Yates's correction for continuity (the default is \code{TRUE}).}
   \item{data}{optional data frame containing the variables given to the arguments above.}
   \item{include}{optional (logical or numeric) vector to specify the subset of studies for which the cell frequencies should be reconstructed.}
   \item{var.names}{character vector with four elements to specify the names of the variables for the reconstructed cell frequencies (the default is \code{c("ai","bi","ci","di")}).}
   \item{append}{logical to specify whether the data frame provided via the \code{data} argument should be returned together with the reconstructed values (the default is \code{TRUE}).}
   \item{replace}{character string to specify how values in \code{var.names} should be replaced (only relevant when using the \code{data} argument and if variables in \code{var.names} already exist in the data frame). See the \sQuote{Value} section for more details.}
}
\details{
   For meta-analyses based on \mjeqn{2 \times 2}{2x2} table data, the problem often arises that some studies do not directly report the cell frequencies. The present function allows the reconstruction of such tables based on other summary statistics.

   In particular, assume that the data of interest for a particular study are of the form:
   \tabular{lcccccc}{
                         \tab \ics \tab variable 2, outcome + \tab \ics \tab variable 2, outcome - \tab \ics \tab total      \cr
   variable 1, outcome + \tab \ics \tab \code{ai}             \tab \ics \tab \code{bi}             \tab \ics \tab \code{n1i} \cr
   variable 1, outcome - \tab \ics \tab \code{ci}             \tab \ics \tab \code{di}             \tab \ics \tab            \cr
   total                 \tab \ics \tab \code{n2i}            \tab \ics \tab                       \tab \ics \tab \code{ni}}
   where \code{ai}, \code{bi}, \code{ci}, and \code{di} denote the cell frequencies (i.e., the number of individuals falling into a particular category), \code{n1i} (i.e., \code{ai+bi}) and \code{n2i} (i.e., \code{ai+ci}) are the marginal totals for the outcome of interest on the first and second variable, respectively, and \code{ni} is the total sample size (i.e., \code{ai+bi+ci+di}) of the study.

   For example, if variable 1 denotes two different groups (e.g., treated versus control) and variable 2 indicates whether a particular outcome of interest has occurred or not (e.g., death, complications, failure to improve under the treatment), then \code{n1i} denotes the number of individuals in the treatment group, but \code{n2i} is \emph{not} the number of individuals in the control group, but the total number of individuals who experienced the outcome of interest on variable 2. Note that the meaning of \code{n2i} is therefore different here compared to the \code{\link{escalc}} function (where \code{n2i} denotes \code{ci+di}).

   If a study does not report the cell frequencies, but it reports the total sample size (which can be specified via the \code{ni} argument), the two marginal counts (which can be specified via the \code{n1i} and \code{n2i} arguments), and some other statistic corresponding to the table, then it may be possible to reconstruct the cell frequencies. The present function currently allows this for three different cases:

   \enumerate{
   \item If the odds ratio \mjdeqn{OR = \frac{a_i d_i}{b_i c_i}}{ai*di/(bi*ci)} is known, then the cell frequencies can be reconstructed (Bonett, 2007). Odds ratios can be specified via the \code{ori} argument.
   \item If the phi coefficient \mjdeqn{\phi = \frac{a_i d_i - b_i c_i}{\sqrt{n_{1i}(n_i-n_{1i})n_{2i}(n_i-n_{2i})}}}{\phi = (ai*di-bi*ci) / \sqrt{n1i*(ni-n1i)*n2i*(ni-n2i)}} is known, then the cell frequencies can again be reconstructed (own derivation). Phi coefficients can be specified via the \code{ri} argument.
   \item If the chi-square statistic from Pearson's chi-square test of independence is known (which can be specified via the \code{x2i} argument), then it can be used to recalculate the phi coefficient and hence again the cell frequencies can be reconstructed. However, the chi-square statistic does not carry information about the sign of the phi coefficient. Therefore, values specified via the \code{x2i} argument can be positive or negative, which allows the specification of the correct sign. Also, when using a chi-square statistic as input, it is assumed that it was computed using Yates's correction for continuity (unless \code{correct=FALSE}). If the chi-square statistic is not known, but its p-value, one can first back-calculate the chi-square statistic using \code{qchisq(<p-value>, df=1, lower.tail=FALSE)}.
   }

   Typically, the odds ratio, phi coefficient, or chi-square statistic (or its p-value) that can be extracted from a study will be rounded to a certain degree. The calculations underlying the function are exact only for unrounded values. Rounding can therefore introduce some discrepancies.

   If a marginal total is unknown, then external information needs to be used to \sQuote{guestimate} the number of individuals that experienced the outcome of interest on this variable. Depending on the accuracy of such an estimate, the reconstructed cell frequencies will be more or less accurate and need to be treated with due caution.

   The true marginal counts also put constraints on the possible values for the odds ratio, phi coefficient, and chi-square statistic. If a marginal count is replaced by a guestimate which is not compatible with the given statistic, one or more reconstructed cell frequencies may be negative. The function issues a warning if this happens and sets the cell frequencies to \code{NA} for such a study.

   If only one of the two marginal counts is unknown but a 95\% CI for the odds ratio is also available, then the \href{https://cran.r-project.org/package=estimraw}{estimraw} package can also be used to reconstruct the corresponding cell frequencies (Di Pietrantonj, 2006; but see Veroniki et al., 2013, for some cautions).
}
\value{
   If the \code{data} argument was not specified or \code{append=FALSE}, a data frame with four variables called \code{var.names} with the reconstructed cell frequencies.

   If \code{data} was specified and \code{append=TRUE}, then the original data frame is returned. If \code{var.names[j]} (for \mjeqn{\textrm{j} \in \\\\{1, \ldots, 4\\\\}}{for j in \{1, ..., 4\}}) is a variable in \code{data} and \code{replace="ifna"}, then missing values in this variable are replaced with the estimated frequencies (where possible) and otherwise a new variable called \code{var.names[j]} is added to the data frame.

   If \code{replace="all"}, then all values in \code{var.names[j]} where a reconstructed cell frequency can be computed are replaced, even for cases where the value in \code{var.names[j]} is not missing.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Bonett, D. G. (2007). Transforming odds ratios into correlations for meta-analytic research. \emph{American Psychologist}, \bold{62}(3), 254--255. \verb{https://doi.org/10.1037/0003-066x.62.3.254}

   Di Pietrantonj, C. (2006). Four-fold table cell frequencies imputation in meta analysis. \emph{Statistics in Medicine}, \bold{25}(13), 2299--2322. \verb{https://doi.org/10.1002/sim.2287}

   Veroniki, A. A., Pavlides, M., Patsopoulos, N. A., & Salanti, G. (2013). Reconstructing 2 x 2 contingency tables from odds ratios using the Di Pietrantonj method: Difficulties, constraints and impact in meta-analysis results. \emph{Research Synthesis Methods}, \bold{4}(1), 78--94. \verb{https://doi.org/10.1002/jrsm.1061}

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{escalc}} for a function to compute various effect size measures based on \mjeqn{2 \times 2}{2x2} table data.
}
\examples{
### demonstration that the reconstruction of the 2x2 table works
### (note: the values in rows 2, 3, and 4 correspond to those in row 1)
dat <- data.frame(ai=c(36,NA,NA,NA), bi=c(86,NA,NA,NA), ci=c(20,NA,NA,NA), di=c(98,NA,NA,NA),
                  oddsratio=NA, phi=NA, chisq=NA, ni=NA, n1i=NA, n2i=NA)
dat$oddsratio[2] <- round(exp(escalc(measure="OR", ai=ai, bi=bi, ci=ci, di=di, data=dat)$yi[1]), 2)
dat$phi[3] <- round(escalc(measure="PHI", ai=ai, bi=bi, ci=ci, di=di, data=dat)$yi[1], 2)
dat$chisq[4] <- round(chisq.test(matrix(c(t(dat[1,1:4])), nrow=2, byrow=TRUE))$statistic, 2)
dat$ni[2:4]  <- with(dat, ai[1] + bi[1] + ci[1] + di[1])
dat$n1i[2:4] <- with(dat, ai[1] + bi[1])
dat$n2i[2:4] <- with(dat, ai[1] + ci[1])
dat

### reconstruct cell frequencies for rows 2, 3, and 4
dat <- conv.2x2(ri=phi, ori=oddsratio, x2i=chisq, ni=ni, n1i=n1i, n2i=n2i, data=dat)
dat

### same example but with cell frequencies that are 10 times as large
dat <- data.frame(ai=c(360,NA,NA,NA), bi=c(860,NA,NA,NA), ci=c(200,NA,NA,NA), di=c(980,NA,NA,NA),
                  oddsratio=NA, phi=NA, chisq=NA, ni=NA, n1i=NA, n2i=NA)
dat$oddsratio[2] <- round(exp(escalc(measure="OR", ai=ai, bi=bi, ci=ci, di=di, data=dat)$yi[1]), 2)
dat$phi[3] <- round(escalc(measure="PHI", ai=ai, bi=bi, ci=ci, di=di, data=dat)$yi[1], 2)
dat$chisq[4] <- round(chisq.test(matrix(c(t(dat[1,1:4])), nrow=2, byrow=TRUE))$statistic, 2)
dat$ni[2:4]  <- with(dat, ai[1] + bi[1] + ci[1] + di[1])
dat$n1i[2:4] <- with(dat, ai[1] + bi[1])
dat$n2i[2:4] <- with(dat, ai[1] + ci[1])
dat <- conv.2x2(ri=phi, ori=oddsratio, x2i=chisq, ni=ni, n1i=n1i, n2i=n2i, data=dat)
dat # slight inaccuracy in row 3 due to rounding

### demonstrate what happens when a true marginal count is guestimated
escalc(measure="PHI", ai=176, bi=24, ci=72, di=128)
conv.2x2(ri=0.54, ni=400, n1i=200, n2i=248) # using the true marginal counts
conv.2x2(ri=0.54, ni=400, n1i=200, n2i=200) # marginal count for variable 2 is guestimated
conv.2x2(ri=0.54, ni=400, n1i=200, n2i=50)  # marginal count for variable 2 is incompatible

### demonstrate that using the correct sign for the chi-square statistic is important
chisq <- round(chisq.test(matrix(c(40,60,60,40), nrow=2, byrow=TRUE))$statistic, 2)
conv.2x2(x2i=-chisq, ni=200, n1i=100, n2i=100) # correct reconstruction
conv.2x2(x2i=chisq, ni=200, n1i=100, n2i=100) # incorrect reconstruction

### demonstrate use of the 'correct' argument
tab <- matrix(c(28,14,12,18), nrow=2, byrow=TRUE)
chisq <- round(chisq.test(tab)$statistic, 2) # chi-square test with Yates' continuity correction
conv.2x2(x2i=chisq, ni=72, n1i=42, n2i=40) # correct reconstruction
chisq <- round(chisq.test(tab, correct=FALSE)$statistic, 2) # without Yates' continuity correction
conv.2x2(x2i=chisq, ni=72, n1i=42, n2i=40) # incorrect reconstruction
conv.2x2(x2i=chisq, ni=72, n1i=42, n2i=40, correct=FALSE) # correct reconstruction

### recalculate chi-square statistic based on p-value
pval <- round(chisq.test(tab)$p.value, 2)
chisq <- qchisq(pval, df=1, lower.tail=FALSE)
conv.2x2(x2i=chisq, ni=72, n1i=42, n2i=40)
}
\keyword{manip}
