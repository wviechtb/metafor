\name{calc2x2}
\alias{calc2x2}
\title{Reconstruct Cell Frequencies of \mjeqn{2 \times 2}{2x2} Tables}
\description{
   Function to reconstruct the cell frequencies of \mjeqn{2 \times 2}{2x2} tables based on other summary statistics. \loadmathjax
}
\usage{
calc2x2(ri, ori, ni, n1i, n2i, data, include,
        var.names=c("ai","bi","ci","di"), append=TRUE, replace="ifna")
}
\arguments{
   \item{ri}{vector with the phi coefficients corresponding to the tables.}
   \item{ori}{vector with the odds ratios corresponding to the tables.}
   \item{ni}{vector with the total sample sizes.}
   \item{n1i}{vector with the marginal counts for the outcome of interest on the first variable.}
   \item{n2i}{vector with the marginal counts for the outcome of interest on the second variable.}
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

   If a study does not report the cell frequencies, but it reports the odds ratio for the table, the total sample size, and the two marginal counts, then the cell frequencies can be reconstructed (Bonett, 2007). If instead of the odds ratio, the phi coefficient is available (in addition to the other information), then again the cell frequencies can be reconstructed (own derivation). The required information for the reconstruction needs to be specified via the \code{ri}, \code{ori}, \code{ni}, \code{n1i}, and \code{n2i} arguments.

   Typically, the odds ratio or phi coefficient that can be extracted from a study will be rounded to a certain degree. The calculations underlying the function are exact only for unrounded values. Discrepancies introduced due to rounding should be relatively minor though.

   If a marginal total is unknown, then external information needs to be used to \sQuote{guestimate} the number of individuals that experienced the outcome of interest on this variable. Depending on the accuracy of such an estimate, the reconstructed cell frequencies will be more or less accurate and need to be treated with due caution. If only one of the two marginal counts is unknown but a 95\% CI for the odds ratio is also available, then the \href{https://cran.r-project.org/package=estimraw}{estimraw} package can also be used to reconstruct the corresponding cell frequencies (Di Pietrantonj, 2006, but see Veroniki et al., 2013, for some cautions).
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
### (note: the values in rows 2 and 3 correspond to those in row 1)
dat <- data.frame(ai=c(49,NA,NA), bi=c(23,NA,NA), ci=c(12,NA,NA), di=c(19,NA,NA),
                  phi=NA, oddsratio=NA, ni=NA, n1i=NA, n2i=NA)
dat$phi[2] <- round(escalc(measure="PHI", ai=ai, bi=bi, ci=ci, di=di, data=dat)$yi[1], 2)
dat$oddsratio[3] <- round(exp(escalc(measure="OR", ai=ai, bi=bi, ci=ci, di=di, data=dat)$yi[1]), 2)
dat$ni[2:3]  <- with(dat, ai[1] + bi[1] + ci[1] + di[1])
dat$n1i[2:3] <- with(dat, ai[1] + bi[1])
dat$n2i[2:3] <- with(dat, ai[1] + ci[1])
dat

### reconstruct cell frequencies for rows 2 and 3
dat <- calc2x2(ri=phi, ori=oddsratio, ni=ni, n1i=n1i, n2i=n2i, data=dat)
dat

### same example but with cell frequencies that are 10 times as large
dat <- data.frame(ai=c(490,NA,NA), bi=c(230,NA,NA), ci=c(120,NA,NA), di=c(190,NA,NA),
                  phi=NA, oddsratio=NA, ni=NA, n1i=NA, n2i=NA)
dat$phi[2] <- round(escalc(measure="PHI", ai=ai, bi=bi, ci=ci, di=di, data=dat)$yi[1], 2)
dat$oddsratio[3] <- round(exp(escalc(measure="OR", ai=ai, bi=bi, ci=ci, di=di, data=dat)$yi[1]), 2)
dat$ni[2:3]  <- with(dat, ai[1] + bi[1] + ci[1] + di[1])
dat$n1i[2:3] <- with(dat, ai[1] + bi[1])
dat$n2i[2:3] <- with(dat, ai[1] + ci[1])
dat <- calc2x2(ri=phi, ori=oddsratio, ni=ni, n1i=n1i, n2i=n2i, data=dat)
dat # slight inaccuracy in row 2 due to rounding
}
\keyword{manip}
