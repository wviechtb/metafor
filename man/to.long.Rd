\name{to.long}
\alias{to.long}
\title{Convert Data from Vector to Long Format}
\description{
   The function converts summary data in vector format to the corresponding long format. \loadmathjax
}
\usage{
to.long(measure, ai, bi, ci, di, n1i, n2i, x1i, x2i, t1i, t2i,
        m1i, m2i, sd1i, sd2i, xi, mi, ri, ti, sdi, ni, data, slab, subset,
        add=1/2, to="none", drop00=FALSE, vlong=FALSE, append=TRUE, var.names)
}
\arguments{
   \item{measure}{a character string to specify the effect size or outcome measure corresponding to the summary data supplied. See \sQuote{Details} and the documentation of the \code{\link{escalc}} function for possible options.}
   \item{ai}{vector to specify the \mjeqn{2 \times 2}{2x2} table frequencies (upper left cell).}
   \item{bi}{vector to specify the \mjeqn{2 \times 2}{2x2} table frequencies (upper right cell).}
   \item{ci}{vector to specify the \mjeqn{2 \times 2}{2x2} table frequencies (lower left cell).}
   \item{di}{vector to specify the \mjeqn{2 \times 2}{2x2} table frequencies (lower right cell).}
   \item{n1i}{vector to specify the group sizes or row totals (first group/row).}
   \item{n2i}{vector to specify the group sizes or row totals (second group/row).}
   \item{x1i}{vector to specify the number of events (first group).}
   \item{x2i}{vector to specify the number of events (second group).}
   \item{t1i}{vector to specify the total person-times (first group).}
   \item{t2i}{vector to specify the total person-times (second group).}
   \item{m1i}{vector to specify the means (first group or time point).}
   \item{m2i}{vector to specify the means (second group or time point).}
   \item{sd1i}{vector to specify the standard deviations (first group or time point).}
   \item{sd2i}{vector to specify the standard deviations (second group or time point).}
   \item{xi}{vector to specify the frequencies of the event of interest.}
   \item{mi}{vector to specify the frequencies of the complement of the event of interest or the group means.}
   \item{ri}{vector to specify the raw correlation coefficients.}
   \item{ti}{vector to specify the total person-times.}
   \item{sdi}{vector to specify the standard deviations.}
   \item{ni}{vector to specify the sample/group sizes.}
   \item{data}{optional data frame containing the variables given to the arguments above.}
   \item{slab}{optional vector with labels for the studies.}
   \item{subset}{optional (logical or numeric) vector to specify the subset of studies that should included in the data frame returned by the function.}
   \item{add}{see the documentation of the \code{\link{escalc}} function.}
   \item{to}{see the documentation of the \code{\link{escalc}} function.}
   \item{drop00}{see the documentation of the \code{\link{escalc}} function.}
   \item{vlong}{optional logical whether a very long format should be used (only relevant for \mjeqn{2 \times 2}{2x2} or \mjeqn{1 \times 2}{1x2} table data).}
   \item{append}{logical to specify whether the data frame specified via the \code{data} argument (if one has been specified) should be returned together with the long format data (the default is \code{TRUE}).}
   \item{var.names}{optional vector with variable names (length depends on the data type). If unspecified, the function sets appropriate variable names by default.}
}
\details{
   The \code{\link{escalc}} function describes a wide variety of effect size or outcome measures that can be computed for a meta-analysis. The summary data used to compute those measures are typically contained in vectors, each element corresponding to a study. The \code{to.long} function takes this information and constructs a long format dataset from these data.

   For example, in various fields (such as the health and medical sciences), the response variable measured is often dichotomous (binary), so that the data from a study comparing two different groups can be expressed in terms of a \mjeqn{2 \times 2}{2x2} table, such as:
   \tabular{lccc}{
           \tab outcome 1 \tab outcome 2 \tab total      \cr
   group 1 \tab \code{ai} \tab \code{bi} \tab \code{n1i} \cr
   group 2 \tab \code{ci} \tab \code{di} \tab \code{n2i}
   } where \code{ai}, \code{bi}, \code{ci}, and \code{di} denote the cell frequencies (i.e., the number of people falling into a particular category) and \code{n1i} and \code{n2i} the row totals (i.e., the group sizes).

   The cell frequencies in \mjseqn{k} such \mjeqn{2 \times 2}{2x2} tables can be specified via the \code{ai}, \code{bi}, \code{ci}, and \code{di} arguments (or alternatively, via the \code{ai}, \code{ci}, \code{n1i}, and \code{n2i} arguments). The function then creates the corresponding long format dataset. The \code{measure} argument should then be set equal to one of the outcome measures that can be computed based on this type of data, such as \code{"RR"}, \code{"OR"}, \code{"RD"} (it is not relevant which specific measure is chosen, as long as it corresponds to the specified summary data). See the documentation of the \code{\link{escalc}} function for more details on the types of data formats available.

   The long format for data of this type consists of two rows per study, a factor indicating the study (default name \code{study}), a dummy variable indicating the group (default name \code{group}, coded as 1 and 2), and two variables indicating the number of individuals experiencing outcome 1 or outcome 2 (default names \code{out1} and \code{out2}). Alternatively, if \code{vlong=TRUE}, then the long format consists of four rows per study, a factor indicating the study (default name \code{study}), a dummy variable indicating the group (default name \code{group}, coded as 1 and 2), a dummy variable indicating the outcome (default name \code{outcome}, coded as 1 and 2), and a variable indicating the frequency of the respective outcome (default name \code{freq}).

   The default variable names can be changed via the \code{var.names} argument (must be of the appropriate length, depending on the data type).

   The examples below illustrate the use of this function.
}
\value{
   A data frame with either \mjseqn{k}, \mjeqn{2 \times k}{2*k}, or \mjeqn{4 \times k}{4*k} rows and an appropriate number of columns (depending on the data type) with the data in long format. If \code{append=TRUE} and a data frame was specified via the \code{data} argument, then the data in long format are appended to the original data frame (with rows repeated an appropriate number of times).
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{escalc}} for a function to compute observed effect sizes or outcomes (and corresponding sampling variances) based on similar inputs.

   \code{\link{to.table}} for a function to turn similar inputs into tabular form.
}
\examples{
### convert data to long format
dat <- to.long(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
dat

### extra long format
dat <- to.long(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, vlong=TRUE)
dat

### convert data to long format
dat <- to.long(measure="IRR", x1i=x1i, x2i=x2i, t1i=t1i, t2i=t2i,
               data=dat.hart1999, var.names=c("id", "group", "events", "ptime"))
dat

### convert data to long format
dat <- to.long(measure="MD", m1i=m1i, sd1i=sd1i, n1i=n1i,
               m2i=m2i, sd2i=sd2i, n2i=n2i, data=dat.normand1999,
               var.names=c("id", "group", "mean", "sd", "n"))
dat
}
\keyword{manip}
