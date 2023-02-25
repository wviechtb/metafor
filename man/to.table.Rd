\name{to.table}
\alias{to.table}
\title{Convert Data from Vector to Table Format}
\description{
   The function converts summary data in vector format to the corresponding table format. \loadmathjax
}
\usage{
to.table(measure, ai, bi, ci, di, n1i, n2i, x1i, x2i, t1i, t2i,
         m1i, m2i, sd1i, sd2i, xi, mi, ri, ti, sdi, ni, data, slab, subset,
         add=1/2, to="none", drop00=FALSE, rows, cols)
}
\arguments{
   \item{measure}{a character string to specify the effect size or outcome measure corresponding to the summary data supplied. See \sQuote{Details} and the documentation of the \code{\link{escalc}} function for possible options.}
   \item{ai}{vector with the \mjeqn{2 \times 2}{2x2} table frequencies (upper left cell).}
   \item{bi}{vector with the \mjeqn{2 \times 2}{2x2} table frequencies (upper right cell).}
   \item{ci}{vector with the \mjeqn{2 \times 2}{2x2} table frequencies (lower left cell).}
   \item{di}{vector with the \mjeqn{2 \times 2}{2x2} table frequencies (lower right cell).}
   \item{n1i}{vector with the group sizes or row totals (first group/row).}
   \item{n2i}{vector with the group sizes or row totals (second group/row).}
   \item{x1i}{vector with the number of events (first group).}
   \item{x2i}{vector with the number of events (second group).}
   \item{t1i}{vector with the total person-times (first group).}
   \item{t2i}{vector with the total person-times (second group).}
   \item{m1i}{vector with the means (first group or time point).}
   \item{m2i}{vector with the means (second group or time point).}
   \item{sd1i}{vector with the standard deviations (first group or time point).}
   \item{sd2i}{vector with the standard deviations (second group or time point).}
   \item{xi}{vector with the frequencies of the event of interest.}
   \item{mi}{vector with the frequencies of the complement of the event of interest or the group means.}
   \item{ri}{vector with the raw correlation coefficients.}
   \item{ti}{vector with the total person-times.}
   \item{sdi}{vector with the standard deviations.}
   \item{ni}{vector with the sample/group sizes.}
   \item{data}{optional data frame containing the variables given to the arguments above.}
   \item{slab}{optional vector with labels for the studies.}
   \item{subset}{optional (logical or numeric) vector to specify the subset of studies that should be included in the array returned by the function.}
   \item{add}{see the documentation of the \code{\link{escalc}} function.}
   \item{to}{see the documentation of the \code{\link{escalc}} function.}
   \item{drop00}{see the documentation of the \code{\link{escalc}} function.}
   \item{rows}{optional vector with row/group names.}
   \item{cols}{optional vector with column/outcome names.}
}
\details{
   The \code{\link{escalc}} function describes a wide variety of effect sizes or outcome measures that can be computed for a meta-analysis. The summary data used to compute those measures are typically contained in vectors, each element corresponding to a study. The \code{to.table} function takes this information and constructs an array of \mjseqn{k} tables from these data.

   For example, in various fields (such as the health and medical sciences), the response variable measured is often dichotomous (binary), so that the data from a study comparing two different groups can be expressed in terms of a \mjeqn{2 \times 2}{2x2} table, such as:
   \tabular{lcccccc}{
           \tab \ics \tab outcome 1 \tab \ics \tab outcome 2 \tab \ics \tab total      \cr
   group 1 \tab \ics \tab \code{ai} \tab \ics \tab \code{bi} \tab \ics \tab \code{n1i} \cr
   group 2 \tab \ics \tab \code{ci} \tab \ics \tab \code{di} \tab \ics \tab \code{n2i}
   } where \code{ai}, \code{bi}, \code{ci}, and \code{di} denote the cell frequencies (i.e., the number of individuals falling into a particular category) and \code{n1i} and \code{n2i} the row totals (i.e., the group sizes).

   The cell frequencies in \mjseqn{k} such \mjeqn{2 \times 2}{2x2} tables can be specified via the \code{ai}, \code{bi}, \code{ci}, and \code{di} arguments (or alternatively, via the \code{ai}, \code{ci}, \code{n1i}, and \code{n2i} arguments). The function then creates the corresponding \mjeqn{2 \times 2 \times k}{2*2*k} array of tables. The \code{measure} argument should then be set equal to one of the outcome measures that can be computed based on this type of data, such as \code{"RR"}, \code{"OR"}, \code{"RD"} (it is not relevant which specific measure is chosen, as long as it corresponds to the specified summary data). See the documentation of the \code{\link{escalc}} function for more details on the types of data formats available.

   The examples below illustrate the use of this function.
}
\value{
   An array with \mjseqn{k} elements each consisting of either 1 or 2 rows and an appropriate number of columns.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{escalc}} for a function to compute observed effect sizes or outcomes (and corresponding sampling variances) based on similar inputs.

   \code{\link{to.long}} for a function to turn similar inputs into a long format dataset.
}
\examples{
### create tables
dat <- to.table(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg,
                data=dat.bcg, slab=paste(author, year, sep=", "),
                rows=c("Vaccinated", "Not Vaccinated"), cols=c("TB+", "TB-"))
dat

### create tables
dat <- to.table(measure="IRR", x1i=x1i, x2i=x2i, t1i=t1i, t2i=t2i,
                data=dat.hart1999, slab=paste(study, year, sep=", "),
                rows=c("Warfarin Group", "Placebo/Control Group"))
dat

### create tables
dat <- to.table(measure="MD", m1i=m1i, sd1i=sd1i, n1i=n1i,
                m2i=m2i, sd2i=sd2i, n2i=n2i, data=dat.normand1999,
                slab=source, rows=c("Specialized Care", "Routine Care"))
dat
}
\keyword{manip}
