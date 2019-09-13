\name{to.wide}
\alias{to.wide}
\title{Convert Data from a Long to a Wide Format}
\description{
   The function can be used to convert data given in long format to a wide format.
}
\usage{
to.wide(data, study, group, ref, grpvars, postfix=c(".1",".2"),
        addid=TRUE, addcomp=TRUE, adddesign=TRUE, minlen=2)
}
\arguments{
   \item{data}{a data frame in long format.}
   \item{study}{either the name (given as a character string) or the position (given as a single number) of the study variable in the data frame.}
   \item{group}{either the name (given as a character string) or the position (given as a single number) of the group variable in the data frame.}
   \item{ref}{optional character string to specify the reference group (must be one of the levels of the group variable). If not given, the first group/level of the group variable is used as the reference group.}
   \item{grpvars}{either the names (given as a character vector) or the positions (given as a numeric vector) of group-level variables.}
   \item{postfix}{a character string of length 2 giving the affix that is placed after the names of the group-level variables for the first and second group.}
   \item{addid}{logical indicating whether a row id variable (called \code{id}) should be added to the data frame (\code{TRUE} by default).}
   \item{addcomp}{logical indicating whether a comparison id variable (called \code{comp}) should be added to the data frame (\code{TRUE} by default).}
   \item{adddesign}{logical indicating whether a design id variable (called \code{design}) should be added to the data frame (\code{TRUE} by default).}
   \item{minlen}{integer indicating the minimum length of the abbreviated group names for the comparison and design id variables.}
}
\details{
   A meta-analytic dataset may be structured in a \sQuote{long} format, where each row in the dataset corresponds to a particular study group (e.g., treatment arm). Using this function, such a dataset can be restructured into a \sQuote{wide} format, where each group within a study is contrasted against a particular reference group (if a study does not include the chosen reference group, then another group is used as the reference group).

   The \code{study} and \code{group} arguments are used to specify the study and group variables in the dataset (either as character strings or as numbers indicating the column positions of these variables in the dataset). Argument \code{ref} is used to specify the reference group (this must be one of the possible levels of the \code{group} variable). Argument \code{grpvars} is used to specify (either as a character vector or by giving the column positions) of those variables in the dataset that correspond to group-level outcomes (the remaining variables are treated as study-level outcomes).

   The dataset is restructured so that a two-group study will yield a single row in the restructured dataset, contrasting the first group against the second group (either the reference group or whatever second group was included in the study). For studies with more than two groups (often called \sQuote{multiarm} studies in the medical literature), the reference group is repeated as many times as needed (so a four-group study would yield three rows in the restructured dataset).

   To distinguish the names of the group-level outcome variables for the two first and second group in the restructured dataset, the strings given for the \code{postfix} argument are placed after the respective variable names.

   The examples below illustrate the use of this function.
}
\value{
   A data frame with rows contrasting groups against a reference group and an appropriate number of columns (depending on the number of group-level outcome variables).
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{http://www.metafor-project.org/}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
   \code{\link{escalc}}, \code{\link{to.long}}, \code{\link{to.table}}
}
\examples{
### data in long format
dat <- dat.senn2013
dat <- dat[c(1,4,3,2,5,6)]
dat

### restructure to wide format
dat.c <- to.wide(dat, study="study", group="treatment", ref="placebo", grpvars=4:6)
dat.c

### data in long format
dat <- dat.hasselblad1998
dat

### restructure to wide format
dat.c <- to.wide(dat, study="study", group="trt", ref="no_contact", grpvars=6:7)
dat.c
}
\keyword{manip}
