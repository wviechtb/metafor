\name{formatters}
\alias{formatters}
\alias{fmtp}
\alias{fmtx}
\alias{fmtt}
\title{Formatter Functions}
\description{
   Functions to format various types of outputs. \loadmathjax
}
\usage{
fmtp(p, digits=4, pname="", equal=FALSE, sep=FALSE, add0=FALSE, quote=FALSE)
fmtx(x, digits=4, flag="", quote=FALSE, \dots)
fmtt(val, tname, df, df1, df2, pval, digits=4,
     pname="p-val", format=1, sep=TRUE, quote=FALSE, call=FALSE, \dots)
}
\arguments{
   \emph{Arguments for \code{fmtp}:}

   \item{p}{vector of p-values to be formatted.}
   \item{digits}{integer to specify the number of decimal places to which the values should be rounded. For \code{fmmt}, can be a vector of length 2, to specify the number of digits for the test statistic and the p-value, respectively.}
   \item{pname}{string to add as a prefix to the p-value (e.g., something like \code{"p-val"} or just \code{"p"}).}
   \item{equal}{logical to specify whether an equal symbol should be shown before the p-value (when it is larger than the rounding cutoff).}
   \item{sep}{logical to specify whether a space should be added between \code{pname}, the equal/lesser symbol, and the p-value.}
   \item{add0}{logical to specify whether a 0 should be shown before the decimal point when the p-value is below the rounding cutoff.}
   \item{quote}{logical to specify whether formatted strings should be quoted when printed.}

   \emph{Arguments specific for \code{fmtx}:}

   \item{x}{vector of numeric values to be formatted.}
   \item{flag}{a character string giving a format modifier as defined for \code{\link{formatC}}.}

   \emph{Arguments specific for \code{fmtt}:}

   \item{val}{test statistic value to be formatted.}
   \item{tname}{character string for the name of the test statistic.}
   \item{df}{optional value for the degrees of freedom of the test statistic.}
   \item{df1}{optional value for the numerator degrees of freedom of the test statistic.}
   \item{df2}{optional value for the denominator degrees of freedom of the test statistic.}
   \item{pval}{the p-value corresponding to the test statistic.}
   \item{format}{either \code{1} or \code{2} to denote whether the degrees of freedom should be given before the test statistic (in parentheses) or after the test statistic.}
   \item{call}{logical to specify whether the formatted test result should be returned as a call or not.}
   \item{\dots}{other arguments.}
}
\details{
   The \code{fmtp} function takes one or multiple p-values as input and rounds them to the chosen number of digits. For p-values that are smaller than \code{10^(-digits)} (e.g., \code{0.0001} for \code{digits=4}), the value is shown to fall below this bound (e.g., \code{<.0001}).

   The \code{fmtx} function takes one or multiple numeric values as input and rounds them to the chosen number of digits, without using scientific notation and without dropping trailing zeros (using \code{\link{formatC}}).

   The \code{fmtt} function takes a single test statistic value as input (and, if applicable, its degrees of freedom via argument \code{df} or its numerator and denominator degrees of freedom via arguments \code{df1} and \code{df2}) and the corresponding p-value and formats it for printing. Two different formats are available (chosen via the \code{format} argument), one giving the degrees of freedom before the test statistic (in parentheses) and one after the test statistic.
}
\value{
   A character vector with the formatted values. By default (i.e., when \code{quote=FALSE}), formatted strings are not quoted when printed.
}
\note{
   The option in \code{fmtt} to return the formatted test result as a call can be useful when adding the output to a plot with \code{\link{text}} and one would like to use \code{\link{plotmath}} formatting for \code{tname}.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\examples{
# examples for fmtp()
fmtp(c(.0002, .00008), quote=TRUE, digits=4)
fmtp(c(.0002, .00008), quote=TRUE, digits=4, equal=TRUE)
fmtp(c(.0002, .00008), quote=TRUE, digits=4, equal=TRUE, sep=TRUE)
fmtp(c(.0002, .00008), quote=TRUE, digits=4, equal=TRUE, sep=TRUE, add0=TRUE)

# examples for fmtx()
fmtx(c(1.0002, 2.00008, 3.00004), digits=4)
fmtx(c(-1, 1), digits=4)
fmtx(c(-1, 1), digits=4, flag=" ")

# examples for fmtt()
fmtt(2.45, "z", pval=0.01429, digits=2)
fmtt(3.45, "z", pval=0.00056, digits=2)
fmtt(2.45, "t", df=23, pval=0.02232, digits=2)
fmtt(3.45, "t", df=23, pval=0.00218, digits=2)
fmtt(3.45, "t", df=23, pval=0.00218, digits=2, format=2)
fmtt(46.23, "Q", df=29, pval=0.0226, digits=2)
fmtt(46.23, "Q", df=29, pval=0.0226, digits=2, format=2)
fmtt(8.75, "F", df1=2, df2=35, pval=0.00083, digits=c(2,3))
fmtt(8.75, "F", df1=2, df2=35, pval=0.00083, digits=c(2,3), format=2, pname="p")
fmtt(8.75, "F", df1=2, df2=35, pval=0.00083, digits=c(2,3), format=2, pname="p", sep=FALSE)
}
\keyword{manip}
