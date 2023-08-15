\name{mfopt}
\alias{mfopt}
\alias{getmfopt}
\alias{setmfopt}
\title{Getting and Setting Package Options}
\description{
   Functions for getting and setting \pkg{metafor} package options. \loadmathjax
}
\usage{
getmfopt(x, default=NULL)
setmfopt(...)
}
\arguments{
   \item{x}{The name of an option. If unspecified, all options are returned.}
   \item{default}{value to return if the option name does not exist.}
   \item{\dots}{one or more option names and the corresponding values to which they should be set.}
}
\details{
   The \pkg{metafor} package stores some of its options as a list element called \code{"metafor"} in the system options (see \code{\link{options}}). Hence, \code{getmfopt()} is the same as \code{getOption("metafor")}. One can also set \code{x} to the name of an option to return. With \code{setmfopt()}, one can set one or more options to their desired values.

   Currently, the following options are supported:

   \describe{
      \item{\code{check}}{logical to specify whether a version check should be carried out when loading the package (the default is \code{TRUE}). See \link[=misc-options]{here} for details. Obviously, this option must be set before loading the package (e.g., with \code{options(metafor=list(check=FALSE))}).}
      \item{\code{silent}}{logical to specify whether a startup message should be issued when loading the package (the default is \code{FALSE}). Obviously, this option must be set before loading the package (e.g., with \code{options(metafor=list(silent=TRUE))}). Note that messages about required packages that are automatically loaded are not suppressed by this. To fully suppress all startup messages, load the package with \code{\link{suppressPackageStartupMessages}}.}
      \item{\code{space}}{logical to specify whether an empty line should be added before and after the output (the default is \code{TRUE}). See \link[=misc-options]{here} for details.}
      \item{\code{digits}}{a named vector to specify how various aspects of the output should be rounded (unset by default). See \link[=misc-options]{here} for details.}
      \item{\code{style}}{a list whose elements specify the styles for various parts of the output when the \href{https://cran.r-project.org/package=crayon}{crayon} package is loaded and a terminal is used that support \sQuote{ANSI} color/highlight codes (unset by default). See \link[=misc-options]{here} for details.}
      \item{\code{theme}}{character string to specify how plots created by the package should be themed. The default is \code{"default"}, which means that the default foreground and background colors of plotting devices are used. Alternative options are \code{"light"} and \code{"dark"}, which forces plots to be drawn with a light or dark background, respectively. See \link[=misc-options]{here} for further details. RStudio users can also set this to \code{"auto"}, in which case plotting colors are chosen depending on the RStudio theme used (for some themes, using \code{"auto2"} might be aesthetically more pleasing). One can also use \code{setmfopt(theme="custom", fg=<color>, bg=<color>)} to set the foreground and background colors to custom choices (depending on the colors chosen, using \code{"custom2"} might be aesthetically more pleasing).}
   }

}
\value{
   Either a vector with the value for the chosen option or a list with all options.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\examples{
getmfopt()
getmfopt(space)
setmfopt(space=FALSE)
getmfopt()
setmfopt(space=TRUE)
getmfopt()
}
\keyword{manip}
