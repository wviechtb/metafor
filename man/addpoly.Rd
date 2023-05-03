\name{addpoly}
\alias{addpoly}
\title{Add Polygons to Forest Plots}
\description{
   The \code{addpoly} function can be used to add polygons, sometimes called \sQuote{diamonds}, to a forest plot, for example to indicate summary estimates for subgroups of studies or to indicate fitted/predicted values based on models involving moderators.
}
\usage{
addpoly(x, \dots)
}
\arguments{
   \item{x}{either an object of class \code{"rma"}, an object of class \code{"predict.rma"}, or the values at which polygons should be drawn. See \sQuote{Details}.}
   \item{\dots}{other arguments.}
}
\details{
   Currently, methods exist for three types of situations.

   In the first case, object \code{x} is a fitted model coming from the \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}, \code{\link{rma.glmm}}, or \code{\link{rma.mv}} functions. The model must either be an equal- or a random-effects model, that is, the model should not contain any moderators. The corresponding method is \code{\link{addpoly.rma}}. It can be used to add a polygon to an existing forest plot (usually at the bottom), showing the summary estimate (with its confidence interval) based on the fitted model.

   Alternatively, \code{x} can be an object of class \code{"predict.rma"} obtained with the \code{\link[=predict.rma]{predict}} function. In this case, polygons based on the predicted values are drawn. The corresponding method is \code{\link{addpoly.predict.rma}}.

   Alternatively, object \code{x} can be a vector with values at which one or more polygons should be drawn. The corresponding method is \code{\link{addpoly.default}}.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{addpoly.rma}}, \code{\link{addpoly.predict.rma}}, and \code{\link{addpoly.default}} for the specific method functions.

   \code{\link{forest}} for functions to draw forest plots to which polygons can be added.
}
\keyword{aplot}
