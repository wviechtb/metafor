\name{dat.crede2010}
\docType{data}
\alias{dat.crede2010}
\title{Studies on the Relationship between Class Attendance and Grades in College Students}
\description{Results from 68 studies on the relationship between class attendence and class performance and/or grade point average in college students.}
\usage{dat.crede2010}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{studyid}   \tab \code{numeric}   \tab study number \cr
\bold{year}      \tab \code{numeric}   \tab publication year \cr
\bold{source}    \tab \code{character} \tab study source (journal, dissertation, other) \cr
\bold{sampleid}  \tab \code{numeric}   \tab sample within study number \cr
\bold{criterion} \tab \code{character} \tab criterion variable (grade, gpa) \cr
\bold{class}     \tab \code{character} \tab class type (science, nonscience) \cr
\bold{ni}        \tab \code{numeric}   \tab sample size \cr
\bold{ri}        \tab \code{numeric}   \tab observed correlation
}
}
\details{
   The 68 studies included in this dataset provide information about the relationship between class attendance of college students and their performance (i.e., grade) in the class and/or their overall grade point average. Some studies included multiple samples and hence the dataset actually contains 97 correlation coefficients.

   The dataset was obtained via personal communication. Note that this dataset differs just slightly from the one used by Credé et al. (2010).
}
\source{
   Personal communication.
}
\references{
   Credé, M., Roch, S. G., & Kieszczynka, U. M. (2010). Class attendance in college: A meta-analytic review of the relationship of class attendance with grades and student characteristics. \emph{Review of Educational Research}, \bold{80}(2), 272--295. \verb{https://doi.org/10.3102/0034654310362998}
}
\examples{
### copy data into 'dat'
dat <- dat.crede2010

### calculate r-to-z transformed correlations and corresponding sampling variances
dat <- escalc(measure="ZCOR", ri=ri, ni=ni, data=dat)

############################################################################

### meta-analysis for the relationship between attendance and grades
res <- rma(yi, vi, data=dat, subset=criterion=="grade")
res

### estimated average correlation with 95\% CI/PI
predict(res, transf=transf.ztor, digits=2)

### examine if relationship between attendance and grades differs for nonscience/science classes
res <- rma(yi, vi, mods = ~ class, data=dat, subset=criterion=="grade")
res

### estimated average correlations for nonscience and science classes
predict(res, newmods=c(0,1), transf=transf.ztor, digits=2)

### examine if relationship between attendance and grades has changed over time
res <- rma(yi, vi, mods = ~ year, data=dat, subset=criterion=="grade")
res

############################################################################

### meta-analysis for the relationship between attendance and GPA
res <- rma(yi, vi, data=dat, subset=criterion=="gpa")
res

### estimated average correlation with 95\% CI/PI
predict(res, transf=transf.ztor, digits=2)

### examine if relationship between attendance and GPA has changed over time
res <- rma(yi, vi, mods = ~ year, data=dat, subset=criterion=="gpa")
res

############################################################################

### use a multilevel model to examine the relationship between attendance and grades
res <- rma.mv(yi, vi, random = ~ 1 | studyid/sampleid, data=dat, subset=criterion=="grade")
res
predict(res, transf=transf.ztor, digits=2)

### use a multilevel model to examine the relationship between attendance and gpa
res <- rma.mv(yi, vi, random = ~ 1 | studyid/sampleid, data=dat, subset=criterion=="gpa")
res
predict(res, transf=transf.ztor, digits=2)
}
\keyword{datasets}
