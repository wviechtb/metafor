\name{dat.mcdaniel1994}
\docType{data}
\alias{dat.mcdaniel1994}
\title{Studies on the Validity of Employment Interviews}
\description{Results from 160 studies on the correlation between employment interview assessments and job performance.}
\usage{dat.mcdaniel1994}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{study}    \tab \code{numeric}   \tab study number \cr
\bold{ni}       \tab \code{numeric}   \tab sample size of the study \cr
\bold{ri}       \tab \code{numeric}   \tab observed correlation \cr
\bold{type}     \tab \code{character} \tab interview type (j = job-related, s = situational, p = psychological) \cr
\bold{struct}   \tab \code{character} \tab interview structure (u = unstructured, s = structured)
}
}
\details{
   The 160 studies provide data in terms of the correlation between employment interview performance and actual job performance. In addition, the interview type and the interview structure are indicated. 
   
   McDaniel et al. (1994) describe the interview type and structure variables as follows. "Questions in situational interviews [...] focus on the individual's ability to project what his or her behavior would be in a given situation. [...] Job-related interviews are those in which the interviewer is a personnel officer or hiring authority and the questions attempt to assess past behaviors and job-related information, but most questions are not considered situational. Psychological interviews are conducted by a psychologist, and the questions are intended to assess personal traits, such as dependability." In structured interviews, "the questions and acceptable responses were specified in advance and the responses were rated for appropriateness of content. [...] Unstructured interviews gather applicant information in a less systematic manner than do structured interviews. Although the questions may be specified in advance, they usually are not, and there is seldom a formalized scoring guide. Also, all persons being interviewed are not typically asked the same questions."
   
   The goal of the meta-analysis was to examine the overall criterion-related validity of employment interviews and to examine whether the validity depends on the type and structure of the interview.

   The data in this dataset were obtained from Table A.2 in Rothstein, Sutton, and Borenstein (2005, p. 325-329). Note that the \code{type} and \code{struct} variables contain some \code{NA}s.
}
\source{
   Rothstein, H. R., Sutton, A. J., & Borenstein, M. (Eds.). (2005). \emph{Publication bias in meta-analysis: Prevention, assessment, and adjustments}. Chichester, England: Wiley.
}
\references{
   McDaniel, M. A., Whetzel, D. L., Schmidt, F. L., & Maurer, S. D. (1994). The validity of employment interviews: A comprehensive review and meta-analysis. \emph{Journal of Applied Psychology}, \bold{79}, 599--616.
}
\examples{
### load data
dat <- get(data(dat.mcdaniel1994))

### calculate r-to-z transformed correlations and corresponding sampling variances
dat <- escalc(measure="ZCOR", ri=ri, ni=ni, data=dat)
dat

### meta-analysis of the transformed correlations using a random-effects model
res <- rma(yi, vi, data=dat)
res

### average correlation with 95\% CI
predict(res, transf=transf.ztor)

### mixed-effects model with interview type as factor
### note: job-related interviews is the reference level
rma(yi, vi, mods = ~ factor(type), data=dat)

### mixed-effects model with interview structure as factor
### note: structured interviews is the reference level
rma(yi, vi, mods = ~ factor(struct), data=dat)

### note: the interpretation of the results is difficult since all 
### situational interviews were structured, almost all psychological 
### interviews were unstructured, and actually for the majority of 
### the psychological interviews it was unknown whether the interview 
### was structured or unstructured
table(dat$type, dat$struct, useNA="always")

### meta-analysis of raw correlations using a random-effects model
res <- rma(measure="COR", ri=ri, ni=ni, data=dat.mcdaniel1994)
res
}
\keyword{datasets}
