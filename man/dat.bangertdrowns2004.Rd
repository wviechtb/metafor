\name{dat.bangertdrowns2004}
\docType{data}
\alias{dat.bangertdrowns2004}
\title{Studies on the Effectiveness of Writing-to-Learn Interventions}
\description{Results from 48 studies on the effectiveness of school-based writing-to-learn interventions on academic achievement. \loadmathjax}
\usage{dat.bangertdrowns2004}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{id}       \tab \code{numeric}   \tab study number \cr
\bold{author}   \tab \code{character} \tab study author(s) \cr
\bold{year}     \tab \code{numeric}   \tab publication year \cr
\bold{grade}    \tab \code{numeric}   \tab grade level (1 = elementary; 2 = middle; 3 = high-school; 4 = college) \cr
\bold{length}   \tab \code{numeric}   \tab treatment length (in weeks) \cr
\bold{minutes}  \tab \code{numeric}   \tab minutes per assignment \cr
\bold{wic}      \tab \code{numeric}   \tab writing tasks were completed in class (0 = no; 1 = yes) \cr
\bold{feedback} \tab \code{numeric}   \tab feedback on writing was provided (0 = no; 1 = yes) \cr
\bold{info}     \tab \code{numeric}   \tab writing contained informational components (0 = no; 1 = yes) \cr
\bold{pers}     \tab \code{numeric}   \tab writing contained personal components (0 = no; 1 = yes) \cr
\bold{imag}     \tab \code{numeric}   \tab writing contained imaginative components (0 = no; 1 = yes) \cr
\bold{meta}     \tab \code{numeric}   \tab prompts for metacognitive reflection (0 = no; 1 = yes) \cr
\bold{subject}  \tab \code{character} \tab subject matter \cr
\bold{ni}       \tab \code{numeric}   \tab total sample size of the study \cr
\bold{yi}       \tab \code{numeric}   \tab standardized mean difference \cr
\bold{vi}       \tab \code{numeric}   \tab corresponding sampling variance
}
}
\details{
   In each of the studies included in this meta-analysis, an experimental group (i.e., a group of students that received instruction with increased emphasis on writing tasks) was compared against a control group (i.e., a group of students that received conventional instruction) with respect to some content-related measure of academic achievement (e.g., final grade, an exam/quiz/test score). The outcome measure for this meta-analysis was the standardized mean difference (with positive values indicating a higher mean level of academic achievement in the intervention group).

   The standardized mean differences given here are bias-corrected and therefore differ slighty from the values reported in the article. Also, since only the total sample size is given in the article, the sampling variances were computed under the assumption that \mjeqn{n_{i1} = n_{i2} = n_i / 2}{n_i1 = n_i2 = n_i / 2}.
}
\source{
   Bangert-Drowns, R. L., Hurley, M. M., & Wilkinson, B. (2004). The effects of school-based writing-to-learn interventions on academic achievement: A meta-analysis. \emph{Review of Educational Research}, \bold{74}(1), 29--58. \verb{https://doi.org/10.3102/00346543074001029}
}
\examples{
### copy data into 'dat' and examine data
dat <- dat.bangertdrowns2004
dat

### random-effects model
res <- rma(yi, vi, data=dat)
res
}
\keyword{datasets}
