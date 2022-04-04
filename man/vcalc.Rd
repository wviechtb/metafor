\name{vcalc}
\alias{vcalc}
\title{Construct or Approximate the Variance-Covariance Matrix of Dependent Effect Sizes or Outcomes}
\description{
   The function can be used to construct or approximate the variance-covariance matrix of dependent effect sizes or outcomes, or more precisely, of their sampling errors (i.e., the \code{V} matrix in \code{\link{rma.mv}}). \loadmathjax
}
\usage{
vcalc(vi, cluster, subgroup, obs, type, time1, time2, grp1, grp2, w1, w2,
      data, rho, phi, rvars, checkpd=TRUE, nearpd=FALSE, \dots)
}
\arguments{
   \item{vi}{numeric vector to specify the sampling variances of the observed effect sizes or outcomes.}
   \item{cluster}{vector to specify the clustering variable (e.g., study).}
   \item{subgroup}{optional vector to specify different (independent) subgroups within clusters.}
   \item{obs}{optional vector to distinguish different observed effect sizes or outcomes corresponding to the same construct or response/dependent variable.}
   \item{type}{optional vector to distinguish different types of constructs or response/dependent variables underlying the observed effect sizes or outcomes.}
   \item{time1}{optional numeric vector to specify the time points when the observed effect sizes or outcomes were obtained (in the first condition if the observed effect sizes or outcomes represent contrasts between two conditions).}
   \item{time2}{optional numeric vector to specify the time points when the observed effect sizes or outcomes were obtained in the second condition (only relevant when the observed effect sizes or outcomes represent contrasts between two conditions).}
   \item{grp1}{optional vector to specify the group of the first condition when the observed effect sizes or outcomes represent contrasts between two conditions.}
   \item{grp2}{optional vector to specify the group of the second condition when the observed effect sizes or outcomes represent contrasts between two conditions.}
   \item{w1}{optional numeric vector to specify the size of the group (or more generally, the inverse-sampling variance weight) of the first condition when the observed effect sizes or outcomes represent contrasts between two conditions.}
   \item{w2}{optional numeric vector to specify the size of the group (or more generally, the inverse-sampling variance weight) of the second condition when the observed effect sizes or outcomes represent contrasts between two conditions.}
   \item{data}{optional data frame containing the variables given to the arguments above.}
   \item{rho}{argument to specify the correlation(s) of observed effect sizes or outcomes measured concurrently. See \sQuote{Details}.}
   \item{phi}{argument to specify the autocorrelation of observed effect sizes or outcomes measured at different time points. See \sQuote{Details}.}
   \item{rvars}{optional argument for specifying the variables that correspond to the correlation matrices of the studies (if this is specified, all arguments above except for \code{cluster} and \code{subgroup} are ignored). See \sQuote{Details}.}
   \item{checkpd}{logical to specify whether to check that the variance-covariance matrices within clusters are positive definite (the default is \code{TRUE}). See \sQuote{Note}.}
   \item{nearpd}{logical to specify whether the \code{\link[Matrix]{nearPD}} function from the \href{https://cran.r-project.org/package=Matrix}{Matrix} package should be used on variance-covariance matrices that are not positive definite. See \sQuote{Note}.}
   \item{\dots}{other arguments.}
}
\details{
   Standard meta-analytic models (such as those that can be fitted with the \code{\link{rma.uni}} function) assume that the observed effect sizes or outcomes (or more precisely, their sampling errors) are independent. This assumption is typically violated whenever multiple observed effect sizes or outcomes are computed based on the same sample of subjects (or whatever the study units are) or if there is at least partial overlap of subjects that contribute information to the computation of multiple effect sizes or outcomes.

   The present function can be used to construct or approximate the variance-covariance matrix of the sampling errors of dependent effect sizes or outcomes for a wide variety of circumstances (this variance-covariance matrix is the so-called \code{V} matrix that may be needed as input for multilevel/multivariate meta-analytic models as can be fitted with the \code{\link{rma.mv}} function; see also \link[=misc-recs]{here} for some recommendations on a general workflow for meta-analyses involving complex dependency structures).

   Argument \code{cluster} is used to specify the clustering variable. Rows with the same value of this variable are allowed to be dependent, while rows with different values are assumed to be independent. Typically, \code{cluster} will be a study identifier.

   Within the same cluster, there may be different subgroups with no overlap of subjects across subgroups. Argument \code{subgroup} can be used to distinguish such subgroups. Rows with the same value of this variable within a cluster are allowed to be dependent, while rows with different values are assumed to be independent even if they come from the same cluster. Therefore, from hereon, \sQuote{cluster} really refers to the combination of \code{cluster} and \code{subgroup}.

   Multiple effect sizes or outcomes belonging to the same cluster may be dependent due to a variety of reasons:
   \enumerate{
   \item The same construct of interest (e.g., severity of depression) may have been measured using different scales or instruments within a study (e.g., using the Beck Depression Inventory (BDI) and the Hamilton Depression Rating Scale (HDRS)) based on which multiple effect sizes can be computed for the same group of subjects (e.g., contrasting a treatment versus a control group with respect to each scale). In this case, we have multiple effect sizes that are different \sQuote{observations} of the effect with respect to the same type of construct.

   Argument \code{obs} is then used to distinguish different effect sizes corresponding to the same construct. If \code{obs} is specified, then argument \code{rho} must also be specified to indicate the degree of correlation among the sampling errors of the different effect sizes. Since this correlation is typically not known, the correlation among the various scales (or a rough \sQuote{guestimate} thereof) can be used as a proxy (i.e., the (typical) correlation between BDI and HDRS measurements).

   One can also specify an entire correlation matrix via \code{rho} to indicate, for each possible pair of \code{obs} values, the corresponding correlation. The row/column names of the matrix must then correspond to the unique values of the \code{obs} variable.

   \item Multiple types of constructs (or more generally, types of response/dependent variables) may have been measured in the same group of subjects (e.g., severity of depression as measured with the Beck Depression Inventory (BDI) and severity of anxiety as measured with the State-Trait Anxiety Inventory (STAI)). If this is of interest for a meta-analysis, effect sizes can then be computed with respect to each \sQuote{type} of construct.

   Argument \code{type} is then used to distinguish effect sizes corresponding to these different types of constructs. If \code{type} is specified, then argument \code{rho} must also be specified to indicate the degree of correlation among the sampling errors of effect sizes belonging to these different types. As above, the correlation among the various scales is typically used here as a proxy (i.e., the (typical) correlation between BDI and STAI measurements).

   One can also specify an entire correlation matrix via \code{rho} to indicate, for each possible pair of \code{type} values, the corresponding correlation. The row/column names of the matrix must then correspond to the unique values of the \code{type} variable.

   \item If there are multiple types of constructs, multiple scales or instruments may also have been used (in at least some of the studies) to measure the same construct and hence there may again be multiple effect sizes that are \sQuote{observations} of the same type of construct. Arguments \code{type} and \code{obs} should then be used together to indicate the various construct types and observations thereof. In this case, argument \code{rho} must be a vector of two values, the first to specify the within-construct correlation and the second to specify the between-construct correlation.

   One can also specify a list with two elements for \code{rho}, the first element being either a scalar or an entire correlation matrix for the within-construct correlation(s) and the second element being a scalar or an entire correlation matrix for the between-construct correlation(s). As above, any matrices specified must have row/column names corresponding to the unique values of the \code{obs} and/or \code{type} variables.

   \item The same construct and scale may have been assessed/used multiple times, allowing the computation of multiple effect sizes for the same group of subjects at different time points (e.g., right after the end of a treatment, at a short-term follow-up, and at a long-term follow-up). Argument \code{time1} is then used to specify the time points when the observed effect sizes were obtained. Argument \code{phi} must then also be specified to indicate the autocorrelation among the sampling errors of two effect sizes that differ by one unit on the \code{time1} variable. As above, the autocorrelation of the measurements themselves can be used here as a proxy.

   If multiple constructs and/or multiple scales have also been assessed at the various time points, then arguments \code{type} and/or \code{obs} (together with argument \code{rho}) should be used as needed to differentiate effect sizes corresponding to the different constructs and/or scales.

   \item Many effect size or outcome measures (e.g., raw or standardized mean differences, log-transformed ratios of means, log risk/odds ratios and risk differences) reflect the difference between two conditions (i.e., a contrast). Within a study, there may be more than two conditions, allowing the computation of multiple such contrasts (e.g., treatment A versus a control condition and treatment B versus the same control condition) and hence corresponding effect sizes. The reuse of information from the \sQuote{shared} condition (in this example, the control condition) then induces correlation among the effect sizes.

   To account for this, arguments \code{grp1} and \code{grp2} should be specified to indicate (within each cluster) which two groups were compared in the computation of each effect size (e.g., in the example above, the coding could be \code{grp1=c(1,2)} and \code{grp2=c(3,3)}; whether numbers or strings are used as identifiers is irrelevant).

   The degree of correlation between two contrast-type effect sizes that is induced by the use of a shared condition is a function of the size of the groups involved in the computation of the two effect sizes (or, more generally, the inverse-sampling variance weights of the condition-specific outcomes). By default, the group sizes (weights) are assumed to be identical across conditions, which implies a correlation of 0.5. If the group sizes (weights) are known, they can be specified via arguments \code{w1} and \code{w2} (in which case this information is used by the function to calculate a more accurate estimate of the correlation induced by the shared condition).

   In fact, a contrast-type effect size can be based on a between- or a within-subjects design. When at least one or more of the contrast-type effect sizes are based on a within-subjects design, then \code{time1} and \code{time2} should be used in combination with \code{grp1} and \code{grp2} to indicate for each effect size the group(s) and time point(s) involved.

   For example, \code{grp1=c(1,2)} and \code{grp2=c(3,3)} as above in combination with \code{time1=c(1,1)} and \code{time2=c(1,1)} would imply a between-subjects design involving three groups where two effect size were computed contrasting groups 1 and 2 versus group 3 at a single time point. On the other hand, \code{grp1=c(1,1)} and \code{grp2=c(1,1)} in combination with \code{time1=c(2,1)} and \code{time2=c(3,1)} would imply a within-subjects design where two effect size were computed contrasting time points 2 and 3 versus time point 1 in a single group. Argument \code{phi} is then used as above to indicate the autocorrelation of the measurements within groups (i.e., for the within-subjects design above, it would be the autocorrelation between time points 2 and 1 or equivalently, between time points 3 and 2).
   }

   All of the arguments above can be specified together to account for a fairly wide variety of dependency types.

   \subsection{Using the \code{rvars} Argument}{

      The function also provides an alternative approach for constructing the variance-covariance matrix using the \code{rvars} argument. Here, one must specify the names of the variables in the dataset that correspond to the correlation matrices of the studies (the variables should be specified as a vector; e.g., \code{c(var1, var2, var3)}).

      In particular, let \mjseqn{k_i} denote the number of rows corresponding to the \mjeqn{i\textrm{th}}{ith} cluster. Then the values of the first \mjseqn{k_i} variables from \code{rvars} are used to construct the correlation matrix and, together with the sampling variances (specified via \code{vi}), the variance-covariance matrix. Say there are three studies, the first with two correlated estimates, the second with a single estimate, and the third with four correlated estimates. Then the data structure should look like this:

\preformatted{study  yi  vi  r1  r2  r3  r4
=============================
    1   .   .   1  NA  NA  NA
    1   .   .  .6   1  NA  NA
-----------------------------
    2   .   .   1  NA  NA  NA
-----------------------------
    3   .   .   1  NA  NA  NA
    3   .   .  .8   1  NA  NA
    3   .   .  .5  .5   1  NA
    3   .   .  .5  .5  .8   1
=============================}

      with \code{rvars = c(r1, r2, r3, r4)}. If the \code{rvars} variables are a consecutive set in the data frame (as above), then one can use the shorthand notation \code{rvars = c(r1:r4)}, so \code{r1} denotes the first and \code{r4} the last variable in the set. Note that only the lower triangular part of the submatrices defined by the \code{rvars} variables is used.

      There must be as many variables specified via \code{rvars} as the number of rows in the \emph{largest} cluster (in smaller clusters, the non-relevant variables can just be set to \code{NA}; see above).

   }

}
\value{
   A \mjseqn{k \times k} variance-covariance matrix, where \mjseqn{k} denotes the length of the \code{vi} variable (i.e., the number of rows in the dataset).
}
\note{
   Depending on the data structure, the specified variables, and the specified values for \code{rho} and/or \code{phi}, it is possible that the constructed variance-covariance matrix is not positive definite within one or more clusters (this is checked when \code{checkpd=TRUE}, which is the default). If such non-positive definite submatrices are found, the reasons for this should be carefully checked since this might indicate misapplication of the function and/or the specification of implausible values for \code{rho} and/or \code{phi}.

   When setting \code{nearpd=TRUE}, the \code{\link[Matrix]{nearPD}} function from the \href{https://cran.r-project.org/package=Matrix}{Matrix} package is used on variance-covariance submatrices that are not positive definite. This should only be used cautiously and after understanding why these matrices are not positive definite.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{escalc}} for a function to compute the observed effect sizes or outcomes (and corresponding sampling variances) for which a variance-covariance matrix could be constructed.

   \code{\link{rcalc}} for a function to construct the variance-covariance matrix of dependent correlation coefficients.

   \code{\link{rma.mv}} for a model fitting function that can be used to meta-analyze dependent effect sizes or outcomes.
}
\examples{
############################################################################

### see help(dat.assink2016) for further details on this dataset

dat <- dat.assink2016
head(dat, 9)

### assume that the effect sizes within studies are correlated with rho=0.6
V <- vcalc(vi, cluster=study, obs=esid, data=dat, rho=0.6)

### show part of V matrix for studies 1 and 2
round(V[dat$study \%in\% c(1,2), dat$study \%in\% c(1,2)], 4)

### or show as list of matrices
blsplit(V, dat$study, round, 4)[1:2]

### use a correlation of 0.7 for effect sizes corresponding to the same type of
### delinquent behavior and a correlation of 0.5 for effect sizes corresponding
### to different types of delinquent behavior
V <- vcalc(vi, cluster=study, type=deltype, obs=esid, data=dat, rho=c(0.7, 0.5))
blsplit(V, dat$study, round, 3)[16]

### examine the correlation matrix for study 16
blsplit(V, dat$study, cov2cor)[16]

############################################################################

### see help(dat.ishak2007) for further details on this dataset

dat <- dat.ishak2007
head(dat, 5)

### create long format dataset
dat <- reshape(dat, direction="long", idvar="study", v.names=c("yi","vi"),
               varying=list(c(2,4,6,8), c(3,5,7,9)))
dat <- dat[order(study, time),]

### remove missing measurement occasions from dat
dat <- dat[!is.na(yi),]
rownames(dat) <- NULL

### show the data for the first 5 studies
head(dat, 8)

### construct the full (block diagonal) V matrix with an AR(1) structure
### assuming an autocorrelation of 0.97 as estimated by Ishak et al. (2007)
V <- vcalc(vi, cluster=study, time1=time, phi=0.97, data=dat)
V[1:8, 1:8]
cov2cor(V[1:8, 1:8])

### or show as a list of matrices
blsplit(V, dat$study)[1:5]
blsplit(V, dat$study, cov2cor)[1:5]

############################################################################

### see help(dat.kalaian1996) for further details on this dataset

dat <- dat.kalaian1996
head(dat, 12)

### construct the variance-covariance matrix assuming rho = 0.66 for effect sizes
### corresponding to the 'verbal' and 'math' outcome types
V <- vcalc(vi, cluster=study, type=outcome, data=dat, rho=0.66)
round(V[1:12,1:12], 4)

############################################################################

### see help(dat.berkey1998) for further details on this dataset

dat <- dat.berkey1998

### variables v1i and v2i correspond to the 2x2 var-cov matrices of the studies;
### so use these variables to construct the V matrix (note: since v1i and v2i are
### var-cov matrices and not correlation matrices, set vi=1 for all rows)
V <- vcalc(vi=1, cluster=author, rvars=c(v1i, v2i), data=dat)
V
round(cov2cor(V), 2)

### or show as a list of matrices
blsplit(V, dat$author, function(x) round(cov2cor(x), 2))

### construct the variance-covariance matrix assuming rho = 0.4 for effect sizes
### corresponding to the 'PD' and 'AL' outcome types
V <- vcalc(vi=vi, cluster=trial, type=outcome, data=dat, rho=0.4)
round(V,4)
cov2cor(V)

############################################################################

### see help(dat.knapp2017) for further details on this dataset

dat <- dat.knapp2017
dat[-c(1:2)]

### create variable that indicates the task and difficulty combination as increasing integers
dat$task.diff <- unlist(lapply(split(dat, dat$study), function(x) {
   task.int <- as.integer(factor(x$task))
   diff.int <- as.integer(factor(x$difficulty))
   diff.int[is.na(diff.int)] <- 1
   paste0(task.int, ".", diff.int)}))

### construct correlation matrix for two tasks with four different difficulties where the
### correlation is 0.4 for different difficulties of the same task, 0.7 for the same
### difficulty of different tasks, and 0.28 for different difficulties of different tasks
R <- matrix(0.4, nrow=8, ncol=8)
R[5:8,1:4] <- R[1:4,5:8] <- 0.28
diag(R[1:4,5:8]) <- 0.7
diag(R[5:8,1:4]) <- 0.7
diag(R) <- 1
rownames(R) <- colnames(R) <- paste0(rep(1:2, each=4), ".", 1:4)
R

### construct an approximate V matrix accounting for the use of shared groups and
### for correlations among tasks/difficulties as specified in the R matrix above
V <- vcalc(vi, cluster=study, grp1=group1, grp2=group2, w1=n_sz, w2=n_hc,
           obs=task.diff, rho=R, data=dat)
Vs <- blsplit(V, dat$study)
cov2cor(Vs[[3]])  # study with multiple SZ groups and a single HC group
cov2cor(Vs[[6]])  # study with two task types and multiple difficulties
cov2cor(Vs[[12]]) # study with multiple difficulties for the same task
cov2cor(Vs[[24]]) # study with separate rows for males and females
cov2cor(Vs[[29]]) # study with separate rows for three genotypes

############################################################################
}
\keyword{datagen}
