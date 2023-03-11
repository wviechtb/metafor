\name{escalc}
\alias{escalc}
\title{Calculate Effect Sizes and Outcome Measures}
\description{
   The function can be used to calculate various effect sizes or outcome measures (and the corresponding sampling variances) that are commonly used in meta-analyses. \loadmathjax
}
\usage{
escalc(measure, ai, bi, ci, di, n1i, n2i, x1i, x2i, t1i, t2i,
       m1i, m2i, sd1i, sd2i, xi, mi, ri, ti, pi, sdi, r2i, ni, yi, vi, sei,
       data, slab, subset, include,
       add=1/2, to="only0", drop00=FALSE, vtype="LS",
       var.names=c("yi","vi"), add.measure=FALSE,
       append=TRUE, replace=TRUE, digits, \dots)
}
\arguments{
   \item{measure}{a character string to specify which effect size or outcome measure should be calculated. See \sQuote{Details} for possible options and how the data needed to compute the selected effect size or outcome measure should then be specified.}
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
   \item{ti}{vector with the total person-times or t-test statistics.}
   \item{pi}{vector with the (signed) p-values.}
   \item{sdi}{vector with the standard deviations.}
   \item{r2i}{vector with the \mjseqn{R^2} values.}
   \item{ni}{vector with the sample/group sizes.}
   \item{yi}{vector with the observed effect sizes or outcomes.}
   \item{vi}{vector with the corresponding sampling variances.}
   \item{sei}{vector with the corresponding standard errors.}
   \item{data}{optional data frame containing the variables given to the arguments above.}
   \item{slab}{optional vector with labels for the studies.}
   \item{subset}{optional (logical or numeric) vector to specify the subset of studies that will be included in the data frame returned by the function.}
   \item{include}{optional (logical or numeric) vector to specify the subset of studies for which the measure should be calculated. See the \sQuote{Value} section for more details.}
   \item{add}{a non-negative number to specify the amount to add to zero cells, counts, or frequencies. See \sQuote{Details}.}
   \item{to}{a character string to specify when the values under \code{add} should be added (either \code{"all"}, \code{"only0"}, \code{"if0all"}, or \code{"none"}). See \sQuote{Details}.}
   \item{drop00}{logical to specify whether studies with no cases/events (or only cases) in both groups should be dropped when calculating the observed effect sizes or outcomes. See \sQuote{Details}.}
   \item{vtype}{a character string to specify the type of sampling variances to calculate. See \sQuote{Details}.}
   \item{var.names}{character vector with two elements to specify the name of the variable for the observed effect sizes or outcomes and the name of the variable for the corresponding sampling variances (the defaults are \code{"yi"} and \code{"vi"}).}
   \item{add.measure}{logical to specify whether a variable should be added to the data frame (with default name \code{"measure"}) that indicates the type of outcome measure computed. When using this option, \code{var.names} can have a third element to change this variable name.}
   \item{append}{logical to specify whether the data frame provided via the \code{data} argument should be returned together with the observed effect sizes or outcomes and corresponding sampling variances (the default is \code{TRUE}).}
   \item{replace}{logical to specify whether existing values for \code{yi} and \code{vi} in the data frame should be replaced. Only relevant when \code{append=TRUE} and the data frame already contains the \code{yi} and \code{vi} variables. If \code{replace=TRUE} (the default), all of the existing values will be overwritten. If \code{replace=FALSE}, only \code{NA} values will be replaced. See the \sQuote{Value} section for more details.}
   \item{digits}{optional integer to specify the number of decimal places to which the printed results should be rounded. If unspecified, the default is 4. Note that the values are stored without rounding in the returned object. See also \link[=misc-options]{here} for further details on how to control the number of digits in the output.}
   \item{\dots}{other arguments.}
}
\details{
   Before a meta-analysis can be conducted, the relevant results from each study must be quantified in such a way that the resulting values can be further aggregated and compared. Depending on (a) the goals of the meta-analysis, (b) the design and types of studies included, and (c) the information provided therein, one of the various effect sizes or outcome measures described below may be appropriate for the meta-analysis and can be computed with the \code{escalc} function.

   The \code{measure} argument is a character string to specify the outcome measure that should be calculated (see below for the various options), arguments \code{ai} through \code{ni} are then used to specify the information needed to calculate the various measures (depending on the chosen outcome measure, different arguments need to be specified), and \code{data} can be used to specify a data frame containing the variables given to the previous arguments. The \code{add}, \code{to}, and \code{drop00} arguments may be needed when dealing with frequency or count data that may need special handling when some of the frequencies or counts are equal to zero (see below for details). Finally, the \code{vtype} argument is used to specify how the sampling variances should be estimated (again, see below for details).

   To provide a structure to the various effect sizes or outcome measures that can be calculated with the \code{escalc} function, we can distinguish between measures that are used to:

   \tabular{lll}{
   \ics \tab (1) \tab contrast two independent (either experimentally created or naturally occurring) groups, \cr
   \ics \tab (2) \tab describe the direction and strength of the association between two variables, \cr
   \ics \tab (3) \tab summarize some characteristic or attribute of individual groups, or \cr
   \ics \tab (4) \tab quantify change within a single group or the difference between two matched/paired samples.}

   Furthermore, where appropriate, we can further distinguish between measures that are applicable when the characteristic, response, or dependent variable assessed within the individual studies is:

   \tabular{lll}{
   \ics \tab (a) \tab a quantitative variable (e.g., amount of depression as assessed by a rating scale), \cr
   \ics \tab (b) \tab a dichotomous (binary) variable (e.g., remission versus no remission), \cr
   \ics \tab (c) \tab a count of events per time unit (e.g., number of migraines per year), or \cr
   \ics \tab (d) \tab a mix of the types above.}

   Below, these number and letter codes are used (also in combination) to make it easier to quickly find a measure suitable for a particular meta-analysis (e.g., search for \code{(1b)} to find measures that describe the difference between two groups with respect to a dichotomous variable or \code{(2a)} for measures that quantify the association between two quantitative variables).

   \subsection{(1) Outcome Measures for Two-Group Comparisons}{

      In many meta-analyses, the goal is to synthesize the results from studies that compare or contrast two groups. The groups may be experimentally defined (e.g., a treatment and a control group created via random assignment) or may occur naturally (e.g., men and women, employees working under high- versus low-stress conditions, people/animals/plants exposed to some environmental risk factor versus those not exposed).

      \subsection{(1a) Measures for Quantitative Variables}{

         When the response or dependent variable assessed within the individual studies is measured on some quantitative scale, it is customary to report certain summary statistics, such as the mean and standard deviation of the observations. The data layout for a study comparing two groups with respect to such a variable is then of the form:
         \tabular{lcccccc}{
                 \tab \ics \tab mean       \tab \ics \tab standard deviation \tab \ics \tab group size \cr
         group 1 \tab \ics \tab \code{m1i} \tab \ics \tab \code{sd1i}        \tab \ics \tab \code{n1i} \cr
         group 2 \tab \ics \tab \code{m2i} \tab \ics \tab \code{sd2i}        \tab \ics \tab \code{n2i}}
         where \code{m1i} and \code{m2i} are the observed means of the two groups, \code{sd1i} and \code{sd2i} are the observed standard deviations, and \code{n1i} and \code{n2i} denote the number of individuals in each group. The two groups may be experimentally created (e.g., a treatment and control group based on random assignment) or naturally occurring (e.g., men and women). In either case, the raw mean difference, the standardized mean difference, and the (log transformed) ratio of means (also called the log \sQuote{response ratio}) are useful outcome measures when meta-analyzing studies of this type.

         The options for the \code{measure} argument are then:
         \itemize{
         \item \code{"MD"} for the \emph{raw mean difference} (e.g., Borenstein, 2009),
         \item \code{"SMD"} for the \emph{standardized mean difference} (Hedges, 1981),
         \item \code{"SMDH"} for the \emph{standardized mean difference} with heteroscedastic population variances in the two groups (Bonett, 2008, 2009),
         \item \code{"SMD1"} for the \emph{standardized mean difference} where the mean difference is divided by the standard deviation of the second group (and \code{"SMD1H"} for the same but with heteroscedastic population variances),
         \item \code{"ROM"} for the \emph{log transformed ratio of means} (Hedges et al., 1999; Lajeunesse, 2011).
         }

         The raw mean difference is simply \code{(m1i - m2i)}, while the standardized mean difference is computed with \code{(m1i - m2i) / sdpi}. For \code{measure="SMD"}, \code{sdpi = sqrt(((n1i-1)*sd1i^2 + (n2i-1)*sd2i^2) / (n1i+n2i-2))} is the pooled standard deviation of the two groups. For \code{measure="SMDH"}, \code{sdpi = sqrt((sd1i^2 + sd2i^2) / 2)} is the square root of the average variance. Finally, for \code{measure="SMD1"} and \code{measure="SMD1H"}, we simply use \code{sdpi = sd2i} (note: for \code{measure="SMD1"}, only \code{sd2i} needs to be specified and \code{sd1i} is ignored).

         For \code{measure="SMD"}, the positive bias in the standardized mean difference (i.e., in a Cohen's d value) is automatically corrected for within the function, yielding Hedges' g (Hedges, 1981). Similarly, the same bias correction is applied for \code{measure="SMDH"} (Bonett, 2009), \code{measure="SMD1"} (Hedges, 1981), and \code{measure="SMD1H"}.

         For \code{measure="ROM"}, the log is taken of the ratio of means (i.e., \code{log(m1i/m2i)}), which makes this outcome measure symmetric around 0 and results in a sampling distribution that is closer to normality. Hence, this measure cannot be computed when \code{m1i} and \code{m2i} have opposite signs (in fact, this measure is only meant to be used for ratio scale measurements, where both means should be positive anyway).

         For \code{measure="SMD"}, if the means and standard deviations are unknown for some studies, but the standardized mean differences (Cohen's d values) are directly available (e.g., if they are reported in those studies), then these can be specified via argument \code{di}. Also, if the t-statistics from an independent samples (Student's) t-test are available for some studies, one can specify those values via argument \code{ti}, which are then transformed into the corresponding standardized mean differences within the function (the sign of the t-statistics is then taken to be the sign of the standardized mean differences). If only the (two-sided) p-values corresponding to the t-tests are known, one can specify those values via argument \code{pi} (which are then transformed into the t-statistics and then further into the standardized mean differences). However, since a two-sided p-value does not carry information about the sign of the test statistic (and hence neither about the standardized mean difference), the sign of the p-values (which can be negative) is used as the sign of the standardized mean differences (e.g., \code{escalc(measure="SMD", pi=-0.018, n1i=20, n2i=20)} yields a negative standardized mean difference of \code{-0.7664}). See \href{https://www.metafor-project.org/doku.php/tips:assembling_data_smd}{here} for a more detailed illustration of using the \code{ti} and \code{pi} arguments.

         For \code{measure="MD"}, one can choose between \code{vtype="LS"} (the default) and \code{vtype="HO"}. The former computes the sampling variances without assuming homoscedasticity (i.e., that the true variances of the measurements are the same in group 1 and group 2 within each study), while the latter assumes homoscedasticity (equations 12.5 and 12.3 in Borenstein, 2009, respectively). For \code{measure="SMD"}, one can choose between \code{vtype="LS"} (the default) for the usual large-sample approximation to compute the sampling variances (equation 8 in Hedges, 1982), \code{vtype="UB"} to compute unbiased estimates of the sampling variances (equation 9 in Hedges, 1983), \code{vtype="LS2"} to compute the sampling variances as described in Borenstein (2009) (i.e., equation 12.17), and \code{vtype="AV"} to compute the sampling variances with the usual large-sample approximation but plugging the sample-size weighted average of the Hedges' g values into the equation. The same choices also apply to \code{measure="SMD1"}. For \code{measure="ROM"}, one can choose between \code{vtype="LS"} (the default) for the usual large-sample approximation to compute the sampling variances (equation 1 in Hedges et al., 1999), \code{vtype="HO"} to compute the sampling variances assuming homoscedasticity (the unnumbered equation after equation 1 in Hedges et al., 1999), \code{vtype="AV"} to compute the sampling variances assuming homoscedasticity of the coefficient of variation within each group across studies, and \code{vtype="AVHO"} to compute the sampling variances assuming homoscedasticity of the coefficient of variation for both groups across studies.

         Datasets corresponding to data of this type are provided in \code{\link[metadat]{dat.normand1999}} and \code{\link[metadat]{dat.curtis1998}}.

         Interest may also be focused on differences between the two groups with respect to their variability. Here, the (log transformed) ratio of the coefficient of variation of the two groups (also called the coefficient of variation ratio) can be a useful measure (Nakagawa et al., 2015). If focus is solely on the variability of the measurements within the two groups, then the (log transformed) ratio of the standard deviations (also called the variability ratio) can be used (Nakagawa et al., 2015). For the latter, one only needs to specify \code{sd1i}, \code{sd2i}, \code{n1i}, and \code{n2i}. The options for the \code{measure} argument are:
         \itemize{
         \item \code{"CVR"} for the \emph{log transformed coefficient of variation ratio},
         \item \code{"VR"} for the \emph{log transformed variability ratio}.
         }

         Measure \code{"CVR"} is computed with \code{log(sd1i/m1i)-log(sd2i/m2i)}, while \code{"VR"} is simply \code{log(sd1i/sd2i)}, but note that a slight bias correction is applied for both of these measures (Nakagawa et al., 2015). Also, the sampling variance for \code{measure="CVR"} is computed as given by equation 12 in Nakagawa et al. (2015), but without the \sQuote{\mjseqn{-2 \rho \ldots}} terms, since for normally distributed data (which we assume here) the mean and variance (and transformations thereof) are independent.

      }

      \subsection{(1b) Measures for Dichotomous Variables}{

         In various fields of research (such as the health and medical sciences), the response variable measured within the individual studies is often dichotomous (binary), so that the data from a study comparing two different groups can be expressed in terms of a \mjeqn{2 \times 2}{2x2} table, such as:
         \tabular{lcccccc}{
                 \tab \ics \tab outcome 1 \tab \ics \tab outcome 2 \tab \ics \tab total      \cr
         group 1 \tab \ics \tab \code{ai} \tab \ics \tab \code{bi} \tab \ics \tab \code{n1i} \cr
         group 2 \tab \ics \tab \code{ci} \tab \ics \tab \code{di} \tab \ics \tab \code{n2i}}
         where \code{ai}, \code{bi}, \code{ci}, and \code{di} denote the cell frequencies (i.e., the number of individuals falling into a particular category) and \code{n1i} and \code{n2i} are the row totals (i.e., the group sizes).

         For example, in a set of randomized clinical trials, group 1 and group 2 may refer to the treatment and placebo/control group, respectively, with outcome 1 denoting some event of interest (e.g., death, complications, failure to improve under the treatment) and outcome 2 its complement. Similarly, in a set of cohort studies, group 1 and group 2 may denote those who engage in and those who do not engage in a potentially harmful behavior (e.g., smoking), with outcome 1 denoting the development of a particular disease (e.g., lung cancer) during the follow-up period. Finally, in a set of case-control studies, group 1 and group 2 may refer to those with the disease (i.e., cases) and those free of the disease (i.e., controls), with outcome 1 denoting, for example, exposure to some environmental risk factor in the past and outcome 2 non-exposure. Note that in all of these examples, the stratified sampling scheme fixes the row totals (i.e., the group sizes) by design.

         A meta-analysis of studies reporting results in terms of \mjeqn{2 \times 2}{2x2} tables can be based on one of several different outcome measures, including the risk ratio (also called the relative risk), the odds ratio, the risk difference, and the arcsine square root transformed risk difference (e.g., Fleiss & Berlin, 2009, \enc{Rücker}{Ruecker} et al., 2009). For any of these outcome measures, one needs to specify the cell frequencies via the \code{ai}, \code{bi}, \code{ci}, and \code{di} arguments (or alternatively, one can use the \code{ai}, \code{ci}, \code{n1i}, and \code{n2i} arguments).

         The options for the \code{measure} argument are then:
         \itemize{
         \item \code{"RR"} for the \emph{log risk ratio},
         \item \code{"OR"} for the \emph{log odds ratio},
         \item \code{"RD"} for the \emph{risk difference},
         \item \code{"AS"} for the \emph{arcsine square root transformed risk difference} (\enc{Rücker}{Ruecker} et al., 2009),
         \item \code{"PETO"} for the \emph{log odds ratio} estimated with Peto's method (Yusuf et al., 1985).
         }

         Note that the log is taken of the risk ratio and the odds ratio, which makes these outcome measures symmetric around 0 and results in corresponding sampling distributions that are closer to normality. Also, when multiplied by 2, the arcsine square root transformed risk difference is actually identical to Cohen's h (Cohen, 1988).

         Cell entries with a zero count can be problematic, especially for the risk ratio and the odds ratio. Adding a small constant to the cells of the \mjeqn{2 \times 2}{2x2} tables is a common solution to this problem. When \code{to="only0"} (the default), the value of \code{add} (the default is \code{1/2}; but see \sQuote{Note}) is added to each cell of those \mjeqn{2 \times 2}{2x2} tables with at least one cell equal to 0. When \code{to="all"}, the value of \code{add} is added to each cell of all \mjeqn{2 \times 2}{2x2} tables. When \code{to="if0all"}, the value of \code{add} is added to each cell of all \mjeqn{2 \times 2}{2x2} tables, but only when there is at least one \mjeqn{2 \times 2}{2x2} table with a zero cell. Setting \code{to="none"} or \code{add=0} has the same effect: No adjustment to the observed table frequencies is made. Depending on the outcome measure and the data, this may lead to division by zero (when this occurs, the resulting value is recoded to \code{NA}). Also, studies where \code{ai=ci=0} or \code{bi=di=0} may be considered to be uninformative about the size of the effect and dropping such studies has sometimes been recommended (Higgins et al., 2019). This can be done by setting \code{drop00=TRUE}. The values for such studies will then be set to \code{NA} (i.e., missing).

         Datasets corresponding to data of this type are provided in \code{\link[metadat]{dat.bcg}}, \code{\link[metadat]{dat.collins1985a}}, \code{\link[metadat]{dat.collins1985b}}, \code{\link[metadat]{dat.egger2001}}, \code{\link[metadat]{dat.hine1989}}, \code{\link[metadat]{dat.laopaiboon2015}}, \code{\link[metadat]{dat.lee2004}}, \code{\link[metadat]{dat.li2007}}, \code{\link[metadat]{dat.linde2005}}, \code{\link[metadat]{dat.nielweise2007}}, and \code{\link[metadat]{dat.yusuf1985}}.

         If the \mjeqn{2 \times 2}{2x2} table is not available (or cannot be reconstructed) for a study, but an odds ratio and a corresponding confidence interval is reported, one can easily transform these values into the corresponding log odds ratio and sampling variance (and combine such a study with those that do report \mjeqn{2 \times 2}{2x2} table data). See \href{https://www.metafor-project.org/doku.php/tips:assembling_data_or}{here} for an illustration/discussion of this.

      }

      \subsection{(1c) Measures for Event Counts}{

         In medical and epidemiological studies comparing two different groups (e.g., treated versus untreated patients, exposed versus unexposed individuals), results are sometimes reported in terms of event counts (i.e., the number of events, such as strokes or myocardial infarctions) over a certain period of time. Data of this type are also referred to as \sQuote{person-time data}. Assume that the studies report data in the form:
         \tabular{lcccc}{
                 \tab \ics \tab number of events \tab \ics \tab total person-time \cr
         group 1 \tab \ics \tab \code{x1i}       \tab \ics \tab \code{t1i} \cr
         group 2 \tab \ics \tab \code{x2i}       \tab \ics \tab \code{t2i}}
         where \code{x1i} and \code{x2i} denote the number of events in the first and the second group, respectively, and \code{t1i} and \code{t2i} the corresponding total person-times at risk. Often, the person-time is measured in years, so that \code{t1i} and \code{t2i} denote the total number of follow-up years in the two groups.

         This form of data is fundamentally different from what was described in the previous section, since the total follow-up time may differ even for groups of the same size and the individuals studied may experience the event of interest multiple times. Hence, different outcome measures than the ones described in the previous section need to be considered when data are reported in this format. These include the incidence rate ratio, the incidence rate difference, and the square root transformed incidence rate difference (Bagos & Nikolopoulos, 2009; Rothman et al., 2008). For any of these outcome measures, one needs to specify the total number of events via the \code{x1i} and \code{x2i} arguments and the corresponding total person-time values via the \code{t1i} and \code{t2i} arguments.

         The options for the \code{measure} argument are then:
         \itemize{
         \item \code{"IRR"} for the \emph{log incidence rate ratio},
         \item \code{"IRD"} for the \emph{incidence rate difference},
         \item \code{"IRSD"} for the \emph{square root transformed incidence rate difference}.
         }

         Note that the log is taken of the incidence rate ratio, which makes this outcome measure symmetric around 0 and results in a sampling distribution that is closer to normality.

         Studies with zero events in one or both groups can be problematic, especially for the incidence rate ratio. Adding a small constant to the number of events is a common solution to this problem. When \code{to="only0"} (the default), the value of \code{add} (the default is \code{1/2}; but see \sQuote{Note}) is added to \code{x1i} and \code{x2i} only in the studies that have zero events in one or both groups. When \code{to="all"}, the value of \code{add} is added to \code{x1i} and \code{x2i} in all studies. When \code{to="if0all"}, the value of \code{add} is added to \code{x1i} and \code{x2i} in all studies, but only when there is at least one study with zero events in one or both groups. Setting \code{to="none"} or \code{add=0} has the same effect: No adjustment to the observed number of events is made. Depending on the outcome measure and the data, this may lead to division by zero (when this occurs, the resulting value is recoded to \code{NA}). Like for \mjeqn{2 \times 2}{2x2} table data, studies where \code{x1i=x2i=0} may be considered to be uninformative about the size of the effect and dropping such studies has sometimes been recommended. This can be done by setting \code{drop00=TRUE}. The values for such studies will then be set to \code{NA}.

         Datasets corresponding to data of this type are provided in \code{\link[metadat]{dat.hart1999}} and \code{\link[metadat]{dat.nielweise2008}}.

      }

      \subsection{(1d) Transforming SMDs to ORs and Vice-Versa}{

         In some meta-analyses, one may encounter studies that contrast two groups with respect to a quantitative response variable (case 1a above) and other studies that contrast the same two groups with respect to a dichotomous variable (case 2b above). If both types of studies are to be combined in the same analysis, one needs to compute the same outcome measure across all studies.

         For this, one may need to transform standardized mean differences into log odds ratios (e.g., Cox & Snell, 1989; Chinn, 2000; Hasselblad & Hedges, 1995; \enc{Sánchez-Meca}{Sanchez-Meca} et al., 2003). Here, the data need to be specified as described under (1a) and the options for the \code{measure} argument are then:
         \itemize{
         \item \code{"D2ORN"} for the \emph{transformed standardized mean difference} assuming normal distributions,
         \item \code{"D2ORL"} for the \emph{transformed standardized mean difference} assuming logistic distributions.
         }

         Both of these transformations provide an estimate of the log odds ratio, the first assuming that the responses within the two groups are normally distributed, while the second assumes that the responses follow logistic distributions.

         Alternatively, assuming that the dichotomous outcome in a \mjeqn{2 \times 2}{2x2} table is actually a dichotomized version of the responses on an underlying quantitative scale, it is also possible to estimate the standardized mean difference based on \mjeqn{2 \times 2}{2x2} table data, using either the probit transformed risk difference or a transformation of the odds ratio (e.g., Cox & Snell, 1989; Chinn, 2000; Hasselblad & Hedges, 1995; \enc{Sánchez-Meca}{Sanchez-Meca} et al., 2003). Here, the data need to be specified as described under (1b) and the options for the \code{measure} argument are then:
         \itemize{
         \item \code{"PBIT"} for the \emph{probit transformed risk difference},
         \item \code{"OR2DN"} for the \emph{transformed odds ratio} assuming normal distributions,
         \item \code{"OR2DL"} for the \emph{transformed odds ratio} assuming logistic distributions.
         }

         All of these transformations provide an estimate of the standardized mean difference, the first two assuming that the responses on the underlying quantitative scale are normally distributed, while the third assumes that the responses follow logistic distributions.

         A dataset illustrating the combined analysis of standardized mean differences and probit transformed risk differences is provided in \code{\link[metadat]{dat.gibson2002}}.

      }

   }

   \subsection{(2) Outcome Measures for Variable Association}{

      Meta-analyses are often used to synthesize studies that examine the direction and strength of the association between two variables measured concurrently and/or without manipulation by experimenters. In this section, a variety of outcome measures will be discussed that may be suitable for a meta-analyses with this purpose. We can distinguish between measures that are applicable when both variables are measured on quantitative scales, when both variables measured are dichotomous, and when the two variables are of mixed types.

      \subsection{(2a) Measures for Two Quantitative Variables}{

         The (Pearson or product-moment) correlation coefficient quantifies the direction and strength of the (linear) relationship between two quantitative variables and is therefore frequently used as the outcome measure for meta-analyses. Two alternative measures are a bias-corrected version of the correlation coefficient and Fisher's r-to-z transformed correlation coefficient.

         For these measures, one needs to specify \code{ri}, the vector with the raw correlation coefficients, and \code{ni}, the corresponding sample sizes. The options for the \code{measure} argument are then:
         \itemize{
         \item \code{"COR"} for the \emph{raw correlation coefficient},
         \item \code{"UCOR"} for the \emph{raw correlation coefficient} corrected for its slight negative bias (based on equation 2.3 in Olkin & Pratt, 1958),
         \item \code{"ZCOR"} for \emph{Fisher's r-to-z transformed correlation coefficient} (Fisher, 1921).
         }

         If the correlation coefficient is unknown for some studies, but the t-statistics \mjseqn{t_i = r_i \sqrt{n_i - 2} / \sqrt{1 - r^2}} are available for those studies (for the standard t-test of \mjeqn{\mbox{H}_0{:}\; \rho = 0}{H_0: \rho = 0}), one can specify those values via argument \code{ti}, which are then transformed into the corresponding correlation coefficients within the function (the sign of the t-statistics is then taken to be the sign of the correlations). If only the (two-sided) p-values corresponding to the t-tests are known, one can specify those values via argument \code{pi}. However, since a two-sided p-value does not carry information about the sign of the test statistic (and hence neither about the correlation), the sign of the p-values (which can be negative) is used as the sign of the correlation coefficients (e.g., \code{escalc(measure="COR", pi=-.07, ni=30)} yields a negative correlation of \code{-0.3354}).

         For \code{measure="COR"} and \code{measure="UCOR"}, one can choose between \code{vtype="LS"} (the default) for the usual large-sample approximation to compute the sampling variances (i.e., plugging the (biased-corrected) correlation coefficients into equation 12.27 in Borenstein, 2009), \code{vtype="UB"} to compute unbiased estimates of the sampling variances (see Hedges, 1989, but using the exact equation instead of the approximation), and \code{vtype="AV"} to compute the sampling variances with the usual large-sample approximation but plugging the sample-size weighted average of the (bias-corrected) correlation coefficients into the equation.

         Datasets corresponding to data of this type are provided in \code{\link[metadat]{dat.mcdaniel1994}} and \code{\link[metadat]{dat.molloy2014}}.

         For meta-analyses involving multiple (dependent) correlations extracted from the same sample, see also the \code{\link{rcalc}} function.

      }

      \subsection{(2b) Measures for Two Dichotomous Variables}{

         When the goal of a meta-analysis is to examine the relationship between two dichotomous variables, the data for each study can again be presented in the form of a \mjeqn{2 \times 2}{2x2} table, except that there may not be a clear distinction between the grouping variable and the outcome variable. Moreover, the table may be a result of cross-sectional (i.e., multinomial) sampling, where none of the table margins (except the total sample size) are fixed by the study design.

         In particular, assume that the data of interest for a particular study are of the form:
         \tabular{lcccccc}{
                               \tab \ics \tab variable 2, outcome + \tab \ics \tab variable 2, outcome - \tab \ics \tab total      \cr
         variable 1, outcome + \tab \ics \tab \code{ai}             \tab \ics \tab \code{bi}             \tab \ics \tab \code{n1i} \cr
         variable 1, outcome - \tab \ics \tab \code{ci}             \tab \ics \tab \code{di}             \tab \ics \tab \code{n2i}}
         where \code{ai}, \code{bi}, \code{ci}, and \code{di} denote the cell frequencies (i.e., the number of individuals falling into a particular category) and \code{n1i} and \code{n2i} are the row totals.

         The phi coefficient and the odds ratio are commonly used measures of association for \mjeqn{2 \times 2}{2x2} table data (e.g., Fleiss & Berlin, 2009). The latter is particularly advantageous, as it is directly comparable to values obtained from stratified sampling (as described earlier). Yule's Q and Yule's Y (Yule, 1912) are additional measures of association for \mjeqn{2 \times 2}{2x2} table data (although they are not typically used in meta-analyses). Finally, assuming that the two dichotomous variables are actually dichotomized versions of the responses on two underlying quantitative scales (and assuming that the two variables follow a bivariate normal distribution), it is also possible to estimate the correlation between the two quantitative variables using the tetrachoric correlation coefficient (Pearson, 1900; Kirk, 1973).

         For any of these outcome measures, one needs to specify the cell frequencies via the \code{ai}, \code{bi}, \code{ci}, and \code{di} arguments (or alternatively, one can use the \code{ai}, \code{ci}, \code{n1i}, and \code{n2i} arguments). The options for the \code{measure} argument are then:
         \itemize{
         \item \code{"OR"} for the \emph{log odds ratio},
         \item \code{"PHI"} for the \emph{phi coefficient},
         \item \code{"YUQ"} for \emph{Yule's Q} (Yule, 1912),
         \item \code{"YUY"} for \emph{Yule's Y} (Yule, 1912),
         \item \code{"RTET"} for the \emph{tetrachoric correlation coefficient}.
         }

         There are also measures \code{"ZPHI"} and \code{"ZTET"} for applying Fisher's r-to-z transformation to these measures. This may be useful when combining these with other types of correlation coefficients that were r-to-z transformed. However, note that the r-to-z transformation is \emph{not} a variance-stabilizing transformation for these measures.

         Tables with one or more zero counts are handled as described earlier. For \code{measure="PHI"}, one must indicate via \code{vtype="ST"} or \code{vtype="CS"} whether the data for the studies were obtained using stratified or cross-sectional (i.e., multinomial) sampling, respectively (it is also possible to specify an entire vector for the \code{vtype} argument in case the sampling scheme differed for the various studies).

         A dataset corresponding to data of this type is provided in \code{\link[metadat]{dat.bourassa1996}}.

      }

      \subsection{(2d) Measures for Mixed Variable Types}{

         Finally, we can consider outcome measures that can be used to describe the relationship between two variables, where one variable is dichotomous and the other variable measures some quantitative characteristic. In that case, it is likely that study authors again report summary statistics, such as the mean and standard deviation of the measurements within the two groups (defined by the dichotomous variable). Based on this information, one can compute the point-biserial correlation coefficient (Tate, 1954) as a measure of association between the two variables. If the dichotomous variable is actually a dichotomized version of the responses on an underlying quantitative scale (and assuming that the two variables follow a bivariate normal distribution), it is also possible to estimate the correlation between the two variables using the biserial correlation coefficient (Pearson, 1909; Soper, 1914; Jacobs & Viechtbauer, 2017).

         Here, one again needs to specify \code{m1i} and \code{m2i} for the observed means of the two groups, \code{sd1i} and \code{sd2i} for the observed standard deviations, and \code{n1i} and \code{n2i} for the number of individuals in each group. The options for the \code{measure} argument are then:
         \itemize{
         \item \code{"RPB"} for the \emph{point-biserial correlation coefficient},
         \item \code{"RBIS"} for the \emph{biserial correlation coefficient}.
         }

         There are also measures \code{"ZPB"} and \code{"ZBIS"} for applying Fisher's r-to-z transformation to these measures. This may be useful when combining these with other types of correlation coefficients that were r-to-z transformed. However, note that the r-to-z transformation is \emph{not} a variance-stabilizing transformation for these measures.

         If the means and standard deviations are unknown for some studies, one can also use argument \code{di}, \code{ti}, or \code{pi} to specify standardized mean differences (Cohen's d values), t-statistics from an independent samples t-test, or (signed) p-values for the t-test, respectively, as described earlier under (1a) (together with the group sizes, these are sufficient statistics for computing the (point-)biserial correlation coefficients).

         For \code{measure="RPB"}, one must indicate via \code{vtype="ST"} or \code{vtype="CS"} whether the data for the studies were obtained using stratified or cross-sectional (i.e., multinomial) sampling, respectively (it is also possible to specify an entire vector for the \code{vtype} argument in case the sampling scheme differed for the various studies).

      }

   }

   \subsection{(3) Outcome Measures for Individual Groups}{

      In this section, outcome measures will be described which may be useful when the goal of a meta-analysis is to synthesize studies that characterize some property of individual groups. We will again distinguish between measures that are applicable when the characteristic assessed is a quantitative variable, a dichotomous variable, or when the characteristic represents an event count.

      \subsection{(3a) Measures for Quantitative Variables}{

         The goal of a meta-analysis may be to characterize individual groups, where the response, characteristic, or dependent variable assessed in the individual studies is measured on some quantitative scale. In the simplest case, the raw mean for the quantitative variable is reported for each group, which then becomes the observed outcome for the meta-analysis. Here, one needs to specify \code{mi}, \code{sdi}, and \code{ni} for the observed means, the observed standard deviations, and the sample sizes, respectively. For ratio scale measurements, the log transformed mean or the log transformed coefficient of variation (with bias correction) may also be of interest (Nakagawa et al., 2015). If focus is solely on the variability of the measurements, then the log transformed standard deviation (with bias correction) is a useful measure (Nakagawa et al., 2015; Raudenbush & Bryk, 1987). Here, one only needs to specify \code{sdi} and \code{ni}.

         The options for the \code{measure} argument are:
         \itemize{
         \item \code{"MN"} for the \emph{raw mean},
         \item \code{"MNLN"} for the \emph{log transformed mean},
         \item \code{"CVLN"} for the \emph{log transformed coefficient of variation},
         \item \code{"SDLN"} for the \emph{log transformed standard deviation}.
         }

         Note that \code{sdi} is used to specify the standard deviations of the observed values of the response, characteristic, or dependent variable and not the standard errors of the means. Also, the sampling variance for \code{measure="CVLN"} is computed as given by equation 27 in Nakagawa et al. (2015), but without the \sQuote{\mjseqn{-2 \rho \ldots}} term, since for normally distributed data (which we assume here) the mean and variance (and transformations thereof) are independent.

      }

      \subsection{(3b) Measures for Dichotomous Variables}{

         A meta-analysis may also be conducted to aggregate studies that provide data about individual groups with respect to a dichotomous dependent variable. Here, one needs to specify \code{xi} and \code{ni}, denoting the number of individuals experiencing the event of interest and the total number of individuals within each study, respectively. Instead of specifying \code{ni}, one can use \code{mi} to specify the number of individuals that do not experience the event of interest. The options for the \code{measure} argument are then:
         \itemize{
         \item \code{"PR"} for the \emph{raw proportion},
         \item \code{"PLN"} for the \emph{log transformed proportion},
         \item \code{"PLO"} for the \emph{logit transformed proportion} (i.e., log odds),
         \item \code{"PAS"} for the \emph{arcsine square root transformed proportion} (i.e., the angular transformation),
         \item \code{"PFT"} for the \emph{Freeman-Tukey double arcsine transformed proportion} (Freeman & Tukey, 1950).
         }

         Zero cell entries can be problematic for certain outcome measures. When \code{to="only0"} (the default), the value of \code{add} (the default is \code{1/2}; but see \sQuote{Note}) is added to \code{xi} and \code{mi} only for studies where \code{xi} or \code{mi} is equal to 0. When \code{to="all"}, the value of \code{add} is added to \code{xi} and \code{mi} in all studies. When \code{to="if0all"}, the value of \code{add} is added in all studies, but only when there is at least one study with a zero value for \code{xi} or \code{mi}. Setting \code{to="none"} or \code{add=0} has the same effect: No adjustment to the observed values is made. Depending on the outcome measure and the data, this may lead to division by zero (when this occurs, the resulting value is recoded to \code{NA}).

         Datasets corresponding to data of this type are provided in \code{\link[metadat]{dat.pritz1997}} and \code{\link[metadat]{dat.debruin2009}}.

      }

      \subsection{(3c) Measures for Event Counts}{

         Various measures can be used to characterize individual groups when the dependent variable assessed is an event count. Here, one needs to specify \code{xi} and \code{ti}, denoting the number of events that occurred and the total person-times at risk, respectively. The options for the \code{measure} argument are then:
         \itemize{
         \item \code{"IR"} for the \emph{raw incidence rate},
         \item \code{"IRLN"} for the \emph{log transformed incidence rate},
         \item \code{"IRS"} for the \emph{square root transformed incidence rate},
         \item \code{"IRFT"} for the \emph{Freeman-Tukey transformed incidence rate} (Freeman & Tukey, 1950).
         }

         Measures \code{"IR"} and \code{"IRLN"} can also be used when meta-analyzing standardized incidence ratios (SIRs), where the observed number of events is divided by the expected number of events. In this case, arguments \code{xi} and \code{ti} are used to specify the observed and expected number of events in the studies. Since SIRs are not symmetric around 1, it is usually more appropriate to meta-analyze the log transformed SIRs (i.e., using measure \code{"IRLN"}), which are symmetric around 0.

         Studies with zero events can be problematic, especially for the log transformed incidence rate. Adding a small constant to the number of events is a common solution to this problem. When \code{to="only0"} (the default), the value of \code{add} (the default is \code{1/2}; but see \sQuote{Note}) is added to \code{xi} only in the studies that have zero events. When \code{to="all"}, the value of \code{add} is added to \code{xi} in all studies. When \code{to="if0all"}, the value of \code{add} is added to \code{xi} in all studies, but only when there is at least one study with zero events. Setting \code{to="none"} or \code{add=0} has the same effect: No adjustment to the observed number of events is made. Depending on the outcome measure and the data, this may lead to division by zero (when this occurs, the resulting value is recoded to \code{NA}).

      }

   }

   \subsection{(4) Outcome Measures for Change or Matched Pairs}{

      The purpose of a meta-analysis may be to assess the amount of change within individual groups (e.g., before and after a treatment or under two different treatments) or when dealing with matched pairs designs.

      \subsection{(4a) Measures for Quantitative Variables}{

         When the response or dependent variable assessed in the individual studies is measured on some quantitative scale, the raw mean change, standardized versions thereof, or the (log transformed) ratio of means (log response ratio) can be used as outcome measures (Becker, 1988; Gibbons et al., 1993; Lajeunesse, 2011; Morris, 2000). Here, one needs to specify \code{m1i} and \code{m2i}, the observed means at the two measurement occasions, \code{sd1i} and \code{sd2i} for the corresponding observed standard deviations, \code{ri} for the correlation between the measurements at the two measurement occasions, and \code{ni} for the sample size. The options for the \code{measure} argument are then:
         \itemize{
         \item \code{"MC"} for the \emph{raw mean change},
         \item \code{"SMCC"} for the \emph{standardized mean change} using change score standardization (Gibbons et al., 1993),
         \item \code{"SMCR"} for the \emph{standardized mean change} using raw score standardization (Becker, 1988),
         \item \code{"SMCRH"} for the \emph{standardized mean change} using raw score standardization with heteroscedastic population variances at the two measurement occasions (Bonett, 2008),
         \item \code{"ROMC"} for the \emph{log transformed ratio of means} (Lajeunesse, 2011).
         }

         Note that \code{"MC"} is simply \code{m1i-m2i} and \code{"SMCC"}, \code{"SMCR"}, and \code{"SMCRH"} are standardized versions thereof (e.g., \code{"SMCR"} is computed as \code{(m1i-m2i)/sd1i}). See also Morris and DeShon (2002) for a thorough discussion of the difference between the \code{"SMCC"} and \code{"SMCR"} change score measures.

         A few notes about the change score measures. In practice, one often has a mix of information available from the individual studies to compute these measures. In particular, if \code{m1i} and \code{m2i} are unknown, but the raw mean change is directly reported in a particular study, then one can set \code{m1i} to that value and \code{m2i} to 0 (making sure that the raw mean change was computed as \code{m1i-m2i} within that study and not the other way around). Also, for the raw mean change (\code{"MC"}) or the standardized mean change using change score standardization (\code{"SMCC"}), if \code{sd1i}, \code{sd2i}, and \code{ri} are unknown, but the standard deviation of the change scores is directly reported, then one can set \code{sd1i} to that value and both \code{sd2i} and \code{ri} to 0. Finally, for the standardized mean change using raw score standardization (\code{"SMCR"}), argument \code{sd2i} is actually not needed, as the standardization is only based on \code{sd1i} (Becker, 1988; Morris, 2000), which is usually the pre-test standard deviation (if the post-test standard deviation should be used, then set \code{sd1i} to that). Note that all of these measures are also applicable for matched pairs designs (subscripts 1 and 2 then simply denote the first and second group that are formed by the matching).

         Finally, interest may also be focused on differences in the variability of the measurements at the two measurement occasions (or between the two matched groups). Here, the (log transformed) ratio of the coefficient of variation (also called the coefficient of variation ratio) can be a useful measure (Nakagawa et al., 2015). If focus is solely on the variability of the measurements, then the (log transformed) ratio of the standard deviations (also called the variability ratio) can be used (Nakagawa et al., 2015). For the latter, one only needs to specify \code{sd1i}, \code{sd2i}, \code{ni}, and \code{ri}. The options for the \code{measure} argument are:
         \itemize{
         \item \code{"CVRC"} for the \emph{log transformed coefficient of variation ratio},
         \item \code{"VRC"} for the \emph{log transformed variability ratio}.
         }

         The definitions of these measures are the same as given in Nakagawa et al. (2015) but are here computed for two sets of dependent measurements. Hence, the computation of the sampling variances are adjusted to take the correlation between the measurements into consideration.

      }

      \subsection{(4b) Measures for Dichotomous Variables}{

         The data for a study examining change in a dichotomous variable gives rise to a paired \mjeqn{2 \times 2}{2x2} table, which is of the form:
         \tabular{lcccc}{
                         \ics \tab \tab trt 2 outcome 1 \tab \ics \tab trt 2 outcome 2 \cr
         trt 1 outcome 1 \ics \tab \tab \code{ai}       \tab \ics \tab \code{bi} \cr
         trt 1 outcome 2 \ics \tab \tab \code{ci}       \tab \ics \tab \code{di}}
         where \code{ai}, \code{bi}, \code{ci}, and \code{di} denote the cell frequencies. Note that \sQuote{trt1} and \sQuote{trt2} may be applied to a single group of subjects or to matched pairs of subjects. Also, \sQuote{trt1} and \sQuote{trt2} might refer to two different time points (e.g., before and after a treatment). In any case, the data from such a study can be rearranged into a marginal table of the form:
         \tabular{lcccc}{
               \tab \ics \tab outcome 1    \tab \ics \tab outcome 2 \cr
         trt 1 \tab \ics \tab \code{ai+bi} \tab \ics \tab \code{ci+di} \cr
         trt 2 \tab \ics \tab \code{ai+ci} \tab \ics \tab \code{bi+di}}
         which is of the same form as a \mjeqn{2 \times 2}{2x2} table that would arise in a study comparing/contrasting two independent groups.

         The options for the \code{measure} argument that will compute outcome measures based on the marginal table are:
         \itemize{
         \item \code{"MPRR"} for the matched pairs \emph{marginal log risk ratio},
         \item \code{"MPOR"} for the matched pairs \emph{marginal log odds ratio},
         \item \code{"MPRD"} for the matched pairs \emph{marginal risk difference}.
         }

         See Becker and Balagtas (1993), Curtin et al. (2002), Elbourne et al. (2002), Fagerland et al. (2014), May and Johnson (1997), Newcombe (1998), Stedman et al. (2011), and Zou (2007) for discussions of these measures.

         The options for the \code{measure} argument that will compute outcome measures based on the paired table are:
         \itemize{
         \item \code{"MPORC"} for the \emph{conditional log odds ratio},
         \item \code{"MPPETO"} for the \emph{conditional log odds ratio} estimated with Peto's method.
         }

         See Curtin et al. (2002) and Zou (2007) for discussions of these measures.

         If only marginal tables are available, then another possibility is to compute the marginal log odds ratios based on these table directly. However, for the correct computation of the sampling variances, the correlations (phi coefficients) from the paired tables must be known (or \sQuote{guestimated}). To use this approach, set \code{measure="MPORM"} and use argument \code{ri} to specify the correlation coefficients. Instead of specifying \code{ri}, one can use argument \code{pi} to specify the proportions (or \sQuote{guestimates} thereof) of individuals (or pairs) that experienced the outcome of interest (i.e., \sQuote{outcome1} in the paired \mjeqn{2 \times 2}{2x2} table) under both treatments (i.e., \code{ai / (ai + bi + ci + di)}). Based on these proportions, the correlation coefficients are then back-calculated and used to calculate the correct sampling variances. Note that the values in the marginal tables put constraints on the possible values for \code{ri} and \code{pi}. If a specified value for \code{ri} or \code{pi} is not feasible under a given table, the corresponding sampling variance will be \code{NA}.

      }

   }

   \subsection{(5) Other Outcome Measures for Meta-Analyses}{

      Other outcome measures are sometimes used for meta-analyses that do not directly fall into the categories above. These are described in this section.

      \subsection{Cronbach's alpha and Transformations Thereof}{

         Meta-analytic methods can also be used to aggregate Cronbach's alpha values from multiple studies. This is usually referred to as a \sQuote{reliability generalization meta-analysis} (Vacha-Haase, 1998). Here, one needs to specify \code{ai}, \code{mi}, and \code{ni} for the observed alpha values, the number of items/replications/parts of the measurement instrument, and the sample sizes, respectively. One can either directly analyze the raw Cronbach's alpha values or transformations thereof (Bonett, 2002, 2010; Hakstian & Whalen, 1976). The options for the \code{measure} argument are then:
         \itemize{
         \item \code{"ARAW"} for \emph{raw alpha} values,
         \item \code{"AHW"} for \emph{transformed alpha values} (Hakstian & Whalen, 1976),
         \item \code{"ABT"} for \emph{transformed alpha values} (Bonett, 2002).
         }

         Note that the transformations implemented here are slightly different from the ones described by Hakstian and Whalen (1976) and Bonett (2002). In particular, for \code{"AHW"}, the transformation \mjeqn{1-(1-\alpha)^{1/3}}{1-(1-\alpha)^(1/3)} is used, while for \code{"ABT"}, the transformation \mjeqn{-\ln(1-\alpha)}{-ln(1-\alpha)} is used. This ensures that the transformed values are monotonically increasing functions of \mjseqn{\alpha}.

         A dataset corresponding to data of this type is provided in \code{\link[metadat]{dat.bonett2010}}.

      }

      \subsection{Partial and Semi-Partial Correlations}{

         Aloe and Becker (2012), Aloe and Thompson (2013), and Aloe (2014) describe the use of partial and semi-partial correlation coefficients for meta-analyzing the results from regression models (when the focus is on a common regression coefficient of interest across studies). To compute these measures, one needs to specify \code{ti} for the test statistics (i.e., t-tests) of the regression coefficient of interest, \code{ni} for the sample sizes of the studies, \code{mi} for the number of predictors in the regression models (not counting the intercept term), and \code{r2i} for the \mjseqn{R^2} value of the regression models (the latter is only needed when \code{measure="SPCOR"} or \code{measure="ZSPCOR"}). The options for the \code{measure} argument are then:
         \itemize{
         \item \code{"PCOR"} for the \emph{partial correlation coefficient},
         \item \code{"ZPCOR"} for \emph{Fisher's r-to-z transformed partial correlation coefficient},
         \item \code{"SPCOR"} for the \emph{semi-partial correlation coefficient},
         \item \code{"ZSPCOR"} for \emph{Fisher's r-to-z transformed semi-partial correlation coefficient}.
         }

         Note that the signs of the (semi-)partial correlation coefficients is determined based on the signs of the values specified via the \code{ti} argument. Also, while the Fisher transformation can be applied to both measures, it is only a variance-stabilizing transformation for partial correlation coefficients.

         If the test statistic (i.e., t-test) of the regression coefficient of interest is unknown for some studies, but the (two-sided) p-values corresponding to the t-tests are known, one can specify those values via argument \code{pi}. However, since a two-sided p-value does not carry information about the sign of the test statistic (and hence neither about the correlation), the sign of the p-values (which can be negative) is used as the sign of the correlation coefficients (e.g., \code{escalc(measure="PCOR", pi=-.07, mi=5, ni=30)} yields a negative correlation of \code{-0.3610}).

         In the rare case that the (semi-)partial correlations are known for some of the studies, then these can be directly specified via the \code{ri} argument. This can be useful, for example, when \mjseqn{\eta^2_p} (i.e., partial eta squared) is known for the regression coefficient of interest, since the square root thereof is identical to the absolute value of the partial correlation (although the correct sign then still needs to be reconstructed based on other information).

         A dataset corresponding to data of this type is provided in \code{\link[metadat]{dat.aloe2013}}.

      }

      \subsection{Relative Excess Heterozygosity}{

         Ziegler et al. (2011) describe the use of meta-analytic methods to examine deviations from the Hardy-Weinberg equilibrium across multiple studies. The relative excess heterozygosity (REH) is the proposed measure for such a meta-analysis, which can be computed by setting \code{measure="REH"}. Here, one needs to specify \code{ai} for the number of individuals with homozygous dominant alleles, \code{bi} for the number of individuals with heterozygous alleles, and \code{ci} for the number of individuals with homozygous recessives alleles.

         Note that the log is taken of the REH values, which makes this outcome measure symmetric around 0 and results in a sampling distribution that is closer to normality.

         A dataset corresponding to data of this type is provided in \code{\link[metadat]{dat.frank2008}}.

      }

   }

   \subsection{(6) Converting a Data Frame to an 'escalc' Object}{

      The function can also be used to convert a regular data frame to an \sQuote{escalc} object. One simply sets the \code{measure} argument to one of the options described above (or to \code{measure="GEN"} for a generic outcome measure not further specified) and passes the observed effect sizes or outcomes via the \code{yi} argument and the corresponding sampling variances via the \code{vi} argument (or the standard errors via the \code{sei} argument) to the function.

   }

}
\value{
   An object of class \code{c("escalc","data.frame")}. The object is a data frame containing the following components:
   \item{yi}{vector with the observed effect sizes or outcomes.}
   \item{vi}{vector with the corresponding sampling variances.}

   If a data frame was specified via the \code{data} argument and \code{append=TRUE}, then variables \code{yi} and \code{vi} are appended to this data frame. Note that the \code{var.names} argument actually specifies the names of these two variables (\code{"yi"} and \code{"vi"} are the defaults).

   If the data frame already contains two variables with names as specified by the \code{var.names} argument, the values for these two variables will be overwritten when \code{replace=TRUE} (which is the default). By setting \code{replace=FALSE}, only values that are \code{NA} will be replaced.

   The \code{subset} argument can be used to select the studies that will be included in the data frame returned by the function. On the other hand, the \code{include} argument simply selects for which studies the measure will be computed (if it shouldn't be computed for all of them).

   The object is formatted and printed with the \code{\link[=print.escalc]{print}} function. The \code{\link[=summary.escalc]{summary}} function can be used to obtain confidence intervals for the individual outcomes. See \code{\link{methods.escalc}} for some additional method functions for \code{"escalc"} objects.

   With the \code{\link[=aggregate.escalc]{aggregate}} function, one can aggregate multiple effect sizes or outcomes belonging to the same study (or some other clustering variable) into a single combined effect size or outcome.
}
\note{
   The variable names specified under \code{var.names} should be syntactically valid variable names. If necessary, they are adjusted so that they are.

   Although the default value for \code{add} is \code{1/2}, for certain measures the use of such a bias correction makes little sense and for these measures, the function internally sets \code{add=0}. This applies to the following measures: \code{"AS"}, \code{"PHI"}, \code{"ZPHI"}, \code{"RTET"}, \code{"ZTET"}, \code{"IRSD"}, \code{"PAS"}, \code{"PFT"}, \code{"IRS"}, and \code{"IRFT"}. One can still force the use of the bias correction by explicitly setting the \code{add} argument to some non-zero value when calling the function.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Aloe, A. M. (2014). An empirical investigation of partial effect sizes in meta-analysis of correlational data. \emph{Journal of General Psychology}, \bold{141}(1), 47--64. \verb{https://doi.org/10.1080/00221309.2013.853021}

   Aloe, A. M., & Becker, B. J. (2012). An effect size for regression predictors in meta-analysis. \emph{Journal of Educational and Behavioral Statistics}, \bold{37}(2), 278--297. \verb{https://doi.org/10.3102/1076998610396901}

   Aloe, A. M., & Thompson, C. G. (2013). The synthesis of partial effect sizes. \emph{Journal of the Society for Social Work and Research}, \bold{4}(4), 390--405. \verb{https://doi.org/10.5243/jsswr.2013.24}

   Bagos, P. G., & Nikolopoulos, G. K. (2009). Mixed-effects Poisson regression models for meta-analysis of follow-up studies with constant or varying durations. \emph{The International Journal of Biostatistics}, \bold{5}(1). \verb{https://doi.org/10.2202/1557-4679.1168}

   Becker, B. J. (1988). Synthesizing standardized mean-change measures. \emph{British Journal of Mathematical and Statistical Psychology}, \bold{41}(2), 257--278.  \verb{https://doi.org/10.1111/j.2044-8317.1988.tb00901.x}

   Becker, M. P., & Balagtas, C. C. (1993). Marginal modeling of binary cross-over data. \emph{Biometrics}, \bold{49}(4), 997--1009. \verb{https://doi.org/10.2307/2532242}

   Bonett, D. G. (2002). Sample size requirements for testing and estimating coefficient alpha. \emph{Journal of Educational and Behavioral Statistics}, \bold{27}(4), 335--340. \verb{https://doi.org/10.3102/10769986027004335}

   Bonett, D. G. (2008). Confidence intervals for standardized linear contrasts of means. \emph{Psychological Methods}, \bold{13}(2), 99--109. \verb{https://doi.org/10.1037/1082-989X.13.2.99}

   Bonett, D. G. (2009). Meta-analytic interval estimation for standardized and unstandardized mean differences. \emph{Psychological Methods}, \bold{14}(3), 225--238. \verb{https://doi.org/10.1037/a0016619}

   Bonett, D. G. (2010). Varying coefficient meta-analytic methods for alpha reliability. \emph{Psychological Methods}, \bold{15}(4), 368--385. \verb{https://doi.org/10.1037/a0020142}

   Borenstein, M. (2009). Effect sizes for continuous data. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta-analysis} (2nd ed., pp. 221--235). New York: Russell Sage Foundation.

   Chinn, S. (2000). A simple method for converting an odds ratio to effect size for use in meta-analysis. \emph{Statistics in Medicine}, \bold{19}(22), 3127--3131. \verb{https://doi.org/10.1002/1097-0258(20001130)19:22<3127::aid-sim784>3.0.co;2-m}

   Cohen, J. (1988). \emph{Statistical power analysis for the behavioral sciences} (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum Associates.

   Cox, D. R., & Snell, E. J. (1989). \emph{Analysis of binary data} (2nd ed.). London: Chapman & Hall.

   Curtin, F., Elbourne, D., & Altman, D. G. (2002). Meta-analysis combining parallel and cross-over clinical trials. II: Binary outcomes. \emph{Statistics in Medicine}, \bold{21}(15), 2145--2159. \verb{https://doi.org/10.1002/sim.1206}

   Elbourne, D. R., Altman, D. G., Higgins, J. P. T., Curtin, F., Worthington, H. V., & Vail, A. (2002). Meta-analyses involving cross-over trials: Methodological issues. \emph{International Journal of Epidemiology}, \bold{31}(1), 140--149. \verb{https://doi.org/10.1093/ije/31.1.140}

   Fagerland, M. W., Lydersen, S., & Laake, P. (2014). Recommended tests and confidence intervals for paired binomial proportions. \emph{Statistics in Medicine}, \bold{33}(16), 2850--2875. \verb{https://doi.org/10.1002/sim.6148}

   Fisher, R. A. (1921). On the \dQuote{probable error} of a coefficient of correlation deduced from a small sample. \emph{Metron}, \bold{1}, 1--32. \verb{http://hdl.handle.net/2440/15169}

   Fleiss, J. L., & Berlin, J. (2009). Effect sizes for dichotomous data. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta-analysis} (2nd ed., pp. 237--253). New York: Russell Sage Foundation.

   Freeman, M. F., & Tukey, J. W. (1950). Transformations related to the angular and the square root. \emph{Annals of Mathematical Statistics}, \bold{21}(4), 607--611. \verb{https://doi.org/10.1214/aoms/1177729756}

   Gibbons, R. D., Hedeker, D. R., & Davis, J. M. (1993). Estimation of effect size from a series of experiments involving paired comparisons. \emph{Journal of Educational Statistics}, \bold{18}(3), 271--279. \verb{https://doi.org/10.3102/10769986018003271}

   Hakstian, A. R., & Whalen, T. E. (1976). A k-sample significance test for independent alpha coefficients. \emph{Psychometrika}, \bold{41}(2), 219--231. \verb{https://doi.org/10.1007/BF02291840}

   Hasselblad, V., & Hedges, L. V. (1995). Meta-analysis of screening and diagnostic tests. Psychological Bulletin, 117(1), 167-178. \verb{https://doi.org/10.1037/0033-2909.117.1.167}

   Hedges, L. V. (1981). Distribution theory for Glass's estimator of effect size and related estimators. \emph{Journal of Educational Statistics}, \bold{6}(2), 107--128. \verb{https://doi.org/10.3102/10769986006002107}

   Hedges, L. V. (1982). Estimation of effect size from a series of independent experiments. \emph{Psychological Bulletin}, \bold{92}(2), 490--499. \verb{https://doi.org/10.1037/0033-2909.92.2.490}

   Hedges, L. V. (1983). A random effects model for effect sizes. \emph{Psychological Bulletin}, \bold{93}(2), 388--395. \verb{https://doi.org/10.1037/0033-2909.93.2.388}

   Hedges, L. V. (1989). An unbiased correction for sampling error in validity generalization studies. \emph{Journal of Applied Psychology}, \bold{74}(3), 469--477. \verb{https://doi.org/10.1037/0021-9010.74.3.469}

   Hedges, L. V., Gurevitch, J., & Curtis, P. S. (1999). The meta-analysis of response ratios in experimental ecology. \emph{Ecology}, \bold{80}(4), 1150--1156. \verb{https://doi.org/10.1890/0012-9658(1999)080[1150:TMAORR]2.0.CO;2}

   Higgins, J. P. T., Thomas, J., Chandler, J., Cumpston, M., Li, T., Page, M. J., & Welch, V. A. (Eds.) (2019). \emph{Cochrane handbook for systematic reviews of interventions} (2nd ed.). Chichester, UK: Wiley. \verb{https://training.cochrane.org/handbook}

   Jacobs, P., & Viechtbauer, W. (2017). Estimation of the biserial correlation and its sampling variance for use in meta-analysis. \emph{Research Synthesis Methods}, \bold{8}(2), 161--180. \verb{https://doi.org/10.1002/jrsm.1218}

   Kirk, D. B. (1973). On the numerical approximation of the bivariate normal (tetrachoric) correlation coefficient. \emph{Psychometrika}, \bold{38}(2), 259--268. \verb{https://doi.org/10.1007/BF02291118}

   Lajeunesse, M. J. (2011). On the meta-analysis of response ratios for studies with correlated and multi-group designs. \emph{Ecology}, \bold{92}(11), 2049--2055. \verb{https://doi.org/10.1890/11-0423.1}

   May, W. L., & Johnson, W. D. (1997). Confidence intervals for differences in correlated binary proportions. \emph{Statistics in Medicine}, \bold{16}(18), 2127--2136. \verb{https://doi.org/10.1002/(SICI)1097-0258(19970930)16:18<2127::AID-SIM633>3.0.CO;2-W}

   Morris, S. B. (2000). Distribution of the standardized mean change effect size for meta-analysis on repeated measures. \emph{British Journal of Mathematical and Statistical Psychology}, \bold{53}(1), 17--29. \verb{https://doi.org/10.1348/000711000159150}

   Morris, S. B., & DeShon, R. P. (2002). Combining effect size estimates in meta-analysis with repeated measures and independent-groups designs. \emph{Psychological Methods}, \bold{7}(1), 105--125. \verb{https://doi.org/10.1037/1082-989x.7.1.105}

   Nakagawa, S., Poulin, R., Mengersen, K., Reinhold, K., Engqvist, L., Lagisz, M., & Senior, A. M. (2015). Meta-analysis of variation: Ecological and evolutionary applications and beyond. \emph{Methods in Ecology and Evolution}, \bold{6}(2), 143--152. \verb{https://doi.org/10.1111/2041-210x.12309}

   Newcombe, R. G. (1998). Improved confidence intervals for the difference between binomial proportions based on paired data. \emph{Statistics in Medicine}, \bold{17}(22), 2635--2650. \verb{https://doi.org/10.1002/(SICI)1097-0258(19981130)17:22<2635::AID-SIM954>3.0.CO;2-C}

   Olkin, I., & Pratt, J. W. (1958). Unbiased estimation of certain correlation coefficients. \emph{Annals of Mathematical Statistics}, \bold{29}(1), 201--211. \verb{https://doi.org/10.1214/aoms/1177706717}

   Pearson, K. (1900). Mathematical contributions to the theory of evolution. VII. On the correlation of characters not quantitatively measurable. \emph{Philosophical Transactions of the Royal Society of London, Series A}, \bold{195}, 1--47. \verb{https://doi.org/10.1098/rsta.1900.0022}

   Pearson, K. (1909). On a new method of determining correlation between a measured character A, and a character B, of which only the percentage of cases wherein B exceeds (or falls short of) a given intensity is recorded for each grade of A. \emph{Biometrika}, \bold{7}(1/2), 96--105. \verb{https://doi.org/10.1093/biomet/7.1-2.96}

   Raudenbush, S. W., & Bryk, A. S. (1987). Examining correlates of diversity. \emph{Journal of Educational Statistics}, \bold{12}(3), 241--269. \verb{https://doi.org/10.3102/10769986012003241}

   Rothman, K. J., Greenland, S., & Lash, T. L. (2008). \emph{Modern epidemiology} (3rd ed.). Philadelphia: Lippincott Williams & Wilkins.

   \enc{Rücker}{Ruecker}, G., Schwarzer, G., Carpenter, J., & Olkin, I. (2009). Why add anything to nothing? The arcsine difference as a measure of treatment effect in meta-analysis with zero cells. \emph{Statistics in Medicine}, \bold{28}(5), 721--738. \verb{https://doi.org/10.1002/sim.3511}

   \enc{Sánchez-Meca}{Sanchez-Meca}, J., \enc{Marín-Martínez}{Marin-Martinez}, F., & \enc{Chacón-Moscoso}{Chacon-Moscoso}, S. (2003). Effect-size indices for dichotomized outcomes in meta-analysis. \emph{Psychological Methods}, \bold{8}(4), 448--467. \verb{https://doi.org/10.1037/1082-989X.8.4.448}

   Soper, H. E. (1914). On the probable error of the bi-serial expression for the correlation coefficient. \emph{Biometrika}, \bold{10}(2/3), 384--390. \verb{https://doi.org/10.1093/biomet/10.2-3.384}

   Stedman, M. R., Curtin, F., Elbourne, D. R., Kesselheim, A. S., & Brookhart, M. A. (2011). Meta-analyses involving cross-over trials: Methodological issues. \emph{International Journal of Epidemiology}, \bold{40}(6), 1732--1734. \verb{https://doi.org/10.1093/ije/dyp345}

   Tate, R. F. (1954). Correlation between a discrete and a continuous variable: Point-biserial correlation. \emph{Annals of Mathematical Statistics}, \bold{25}(3), 603--607. \verb{https://doi.org/10.1214/aoms/1177728730}

   Vacha-Haase, T. (1998). Reliability generalization: Exploring variance in measurement error affecting score reliability across studies. \emph{Educational and Psychological Measurement}, \bold{58}(1), 6--20. \verb{https://doi.org/10.1177/0013164498058001002}

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}

   Yule, G. U. (1912). On the methods of measuring association between two attributes. \emph{Journal of the Royal Statistical Society}, \bold{75}(6), 579--652. \verb{https://doi.org/10.2307/2340126}

   Yusuf, S., Peto, R., Lewis, J., Collins, R., & Sleight, P. (1985). Beta blockade during and after myocardial infarction: An overview of the randomized trials. \emph{Progress in Cardiovascular Disease}, \bold{27}(5), 335--371. \verb{https://doi.org/10.1016/s0033-0620(85)80003-7}

   Ziegler, A., Steen, K. V. & Wellek, S. (2011). Investigating Hardy-Weinberg equilibrium in case-control or cohort studies or meta-analysis. \emph{Breast Cancer Research and Treatment}, \bold{128}(1), 197--201. \verb{https://doi.org/10.1007/s10549-010-1295-z}

   Zou, G. Y. (2007). One relative risk versus two odds ratios: Implications for meta-analyses involving paired and unpaired binary data. \emph{Clinical Trials}, \bold{4}(1), 25--31. \verb{https://doi.org/10.1177/1740774506075667}
}
\seealso{
   \code{\link{print.escalc}} and \code{\link{summary.escalc}} for the print and summary methods.

   \code{\link{conv.2x2}} for a function to reconstruct the cell frequencies of \mjeqn{2 \times 2}{2x2} tables based on other summary statistics.

   \code{\link{conv.fivenum}} for a function to convert five-number summary values to means and standard deviations (needed to compute various effect size measures, such as raw or standardized mean differences and ratios of means / response ratios).

   \code{\link{conv.wald}} for a function to convert Wald-type confidence intervals and test statistics to sampling variances.

   \code{\link{conv.delta}} for a function to transform observed effect sizes or outcomes and their sampling variances using the delta method.

   \code{\link{vcalc}} for a function to construct or approximate the variance-covariance matrix of dependent effect sizes or outcomes.

   \code{\link{rcalc}} for a function to construct the variance-covariance matrix of dependent correlation coefficients.

   \code{\link{rma.uni}} and \code{\link{rma.mv}} for model fitting functions that can take the calculated effect sizes or outcomes (and the corresponding sampling variances) as input.

   \code{\link{rma.mh}}, \code{\link{rma.peto}}, and \code{\link{rma.glmm}} for model fitting functions that take similar inputs.
}
\examples{
############################################################################

### data from the meta-analysis by Coliditz et al. (1994) on the efficacy of
### BCG vaccine in the prevention of tuberculosis dat.bcg
dat.bcg

### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
dat

### suppose that for a particular study, yi and vi are known (i.e., have
### already been calculated) but the 2x2 table counts are not known; with
### replace=FALSE, the yi and vi values for that study are not replaced
dat[1:12,10:11] <- NA
dat[13,4:7] <- NA
dat
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat, replace=FALSE)
dat

### illustrate difference between 'subset' and 'include' arguments
escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, subset=1:6)
escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, include=1:6)

### illustrate the 'var.names' argument
escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, var.names=c("lnrr","var"))

### illustrate the 'slab' argument
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg,
              data=dat.bcg, slab=paste0(author, ", ", year))
dat

### note: the output looks the same but the study labels are stored as an attribute with the
### effect size estimates (together with the total sample size of the studies and the chosen
### effect size measure)
dat$yi

### this information can then be used by other functions; for example in a forest plot
forest(dat$yi, dat$vi, header=TRUE)

############################################################################

### convert a regular data frame to an 'escalc' object
### dataset from Lipsey & Wilson (2001), Table 7.1, page 130
dat <- data.frame(id = c(100, 308, 1596, 2479, 9021, 9028, 161, 172, 537, 7049),
                  yi = c(-0.33, 0.32, 0.39, 0.31, 0.17, 0.64, -0.33, 0.15, -0.02, 0.00),
                  vi = c(0.084, 0.035, 0.017, 0.034, 0.072, 0.117, 0.102, 0.093, 0.012, 0.067),
                  random = c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1),
                  intensity = c(7, 3, 7, 5, 7, 7, 4, 4, 5, 6))
dat
dat <- escalc(measure="SMD", yi=yi, vi=vi, data=dat, slab=paste("Study ID:", id), digits=3)
dat

############################################################################
}
\keyword{datagen}
