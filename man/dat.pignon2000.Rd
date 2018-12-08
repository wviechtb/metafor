\name{dat.pignon2000}
\docType{data}
\alias{dat.pignon2000}
\title{Studies on the Effectiveness of Locoregional Treatment plus Chemotherapy for Head and Neck Squamous-Cell Carcinoma}
\description{Results from studies examining mortality risk in patients with nonmetastatic head and neck squamous-cell carcinoma receiving either locoregional treatment plus chemotherapy versus locoregional treatment alone.}
\usage{dat.pignon2000}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{id}         \tab \code{numeric}    \tab study id number \cr
\bold{trial}      \tab \code{character}  \tab trial abbreviation \cr
\bold{OmE}        \tab \code{numeric}    \tab observed minus expected number of deaths in the locoregional treatment plus chemotherapy group \cr
\bold{V}          \tab \code{numeric}    \tab corresponding variance \cr
\bold{grp}        \tab \code{numeric}    \tab timing of chemotherapy: 1 = adjuvant, 2 = neoadjuvant, 3 = concomitant
}
}
\details{
   The purpose of this meta-analysis was to examine the mortality risk in patients with nonmetastatic head and neck squamous-cell carcinoma receiving either locoregional treatment plus chemotherapy versus locoregional treatment alone. For 65 trials, the dataset provides the observed minus expected number of deaths and corresponding variances in the locoregional treatment plus chemotherapy group. Based on these values, we can estimate the log hazard ratios with \code{OmE/V} and the corresponding sampling variance with \code{1/V}.

   The trials were also divided according to the timing of the chomotherapy: (1) adjuvant, after the locoregional treatment, (2) neoadjuvant, before the locoregional treatment, and (3) concomitant, chemotherapy given concomitantly or alternating with radiotherapy.
}
\source{
   Pignon, J. P., Bourhis, J., Domenge, C., & Designe, L. (2000). Chemotherapy added to locoregional treatment for head and neck squamous-cell carcinoma: Three meta-analyses of updated individual data. \emph{Lancet}, \bold{355}(9208), 949--955.
}
\examples{
### copy data into 'dat'
dat <- dat.pignon2000

### compute estimated log hazard ratios and sampling variances
dat$yi <- with(dat, OmE/V)
dat$vi <- with(dat, 1/V)
dat

### meta-analysis based on all 65 trials
res <- rma(yi, vi, data=dat, method="FE", digits=2)
res
predict(res, transf=exp)

### only adjuvant trials
res <- rma(yi, vi, data=dat, method="FE", subset=grp==1, digits=2)
res
predict(res, transf=exp)

### only neoadjuvant trials
res <- rma(yi, vi, data=dat, method="FE", subset=grp==2, digits=2)
res
predict(res, transf=exp)

### only concomitant trials
res <- rma(yi, vi, data=dat, method="FE", subset=grp==3, digits=2)
res
predict(res, transf=exp)
}
\keyword{datasets}
