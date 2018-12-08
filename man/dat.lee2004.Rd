\name{dat.lee2004}
\docType{data}
\alias{dat.lee2004}
\title{Studies on Acupoint P6 Stimulation for Preventing Nausea}
\description{Results from studies examining the effectiveness of wrist acupuncture point P6 stimulation for preventing postoperative nausea.}
\usage{dat.lee2004}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{id}     \tab \code{numeric}   \tab trial id number \cr
\bold{study}  \tab \code{character} \tab first author \cr
\bold{year}   \tab \code{numeric}   \tab study year \cr
\bold{ai}     \tab \code{numeric}   \tab number of patients experiencing nausea in the treatment group \cr
\bold{n1i}    \tab \code{numeric}   \tab total number of patients in treatment group \cr
\bold{ci}     \tab \code{numeric}   \tab number of patients experiencing nausea in the sham group \cr
\bold{n2i}    \tab \code{numeric}   \tab total number of patients in the sham group
}
}
\details{
   Postoperative nausea and vomiting are common complications following surgery and anaesthesia. As an alternative to drug therapy, acupuncture has been studied as a potential treatment in several trials. The dataset contains the results from 16 clinical trials examining the effectiveness of wrist acupuncture point P6 stimulation for preventing postoperative nausea.
}
\source{
   Lee, A., & Done, M. L. (2004). Stimulation of the wrist acupuncture point P6 for preventing postoperative nausea and vomiting. \emph{Cochrane Database of Systematic Reviews}, \bold{3}, CD003281.
}
\examples{
### copy data into 'dat' and examine data
dat <- dat.lee2004
dat

### meta-analysis based on log risk ratios
res <- rma(measure="RR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat)
res
predict(res, transf=exp, digits=2)
}
\keyword{datasets}
