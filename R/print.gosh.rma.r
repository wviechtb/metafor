print.gosh.rma <- function(x, digits=x$digits, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(x), must="gosh.rma")

   digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)

   .space()

   cat(mstyle$text("Model fits attempted: "))
   cat(mstyle$result(length(x$fit)))
   cat("\n")
   cat(mstyle$text("Model fits succeeded: "))
   cat(mstyle$result(sum(x$fit)))
   cat("\n\n")

   res.table <- matrix(NA_real_, nrow=ncol(x$res), ncol=6)

   res.table[,1] <- apply(x$res, 2, mean, na.rm=TRUE)
   res.table[,2] <- apply(x$res, 2, min, na.rm=TRUE)
   res.table[,3] <- apply(x$res, 2, quantile, .25, na.rm=TRUE)
   res.table[,4] <- apply(x$res, 2, quantile, .50, na.rm=TRUE)
   res.table[,5] <- apply(x$res, 2, quantile, .75, na.rm=TRUE)
   res.table[,6] <- apply(x$res, 2, max, na.rm=TRUE)

   res.table <- fmtx(res.table, digits[["est"]])

   colnames(res.table) <- c("mean", "min", "q1", "median", "q3", "max")
   rownames(res.table) <- colnames(x$res)

   if (ncol(x$res) == 6)
      rownames(res.table)[2] <- "Q"

   ### add blank row before the model coefficients in meta-regression models

   if (ncol(x$res) > 6)
      res.table <- rbind(res.table[seq_len(5),], "", res.table[6:nrow(res.table),,drop=FALSE])

   ### remove row for tau^2 in FE/EE/CE models

   if (is.element(x$method, c("FE","EE","CE")))
      res.table <- res.table[-5,]

   tmp <- capture.output(print(res.table, quote=FALSE, right=TRUE))
   .print.table(tmp, mstyle)

   .space()

   invisible()

}
