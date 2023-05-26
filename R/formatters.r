############################################################################

fmtp <- function(p, digits=4, pname="", equal=FALSE, sep=FALSE, add0=FALSE, quote=FALSE) {

   p[p < 0] <- 0
   p[p > 1] <- 1

   digits  <- max(digits, 1)
   cutoff  <- paste(c(".", rep(0,digits-1),1), collapse="")
   ncutoff <- as.numeric(cutoff)

   equal <- ifelse(equal, "=", "")

   if (sep) {
      if (pname != "")
         pname <- paste0(pname, " ")
      sep <- " "
   } else {
      sep <- ""
   }

   out <- ifelse(is.na(p),            paste0(pname, equal, sep, "NA"),
                 ifelse(p >= ncutoff, paste0(pname, equal, sep, formatC(p, digits=digits, format="f")),
                                      paste0(pname, "<",   sep, ifelse(add0, "0", ""), cutoff)))

   if (!quote)
      out <- noquote(out)

   return(out)

}

fmtx <- function(x, digits=4, flag="", quote=FALSE, ...) {

   # in case x is a data frame / matrix with two dimensions

   if (length(dim(x)) == 2L) {

      if (length(digits) == 1L)
         digits <- rep(digits, ncol(x))

      out <- matrix("", nrow=nrow(x), ncol=ncol(x))

      rownames(out) <- rownames(x)
      colnames(out) <- colnames(x)

      for (j in seq_len(ncol(x)))
         out[,j] <- fmtx(x[,j], digits=digits[[j]], flag=flag, ...)

      if (!quote)
         out <- noquote(out, right=TRUE)

      return(out)

   }

   ddd <- list(...)

   if (is.null(ddd$addwidth)) {
      width <- NULL
   } else {
      width <- digits + ddd$addwidth
   }

   if (is.null(ddd$drop0ifint)) {
      drop0ifint <- FALSE
   } else {
      drop0ifint <- ddd$drop0ifint
   }

   if (is.null(ddd$add0)) {
      add0 <- TRUE
   } else {
      add0 <- ddd$add0
   }

   if (!is.null(ddd$thresh)) {
      if (length(x) != 1L)
         stop("Can only use 'thresh' when 'x' is a scalar.")
      if (isTRUE(abs(x) <= ddd$thresh))
         digits <- 0
   }

   if (is.null(ddd$postfix)) {
      postfix <- ""
   } else {
      postfix <- ddd$postfix
   }

   out <- sapply(x, function(x) {
      if (is.na(x))
         return(paste0("NA", postfix))
      out <- formatC(x, format="f", digits=digits, flag=flag, width=width, drop0trailing=drop0ifint && is.integer(digits))
      if (!add0)
         out <- gsub("0\\.", ".", out)
      out <- paste0(out, postfix)
      return(out)
   })

   if (!quote)
      out <- noquote(out, right=TRUE)

   return(out)

}

############################################################################

fmtt <- function(val, tname, df, df1, df2, pval, digits=4, pname="p-val", format=1, sep=TRUE, quote=FALSE, call=FALSE, ...) {

   if (length(val) != 1L)
      stop("Argument 'val' must be a scalar.")

   if (!is.element(format, c(1,2)))
      stop("Argument 'format' can only be equal to 1 or 2.")

   if (missing(pval))
      stop("Must specify the 'pval' argument.")

   sepset <- sep

   if (sep) {
      sep <- " "
   } else {
      sep <- ""
   }

   ddd <- list(...)

   if (is.null(ddd$flag)) {
      flag <- ""
   } else {
      flag <- ddd$flag
   }

   if (length(digits) == 1L)
      digits <- c(test = digits, pval = digits)
   if (length(digits) == 2L)
      names(digits) <- c("test", "pval")
   if (any(!is.element(c("test","pval"), names(digits))))
      stop("Argument 'digits' must have a 'test' and a 'pval' element.")

   if (format == 1) {

      if (missing(df)) {
         if (!missing(df1) && !missing(df2)) {
            out <- bquote(paste(.(tname), "(df1", .(sep), "=", .(sep), .(df1), ",", .(sep), "df2", .(sep), "=", .(sep), .(round(df2,2)), ")", .(sep), "=", .(sep), .(fmtx(val, digits[["test"]], flag=flag)), ", ", .(pname), .(sep), .(fmtp(pval, digits[["pval"]], equal=TRUE, sep=sepset)), sep=""))
            #paste0(tname, "(df1 = ", df1, ", df2 = ", round(df2,2), ") = ", fmtx(val, digits[["test"]]), ", ", pname, fmtp(pval, digits[["pval"]], equal=TRUE, sep=TRUE))
         } else {
            out <- bquote(paste(.(tname), .(sep), "=", .(sep), .(fmtx(val, digits[["test"]], flag=flag)), ", ", .(pname), .(sep), .(fmtp(pval, digits[["pval"]], equal=TRUE, sep=sepset)), sep=""))
         }
      } else {
         out <- bquote(paste(.(tname), "(df", .(sep), "=", .(sep), .(df), ")", .(sep), "=", .(sep), .(fmtx(val, digits[["test"]], flag=flag)), ", ", .(pname), .(sep), .(fmtp(pval, digits[["pval"]], equal=TRUE, sep=sepset)), sep=""))
         #paste0(tname, "(df = ", df, ") = ", fmtx(val, digits[["test"]]), ", ", pname, fmtp(pval, digits[["pval"]], equal=TRUE, sep=TRUE))
      }

   }

   if (format[[1]] == 2) {

      if (missing(df)) {
         if (!missing(df1) && !missing(df2)) {
            out <- bquote(paste(.(tname), .(sep), "=", .(sep), .(fmtx(val, digits[["test"]], flag=flag)), ", df1", .(sep), "=", .(sep), .(df1), ", df2", .(sep), "=", .(sep), .(round(df2,2)), ", ", .(pname), .(sep), .(fmtp(pval, digits[["pval"]], equal=TRUE, sep=sepset)), sep=""))
         } else {
            out <- bquote(paste(.(tname), .(sep), "=", .(sep), .(fmtx(val, digits[["test"]], flag=flag)), ", ", .(pname), .(sep), .(fmtp(pval, digits[["pval"]], equal=TRUE, sep=sepset)), sep=""))
         }
      } else {
         out <- bquote(paste(.(tname), .(sep), "=", .(sep), .(fmtx(val, digits[["test"]], flag=flag)), ", df", .(sep), "=", .(sep), .(df), ", ", .(pname), .(sep), .(fmtp(pval, digits[["pval"]], equal=TRUE, sep=sepset)), sep=""))
      }

   }

   if (call) {
      out$sep <- NULL
      return(out)
   } else {
      out <- eval(out)
      if (!quote)
         out <- noquote(out)
      return(out)
   }

}

############################################################################
