conv.2x2 <- function(ori, ri, x2i, ni, n1i, n2i, sens, spec, ppv, npv, correct=TRUE, drop01=TRUE, data,
                     include, var.names=c("ai","bi","ci","di"), append=TRUE, replace="ifna", ...) {

   mstyle <- .get.mstyle()

   if (is.logical(replace)) {
      if (isTRUE(replace)) {
         replace <- "all"
      } else {
         replace <- "ifna"
      }
   }

   replace <- match.arg(replace, c("ifna","all"))

   ### get ... argument and check for extra/superfluous arguments

   ddd <- list(...)

   .chkdots(ddd, "method")

   method <- .chkddd(ddd$method, "optim", match.arg(ddd$method, c("simple","optim", "search")))

   #########################################################################

   if (missing(data))
      data <- NULL

   has.data <- !is.null(data)

   if (is.null(data)) {
      data <- sys.frame(sys.parent())
   } else {
      if (!is.data.frame(data))
         data <- data.frame(data)
   }

   # checks on 'var.names' argument

   if (length(var.names) != 4L)
      stop(mstyle$stop("Argument 'var.names' must be of length 4."))

   if (any(var.names != make.names(var.names, unique=TRUE))) {
      var.names <- make.names(var.names, unique=TRUE)
      warning(mstyle$warning(paste0("Argument 'var.names' does not contain syntactically valid variable names.\nVariable names adjusted to: var.names = c('", var.names[1], "','", var.names[2], "','", var.names[3], "','", var.names[4], "').")), call.=FALSE)
   }

   #########################################################################

   mf <- match.call()

   ori     <- .getx("ori",     mf=mf, data=data, checknumeric=TRUE)
   ri      <- .getx("ri",      mf=mf, data=data, checknumeric=TRUE)
   x2i     <- .getx("x2i",     mf=mf, data=data, checknumeric=TRUE)
   sens    <- .getx("sens",    mf=mf, data=data, checknumeric=TRUE)
   spec    <- .getx("spec",    mf=mf, data=data, checknumeric=TRUE)
   ppv     <- .getx("ppv",     mf=mf, data=data, checknumeric=TRUE)
   npv     <- .getx("npv",     mf=mf, data=data, checknumeric=TRUE)
   ni      <- .getx("ni",      mf=mf, data=data, checknumeric=TRUE)
   n1i     <- .getx("n1i",     mf=mf, data=data, checknumeric=TRUE)
   n2i     <- .getx("n2i",     mf=mf, data=data, checknumeric=TRUE)
   correct <- .getx("correct", mf=mf, data=data, default=TRUE)
   include <- .getx("include", mf=mf, data=data)

   if (!.equal.length(ori, ri, x2i, sens, spec, ppv, npv, ni, n1i, n2i))
      stop(mstyle$stop("Supplied data vectors are not all of the same length."))

   k <- .maxlength(ori, ri, x2i, sens, spec, ppv, npv, ni, n1i, n2i)

   if (is.null(ori))
      ori <- rep(NA_real_, k)
   if (is.null(ri))
      ri <- rep(NA_real_, k)
   if (is.null(x2i))
      x2i <- rep(NA_real_, k)
   if (is.null(sens))
      sens <- rep(NA_real_, k)
   if (is.null(spec))
      spec <- rep(NA_real_, k)
   if (is.null(ppv))
      ppv <- rep(NA_real_, k)
   if (is.null(npv))
      npv <- rep(NA_real_, k)
   if (is.null(ni))
      ni <- rep(NA_real_, k)
   if (is.null(n1i))
      n1i <- rep(NA_real_, k)
   if (is.null(n2i))
      n2i <- rep(NA_real_, k)

   # handle 'correct' argument

   correct <- .expand1(correct, k)

   if (length(correct) != k)
      stop(mstyle$stop(paste0("Length of the 'correct' argument (", length(correct), ") does not match the length of the data (", k, ").")))

   correct[is.na(correct)] <- TRUE

   # if 'include' is NULL, set to TRUE vector

   if (is.null(include))
      include <- rep(TRUE, k)

   # turn numeric 'include' vector into a logical vector

   include <- .chksubset(include, k, stoponk0=FALSE)

   # set inputs to NA for rows not to be included

   ori[!include]  <- NA_real_
   ri[!include]   <- NA_real_
   x2i[!include]  <- NA_real_
   sens[!include] <- NA_real_
   spec[!include] <- NA_real_
   ppv[!include]  <- NA_real_
   npv[!include]  <- NA_real_
   ni[!include]   <- NA_real_
   n1i[!include]  <- NA_real_
   n2i[!include]  <- NA_real_

   # round ni, n1i, and n2i

   ni  <- round(ni)
   n1i <- round(n1i)
   n2i <- round(n2i)

   # checks on values

   if (any(c(ni < 0, n1i < 0, n2i < 0), na.rm=TRUE))
      stop(mstyle$stop("One or more sample sizes or marginal counts are negative."))

   if (any(c(n1i > ni, n2i > ni), na.rm=TRUE))
      stop(mstyle$stop("One or more marginal counts are larger than the sample sizes."))

   if (any(abs(ri) > 1, na.rm=TRUE))
      stop(mstyle$stop("One or more phi coefficients are > 1 or < -1."))

   if (any(sens < 0 | sens > 1, na.rm=TRUE))
      stop(mstyle$stop("One or more sensitivity values are < 0 or > 1."))

   if (any(spec < 0 | spec > 1, na.rm=TRUE))
      stop(mstyle$stop("One or more specificity values are < 0 or > 1."))

   if (any(ppv < 0 | ppv > 1, na.rm=TRUE))
      stop(mstyle$stop("One or more positive predictive values are < 0 or > 1."))

   if (any(npv < 0 | npv > 1, na.rm=TRUE))
      stop(mstyle$stop("One or more negative predictive values are < 0 or > 1."))

   # compute marginal proportions for the two variables

   p1i <- n1i / ni
   p2i <- n2i / ni

   #########################################################################

   p11i <- rep(NA_real_, k)
   ai <- rep(NA_real_, k)
   bi <- rep(NA_real_, k)
   ci <- rep(NA_real_, k)
   di <- rep(NA_real_, k)

   incons <- rep(NA, k)

   for (i in seq_len(k)) {

      if (is.na(ni[i]))
         next

      if (sum(!is.na(sens[i]), !is.na(spec[i]), !is.na(ppv[i]), !is.na(npv[i])) >= 3L) {

         # if drop01=TRUE, drop studies where sens, spec, ppv, and/or npv is equal to 0 or 1

         if (drop01 && (isTRUE(sens[i]==0) || isTRUE(sens[i]==1) || isTRUE(spec[i]==0) || isTRUE(spec[i]==1) || isTRUE(ppv[i]==0) || isTRUE(ppv[i]==1) || isTRUE(npv[i]==0) || isTRUE(npv[i]==1)))
            next

         # if at least three of (sens, spec, ppv, npv) are available, use reconstruction based on these statistics

         if (is.na(sens[i]))
            sens[i] <- (1-spec[i]) * ppv[i] * npv[i] / ((1-spec[i]) * ppv[i] * npv[i] + spec[i] * (1-ppv[i])*(1-npv[i]))
         if (is.na(spec[i]))
            spec[i] <- (1-sens[i]) * ppv[i] * npv[i] / ((1-sens[i]) * ppv[i] * npv[i] + sens[i] * (1-ppv[i])*(1-npv[i]))
         if (is.na(ppv[i]))
            ppv[i] <- sens[i] * spec[i] * (1-npv[i]) / (sens[i] * spec[i] * (1-npv[i]) + (1-sens[i])*(1-spec[i]) * npv[i])
         if (is.na(npv[i]))
            npv[i] <- sens[i] * spec[i] * (1-ppv[i]) / (sens[i] * spec[i] * (1-ppv[i]) + (1-sens[i])*(1-spec[i]) * ppv[i])

         if (FALSE) {

            # check consistency of inputs (skipped for now - discrepancies can be large when some of the diagnostic statistics in the denominator are close to 0)

            sens.imp <- (1-spec[i]) * ppv[i] * npv[i] / ((1-spec[i]) * ppv[i] * npv[i] + spec[i] * (1-ppv[i])*(1-npv[i]))
            spec.imp <- (1-sens[i]) * ppv[i] * npv[i] / ((1-sens[i]) * ppv[i] * npv[i] + sens[i] * (1-ppv[i])*(1-npv[i]))
            ppv.imp  <- sens[i] * spec[i] * (1-npv[i]) / (sens[i] * spec[i] * (1-npv[i]) + (1-sens[i])*(1-spec[i]) * npv[i])
            npv.imp  <- sens[i] * spec[i] * (1-ppv[i]) / (sens[i] * spec[i] * (1-ppv[i]) + (1-sens[i])*(1-spec[i]) * ppv[i])

            cutoff <- 0.1
            incons[i] <- abs(sens[i] - sens.imp) > cutoff || abs(spec[i] - spec.imp) > cutoff || abs(ppv[i] - ppv.imp) > cutoff || abs(npv[i] - npv.imp) > cutoff

            if (incons[i])
               next

         }

         #print(c(sens=sens[i], spec=spec[i], ppv=ppv[i], npv=npv[i]))

         if (isTRUE(sens[i] * (1 - ppv[i]) + ppv[i] * (1 - spec[i]) == 0))
            next

         if (method=="simple")
            tmp <- .rec2x2diag(sens=sens[i], spec=spec[i], ppv=ppv[i], ni=ni[i], round=FALSE)

         if (method=="optim") {

            obs <- c(sens=sens[i], spec=spec[i], ppv=ppv[i], npv=npv[i])
            start <- c(sens=sens[i], spec=spec[i], ppv=ppv[i])
            res <- try(optim(start, .funconv2x2diag, method="BFGS", obs=obs, ni=ni[i]), silent=TRUE)
            #res <- try(optim(start, .funconv2x2diag, method="L-BFGS-B", obs=obs, ni=ni[i], lower=pmax(0.005,start-0.005), upper=pmin(0.995,start+0.005)), silent=TRUE)
            if (inherits(res, "try-error")) {
               tmp <- rep(NA, 4)
            } else {
               tmp <- .rec2x2diag(sens=res$par[1], spec=res$par[2], ppv=res$par[3], ni=ni[i], round=FALSE)
            }

         }

         if (method=="search") {

            obs <- c(sens=sens[i], spec=spec[i], ppv=ppv[i], npv=npv[i])

            sens.int <- seq(max(0, sens[i] - 0.005), min(1, sens[i] + 0.005), length=21)
            spec.int <- seq(max(0, spec[i] - 0.005), min(1, spec[i] + 0.005), length=21)
            ppv.int  <- seq(max(0, ppv[i]  - 0.005), min(1, ppv[i]  + 0.005), length=21)

            loss <- array(NA_real_, dim=c(length(sens.int), length(spec.int), length(ppv.int)))

            for (i1 in 1:length(sens.int)) {
               for (i2 in 1:length(spec.int)) {
                  for (i3 in 1:length(ppv.int)) {
                     loss[i1,i2,i3] <- .funconv2x2diag(c(sens=sens.int[i1], spec=spec.int[i2], ppv=ppv.int[i3]), obs=obs, ni=ni[i], round=FALSE)
                  }
               }
            }

            hits <- which(loss == min(loss, na.rm=TRUE), arr.ind=TRUE)
            rec <- apply(hits, 1, function(x) {
               sens.sel <- sens.int[x[1]]
               spec.sel <- spec.int[x[2]]
               ppv.sel  <- ppv.int[x[3]]
               rec <- .rec2x2diag(sens.sel, spec.sel, ppv.sel, ni=ni[i], round=FALSE)
            })
            #print(table(apply(rec, 2, function(x) x == rec[,1])))
            tmp <- rec[,1]

         }

         tmp <- .largestremaindermethod(tmp, ni[i])
         ai[i] <- tmp[1]
         bi[i] <- tmp[2]
         ci[i] <- tmp[3]
         di[i] <- tmp[4]

      } else {

         if (is.na(n1i[i]) || is.na(n2i[i]))
            next

         if (!is.na(ori[i])) {

            if (ori[i] == 1) {

               p11i[i] <- n1i[i] * n2i[i] / ni[i]^2

            } else {

               p1. <- p1i[i]
               p2. <- 1-p1i[i]
               p.1 <- p2i[i]
               p.2 <- 1-p2i[i]

               x <- ori[i] * (p1. + p.1) + p2. - p.1
               y <- sqrt(x^2 - 4 * p1. * p.1 * ori[i] * (ori[i]-1))

               p11i[i] <- (x - y) / (2 * (ori[i] - 1))

            }

         }

         # note: when x2i=0, then sign(0) = 0 and hence ri is automatically 0, which is correct (i.e., we do not want to use the continuity correction in this case)
         if (is.na(ri[i]) && !is.na(x2i[i])) {
            if (correct[i]) {
               ri[i] <- sign(x2i[i]) * (sqrt(abs(x2i[i])/ni[i]) + ni[i] / (2*sqrt(n1i[i]*(ni[i]-n1i[i])*n2i[i]*(ni[i]-n2i[i]))))
            } else {
               ri[i] <- sign(x2i[i]) * sqrt(abs(x2i[i])/ni[i])
            }
         }

         if (is.na(p11i[i]) && !is.na(ri[i]))
            p11i[i] <- p1i[i]*p2i[i] + ri[i] * sqrt(p1i[i]*(1-p1i[i])*p2i[i]*(1-p2i[i]))

         ai[i] <- round(ni[i] * p11i[i])
         bi[i] <- n1i[i] - ai[i]
         ci[i] <- n2i[i] - ai[i]
         di[i] <- ni[i] - ai[i] - bi[i] - ci[i]

      }

   }

   #print(matrix(c(ai,bi,ci,di), nrow=2, byrow=TRUE))

   if (any(incons, na.rm=TRUE)) {

      warning(mstyle$warning(paste0("There are inconsistency diagnostic statistics in table", ifelse(sum(incons, na.rm=TRUE) > 1, "s ", " "),
                                    paste0(which(incons), collapse=","), ".")), call.=FALSE)

   }

   # check for negative cell frequencies

   hasneg <- (ai < 0) | (bi < 0) | (ci < 0) | (di < 0)

   if (any(hasneg, na.rm=TRUE)) {

      warning(mstyle$warning(paste0("There are negative cell frequencies in table", ifelse(sum(hasneg, na.rm=TRUE) > 1, "s ", " "),
                                    paste0(which(hasneg), collapse=","), ".")), call.=FALSE)

      ai[hasneg] <- NA_real_
      bi[hasneg] <- NA_real_
      ci[hasneg] <- NA_real_
      di[hasneg] <- NA_real_

   }

   #########################################################################

   if (has.data && append) {

      if (is.element(var.names[1], names(data))) {
         if (replace=="ifna") {
            data[[var.names[1]]] <- replmiss(data[[var.names[1]]], ai)
         } else {
            data[[var.names[1]]][!is.na(ai)] <- ai[!is.na(ai)]
         }
      } else {
         data <- cbind(data, ai)
         names(data)[length(names(data))] <- var.names[1]
      }

      if (is.element(var.names[2], names(data))) {
         if (replace=="ifna") {
            data[[var.names[2]]] <- replmiss(data[[var.names[2]]], bi)
         } else {
            data[[var.names[2]]][!is.na(bi)] <- bi[!is.na(bi)]
         }
      } else {
         data <- cbind(data, bi)
         names(data)[length(names(data))] <- var.names[2]
      }

      if (is.element(var.names[3], names(data))) {
         if (replace=="ifna") {
            data[[var.names[3]]] <- replmiss(data[[var.names[3]]], ci)
         } else {
            data[[var.names[3]]][!is.na(ci)] <- ai[!is.na(ci)]
         }
      } else {
         data <- cbind(data, ci)
         names(data)[length(names(data))] <- var.names[3]
      }

      if (is.element(var.names[4], names(data))) {
         if (replace=="ifna") {
            data[[var.names[4]]] <- replmiss(data[[var.names[4]]], di)
         } else {
            data[[var.names[4]]][!is.na(di)] <- ai[!is.na(di)]
         }
      } else {
         data <- cbind(data, di)
         names(data)[length(names(data))] <- var.names[4]
      }

   } else {

      data <- data.frame(ai, bi, ci, di)
      names(data) <- var.names

   }

   return(data)

}
