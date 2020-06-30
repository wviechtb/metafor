rma.peto <- function(ai, bi, ci, di, n1i, n2i,
data, slab, subset,
add=1/2, to="only0", drop00=TRUE, ### for add/to/drop00, 1st element for escalc(), 2nd for Peto's method
level=95, digits, verbose=FALSE, ...) {

   #########################################################################

   ###### setup

   mstyle <- .get.mstyle("crayon" %in% .packages())

   ### check argument specifications

   if (length(add) == 1L)
      add <- c(add, 0)

   if (length(add) != 2L)
      stop(mstyle$stop("Argument 'add' should specify one or two values (see 'help(rma.peto)')."))

   if (length(to) == 1L)
      to <- c(to, "none")

   if (length(to) != 2L)
      stop(mstyle$stop("Argument 'to' should specify one or two values (see 'help(rma.peto)')."))

   if (length(drop00) == 1L)
      drop00 <- c(drop00, FALSE)

   if (length(drop00) != 2L)
      stop(mstyle$stop("Argument 'drop00' should specify one or two values (see 'help(rma.peto)')."))

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   if (!is.element(to[1], c("all","only0","if0all","none")))
      stop(mstyle$stop("Unknown 'to' argument specified."))

   if (!is.element(to[2], c("all","only0","if0all","none")))
      stop(mstyle$stop("Unknown 'to' argument specified."))

   time.start <- proc.time()

   ### get ... argument and check for extra/superfluous arguments

   ddd <- list(...)

   .chkdots(ddd, c("outlist", "time"))

   measure <- "PETO" ### set measure here so that it can be added below

   ### set defaults for digits

   if (missing(digits)) {
      digits <- .set.digits(dmiss=TRUE)
   } else {
      digits <- .set.digits(digits, dmiss=FALSE)
   }

   ### set options(warn=1) if verbose > 2

   if (verbose > 2) {
      opwarn <- options(warn=1)
      on.exit(options(warn=opwarn$warn))
   }

   #########################################################################

   if (verbose)
      message(mstyle$message("\nExtracting data and computing yi/vi values ..."))

   ### check if data argument has been specified

   if (missing(data))
      data <- NULL

   if (is.null(data)) {
      data <- sys.frame(sys.parent())
   } else {
      if (!is.data.frame(data))
         data <- data.frame(data)
   }

   mf <- match.call()

   ### extract slab and subset values, possibly from the data frame specified via data (arguments not specified are NULL)

   mf.slab   <- mf[[match("slab",   names(mf))]]
   mf.subset <- mf[[match("subset", names(mf))]]
   slab   <- eval(mf.slab,   data, enclos=sys.frame(sys.parent()))
   subset <- eval(mf.subset, data, enclos=sys.frame(sys.parent()))

   ### extract/calculate ai,bi,ci,di,n1i,n2i values

   mf.ai  <- mf[[match("ai",  names(mf))]]
   mf.bi  <- mf[[match("bi",  names(mf))]]
   mf.ci  <- mf[[match("ci",  names(mf))]]
   mf.di  <- mf[[match("di",  names(mf))]]
   mf.n1i <- mf[[match("n1i", names(mf))]]
   mf.n2i <- mf[[match("n2i", names(mf))]]
   ai  <- eval(mf.ai,  data, enclos=sys.frame(sys.parent()))
   bi  <- eval(mf.bi,  data, enclos=sys.frame(sys.parent()))
   ci  <- eval(mf.ci,  data, enclos=sys.frame(sys.parent()))
   di  <- eval(mf.di,  data, enclos=sys.frame(sys.parent()))
   n1i <- eval(mf.n1i, data, enclos=sys.frame(sys.parent()))
   n2i <- eval(mf.n2i, data, enclos=sys.frame(sys.parent()))
   if (is.null(bi)) bi <- n1i - ai
   if (is.null(di)) di <- n2i - ci
   ni <- ai + bi + ci + di

   k <- length(ai) ### number of outcomes before subsetting
   k.all <- k

   ids <- seq_len(k)

   ### generate study labels if none are specified

   if (verbose)
      message(mstyle$message("Generating/extracting study labels ..."))

   if (is.null(slab)) {

      slab.null <- TRUE
      slab      <- ids

   } else {

      if (anyNA(slab))
         stop(mstyle$stop("NAs in study labels."))

      if (length(slab) != k)
         stop(mstyle$stop("Study labels not of same length as data."))

      if (is.factor(slab))
         slab <- as.character(slab)

      slab.null <- FALSE

   }

   ### if a subset of studies is specified

   if (!is.null(subset)) {

      if (verbose)
         message(mstyle$message("Subsetting ..."))

      ai   <- ai[subset]
      bi   <- bi[subset]
      ci   <- ci[subset]
      di   <- di[subset]
      ni   <- ni[subset]
      slab <- slab[subset]
      ids  <- ids[subset]
      k    <- length(ai)

   }

   ### check if study labels are unique; if not, make them unique

   if (anyDuplicated(slab))
      slab <- .make.unique(slab)

   ### calculate observed effect estimates and sampling variances

   dat <- escalc(measure="PETO", ai=ai, bi=bi, ci=ci, di=di, add=add[1], to=to[1], drop00=drop00[1])
   yi  <- dat$yi ### one or more yi/vi pairs may be NA/NA
   vi  <- dat$vi ### one or more yi/vi pairs may be NA/NA

   ### if drop00[2]=TRUE, set counts to NA for studies that have no events (or all events) in both arms

   if (drop00[2]) {
      id00 <- c(ai == 0L & ci == 0L) | c(bi == 0L & di == 0L)
      id00[is.na(id00)] <- FALSE
      ai[id00] <- NA
      bi[id00] <- NA
      ci[id00] <- NA
      di[id00] <- NA
   }

   ### save the actual cell frequencies and yi/vi values (including potential NAs)

   ai.f <- ai
   bi.f <- bi
   ci.f <- ci
   di.f <- di
   yi.f <- yi
   vi.f <- vi
   ni.f <- ni

   k.f <- k ### total number of tables including all NAs

   ### check for NAs in table data and act accordingly

   has.na <- is.na(ai) | is.na(bi) | is.na(ci) | is.na(di)
   not.na <- !has.na

   if (any(has.na)) {

      if (verbose)
         message(mstyle$message("Handling NAs in table data ..."))

      if (na.act == "na.omit" || na.act == "na.exclude" || na.act == "na.pass") {
         ai   <- ai[not.na]
         bi   <- bi[not.na]
         ci   <- ci[not.na]
         di   <- di[not.na]
         k    <- length(ai)
         warning(mstyle$warning("Tables with NAs omitted from model fitting."), call.=FALSE)
      }

      if (na.act == "na.fail")
         stop(mstyle$stop("Missing values in tables."))

   }

   ### at least one study left?

   if (k < 1)
      stop(mstyle$stop("Processing terminated since k = 0."))

   ### check for NAs in yi/vi and act accordingly

   yivi.na <- is.na(yi) | is.na(vi)
   not.na.yivi <- !yivi.na

   if (any(yivi.na)) {

      if (verbose)
         message(mstyle$message("Handling NAs in yi/vi ..."))

      if (na.act == "na.omit" || na.act == "na.exclude" || na.act == "na.pass") {

         yi <- yi[not.na.yivi]
         vi <- vi[not.na.yivi]
         ni <- ni[not.na.yivi]
         warning(mstyle$warning("Some yi/vi values are NA."), call.=FALSE)

         attr(yi, "measure") <- measure ### add measure attribute back
         attr(yi, "ni")      <- ni      ### add ni attribute back

      }

      if (na.act == "na.fail")
         stop(mstyle$stop("Missing yi/vi values."))

   }

   k.yi <- length(yi) ### number of yi/vi pairs that are not NA (needed for QE df and fit.stats calculation)

   ### add/to procedures for the 2x2 tables for the actual meta-analysis
   ### note: technically, nothing needs to be added, but Stata/RevMan add 1/2 by default for only0 studies (but drop studies with no/all events)

   if (to[2] == "all") {

      ### always add to all cells in all studies

      ai <- ai + add[2]
      bi <- bi + add[2]
      ci <- ci + add[2]
      di <- di + add[2]

   }

   if (to[2] == "only0") {

      ### add to cells in studies with at least one 0 entry

      id0 <- c(ai == 0L | bi == 0L | ci == 0L | di == 0L)

      ai[id0] <- ai[id0] + add[2]
      bi[id0] <- bi[id0] + add[2]
      ci[id0] <- ci[id0] + add[2]
      di[id0] <- di[id0] + add[2]

   }

   if (to[2] == "if0all") {

      ### add to cells in all studies if there is at least one 0 entry

      id0 <- c(ai == 0L | bi == 0L | ci == 0L | di == 0L)

      if (any(id0)) {

         ai <- ai + add[2]
         bi <- bi + add[2]
         ci <- ci + add[2]
         di <- di + add[2]

      }

   }

   n1i <- ai + bi
   n2i <- ci + di
   Ni  <- ai + bi + ci + di

   #########################################################################

   level <- ifelse(level == 0, 1, ifelse(level >= 1, (100-level)/100, ifelse(level > .5, 1-level, level)))

   ###### model fitting, test statistics, and confidence intervals

   if (verbose)
      message(mstyle$message("Model fitting ..."))

   xt <- ai + ci ### frequency of outcome1 in both groups combined
   yt <- bi + di ### frequency of outcome2 in both groups combined
   Ei <- xt * n1i / Ni
   Vi <- xt * yt * (n1i/Ni) * (n2i/Ni) / (Ni - 1) ### 0 when xt = 0 or yt = 0 in a table

   sumVi <- sum(Vi)

   if (sumVi == 0L) ### sumVi = 0 when xt or yt = 0 in *all* tables
      stop(mstyle$stop("One of the two outcomes never occurred in any of the tables. Peto's method cannot be used."))

   beta  <- sum(ai - Ei) / sumVi
   se    <- sqrt(1/sumVi)
   zval  <- beta / se
   pval  <- 2*pnorm(abs(zval), lower.tail=FALSE)
   ci.lb <- beta - qnorm(level/2, lower.tail=FALSE) * se
   ci.ub <- beta + qnorm(level/2, lower.tail=FALSE) * se

   names(beta) <- "intrcpt"
   vb <- matrix(se^2, dimnames=list("intrcpt", "intrcpt"))

   #########################################################################

   ### heterogeneity test (Peto's method)

   if (verbose)
      message(mstyle$message("Heterogeneity testing ..."))

   k.pos <- sum(Vi > 0) ### number of tables with positive sampling variance
   Vi[Vi == 0] <- NA    ### set 0 sampling variances to NA
   QE <- max(0, sum((ai - Ei)^2 / Vi, na.rm=TRUE) - sum(ai - Ei)^2 / sum(Vi, na.rm=TRUE))

   if (k.pos > 1) {
      QEp <- pchisq(QE, df=k.yi-1, lower.tail=FALSE)
      I2  <- max(0, 100 * (QE - (k.yi-1)) / QE)
      H2  <- QE / (k.yi-1)
   } else {
      QEp <- 1
      I2  <- 0
      H2  <- 1
   }

   wi  <- 1/vi
   RSS <- sum(wi*(yi-beta)^2)

   #########################################################################

   ###### fit statistics

   if (verbose)
      message(mstyle$message("Computing fit statistics and log likelihood ..."))

   ll.ML     <- -1/2 * (k.yi)   * log(2*base::pi)                   - 1/2 * sum(log(vi))                      - 1/2 * RSS
   ll.REML   <- -1/2 * (k.yi-1) * log(2*base::pi) + 1/2 * log(k.yi) - 1/2 * sum(log(vi)) - 1/2 * log(sum(wi)) - 1/2 * RSS
   dev.ML    <- -2 * (ll.ML - sum(dnorm(yi, mean=yi, sd=sqrt(vi), log=TRUE)))
   AIC.ML    <- -2 * ll.ML   + 2
   BIC.ML    <- -2 * ll.ML   + log(k.yi)
   AICc.ML   <- -2 * ll.ML   + 2 * max(k.yi, 3) / (max(k.yi, 3) - 2)
   dev.REML  <- -2 * (ll.REML - 0)
   AIC.REML  <- -2 * ll.REML + 2
   BIC.REML  <- -2 * ll.REML + log(k.yi-1)
   AICc.REML <- -2 * ll.REML + 2 * max(k.yi-1, 3) / (max(k.yi-1, 3) - 2)

   fit.stats <- matrix(c(ll.ML, dev.ML, AIC.ML, BIC.ML, AICc.ML, ll.REML, dev.REML, AIC.REML, BIC.REML, AICc.REML), ncol=2, byrow=FALSE)
   dimnames(fit.stats) <- list(c("ll","dev","AIC","BIC","AICc"), c("ML","REML"))
   fit.stats <- data.frame(fit.stats)

   #########################################################################

   ###### prepare output

   if (verbose)
      message(mstyle$message("Preparing output ..."))

   parms     <- 1
   p         <- 1
   p.eff     <- 1
   k.eff     <- k
   tau2      <- 0
   X.f       <- cbind(rep(1,k.f))
   intercept <- TRUE
   int.only  <- TRUE
   btt       <- 1
   m         <- 1

   method    <- "FE"
   weighted  <- TRUE
   test      <- "z"
   dfs       <- NA

   if (is.null(ddd$outlist)) {

      res <- list(b=beta, beta=beta, se=se, zval=zval, pval=pval, ci.lb=ci.lb, ci.ub=ci.ub, vb=vb,
                  tau2=tau2, tau2.f=tau2,
                  I2=I2, H2=H2,
                  QE=QE, QEp=QEp,
                  k=k, k.f=k.f, k.yi=k.yi, k.pos=k.pos, k.eff=k.eff, k.all=k.all, p=p, p.eff=p.eff, parms=parms,
                  int.only=int.only, intercept=intercept,
                  yi=yi, vi=vi, yi.f=yi.f, vi.f=vi.f, X.f=X.f, ai=ai, bi=bi, ci=ci, di=di, ai.f=ai.f, bi.f=bi.f, ci.f=ci.f, di.f=di.f, ni=ni, ni.f=ni.f,
                  ids=ids, not.na=not.na, subset=subset, not.na.yivi=not.na.yivi, slab=slab, slab.null=slab.null,
                  measure=measure, method=method, weighted=weighted,
                  test=test, dfs=dfs, btt=btt, m=m,
                  digits=digits, level=level,
                  add=add, to=to, drop00=drop00,
                  fit.stats=fit.stats,
                  formula.yi=NULL, formula.mods=NULL, version=packageVersion("metafor"), call=mf)

   }

   if (!is.null(ddd$outlist)) {
      if (ddd$outlist == "minimal") {
         res <- list(b=beta, beta=beta, se=se, zval=zval, pval=pval, ci.lb=ci.lb, ci.ub=ci.ub, vb=vb,
                     tau2=tau2,
                     I2=I2, H2=H2,
                     QE=QE, QEp=QEp,
                     k=k, k.pos=k.pos, k.eff=k.eff, p=p, p.eff=p.eff, parms=parms,
                     int.only=int.only,
                     measure=measure, method=method,
                     test=test, dfs=dfs, btt=btt, m=m,
                     digits=digits,
                     fit.stats=fit.stats)
      } else {
         res <- eval(parse(text=paste0("list(", ddd$outlist, ")")))
      }
   }

   time.end <- proc.time()
   res$time <- unname(time.end - time.start)[3]

   if (.isTRUE(ddd$time))
      .print.time(res$time)

   if (verbose || .isTRUE(ddd$time))
      cat("\n")

   class(res) <- c("rma.peto", "rma")
   return(res)

}
