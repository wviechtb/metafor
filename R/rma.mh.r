rma.mh   <- function(ai, bi, ci, di, n1i, n2i, x1i, x2i, t1i, t2i, measure="OR",
data, slab, subset,
add=1/2, to="only0", drop00=TRUE, ### for add/to/drop00, 1st element for escalc(), 2nd for MH method
correct=TRUE, level=95, verbose=FALSE, digits, ...) {

   #########################################################################

   ###### setup

   mstyle <- .get.mstyle("crayon" %in% .packages())

   ### check argument specifications

   if (!is.element(measure, c("OR","RR","RD","IRR","IRD")))
      stop(mstyle$stop("Mantel-Haenszel method can only be used with measures OR, RR, RD, IRR, and IRD."))

   if (length(add) == 1L)
      add <- c(add, 0)

   if (length(add) != 2L)
      stop(mstyle$stop("Argument 'add' should specify one or two values (see 'help(rma.mh)')."))

   if (length(to) == 1L)
      to <- c(to, "none")

   if (length(to) != 2L)
      stop(mstyle$stop("Argument 'to' should specify one or two values (see 'help(rma.mh)')."))

   if (length(drop00) == 1L)
      drop00 <- c(drop00, FALSE)

   if (length(drop00) != 2L)
      stop(mstyle$stop("Argument 'drop00' should specify one or two values (see 'help(rma.mh)')."))

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

   .chkdots(ddd, c("outlist", "onlyo1", "addyi", "addvi", "time"))

   ### set defaults or get onlyo1, addyi, and addvi arguments

   onlyo1 <- ifelse(is.null(ddd$onlyo1), FALSE, ddd$onlyo1)
   addyi  <- ifelse(is.null(ddd$addyi),  TRUE,  ddd$addyi)
   addvi  <- ifelse(is.null(ddd$addvi),  TRUE,  ddd$addvi)

   ### set defaults for digits

   if (missing(digits)) {
      digits <- .set.digits(dmiss=TRUE)
   } else {
      digits <- .set.digits(digits, dmiss=FALSE)
   }

   ### set options(warn=1) if verbose > 2

   if (verbose > 2) {
      opwarn <- options(warn=1)
      on.exit(options(warn=opwarn$warn), add=TRUE)
   }

   #########################################################################

   if (verbose) .space()

   if (verbose)
      message(mstyle$message("Extracting data and computing yi/vi values ..."))

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

   slab   <- .getx("slab",   mf=mf, data=data)
   subset <- .getx("subset", mf=mf, data=data)

   #########################################################################

   ### for RR, OR, and RD: extract/calculate ai,bi,ci,di,n1i,n2i values

   if (is.element(measure, c("RR","OR","RD"))) {

      x1i <- x2i <- t1i <- t2i <- NA_real_

      ai  <- .getx("ai",  mf=mf, data=data, checknumeric=TRUE)
      bi  <- .getx("bi",  mf=mf, data=data, checknumeric=TRUE)
      ci  <- .getx("ci",  mf=mf, data=data, checknumeric=TRUE)
      di  <- .getx("di",  mf=mf, data=data, checknumeric=TRUE)
      n1i <- .getx("n1i", mf=mf, data=data, checknumeric=TRUE)
      n2i <- .getx("n2i", mf=mf, data=data, checknumeric=TRUE)

      if (is.null(bi)) bi <- n1i - ai
      if (is.null(di)) di <- n2i - ci
      ni <- ai + bi + ci + di

      k <- length(ai) ### number of outcomes before subsetting
      k.all <- k

      if (length(ai)==0L || length(bi)==0L || length(ci)==0L || length(di)==0L)
         stop(mstyle$stop("Must specify arguments ai, bi, ci, di (or ai, ci, n1i, n2i) for this measure."))

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

         subset <- .chksubset(subset, k)

         ai   <- .getsubset(ai,   subset)
         bi   <- .getsubset(bi,   subset)
         ci   <- .getsubset(ci,   subset)
         di   <- .getsubset(di,   subset)
         ni   <- .getsubset(ni,   subset)
         slab <- .getsubset(slab, subset)
         ids  <- .getsubset(ids,  subset)

         k <- length(ai)

      }

      ### check if study labels are unique; if not, make them unique

      if (anyDuplicated(slab))
         slab <- .make.unique(slab)

      ### calculate observed effect estimates and sampling variances

      dat <- .do.call(escalc, measure=measure, ai=ai, bi=bi, ci=ci, di=di, add=add[1], to=to[1], drop00=drop00[1], onlyo1=onlyo1, addyi=addyi, addvi=addvi)
      yi <- dat$yi ### one or more yi/vi pairs may be NA/NA
      vi <- dat$vi ### one or more yi/vi pairs may be NA/NA

      ### if drop00[2]=TRUE, set counts to NA for studies that have no events (or all events) in both arms

      if (drop00[2]) {
         id00 <- c(ai == 0L & ci == 0L) | c(bi == 0L & di == 0L)
         id00[is.na(id00)] <- FALSE
         ai[id00] <- NA_real_
         bi[id00] <- NA_real_
         ci[id00] <- NA_real_
         di[id00] <- NA_real_
      }

      ### save the actual cell frequencies and yi/vi values (including potential NAs)

      outdat.f <- list(ai=ai, bi=bi, ci=ci, di=di)

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
            ai <- ai[not.na]
            bi <- bi[not.na]
            ci <- ci[not.na]
            di <- di[not.na]
            k  <- length(ai)
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

   }

   #########################################################################

   ### for IRR and IRD: extract/calculate x1i,x2i,t1i,t2i values

   if (is.element(measure, c("IRR","IRD"))) {

      ai <- bi <- ci <- di <- NA_real_

      x1i <- .getx("x1i", mf=mf, data=data)
      x2i <- .getx("x2i", mf=mf, data=data)
      t1i <- .getx("t1i", mf=mf, data=data)
      t2i <- .getx("t2i", mf=mf, data=data)

      ni  <- t1i + t2i

      k <- length(x1i) ### number of outcomes before subsetting
      k.all <- k

      if (length(x1i)==0L || length(x2i)==0L || length(t1i)==0L || length(t2i)==0L)
         stop(mstyle$stop("Must specify arguments x1i, x2i, t1i, t2i for this measure."))

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

         slab.null <- FALSE

      }

      ### if a subset of studies is specified

      if (!is.null(subset)) {

         if (verbose)
            message(mstyle$message("Subsetting ..."))

         subset <- .chksubset(subset, k)

         x1i  <- .getsubset(x1i,  subset)
         x2i  <- .getsubset(x2i,  subset)
         t1i  <- .getsubset(t1i,  subset)
         t2i  <- .getsubset(t2i,  subset)
         ni   <- .getsubset(ni,   subset)
         slab <- .getsubset(slab, subset)
         ids  <- .getsubset(ids,  subset)

         k <- length(x1i)

      }

      ### check if study labels are unique; if not, make them unique

      if (anyDuplicated(slab))
         slab <- .make.unique(slab)

      ### calculate observed effect estimates and sampling variances

      dat <- .do.call(escalc, measure=measure, x1i=x1i, x2i=x2i, t1i=t1i, t2i=t2i, add=add[1], to=to[1], drop00=drop00[1], onlyo1=onlyo1, addyi=addyi, addvi=addvi)
      yi  <- dat$yi ### one or more yi/vi pairs may be NA/NA
      vi  <- dat$vi ### one or more yi/vi pairs may be NA/NA

      ### if drop00[2]=TRUE, set counts to NA for studies that have no events in both arms

      if (drop00[2]) {
         id00 <- c(x1i == 0L & x2i == 0L)
         id00[is.na(id00)] <- FALSE
         x1i[id00] <- NA_real_
         x2i[id00] <- NA_real_
      }

      ### save the actual cell frequencies and yi/vi values (including potential NAs)

      outdat.f <- list(x1i=x1i, x2i=x2i, t1i=t1i, t2i=t2i)

      yi.f <- yi
      vi.f <- vi
      ni.f <- ni

      k.f <- k ### total number of tables including all NAs

      ### check for NAs in table data and act accordingly

      has.na <- is.na(x1i) | is.na(x2i) | is.na(t1i) | is.na(t2i)
      not.na <- !has.na

      if (any(has.na)) {

         if (verbose)
            message(mstyle$message("Handling NAs in table data ..."))

         if (na.act == "na.omit" || na.act == "na.exclude" || na.act == "na.pass") {
            x1i  <- x1i[not.na]
            x2i  <- x2i[not.na]
            t1i  <- t1i[not.na]
            t2i  <- t2i[not.na]
            k    <- length(x1i)
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

      k.yi <- length(yi) ### number of yi/vi pairs that are not NA (needed for QE df and fitstats calculation)

      ### add/to procedures for the 2x2 tables for the actual meta-analysis
      ### note: technically, nothing needs to be added

      if (to[2] == "all") {

         ### always add to all cells in all studies

         x1i <- x1i + add[2]
         x2i <- x2i + add[2]

      }

      if (to[2] == "only0") {

         ### add to cells in studies with at least one 0 entry

         id0 <- c(x1i == 0L | x2i == 0L)

         x1i[id0] <- x1i[id0] + add[2]
         x2i[id0] <- x2i[id0] + add[2]

      }

      if (to[2] == "if0all") {

         ### add to cells in all studies if there is at least one 0 entry

         id0 <- c(x1i == 0L | x2i == 0L)

         if (any(id0)) {

            x1i <- x1i + add[2]
            x2i <- x2i + add[2]

         }

      }

      Ti <- t1i + t2i

   }

   #########################################################################

   level <- .level(level)

   CO <- COp <- MH <- MHp <- BD <- BDp <- TA <- TAp <- NA_real_
   k.pos <- NA_integer_

   ###### model fitting, test statistics, and confidence intervals

   if (verbose)
      message(mstyle$message("Model fitting ..."))

   if (measure == "OR") {

      Pi <- ai/Ni + di/Ni
      Qi <- bi/Ni + ci/Ni
      Ri <- (ai/Ni) * di
      Si <- (bi/Ni) * ci
      R  <- sum(Ri)
      S  <- sum(Si)

      if (identical(R,0) || identical(S,0) || identical(R,0L) || identical(S,0L)) {
         beta.exp <- NA_real_
         beta     <- NA_real_
         se       <- NA_real_
         zval     <- NA_real_
         pval     <- NA_real_
         ci.lb    <- NA_real_
         ci.ub    <- NA_real_
      } else {
         beta.exp <- R/S
         beta     <- log(beta.exp)
         se       <- sqrt(1/2 * (sum(Pi*Ri)/R^2 + sum(Pi*Si + Qi*Ri)/(R*S) + sum(Qi*Si)/S^2)) ### based on Robins et al. (1986)
         zval     <- beta / se
         pval     <- 2*pnorm(abs(zval), lower.tail=FALSE)
         ci.lb    <- beta - qnorm(level/2, lower.tail=FALSE) * se
         ci.ub    <- beta + qnorm(level/2, lower.tail=FALSE) * se
      }

      names(beta) <- "intrcpt"
      vb <- matrix(se^2, dimnames=list("intrcpt", "intrcpt"))

      ### Cochran and Cochran-Mantel-Haenszel Statistics

      xt <- ai + ci
      yt <- bi + di

      if (identical(sum(xt),0) || identical(sum(yt),0) || identical(sum(xt),0L) || identical(sum(yt),0L)) {
         CO  <- NA_real_
         COp <- NA_real_
         MH  <- NA_real_
         MHp <- NA_real_
      } else {
         CO  <- (abs(sum(ai - (n1i/Ni)*xt)) - ifelse(correct, 0.5, 0))^2 / sum((n1i/Ni)*(n2i/Ni)*(xt*(yt/Ni)))
         COp <- pchisq(CO, df=1, lower.tail=FALSE)
         MH  <- (abs(sum(ai - (n1i/Ni)*xt)) - ifelse(correct, 0.5, 0))^2 / sum((n1i/Ni)*(n2i/Ni)*(xt*(yt/(Ni-1))))
         MHp <- pchisq(MH, df=1, lower.tail=FALSE)
      }

      ### Breslow-Day and Tarone's Test for Heterogeneity

      if (is.na(beta)) {
         BD    <- NA_real_
         TA    <- NA_real_
         BDp   <- NA_real_
         TAp   <- NA_real_
         k.pos <- 0L
      } else {
         if (identical(beta.exp,1) || identical(beta.exp,1L)) {
            N11 <- (n1i/Ni)*xt
         } else {
            A   <- beta.exp * (n1i + xt) + (n2i - xt)
            B   <- sqrt(A^2 - 4*n1i*xt*beta.exp*(beta.exp-1))
            N11 <- (A-B) / (2*(beta.exp-1))
         }
         pos   <- (N11 > 0) & (xt > 0) & (yt > 0)
         k.pos <- sum(pos)
         N11   <- N11[pos]
         N12   <- n1i[pos] - N11
         N21   <- xt[pos]  - N11
         N22   <- N11 - n1i[pos] - xt[pos] + Ni[pos]
         BD    <- max(0, sum((ai[pos]-N11)^2 / (1/N11 + 1/N12 + 1/N21 + 1/N22)^(-1)))
         TA    <- max(0, BD - sum(ai[pos]-N11)^2 / sum((1/N11 + 1/N12 + 1/N21 + 1/N22)^(-1)))
         if (k.pos > 1L) {
            BDp <- pchisq(BD, df=k.pos-1L, lower.tail=FALSE)
            TAp <- pchisq(TA, df=k.pos-1L, lower.tail=FALSE)
         } else {
            BDp <- NA_real_
            TAp <- NA_real_
         }
      }

   }

   if (measure == "RR") {

      R <- sum(ai * (n2i/Ni))
      S <- sum(ci * (n1i/Ni))

      if (identical(sum(ai),0) || identical(sum(ci),0) || identical(sum(ai),0L) || identical(sum(ci),0L)) {
         beta.exp <- NA_real_
         beta     <- NA_real_
         se       <- NA_real_
         zval     <- NA_real_
         pval     <- NA_real_
         ci.lb    <- NA_real_
         ci.ub    <- NA_real_
      } else {
         beta.exp <- R/S
         beta     <- log(beta.exp)
         se       <- sqrt(sum(((n1i/Ni)*(n2i/Ni)*(ai+ci) - (ai/Ni)*ci)) / (R*S))
         zval     <- beta / se
         pval     <- 2*pnorm(abs(zval), lower.tail=FALSE)
         ci.lb    <- beta - qnorm(level/2, lower.tail=FALSE) * se
         ci.ub    <- beta + qnorm(level/2, lower.tail=FALSE) * se
      }

      names(beta) <- "intrcpt"
      vb <- matrix(se^2, dimnames=list("intrcpt", "intrcpt"))

   }

   if (measure == "RD") {

      beta  <- sum(ai*(n2i/Ni) - ci*(n1i/Ni)) / sum(n1i*(n2i/Ni))
      se    <- sqrt((beta * (sum(ci*(n1i/Ni)^2 - ai*(n2i/Ni)^2 + (n1i/Ni)*(n2i/Ni)*(n2i-n1i)/2)) + sum(ai*(n2i-ci)/Ni + ci*(n1i-ai)/Ni)/2) / sum(n1i*(n2i/Ni))^2) ### equation in: Sato, Greenland, & Robins (1989)
      #se   <- sqrt(sum(((ai/Ni^2)*bi*(n2i^2/n1i) + (ci/Ni^2)*di*(n1i^2/n2i))) / sum(n1i*(n2i/Ni))^2) ### equation in: Greenland & Robins (1985)
      zval  <- beta / se
      pval  <- 2*pnorm(abs(zval), lower.tail=FALSE)
      ci.lb <- beta - qnorm(level/2, lower.tail=FALSE) * se
      ci.ub <- beta + qnorm(level/2, lower.tail=FALSE) * se

      names(beta) <- "intrcpt"
      vb <- matrix(se^2, dimnames=list("intrcpt", "intrcpt"))

   }

   if (measure == "IRR") {

      R <- sum(x1i * (t2i/Ti))
      S <- sum(x2i * (t1i/Ti))

      if (identical(sum(x1i),0) || identical(sum(x2i),0) || identical(sum(x1i),0L) || identical(sum(x2i),0L)) {
         beta.exp <- NA_real_
         beta     <- NA_real_
         se       <- NA_real_
         zval     <- NA_real_
         pval     <- NA_real_
         ci.lb    <- NA_real_
         ci.ub    <- NA_real_
      } else {
         beta.exp <- R/S
         beta     <- log(beta.exp)
         se       <- sqrt(sum((t1i/Ti)*(t2i/Ti)*(x1i+x2i)) / (R*S))
         zval     <- beta / se
         pval     <- 2*pnorm(abs(zval), lower.tail=FALSE)
         ci.lb    <- beta - qnorm(level/2, lower.tail=FALSE) * se
         ci.ub    <- beta + qnorm(level/2, lower.tail=FALSE) * se
      }

      names(beta) <- "intrcpt"
      vb <- matrix(se^2, dimnames=list("intrcpt", "intrcpt"))

      ### Mantel-Haenszel Statistic

      xt <- x1i + x2i
      if (identical(sum(xt),0) || identical(sum(xt),0L)) {
         MH  <- NA_real_
         MHp <- NA_real_
      } else {
         MH  <- (abs(sum(x1i - xt*(t1i/Ti))) - ifelse(correct, 0.5, 0))^2 / sum(xt*(t1i/Ti)*(t2i/Ti))
         MHp <- pchisq(MH, df=1, lower.tail=FALSE)
      }

   }

   if (measure == "IRD") {

      beta  <- sum((x1i*t2i - x2i*t1i)/Ti) / sum((t1i/Ti)*t2i)
      se    <- sqrt(sum(((t1i/Ti)*t2i)^2*(x1i/t1i^2+x2i/t2i^2))) / sum((t1i/Ti)*t2i) ### from Rothland et al. (2008), chapter 15
      zval  <- beta / se
      pval  <- 2*pnorm(abs(zval), lower.tail=FALSE)
      ci.lb <- beta - qnorm(level/2, lower.tail=FALSE) * se
      ci.ub <- beta + qnorm(level/2, lower.tail=FALSE) * se

      names(beta) <- "intrcpt"
      vb <- matrix(se^2, dimnames=list("intrcpt", "intrcpt"))

   }

   #########################################################################

   ### heterogeneity test (inverse variance method)

   if (verbose)
      message(mstyle$message("Heterogeneity testing ..."))

   wi <- 1/vi

   if (k.yi > 1) {
      QE  <- max(0, sum(wi*(yi-beta)^2))
      QEp <- pchisq(QE, df=k.yi-1, lower.tail=FALSE)
      I2  <- max(0, 100 * (QE - (k.yi-1)) / QE)
      H2  <- QE / (k.yi-1)
   } else {
      QE  <- 0
      QEp <- 1
      I2  <- 0
      H2  <- 1
   }

   #########################################################################

   ###### fit statistics

   if (verbose)
      message(mstyle$message("Computing fit statistics and log likelihood ..."))

   if (k.yi >= 1) {

      ll.ML     <- -1/2 * (k.yi)   * log(2*base::pi)                   - 1/2 * sum(log(vi))                      - 1/2 * QE
      ll.REML   <- -1/2 * (k.yi-1) * log(2*base::pi) + 1/2 * log(k.yi) - 1/2 * sum(log(vi)) - 1/2 * log(sum(wi)) - 1/2 * QE
      if (any(vi <= 0)) {
         dev.ML <- -2 * ll.ML
      } else {
         dev.ML <- -2 * (ll.ML - sum(dnorm(yi, mean=yi, sd=sqrt(vi), log=TRUE)))
      }
      AIC.ML    <- -2 * ll.ML   + 2
      BIC.ML    <- -2 * ll.ML   + log(k.yi)
      AICc.ML   <- -2 * ll.ML   + 2 * max(k.yi, 3) / (max(k.yi, 3) - 2)
      dev.REML  <- -2 * (ll.REML - 0)
      AIC.REML  <- -2 * ll.REML + 2
      BIC.REML  <- -2 * ll.REML + log(k.yi-1)
      AICc.REML <- -2 * ll.REML + 2 * max(k.yi-1, 3) / (max(k.yi-1, 3) - 2)

      fit.stats <- matrix(c(ll.ML, dev.ML, AIC.ML, BIC.ML, AICc.ML, ll.REML, dev.REML, AIC.REML, BIC.REML, AICc.REML), ncol=2, byrow=FALSE)

   } else {
      fit.stats <- matrix(NA_real_, nrow=5, ncol=2, byrow=FALSE)
   }

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
   coef.na   <- c(X=FALSE)

   method    <- "FE"
   weighted  <- TRUE
   test      <- "z"
   ddf       <- NA_integer_

   if (is.null(ddd$outlist) || ddd$outlist == "nodata") {

      outdat <- list(ai=ai, bi=bi, ci=ci, di=di, x1i=x1i, x2i=x2i, t1i=t1i, t2i=t2i)

      res <- list(b=beta, beta=beta, se=se, zval=zval, pval=pval, ci.lb=ci.lb, ci.ub=ci.ub, vb=vb,
                  tau2=tau2, tau2.f=tau2,
                  I2=I2, H2=H2,
                  QE=QE, QEp=QEp, CO=CO, COp=COp, MH=MH, MHp=MHp, BD=BD, BDp=BDp, TA=TA, TAp=TAp,
                  k=k, k.f=k.f, k.yi=k.yi, k.pos=k.pos, k.eff=k.eff, k.all=k.all, p=p, p.eff=p.eff, parms=parms,
                  int.only=int.only, intercept=intercept, coef.na=coef.na,
                  yi=yi, vi=vi, yi.f=yi.f, vi.f=vi.f, X.f=X.f,
                  outdat.f=outdat.f, outdat=outdat, ni=ni, ni.f=ni.f,
                  ids=ids, not.na=not.na, subset=subset, not.na.yivi=not.na.yivi, slab=slab, slab.null=slab.null,
                  measure=measure, method=method, weighted=weighted,
                  test=test, ddf=ddf, dfs=ddf, btt=btt, m=m,
                  digits=digits, level=level,
                  add=add, to=to, drop00=drop00, correct=correct,
                  fit.stats=fit.stats,
                  formula.yi=NULL, formula.mods=NULL, version=packageVersion("metafor"), call=mf)

      if (is.null(ddd$outlist))
         res <- append(res, list(data=data), which(names(res) == "fit.stats"))

   } else {

      if (ddd$outlist == "minimal") {
         res <- list(b=beta, beta=beta, se=se, zval=zval, pval=pval, ci.lb=ci.lb, ci.ub=ci.ub, vb=vb,
                     tau2=tau2,
                     I2=I2, H2=H2,
                     QE=QE, QEp=QEp, CO=CO, COp=COp, MH=MH, MHp=MHp, BD=BD, BDp=BDp, TA=TA, TAp=TAp,
                     k=k, k.yi=k.yi, k.pos=k.pos, k.eff=k.eff, p=p, p.eff=p.eff, parms=parms,
                     int.only=int.only,
                     measure=measure, method=method,
                     test=test, ddf=ddf, dfs=ddf, btt=btt, m=m,
                     digits=digits,
                     fit.stats=fit.stats)
      } else {
         res <- eval(str2lang(paste0("list(", ddd$outlist, ")")))
      }

   }

   time.end <- proc.time()
   res$time <- unname(time.end - time.start)[3]

   if (.isTRUE(ddd$time))
      .print.time(res$time)

   if (verbose || .isTRUE(ddd$time))
      cat("\n")

   class(res) <- c("rma.mh", "rma")
   return(res)

}
