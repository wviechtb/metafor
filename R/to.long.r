to.long <- function(measure, ai, bi, ci, di, n1i, n2i, x1i, x2i, t1i, t2i, m1i, m2i, sd1i, sd2i, xi, mi, ri, ti, sdi, ni,
data, slab, subset, add=1/2, to="none", drop00=FALSE, vlong=FALSE, append=TRUE, var.names) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   ### check argument specifications

   if (!is.element(measure, c("RR","OR","PETO","RD","AS","PHI","YUQ","YUY","RTET", ### 2x2 table measures
                              "PBIT","OR2D","OR2DN","OR2DL",                       ### - transformations to SMD
                              "MPRD","MPRR","MPOR","MPORC","MPPETO",               ### - measures for matched pairs data
                              "IRR","IRD","IRSD",                                  ### two-group person-time data measures
                              "MD","SMD","SMDH","ROM",                             ### two-group mean/SD measures
                              "CVR","VR",                                          ### coefficient of variation ratio, variability ratio
                              "RPB","RBIS","D2OR","D2ORN","D2ORL",                 ### - transformations to r_PB, r_BIS, and log(OR)
                              "COR","UCOR","ZCOR",                                 ### correlations (raw and r-to-z transformed)
                              "PCOR","ZPCOR","SPCOR",                              ### partial and semi-partial correlations
                              "PR","PLN","PLO","PAS","PFT",                        ### single proportions (and transformations thereof)
                              "IR","IRLN","IRS","IRFT",                            ### single-group person-time data (and transformations thereof)
                              "MN","MNLN","CVLN","SDLN",                           ### mean, log(mean), log(CV), log(SD)
                              "MC","SMCC","SMCR","SMCRH","ROMC","CVRC","VRC",      ### raw/standardized mean change, log(ROM), CVR, and VR for dependent samples
                              "ARAW","AHW","ABT")))                                ### alpha (and transformations thereof)
      stop(mstyle$stop("Unknown 'measure' specified."))

   if (is.element(measure, c("CVR","VR","PCOR","ZPCOR","SPCOR","CVLN","SDLN","VRC")))
      stop(mstyle$stop("Function not available for this outcome measure."))

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   if (!is.element(to, c("all","only0","if0all","none")))
      stop(mstyle$stop("Unknown 'to' argument specified."))

   ### check if data argument has been specified

   if (missing(data))
      data <- NULL

   ### need this at the end to check if append=TRUE can actually be done

   has.data <- !is.null(data)

   if (is.null(data)) {
      data <- sys.frame(sys.parent())
   } else {
      if (!is.data.frame(data))
         data <- data.frame(data)
   }

   mf <- match.call()

   ### get slab and subset arguments (will be NULL when unspecified)

   mf.slab   <- mf[[match("slab",   names(mf))]]
   mf.subset <- mf[[match("subset", names(mf))]]
   slab   <- eval(mf.slab,   data, enclos=sys.frame(sys.parent()))
   subset <- eval(mf.subset, data, enclos=sys.frame(sys.parent()))

   #########################################################################
   #########################################################################
   #########################################################################

   if (is.element(measure, c("RR","OR","RD","AS","PETO","PHI","YUQ","YUY","RTET","PBIT","OR2D","OR2DN","OR2DL","MPRD","MPRR","MPOR","MPORC","MPPETO"))) {

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

      k <- length(ai) ### number of outcomes before subsetting

      if (!is.null(subset)) {
         ai <- ai[subset]
         bi <- bi[subset]
         ci <- ci[subset]
         di <- di[subset]
      }

      if (length(ai)==0L || length(bi)==0L || length(ci)==0L || length(di)==0L)
         stop(mstyle$stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments."))

      if (!all(length(ai) == c(length(ai),length(bi),length(ci),length(di))))
         stop(mstyle$stop("Supplied data vectors are not all of the same length."))

      if (any(c(ai, bi, ci, di) < 0, na.rm=TRUE))
         stop(mstyle$stop("One or more counts are negative."))

      ni.u <- ai + bi + ci + di ### unadjusted total sample sizes

      ### if drop00=TRUE, set counts to NA for studies that have no events (or all events) in both arms

      if (drop00) {
         id00 <- c(ai == 0L & ci == 0L) | c(bi == 0L & di == 0L)
         id00[is.na(id00)] <- FALSE
         ai[id00] <- NA
         bi[id00] <- NA
         ci[id00] <- NA
         di[id00] <- NA
      }

      if (to == "all") {

         ### always add to all cells in all studies

         ai  <- ai + add
         ci  <- ci + add
         bi  <- bi + add
         di  <- di + add

      }

      if (to == "only0") {

         ### add to cells in studies with at least one 0 entry

         id0 <- c(ai == 0L | ci == 0L | bi == 0L | di == 0L)
         id0[is.na(id0)] <- FALSE

         ai[id0] <- ai[id0] + add
         ci[id0] <- ci[id0] + add
         bi[id0] <- bi[id0] + add
         di[id0] <- di[id0] + add

      }

      if (to == "if0all") {

         ### add to cells in all studies if there is at least one 0 entry

         id0 <- c(ai == 0L | ci == 0L | bi == 0L | di == 0L)
         id0[is.na(id0)] <- FALSE

         if (any(id0)) {

            ai <- ai + add
            ci <- ci + add
            bi <- bi + add
            di <- di + add

         }

      }

   }

   #########################################################################

   if (is.element(measure, c("IRR","IRD","IRSD"))) {

      mf.x1i <- mf[[match("x1i", names(mf))]]
      mf.x2i <- mf[[match("x2i", names(mf))]]
      mf.t1i <- mf[[match("t1i", names(mf))]]
      mf.t2i <- mf[[match("t2i", names(mf))]]
      x1i <- eval(mf.x1i, data, enclos=sys.frame(sys.parent()))
      x2i <- eval(mf.x2i, data, enclos=sys.frame(sys.parent()))
      t1i <- eval(mf.t1i, data, enclos=sys.frame(sys.parent()))
      t2i <- eval(mf.t2i, data, enclos=sys.frame(sys.parent()))

      k <- length(x1i) ### number of outcomes before subsetting

      if (!is.null(subset)) {
         x1i <- x1i[subset]
         x2i <- x2i[subset]
         t1i <- t1i[subset]
         t2i <- t2i[subset]
      }

      if (length(x1i)==0L || length(x2i)==0L || length(t1i)==0L || length(t2i)==0L)
         stop(mstyle$stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments."))

      if (!all(length(x1i) == c(length(x1i),length(x2i),length(t1i),length(t2i))))
         stop(mstyle$stop("Supplied data vectors are not all of the same length."))

      if (any(c(x1i, x2i) < 0, na.rm=TRUE))
         stop(mstyle$stop("One or more counts are negative."))

      if (any(c(t1i, t2i) < 0, na.rm=TRUE))
         stop(mstyle$stop("One or more person-times are negative."))

      ni.u <- t1i + t2i ### unadjusted total sample sizes

      ### if drop00=TRUE, set counts to NA for studies that have no events in both arms

      if (drop00) {
         id00 <- c(x1i == 0L & x2i == 0L)
         id00[is.na(id00)] <- FALSE
         x1i[id00] <- NA
         x2i[id00] <- NA
      }

      if (to == "all") {

         ### always add to all cells in all studies

         x1i <- x1i + add
         x2i <- x2i + add

      }

      if (to == "only0") {

         ### add to cells in studies with at least one 0 entry

         id0 <- c(x1i == 0L | x2i == 0L)
         id0[is.na(id0)] <- FALSE

         x1i[id0] <- x1i[id0] + add
         x2i[id0] <- x2i[id0] + add

      }

      if (to == "if0all") {

         ### add to cells in all studies if there is at least one 0 entry

         id0 <- c(x1i == 0L | x2i == 0L)
         id0[is.na(id0)] <- FALSE

         if (any(id0)) {

            x1i <- x1i + add
            x2i <- x2i + add

         }

      }

   }

   #########################################################################

   if (is.element(measure, c("MD","SMD","SMDH","ROM","RPB","RBIS","D2OR","D2ORN","D2ORL"))) {

      mf.m1i  <- mf[[match("m1i",  names(mf))]]
      mf.m2i  <- mf[[match("m2i",  names(mf))]]
      mf.sd1i <- mf[[match("sd1i", names(mf))]]
      mf.sd2i <- mf[[match("sd2i", names(mf))]]
      mf.n1i  <- mf[[match("n1i",  names(mf))]]
      mf.n2i  <- mf[[match("n2i",  names(mf))]]
      m1i  <- eval(mf.m1i,  data, enclos=sys.frame(sys.parent()))
      m2i  <- eval(mf.m2i,  data, enclos=sys.frame(sys.parent()))
      sd1i <- eval(mf.sd1i, data, enclos=sys.frame(sys.parent()))
      sd2i <- eval(mf.sd2i, data, enclos=sys.frame(sys.parent()))
      n1i  <- eval(mf.n1i,  data, enclos=sys.frame(sys.parent()))
      n2i  <- eval(mf.n2i,  data, enclos=sys.frame(sys.parent()))

      k <- length(n1i) ### number of outcomes before subsetting

      if (!is.null(subset)) {
         m1i  <- m1i[subset]
         m2i  <- m2i[subset]
         sd1i <- sd1i[subset]
         sd2i <- sd2i[subset]
         n1i  <- n1i[subset]
         n2i  <- n2i[subset]
      }

      if (length(m1i)==0L || length(m2i)==0L || length(sd1i)==0L || length(sd2i)==0L || length(n1i)==0L || length(n2i)==0L)
         stop(mstyle$stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments."))

      if (!all(length(m1i) == c(length(m1i),length(m2i),length(sd1i),length(sd2i),length(n1i),length(n2i))))
         stop(mstyle$stop("Supplied data vectors are not all of the same length."))

      if (any(c(sd1i, sd2i) < 0, na.rm=TRUE))
         stop(mstyle$stop("One or more standard deviations are negative."))

      if (any(c(n1i, n2i) < 0, na.rm=TRUE))
         stop(mstyle$stop("One or more sample sizes are negative."))

      ni.u <- n1i + n2i ### unadjusted total sample sizes

   }

   #########################################################################

   if (is.element(measure, c("COR","UCOR","ZCOR"))) {

      mf.ri <- mf[[match("ri", names(mf))]]
      mf.ni <- mf[[match("ni", names(mf))]]
      ri <- eval(mf.ri, data, enclos=sys.frame(sys.parent()))
      ni <- eval(mf.ni, data, enclos=sys.frame(sys.parent()))

      k <- length(ri) ### number of outcomes before subsetting

      if (!is.null(subset)) {
         ri <- ri[subset]
         ni <- ni[subset]
      }

      if (length(ri)==0L || length(ni)==0L)
         stop(mstyle$stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments."))

      if (length(ri) != length(ni))
         stop(mstyle$stop("Supplied data vectors are not of the same length."))

      if (any(abs(ri) > 1, na.rm=TRUE))
         stop(mstyle$stop("One or more correlations are > 1 or < -1."))

      if (any(ni < 0, na.rm=TRUE))
         stop(mstyle$stop("One or more sample sizes are negative."))

      ni.u <- ni ### unadjusted total sample sizes

   }

   #########################################################################

   if (is.element(measure, c("PR","PLN","PLO","PAS","PFT"))) {

      mf.xi <- mf[[match("xi", names(mf))]]
      mf.mi <- mf[[match("mi", names(mf))]]
      mf.ni <- mf[[match("ni", names(mf))]]
      xi <- eval(mf.xi, data, enclos=sys.frame(sys.parent()))
      mi <- eval(mf.mi, data, enclos=sys.frame(sys.parent()))
      ni <- eval(mf.ni, data, enclos=sys.frame(sys.parent()))
      if (is.null(mi)) mi <- ni - xi

      k <- length(xi) ### number of outcomes before subsetting

      if (!is.null(subset)) {
         xi <- xi[subset]
         mi <- mi[subset]
      }

      if (length(xi)==0L || length(mi)==0L)
         stop(mstyle$stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments."))

      if (length(xi) != length(mi))
         stop(mstyle$stop("Supplied data vectors are not all of the same length."))

      if (any(c(xi, mi) < 0, na.rm=TRUE))
         stop(mstyle$stop("One or more counts are negative."))

      ni.u <- xi + mi ### unadjusted total sample sizes

      if (to == "all") {

         ### always add to all cells in all studies

         xi <- xi + add
         mi <- mi + add

      }

      if (to == "only0") {

         ### add to cells in studies with at least one 0 entry

         id0 <- c(xi == 0L | mi == 0L)
         id0[is.na(id0)] <- FALSE

         xi[id0] <- xi[id0] + add
         mi[id0] <- mi[id0] + add

      }

      if (to == "if0all") {

         ### add to cells in all studies if there is at least one 0 entry

         id0 <- c(xi == 0L | mi == 0L)
         id0[is.na(id0)] <- FALSE

         if (any(id0)) {

            xi <- xi + add
            mi <- mi + add

         }

      }

   }

   #########################################################################

   if (is.element(measure, c("IR","IRLN","IRS","IRFT"))) {

      mf.xi <- mf[[match("xi", names(mf))]]
      mf.ti <- mf[[match("ti", names(mf))]]
      xi <- eval(mf.xi, data, enclos=sys.frame(sys.parent()))
      ti <- eval(mf.ti, data, enclos=sys.frame(sys.parent()))

      k <- length(xi) ### number of outcomes before subsetting

      if (!is.null(subset)) {
         xi <- xi[subset]
         ti <- ti[subset]
      }

      if (length(xi)==0L || length(ti)==0L)
         stop(mstyle$stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments."))

      if (length(xi) != length(ti))
         stop(mstyle$stop("Supplied data vectors are not all of the same length."))

      if (any(xi < 0, na.rm=TRUE))
         stop(mstyle$stop("One or more counts are negative."))

      if (any(ti < 0, na.rm=TRUE))
         stop(mstyle$stop("One or more person-times are negative."))

      ni.u <- ti ### unadjusted total sample sizes

      if (to == "all") {

         ### always add to all cells in all studies

         xi <- xi + add

      }

      if (to == "only0") {

         ### add to cells in studies with at least one 0 entry

         id0 <- c(xi == 0L)
         id0[is.na(id0)] <- FALSE

         xi[id0] <- xi[id0] + add

      }

      if (to == "if0all") {

         ### add to cells in all studies if there is at least one 0 entry

         id0 <- c(xi == 0L)
         id0[is.na(id0)] <- FALSE

         if (any(id0)) {

            xi <- xi + add

         }

      }

   }

   #########################################################################

   if (is.element(measure, c("MN","MNLN"))) {

      mf.mi  <- mf[[match("mi",  names(mf))]]
      mf.sdi <- mf[[match("sdi", names(mf))]]
      mf.ni  <- mf[[match("ni",  names(mf))]]
      mi  <- eval(mf.mi,  data, enclos=sys.frame(sys.parent()))
      sdi <- eval(mf.sdi, data, enclos=sys.frame(sys.parent()))
      ni  <- eval(mf.ni,  data, enclos=sys.frame(sys.parent()))

      k <- length(ni) ### number of outcomes before subsetting

      if (!is.null(subset)) {
         mi  <- mi[subset]
         sdi <- sdi[subset]
         ni  <- ni[subset]
      }

      if (length(mi)==0L || length(sdi)==0L || length(ni)==0L)
         stop(mstyle$stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments."))

      if (!all(length(mi) == c(length(mi),length(sdi),length(ni))))
         stop(mstyle$stop("Supplied data vectors are not all of the same length."))

      if (any(sdi < 0, na.rm=TRUE))
         stop(mstyle$stop("One or more standard deviations are negative."))

      if (any(ni < 0, na.rm=TRUE))
         stop(mstyle$stop("One or more sample sizes are negative."))

      if (is.element(measure, c("MNLN","CVLN")) && any(mi < 0, na.rm=TRUE))
         stop(mstyle$stop("One or more means are negative."))

      ni.u <- ni ### unadjusted total sample sizes

  }

   #########################################################################

   if (is.element(measure, c("MC","SMCC","SMCR","SMCRH","ROMC","CVRC"))) {

      mf.m1i  <- mf[[match("m1i",  names(mf))]]
      mf.m2i  <- mf[[match("m2i",  names(mf))]]
      mf.sd1i <- mf[[match("sd1i", names(mf))]]
      mf.sd2i <- mf[[match("sd2i", names(mf))]] ### for SMCR, do not need to supply this
      mf.ni   <- mf[[match("ni",   names(mf))]]
      mf.ri   <- mf[[match("ri",   names(mf))]]
      m1i  <- eval(mf.m1i,  data, enclos=sys.frame(sys.parent()))
      m2i  <- eval(mf.m2i,  data, enclos=sys.frame(sys.parent()))
      sd1i <- eval(mf.sd1i, data, enclos=sys.frame(sys.parent()))
      sd2i <- eval(mf.sd2i, data, enclos=sys.frame(sys.parent()))
      ni   <- eval(mf.ni,   data, enclos=sys.frame(sys.parent()))
      ri   <- eval(mf.ri,   data, enclos=sys.frame(sys.parent()))

      k <- length(m1i) ### number of outcomes before subsetting

      if (!is.null(subset)) {
         m1i  <- m1i[subset]
         m2i  <- m2i[subset]
         sd1i <- sd1i[subset]
         sd2i <- sd2i[subset]
         ni   <- ni[subset]
         ri   <- ri[subset]
      }

      if (is.element(measure, c("MC","SMCC","SMCRH","ROMC","CVRC"))) {

         if (length(m1i)==0L || length(m2i)==0L || length(sd1i)==0L || length(sd2i)==0L || length(ni)==0L || length(ri)==0L)
            stop(mstyle$stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments."))

         if (!all(length(m1i) == c(length(m1i),length(m2i),length(sd1i),length(sd2i),length(ni),length(ri))))
            stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         if (any(c(sd1i, sd2i) < 0, na.rm=TRUE))
            stop(mstyle$stop("One or more standard deviations are negative."))

      } else {

         if (length(m1i)==0L || length(m2i)==0L || length(sd1i)==0L || length(ni)==0L || length(ri)==0L)
            stop(mstyle$stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments."))

         if (!all(length(m1i) == c(length(m1i),length(m2i),length(sd1i),length(ni),length(ri))))
            stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         if (any(sd1i < 0, na.rm=TRUE))
            stop(mstyle$stop("One or more standard deviations are negative."))

      }

      if (any(abs(ri) > 1, na.rm=TRUE))
         stop(mstyle$stop("One or more correlations are > 1 or < -1."))

      if (any(ni < 0, na.rm=TRUE))
         stop(mstyle$stop("One or more sample sizes are negative."))

      ni.u <- ni ### unadjusted total sample sizes

   }

   #########################################################################

   if (is.element(measure, c("ARAW","AHW","ABT"))) {

      mf.ai <- mf[[match("ai", names(mf))]]
      mf.mi <- mf[[match("mi", names(mf))]]
      mf.ni <- mf[[match("ni", names(mf))]]
      ai <- eval(mf.ai, data, enclos=sys.frame(sys.parent()))
      mi <- eval(mf.mi, data, enclos=sys.frame(sys.parent()))
      ni <- eval(mf.ni, data, enclos=sys.frame(sys.parent()))

      k <- length(ai) ### number of outcomes before subsetting

      if (!is.null(subset)) {
         ai <- ai[subset]
         mi <- mi[subset]
         ni <- ni[subset]
      }

      if (length(ai)==0L || length(mi)==0L || length(ni)==0L)
         stop(mstyle$stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments."))

      if (!all(length(ai) == c(length(ai),length(mi),length(ni))))
         stop(mstyle$stop("Supplied data vectors are not all of the same length."))

      if (any(ai > 1, na.rm=TRUE))
         stop(mstyle$stop("One or more alpha values are > 1."))

      if (any(mi < 2, na.rm=TRUE))
         stop(mstyle$stop("One or more mi values are < 2."))

      if (any(ni < 0, na.rm=TRUE))
         stop(mstyle$stop("One or more sample sizes are negative."))

      ni.u <- ni ### unadjusted total sample sizes

   }

   #########################################################################
   #########################################################################
   #########################################################################

   ### generate study labels if none are specified

   if (is.null(slab)) {

      slab  <- seq_len(k)

   } else {

      if (anyNA(slab))
         stop(mstyle$stop("NAs in study labels."))

      if (length(slab) != k)
         stop(mstyle$stop("Study labels not of same length as data."))

      if (is.factor(slab))
         slab <- as.character(slab)

   }

   ### if a subset of studies is specified

   if (!is.null(subset)) {
      slab <- slab[subset]
      if (has.data)
         data <- data[subset,]
   }

   ### check if study labels are unique; if not, make them unique

   if (anyDuplicated(slab))
      slab <- .make.unique(slab)

   #########################################################################
   #########################################################################
   #########################################################################

   if (is.element(measure, c("RR","OR","RD","AS","PETO","PHI","YUQ","YUY","RTET","PBIT","OR2D","OR2DN","OR2DL"))) {

      ### check for NAs in table data and act accordingly

      has.na <- is.na(ai) | is.na(bi) | is.na(ci) | is.na(di)

      if (any(has.na)) {

         not.na <- !has.na

         if (na.act == "na.omit") {
            ai   <- ai[not.na]
            bi   <- bi[not.na]
            ci   <- ci[not.na]
            di   <- di[not.na]
            slab <- slab[not.na]
            if (has.data)
               data <- data[not.na,]
            warning(mstyle$warning("Tables with NAs omitted."), call.=FALSE)
         }

         if (na.act == "na.fail")
            stop(mstyle$stop("Missing values in tables."))

      }

      k <- length(ai)

      ### at least one study left?

      if (k < 1)
         stop(mstyle$stop("Processing terminated since k = 0."))

      ### create long format dataset

      if (vlong) {

         ### create very long format dataset

         dat <- data.frame(rep(slab, each=4), stringsAsFactors=FALSE)

         dat[[2]] <- rep(c(1,1,2,2), k)
         dat[[3]] <- rep(c(1,2,1,2), k)
         dat[[4]] <- c(rbind(ai,bi,ci,di))

         if (missing(var.names)) {
            names(dat) <- c("study", "group", "outcome", "freq")
         } else {
            if (length(var.names) != 4L)
               stop(mstyle$stop("Variable names not of length 4."))
            names(dat) <- var.names
         }

         dat[[1]] <- factor(dat[[1]])
         dat[[2]] <- factor(dat[[2]])
         dat[[3]] <- factor(dat[[3]])

         if (has.data && append)
            dat <- cbind(data[rep(seq_len(k), each=4),], dat)

      } else {

         ### create regular long format dataset

         dat <- data.frame(rep(slab, each=2), stringsAsFactors=FALSE)

         dat[[2]] <- rep(c(1,2), k)
         dat[[3]] <- c(rbind(ai,ci))
         dat[[4]] <- c(rbind(bi,di))

         if (missing(var.names)) {
            names(dat) <- c("study", "group", "out1", "out2")
         } else {
            if (length(var.names) != 4L)
               stop(mstyle$stop("Variable names not of length 4."))
            names(dat) <- var.names
         }

         dat[[1]] <- factor(dat[[1]])
         dat[[2]] <- factor(dat[[2]])

         if (has.data && append)
            dat <- cbind(data[rep(seq_len(k), each=2),], dat)

      }

   }

   #########################################################################

   if (is.element(measure, c("MPRD","MPRR","MPOR"))) {

      ### check for NAs in table data and act accordingly

      has.na <- is.na(ai) | is.na(bi) | is.na(ci) | is.na(di)

      if (any(has.na)) {

         not.na <- !has.na

         if (na.act == "na.omit") {
            ai   <- ai[not.na]
            bi   <- bi[not.na]
            ci   <- ci[not.na]
            di   <- di[not.na]
            slab <- slab[not.na]
            if (has.data)
               data <- data[not.na,]
            warning(mstyle$warning("Tables with NAs omitted."), call.=FALSE)
         }

         if (na.act == "na.fail")
            stop(mstyle$stop("Missing values in tables."))

      }

      k <- length(ai)

      ### at least one study left?

      if (k < 1)
         stop(mstyle$stop("Processing terminated since k = 0."))

      ### create long format dataset

      if (vlong) {

         ### create very long format dataset

         dat <- data.frame(rep(slab, each=4), stringsAsFactors=FALSE)

         dat[[2]] <- rep(c(1,1,2,2), k)
         dat[[3]] <- rep(c(1,2,1,2), k)
         dat[[4]] <- c(rbind(ai+bi,ci+di,ai+ci,bi+di))

         if (missing(var.names)) {
            names(dat) <- c("study", "time", "outcome", "freq")
         } else {
            if (length(var.names) != 4L)
               stop(mstyle$stop("Variable names not of length 4."))
            names(dat) <- var.names
         }

         dat[[1]] <- factor(dat[[1]])
         dat[[2]] <- factor(dat[[2]])
         dat[[3]] <- factor(dat[[3]])

         if (has.data && append)
            dat <- data.frame(data[rep(seq_len(k), each=4),], dat)

      } else {

         ### create regular long format dataset

         dat <- data.frame(rep(slab, each=2), stringsAsFactors=FALSE)

         dat[[2]] <- rep(c(1,2), k)
         dat[[3]] <- c(rbind(ai+bi,ai+ci))
         dat[[4]] <- c(rbind(ci+di,bi+di))

         if (missing(var.names)) {
            names(dat) <- c("study", "time", "out1", "out2")
         } else {
            if (length(var.names) != 4L)
               stop(mstyle$stop("Variable names not of length 4."))
            names(dat) <- var.names
         }

         dat[[1]] <- factor(dat[[1]])
         dat[[2]] <- factor(dat[[2]])

         if (has.data && append)
            dat <- cbind(data[rep(seq_len(k), each=2),], dat)

      }

   }

   #########################################################################

   if (is.element(measure, c("MPORC","MPPETO"))) {

      ### check for NAs in table data and act accordingly

      has.na <- is.na(ai) | is.na(bi) | is.na(ci) | is.na(di)

      if (any(has.na)) {

         not.na <- !has.na

         if (na.act == "na.omit") {
            ai   <- ai[not.na]
            bi   <- bi[not.na]
            ci   <- ci[not.na]
            di   <- di[not.na]
            slab <- slab[not.na]
            if (has.data)
               data <- data[not.na,]
            warning(mstyle$warning("Tables with NAs omitted."), call.=FALSE)
         }

         if (na.act == "na.fail")
            stop(mstyle$stop("Missing values in tables."))

      }

      k <- length(ai)

      ### at least one study left?

      if (k < 1)
         stop(mstyle$stop("Processing terminated since k = 0."))

      ### create long format dataset

      if (vlong) {

         ### create very long format dataset

         dat <- data.frame(rep(slab, each=4), stringsAsFactors=FALSE)

         dat[[2]] <- rep(c(1,1,2,2), k)
         dat[[3]] <- rep(c(1,2,1,2), k)
         dat[[4]] <- c(rbind(ai,bi,ci,di))

         if (missing(var.names)) {
            names(dat) <- c("study", "out.time1", "out.time2", "freq")
         } else {
            if (length(var.names) != 4L)
               stop(mstyle$stop("Variable names not of length 4."))
            names(dat) <- var.names
         }

         dat[[1]] <- factor(dat[[1]])
         dat[[2]] <- factor(dat[[2]])
         dat[[3]] <- factor(dat[[3]])

         if (has.data && append)
            dat <- cbind(data[rep(seq_len(k), each=4),], dat)

      } else {

         ### create regular long format dataset

         dat <- data.frame(rep(slab, each=2), stringsAsFactors=FALSE)

         dat[[2]] <- rep(c(1,2), k)
         dat[[3]] <- c(rbind(ai,ci))
         dat[[4]] <- c(rbind(bi,di))

         if (missing(var.names)) {
            names(dat) <- c("study", "out.time1", "out1.time2", "out2.time2")
         } else {
            if (length(var.names) != 4L)
               stop(mstyle$stop("Variable names not of length 4."))
            names(dat) <- var.names
         }

         dat[[1]] <- factor(dat[[1]])
         dat[[2]] <- factor(dat[[2]])

         if (has.data && append)
            dat <- cbind(data[rep(seq_len(k), each=2),], dat)

      }

   }

   #########################################################################

   if (is.element(measure, c("IRR","IRD","IRSD"))) {

      ### check for NAs in table data and act accordingly

      has.na <- is.na(x1i) | is.na(x2i) | is.na(t1i) | is.na(t2i)

      if (any(has.na)) {

         not.na <- !has.na

         if (na.act == "na.omit") {
            x1i  <- x1i[not.na]
            x2i  <- x2i[not.na]
            t1i  <- t1i[not.na]
            t2i  <- t2i[not.na]
            slab <- slab[not.na]
            if (has.data)
               data <- data[not.na,]
            warning(mstyle$warning("Tables with NAs omitted."), call.=FALSE)
         }

         if (na.act == "na.fail")
            stop(mstyle$stop("Missing values in tables."))

      }

      k <- length(x1i)

      ### at least one study left?

      if (k < 1)
         stop(mstyle$stop("Processing terminated since k = 0."))

      ### create long format dataset

      dat <- data.frame(rep(slab, each=2), stringsAsFactors=FALSE)

      dat[[2]] <- rep(c(1,2), k)
      dat[[3]] <- c(rbind(x1i,x2i))
      dat[[4]] <- c(rbind(t1i,t2i))

      if (missing(var.names)) {
         names(dat) <- c("study", "group", "events", "ptime")
      } else {
         if (length(var.names) != 4L)
            stop(mstyle$stop("Variable names not of length 4."))
         names(dat) <- var.names
      }

      dat[[1]] <- factor(dat[[1]])
      dat[[2]] <- factor(dat[[2]])

      if (has.data && append)
         dat <- cbind(data[rep(seq_len(k), each=2),], dat)

   }

   #########################################################################

   if (is.element(measure, c("MD","SMD","SMDH","ROM","RPB","RBIS","D2OR","D2ORN","D2ORL"))) {

      ### check for NAs in table data and act accordingly

      has.na <- is.na(m1i) | is.na(m2i) | is.na(sd1i) | is.na(sd2i) | is.na(n1i) | is.na(n2i)

      if (any(has.na)) {

         not.na <- !has.na

         if (na.act == "na.omit") {
            m1i  <- m1i[not.na]
            m2i  <- m2i[not.na]
            sd1i <- sd1i[not.na]
            sd2i <- sd2i[not.na]
            n1i  <- n1i[not.na]
            n2i  <- n2i[not.na]
            slab <- slab[not.na]
            if (has.data)
               data <- data[not.na,]
            warning(mstyle$warning("Tables with NAs omitted."), call.=FALSE)
         }

         if (na.act == "na.fail")
            stop(mstyle$stop("Missing values in tables."))

      }

      k <- length(m1i)

      ### at least one study left?

      if (k < 1)
         stop(mstyle$stop("Processing terminated since k = 0."))

      ### create long format dataset

      dat <- data.frame(rep(slab, each=2), stringsAsFactors=FALSE)

      dat[[2]] <- rep(c(1,2), k)
      dat[[3]] <- c(rbind(m1i,m2i))
      dat[[4]] <- c(rbind(sd1i,sd2i))
      dat[[5]] <- c(rbind(n1i,n2i))

      if (missing(var.names)) {
         names(dat) <- c("study", "group", "mean", "sd", "n")
      } else {
         if (length(var.names) != 5L)
            stop(mstyle$stop("Variable names not of length 5."))
         names(dat) <- var.names
      }

      dat[[1]] <- factor(dat[[1]])
      dat[[2]] <- factor(dat[[2]])

      if (has.data && append)
         dat <- cbind(data[rep(seq_len(k), each=2),], dat)

   }

   #########################################################################

   if (is.element(measure, c("COR","UCOR","ZCOR"))) {

      ### check for NAs in table data and act accordingly

      has.na <- is.na(ri) | is.na(ni)

      if (any(has.na)) {

         not.na <- !has.na

         if (na.act == "na.omit") {
            ri   <- ri[not.na]
            ni   <- ni[not.na]
            slab <- slab[not.na]
            if (has.data)
               data <- data[not.na,]
            warning(mstyle$warning("Tables with NAs omitted."), call.=FALSE)
         }

         if (na.act == "na.fail")
            stop(mstyle$stop("Missing values in tables."))

      }

      k <- length(ri)

      ### at least one study left?

      if (k < 1)
         stop(mstyle$stop("Processing terminated since k = 0."))

      ### create long format dataset

      dat <- data.frame(slab, stringsAsFactors=FALSE)

      dat[[2]] <- ri
      dat[[3]] <- ni

      if (missing(var.names)) {
         names(dat) <- c("study", "r", "n")
      } else {
         if (length(var.names) != 3L)
            stop(mstyle$stop("Variable names not of length 3."))
         names(dat) <- var.names
      }

      dat[[1]] <- factor(dat[[1]])

      if (has.data && append)
         dat <- cbind(data, dat)

   }

   #########################################################################

   if (is.element(measure, c("PR","PLN","PLO","PAS","PFT"))) {

      ### check for NAs in table data and act accordingly

      has.na <- is.na(xi) | is.na(mi)

      if (any(has.na)) {

         not.na <- !has.na

         if (na.act == "na.omit") {
            xi   <- xi[not.na]
            mi   <- mi[not.na]
            slab <- slab[not.na]
            if (has.data)
               data <- data[not.na,]
            warning(mstyle$warning("Tables with NAs omitted."), call.=FALSE)
         }

         if (na.act == "na.fail")
            stop(mstyle$stop("Missing values in tables."))

      }

      k <- length(xi)

      ### at least one study left?

      if (k < 1)
         stop(mstyle$stop("Processing terminated since k = 0."))

      ### create long format dataset

      if (vlong) {

         ### create very long format dataset

         dat <- data.frame(rep(slab, each=2), stringsAsFactors=FALSE)

         dat[[2]] <- rep(c(1,2), k)
         dat[[3]] <- c(rbind(xi,mi))

         if (missing(var.names)) {
            names(dat) <- c("study", "outcome", "freq")
         } else {
            if (length(var.names) != 3L)
               stop(mstyle$stop("Variable names not of length 3."))
            names(dat) <- var.names
         }

         dat[[1]] <- factor(dat[[1]])
         dat[[2]] <- factor(dat[[2]])

         if (has.data && append)
            dat <- cbind(data[rep(seq_len(k), each=2),], dat)

      } else {

         ### create regular long format dataset

         dat <- data.frame(slab, stringsAsFactors=FALSE)

         dat[[2]] <- xi
         dat[[3]] <- mi

         if (missing(var.names)) {
            names(dat) <- c("study", "out1", "out2")
         } else {
            if (length(var.names) != 3L)
               stop(mstyle$stop("Variable names not of length 3."))
            names(dat) <- var.names
         }

         dat[[1]] <- factor(dat[[1]])

         if (has.data && append)
            dat <- cbind(data, dat)

      }

   }

   #########################################################################

   if (is.element(measure, c("IR","IRLN","IRS","IRFT"))) {

      ### check for NAs in table data and act accordingly

      has.na <- is.na(xi) | is.na(ti)

      if (any(has.na)) {

         not.na <- !has.na

         if (na.act == "na.omit") {
            xi   <- xi[not.na]
            ti   <- ti[not.na]
            slab <- slab[not.na]
            if (has.data)
               data <- data[not.na,]
            warning(mstyle$warning("Tables with NAs omitted."), call.=FALSE)
         }

         if (na.act == "na.fail")
            stop(mstyle$stop("Missing values in tables."))

      }

      k <- length(xi)

      ### at least one study left?

      if (k < 1)
         stop(mstyle$stop("Processing terminated since k = 0."))

      ### create long format dataset

      dat <- data.frame(slab, stringsAsFactors=FALSE)

      dat[[2]] <- xi
      dat[[3]] <- ti

      if (missing(var.names)) {
         names(dat) <- c("study", "events", "ptime")
      } else {
         if (length(var.names) != 3L)
            stop(mstyle$stop("Variable names not of length 3."))
         names(dat) <- var.names
      }

      dat[[1]] <- factor(dat[[1]])

      if (has.data && append)
         dat <- cbind(data, dat)

   }

   #########################################################################

   if (is.element(measure, c("MN","MNLN"))) {

      ### check for NAs in table data and act accordingly

      has.na <- is.na(mi) | is.na(sdi) | is.na(ni)

      if (any(has.na)) {

         not.na <- !has.na

         if (na.act == "na.omit") {
            mi   <- mi[not.na]
            sdi  <- sdi[not.na]
            ni   <- ni[not.na]
            slab <- slab[not.na]
            if (has.data)
               data <- data[not.na,]
            warning(mstyle$warning("Tables with NAs omitted."), call.=FALSE)
         }

         if (na.act == "na.fail")
            stop(mstyle$stop("Missing values in tables."))

      }

      k <- length(ni)

      ### at least one study left?

      if (k < 1)
         stop(mstyle$stop("Processing terminated since k = 0."))

      ### create long format dataset

      dat <- data.frame(slab, stringsAsFactors=FALSE)

      dat[[2]] <- mi
      dat[[3]] <- sdi
      dat[[4]] <- ni

      if (missing(var.names)) {
         names(dat) <- c("study", "mean", "sd", "n")
      } else {
         if (length(var.names) != 4L)
            stop(mstyle$stop("Variable names not of length 4."))
         names(dat) <- var.names
      }

      dat[[1]] <- factor(dat[[1]])

      if (has.data && append)
         dat <- cbind(data, dat)

   }

   #########################################################################

   if (is.element(measure, c("MC","SMCC","SMCR","SMCRH","ROMC","CVRC"))) {

      ### check for NAs in table data and act accordingly

      if (is.element(measure, c("MC","SMCC","SMCRH","ROMC","CVRC"))) {
         has.na <- is.na(m1i) | is.na(m2i) | is.na(sd1i) | is.na(sd2i) | is.na(ni) | is.na(ri)
      } else {
         has.na <- is.na(m1i) | is.na(m2i) | is.na(sd1i) | is.na(ni) | is.na(ri)
      }

      if (any(has.na)) {

         not.na <- !has.na

         if (na.act == "na.omit") {
            m1i  <- m1i[not.na]
            m2i  <- m2i[not.na]
            sd1i <- sd1i[not.na]
            if (is.element(measure, c("MC","SMCC","SMCRH","ROMC","CVRC")))
               sd2i <- sd2i[not.na]
            ni   <- ni[not.na]
            ri   <- ri[not.na]
            slab <- slab[not.na]
            if (has.data)
               data <- data[not.na,]
            warning(mstyle$warning("Tables with NAs omitted."), call.=FALSE)
         }

         if (na.act == "na.fail")
            stop(mstyle$stop("Missing values in tables."))

      }

      k <- length(m1i)

      ### at least one study left?

      if (k < 1)
         stop(mstyle$stop("Processing terminated since k = 0."))

      ### create long format dataset

      if (is.element(measure, c("MC","SMCC","SMCRH","ROMC","CVRC"))) {

         dat <- data.frame(slab, stringsAsFactors=FALSE)

         dat[[2]] <- m1i
         dat[[3]] <- m2i
         dat[[4]] <- sd1i
         dat[[5]] <- sd2i
         dat[[6]] <- ni
         dat[[7]] <- ri

         if (missing(var.names)) {
            names(dat) <- c("study", "mean1", "mean2", "sd1", "sd2", "n", "r")
         } else {
            if (length(var.names) != 7L)
               stop(mstyle$stop("Variable names not of length 7."))
            names(dat) <- var.names
         }

         dat[[1]] <- factor(dat[[1]])

         if (has.data && append)
            dat <- cbind(data, dat)

      } else {

         dat <- data.frame(slab, stringsAsFactors=FALSE)

         dat[[2]] <- m1i
         dat[[3]] <- m2i
         dat[[4]] <- sd1i
         dat[[5]] <- ni
         dat[[6]] <- ri

         if (missing(var.names)) {
            names(dat) <- c("study", "mean1", "mean2", "sd1", "n", "r")
         } else {
            if (length(var.names) != 6L)
               stop(mstyle$stop("Variable names not of length 6."))
            names(dat) <- var.names
         }

         dat[[1]] <- factor(dat[[1]])

         if (has.data && append)
            dat <- cbind(data, dat)

      }

   }

   #########################################################################

   if (is.element(measure, c("ARAW","AHW","ABT"))) {

      ### check for NAs in table data and act accordingly

      has.na <- is.na(ai) | is.na(mi) | is.na(ni)

      if (any(has.na)) {

         not.na <- !has.na

         if (na.act == "na.omit") {
            ai   <- ai[not.na]
            mi   <- mi[not.na]
            ni   <- ni[not.na]
            slab <- slab[not.na]
            if (has.data)
               data <- data[not.na,]
            warning(mstyle$warning("Tables with NAs omitted."), call.=FALSE)
         }

         if (na.act == "na.fail")
            stop(mstyle$stop("Missing values in tables."))

      }

      k <- length(ai)

      ### at least one study left?

      if (k < 1)
         stop(mstyle$stop("Processing terminated since k = 0."))

      ### create long format dataset

      dat <- data.frame(slab, stringsAsFactors=FALSE)

      dat[[2]] <- ai
      dat[[3]] <- mi
      dat[[4]] <- ni

      if (missing(var.names)) {
         names(dat) <- c("study", "alpha", "m", "n")
      } else {
         if (length(var.names) != 4L)
            stop(mstyle$stop("Variable names not of length 4."))
         names(dat) <- var.names
      }

      dat[[1]] <- factor(dat[[1]])

      if (has.data && append)
         dat <- data.frame(data, dat)

   }

   #########################################################################

   rownames(dat) <- seq_len(nrow(dat))
   return(dat)

}
