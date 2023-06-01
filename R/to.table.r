to.table <- function(measure, ai, bi, ci, di, n1i, n2i, x1i, x2i, t1i, t2i, m1i, m2i, sd1i, sd2i, xi, mi, ri, ti, sdi, ni,
data, slab, subset, add=1/2, to="none", drop00=FALSE, rows, cols) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   ### check argument specifications

   if (missing(measure))
      stop(mstyle$stop("Must specify an effect size or outcome measure via the 'measure' argument."))

   if (!is.character(measure))
      stop(mstyle$stop("The 'measure' argument must be a character string."))

   if (!is.element(measure, c("RR","OR","PETO","RD","AS","PHI","YUQ","YUY","RTET", # 2x2 table measures
                              "PBIT","OR2D","OR2DN","OR2DL",                       # - transformations to SMD
                              "MPRD","MPRR","MPOR","MPORC","MPPETO","MPORM",       # - measures for matched pairs data
                              "IRR","IRD","IRSD",                                  # two-group person-time data measures
                              "MD","SMD","SMDH","ROM",                             # two-group mean/SD measures
                              "CVR","VR",                                          # coefficient of variation ratio, variability ratio
                              "RPB","RBIS","D2OR","D2ORN","D2ORL",                 # - transformations to r_PB, r_BIS, and log(OR)
                              "COR","UCOR","ZCOR",                                 # correlations (raw and r-to-z transformed)
                              "PCOR","ZPCOR","SPCOR",                              # partial and semi-partial correlations
                              "R2","ZR2",                                          # coefficient of determination (raw and r-to-z transformed)
                              "PR","PLN","PLO","PAS","PFT",                        # single proportions (and transformations thereof)
                              "IR","IRLN","IRS","IRFT",                            # single-group person-time data (and transformations thereof)
                              "MN","MNLN","CVLN","SDLN","SMN",                     # mean, log(mean), log(CV), log(SD), standardized mean
                              "MC","SMCC","SMCR","SMCRH","ROMC","CVRC","VRC",      # raw/standardized mean change, log(ROM), CVR, and VR for dependent samples
                              "ARAW","AHW","ABT")))                                # alpha (and transformations thereof)
      stop(mstyle$stop("Unknown 'measure' specified."))

   if (is.element(measure, c("CVR","VR","PCOR","ZPCOR","SPCOR","R2","ZR2","CVLN","SDLN","VRC")))
      stop(mstyle$stop("Function not available for this outcome measure."))

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   if (!is.element(to, c("all","only0","if0all","none")))
      stop(mstyle$stop("Unknown 'to' argument specified."))

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

   ### get slab and subset arguments (will be NULL when unspecified)

   slab   <- .getx("slab",   mf=mf, data=data)
   subset <- .getx("subset", mf=mf, data=data)

   #########################################################################
   #########################################################################
   #########################################################################

   if (is.element(measure, c("RR","OR","RD","AS","PETO","PHI","YUQ","YUY","RTET","PBIT","OR2D","OR2DN","OR2DL","MPRD","MPRR","MPOR","MPORC","MPPETO","MPORM"))) {

      ai  <- .getx("ai",  mf=mf, data=data, checknumeric=TRUE)
      bi  <- .getx("bi",  mf=mf, data=data, checknumeric=TRUE)
      ci  <- .getx("ci",  mf=mf, data=data, checknumeric=TRUE)
      di  <- .getx("di",  mf=mf, data=data, checknumeric=TRUE)
      n1i <- .getx("n1i", mf=mf, data=data, checknumeric=TRUE)
      n2i <- .getx("n2i", mf=mf, data=data, checknumeric=TRUE)

      if (!.equal.length(ai, bi, ci, di, n1i, n2i))
         stop(mstyle$stop("Supplied data vectors are not all of the same length."))

      n1i.inc <- n1i != ai + bi
      n2i.inc <- n2i != ci + di

      if (any(n1i.inc, na.rm=TRUE))
         stop(mstyle$stop("One or more 'n1i' values are not equal to 'ai + bi'."))
      if (any(n2i.inc, na.rm=TRUE))
         stop(mstyle$stop("One or more 'n2i' values are not equal to 'ci + di'."))

      bi <- replmiss(bi, n1i-ai)
      di <- replmiss(di, n2i-ci)

      if (!.all.specified(ai, bi, ci, di))
         stop(mstyle$stop("Cannot compute outcomes. Check that all of the required information is specified\n  via the appropriate arguments (i.e., ai, bi, ci, di or ai, n1i, ci, n2i)."))

      k <- length(ai) ### number of outcomes before subsetting

      if (!is.null(subset)) {
         subset <- .chksubset(subset, k)
         ai <- .getsubset(ai, subset)
         bi <- .getsubset(bi, subset)
         ci <- .getsubset(ci, subset)
         di <- .getsubset(di, subset)
      }

      n1i <- ai + bi
      n2i <- ci + di

      if (any(c(ai > n1i, ci > n2i), na.rm=TRUE))
         stop(mstyle$stop("One or more event counts are larger than the corresponding group sizes."))

      if (any(c(ai, bi, ci, di) < 0, na.rm=TRUE))
         stop(mstyle$stop("One or more counts are negative."))

      if (any(c(n1i < 0, n2i < 0), na.rm=TRUE))
         stop(mstyle$stop("One or more group sizes are negative."))

      ni.u <- ai + bi + ci + di ### unadjusted total sample sizes

      ### if drop00=TRUE, set counts to NA for studies that have no events (or all events) in both arms

      if (drop00) {
         id00 <- c(ai == 0L & ci == 0L) | c(bi == 0L & di == 0L)
         id00[is.na(id00)] <- FALSE
         ai[id00] <- NA_real_
         bi[id00] <- NA_real_
         ci[id00] <- NA_real_
         di[id00] <- NA_real_
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

      x1i <- .getx("x1i", mf=mf, data=data, checknumeric=TRUE)
      x2i <- .getx("x2i", mf=mf, data=data, checknumeric=TRUE)
      t1i <- .getx("t1i", mf=mf, data=data, checknumeric=TRUE)
      t2i <- .getx("t2i", mf=mf, data=data, checknumeric=TRUE)

      if (!.all.specified(x1i, x2i, t1i, t2i))
         stop(mstyle$stop("Cannot compute outcomes. Check that all of the required information is specified\n  via the appropriate arguments (i.e., x1i, x2i, t1i, t2i)."))

      if (!.equal.length(x1i, x2i, t1i, t2i))
         stop(mstyle$stop("Supplied data vectors are not all of the same length."))

      k <- length(x1i) ### number of outcomes before subsetting

      if (!is.null(subset)) {
         subset <- .chksubset(subset, k)
         x1i <- .getsubset(x1i, subset)
         x2i <- .getsubset(x2i, subset)
         t1i <- .getsubset(t1i, subset)
         t2i <- .getsubset(t2i, subset)
      }

      if (any(c(x1i, x2i) < 0, na.rm=TRUE))
         stop(mstyle$stop("One or more counts are negative."))

      if (any(c(t1i, t2i) <= 0, na.rm=TRUE))
         stop(mstyle$stop("One or more person-times are <= 0."))

      ni.u <- t1i + t2i ### unadjusted total sample sizes

      ### if drop00=TRUE, set counts to NA for studies that have no events in both arms

      if (drop00) {
         id00 <- c(x1i == 0L & x2i == 0L)
         id00[is.na(id00)] <- FALSE
         x1i[id00] <- NA_real_
         x2i[id00] <- NA_real_
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

      m1i  <- .getx("m1i",  mf=mf, data=data, checknumeric=TRUE)
      m2i  <- .getx("m2i",  mf=mf, data=data, checknumeric=TRUE)
      sd1i <- .getx("sd1i", mf=mf, data=data, checknumeric=TRUE)
      sd2i <- .getx("sd2i", mf=mf, data=data, checknumeric=TRUE)
      n1i  <- .getx("n1i",  mf=mf, data=data, checknumeric=TRUE)
      n2i  <- .getx("n2i",  mf=mf, data=data, checknumeric=TRUE)

      if (!.all.specified(m1i, m2i, sd1i, sd2i, n1i, n2i))
         stop(mstyle$stop("Cannot compute outcomes. Check that all of the required information is specified\n  via the appropriate arguments (i.e., m1i, m2i, sd1i, sd2i, n1i, n2i)."))

      if (!.equal.length(m1i, m2i, sd1i, sd2i, n1i, n2i))
         stop(mstyle$stop("Supplied data vectors are not all of the same length."))

      k <- length(n1i) ### number of outcomes before subsetting

      if (!is.null(subset)) {
         subset <- .chksubset(subset, k)
         m1i  <- .getsubset(m1i,  subset)
         m2i  <- .getsubset(m2i,  subset)
         sd1i <- .getsubset(sd1i, subset)
         sd2i <- .getsubset(sd2i, subset)
         n1i  <- .getsubset(n1i,  subset)
         n2i  <- .getsubset(n2i,  subset)
      }

      if (any(c(sd1i, sd2i) < 0, na.rm=TRUE))
         stop(mstyle$stop("One or more standard deviations are negative."))

      if (any(c(n1i, n2i) <= 0, na.rm=TRUE))
         stop(mstyle$stop("One or more group sizes are <= 0."))

      ni.u <- n1i + n2i ### unadjusted total sample sizes

   }

   #########################################################################

   if (is.element(measure, c("COR","UCOR","ZCOR"))) {

      ri <- .getx("ri", mf=mf, data=data, checknumeric=TRUE)
      ni <- .getx("ni", mf=mf, data=data, checknumeric=TRUE)
      ti <- .getx("ti", mf=mf, data=data, checknumeric=TRUE)

      if (!.equal.length(ri, ni, ti))
         stop(mstyle$stop("Supplied data vectors are not all of the same length."))

      ri <- replmiss(ri, ti / sqrt(ni - 2 + ti^2))

      if (!.all.specified(ri, ni))
         stop(mstyle$stop("Cannot compute outcomes. Check that all of the required information is specified\n  via the appropriate arguments (i.e., ri, ni)."))

      k <- length(ri) ### number of outcomes before subsetting

      if (!is.null(subset)) {
         subset <- .chksubset(subset, k)
         ri <- .getsubset(ri, subset)
         ni <- .getsubset(ni, subset)
      }

      if (any(abs(ri) > 1, na.rm=TRUE))
         stop(mstyle$stop("One or more correlations are > 1 or < -1."))

      if (any(ni <= 0, na.rm=TRUE))
         stop(mstyle$stop("One or more sample sizes are <= 0."))

      ni.u <- ni ### unadjusted total sample sizes

   }

   #########################################################################

   if (is.element(measure, c("PR","PLN","PLO","PAS","PFT"))) {

      xi <- .getx("xi", mf=mf, data=data, checknumeric=TRUE)
      mi <- .getx("mi", mf=mf, data=data, checknumeric=TRUE)
      ni <- .getx("ni", mf=mf, data=data, checknumeric=TRUE)

      if (!.equal.length(xi, mi, ni))
         stop(mstyle$stop("Supplied data vectors are not all of the same length."))

      ni.inc <- ni != xi + mi

      if (any(ni.inc, na.rm=TRUE))
         stop(mstyle$stop("One or more 'ni' values are not equal to 'xi + mi'."))

      mi <- replmiss(mi, ni-xi)

      if (!.all.specified(xi, mi))
         stop(mstyle$stop("Cannot compute outcomes. Check that all of the required information is specified\n  via the appropriate arguments (i.e., xi, mi or xi, ni)."))

      k <- length(xi) ### number of outcomes before subsetting

      if (!is.null(subset)) {
         subset <- .chksubset(subset, k)
         xi <- .getsubset(xi, subset)
         mi <- .getsubset(mi, subset)
      }

      ni <- xi + mi

      if (any(xi > ni, na.rm=TRUE))
         stop(mstyle$stop("One or more event counts are larger than the corresponding group sizes."))

      if (any(c(xi, mi) < 0, na.rm=TRUE))
         stop(mstyle$stop("One or more counts are negative."))

      if (any(ni <= 0, na.rm=TRUE))
         stop(mstyle$stop("One or more group sizes are <= 0."))

      ni.u <- ni ### unadjusted total sample sizes

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

      xi <- .getx("xi", mf=mf, data=data, checknumeric=TRUE)
      ti <- .getx("ti", mf=mf, data=data, checknumeric=TRUE)

      if (!.all.specified(xi, ti))
         stop(mstyle$stop("Cannot compute outcomes. Check that all of the required information is specified\n  via the appropriate arguments (i.e., xi, ti)."))

      if (!.equal.length(xi, ti))
         stop(mstyle$stop("Supplied data vectors are not all of the same length."))

      k <- length(xi) ### number of outcomes before subsetting

      if (!is.null(subset)) {
         subset <- .chksubset(subset, k)
         xi <- .getsubset(xi, subset)
         ti <- .getsubset(ti, subset)
      }

      if (any(xi < 0, na.rm=TRUE))
         stop(mstyle$stop("One or more counts are negative."))

      if (any(ti <= 0, na.rm=TRUE))
         stop(mstyle$stop("One or more person-times are <= 0."))

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

   if (is.element(measure, c("MN","MNLN","SMN"))) {

      mi  <- .getx("mi",  mf=mf, data=data, checknumeric=TRUE)
      sdi <- .getx("sdi", mf=mf, data=data, checknumeric=TRUE)
      ni  <- .getx("ni",  mf=mf, data=data, checknumeric=TRUE)

      if (!.all.specified(mi, sdi, ni))
         stop(mstyle$stop("Cannot compute outcomes. Check that all of the required information is specified\n  via the appropriate arguments (i.e., mi, sdi, ni)."))

      if (!.equal.length(mi, sdi, ni))
         stop(mstyle$stop("Supplied data vectors are not all of the same length."))

      k <- length(ni) ### number of outcomes before subsetting

      if (!is.null(subset)) {
         subset <- .chksubset(subset, k)
         mi  <- .getsubset(mi,  subset)
         sdi <- .getsubset(sdi, subset)
         ni  <- .getsubset(ni,  subset)
      }

      if (any(sdi < 0, na.rm=TRUE))
         stop(mstyle$stop("One or more standard deviations are negative."))

      if (any(ni <= 0, na.rm=TRUE))
         stop(mstyle$stop("One or more sample sizes are <= 0."))

      if (is.element(measure, c("MNLN","CVLN")) && any(mi < 0, na.rm=TRUE))
         stop(mstyle$stop("One or more means are negative."))

      ni.u <- ni ### unadjusted total sample sizes

  }

   #########################################################################

   if (is.element(measure, c("MC","SMCC","SMCR","SMCRH","ROMC","CVRC"))) {

      m1i  <- .getx("m1i",  mf=mf, data=data, checknumeric=TRUE)
      m2i  <- .getx("m2i",  mf=mf, data=data, checknumeric=TRUE)
      sd1i <- .getx("sd1i", mf=mf, data=data, checknumeric=TRUE)
      sd2i <- .getx("sd2i", mf=mf, data=data, checknumeric=TRUE)
      ri   <- .getx("ri",   mf=mf, data=data, checknumeric=TRUE) ### for SMCR, do not need to supply this
      ni   <- .getx("ni",   mf=mf, data=data, checknumeric=TRUE)

      k <- length(m1i) ### number of outcomes before subsetting

      if (is.element(measure, c("MC","SMCC","SMCRH","ROMC","CVRC"))) {

         if (!.all.specified(m1i, m2i, sd1i, sd2i, ni, ri))
            stop(mstyle$stop("Cannot compute outcomes. Check that all of the required information is specified\n  via the appropriate arguments (i.e., m1i, m2i, sd1i, sd2i, ni, ri)."))

         if (!.equal.length(m1i, m2i, sd1i, sd2i, ni, ri))
            stop(mstyle$stop("Supplied data vectors are not all of the same length."))

      } else {

         if (!.all.specified(m1i, m2i, sd1i, ni, ri))
            stop(mstyle$stop("Cannot compute outcomes. Check that all of the required information is specified\n  via the appropriate arguments (i.e., m1i, m2i, sd1i, ni, ri)."))

         if (!.equal.length(m1i, m2i, sd1i, ni, ri))
            stop(mstyle$stop("Supplied data vectors are not all of the same length."))

      }

      if (!is.null(subset)) {
         subset <- .chksubset(subset, k)
         m1i  <- .getsubset(m1i,  subset)
         m2i  <- .getsubset(m2i,  subset)
         sd1i <- .getsubset(sd1i, subset)
         sd2i <- .getsubset(sd2i, subset)
         ni   <- .getsubset(ni,   subset)
         ri   <- .getsubset(ri,   subset)
      }

      if (is.element(measure, c("MC","SMCC","SMCRH","ROMC","CVRC"))) {
         if (any(c(sd1i, sd2i) < 0, na.rm=TRUE))
            stop(mstyle$stop("One or more standard deviations are negative."))
      } else {
         if (any(sd1i < 0, na.rm=TRUE))
            stop(mstyle$stop("One or more standard deviations are negative."))
      }

      if (any(abs(ri) > 1, na.rm=TRUE))
         stop(mstyle$stop("One or more correlations are > 1 or < -1."))

      if (any(ni <= 0, na.rm=TRUE))
         stop(mstyle$stop("One or more sample sizes are <= 0."))

      ni.u <- ni ### unadjusted total sample sizes

   }

   #########################################################################

   if (is.element(measure, c("ARAW","AHW","ABT"))) {

      ai <- .getx("ai", mf=mf, data=data, checknumeric=TRUE)
      mi <- .getx("mi", mf=mf, data=data, checknumeric=TRUE)
      ni <- .getx("ni", mf=mf, data=data, checknumeric=TRUE)

      if (!.all.specified(ai, mi, ni))
         stop(mstyle$stop("Cannot compute outcomes. Check that all of the required information is specified\n  via the appropriate arguments (i.e., ai, mi, ni)."))

      if (!.equal.length(ai, mi, ni))
         stop(mstyle$stop("Supplied data vectors are not all of the same length."))

      k <- length(ai) ### number of outcomes before subsetting

      if (!is.null(subset)) {
         subset <- .chksubset(subset, k)
         ai <- .getsubset(ai, subset)
         mi <- .getsubset(mi, subset)
         ni <- .getsubset(ni, subset)
      }

      if (any(ai > 1, na.rm=TRUE))
         stop(mstyle$stop("One or more alpha values are > 1."))

      if (any(mi < 2, na.rm=TRUE))
         stop(mstyle$stop("One or more mi values are < 2."))

      if (any(ni <= 0, na.rm=TRUE))
         stop(mstyle$stop("One or more sample sizes are <= 0."))

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

   if (!is.null(subset))
      slab <- .getsubset(slab, subset)

   ### check if study labels are unique; if not, make them unique

   if (anyDuplicated(slab))
      slab <- .make.unique(slab)

   #########################################################################
   #########################################################################
   #########################################################################

   if (is.element(measure, c("RR","OR","RD","AS","PETO","PHI","YUQ","YUY","RTET","PBIT","OR2D","OR2DN","OR2DL","MPORM"))) {

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
            warning(mstyle$warning(paste(sum(has.na), ifelse(sum(has.na) > 1, "tables", "table"), "with NAs omitted.")), call.=FALSE)
         }

         if (na.act == "na.fail")
            stop(mstyle$stop("Missing values in tables."))

      }

      k <- length(ai)

      ### at least one study left?

      if (k < 1)
         stop(mstyle$stop("Processing terminated since k = 0."))

      ### row/group and column/outcome names

      if (missing(rows)) {
         rows <- c("Grp1", "Grp2")
      } else {
         if (length(rows) != 2L)
            stop(mstyle$stop("Group names not of length 2."))
      }

      if (missing(cols)) {
         cols <- c("Out1", "Out2")
      } else {
         if (length(cols) != 2L)
            stop(mstyle$stop("Outcome names not of length 2."))
      }

      dat <- array(NA_real_, dim=c(2,2,k), dimnames=list(rows, cols, slab))

      for (i in seq_len(k)) {
         tab.i <- rbind(c(ai[i],bi[i]), c(ci[i],di[i]))
         dat[,,i] <- tab.i
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
            warning(mstyle$warning(paste(sum(has.na), ifelse(sum(has.na) > 1, "tables", "table"), "with NAs omitted.")), call.=FALSE)
         }

         if (na.act == "na.fail")
            stop(mstyle$stop("Missing values in tables."))

      }

      k <- length(ai)

      ### at least one study left?

      if (k < 1)
         stop(mstyle$stop("Processing terminated since k = 0."))

      ### row/group and column/outcome names

      if (missing(rows)) {
         rows <- c("Time1", "Time2")
      } else {
         if (length(rows) != 2L)
            stop(mstyle$stop("Time names not of length 2."))
      }

      if (missing(cols)) {
         cols <- c("Out1", "Out2")
      } else {
         if (length(cols) != 2L)
            stop(mstyle$stop("Outcome names not of length 2."))
      }

      dat <- array(NA_real_, dim=c(2,2,k), dimnames=list(rows, cols, slab))

      for (i in seq_len(k)) {
         tab.i <- rbind(c(ai[i]+bi[i],ci[i]+di[i]), c(ai[i]+ci[i],bi[i]+di[i]))
         dat[,,i] <- tab.i
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
            warning(mstyle$warning(paste(sum(has.na), ifelse(sum(has.na) > 1, "tables", "table"), "with NAs omitted.")), call.=FALSE)
         }

         if (na.act == "na.fail")
            stop(mstyle$stop("Missing values in tables."))

      }

      k <- length(ai)

      ### at least one study left?

      if (k < 1)
         stop(mstyle$stop("Processing terminated since k = 0."))

      ### row/group and column/outcome names

      if (missing(rows)) {
         rows <- c("Time1.Out1", "Time1.Out2")
      } else {
         if (length(rows) != 2L)
            stop(mstyle$stop("Time1 names not of length 2."))
      }

      if (missing(cols)) {
         cols <- c("Time2.Out1", "Time2.Out2")
      } else {
         if (length(cols) != 2L)
            stop(mstyle$stop("Time2 names not of length 2."))
      }

      dat <- array(NA_real_, dim=c(2,2,k), dimnames=list(rows, cols, slab))

      for (i in seq_len(k)) {
         tab.i <- rbind(c(ai[i],bi[i]), c(ci[i],di[i]))
         dat[,,i] <- tab.i
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
            warning(mstyle$warning(paste(sum(has.na), ifelse(sum(has.na) > 1, "tables", "table"), "with NAs omitted.")), call.=FALSE)
         }

         if (na.act == "na.fail")
            stop(mstyle$stop("Missing values in tables."))

      }

      k <- length(x1i)

      ### at least one study left?

      if (k < 1)
         stop(mstyle$stop("Processing terminated since k = 0."))

      ### row/group and column/outcome names

      if (missing(rows)) {
         rows <- c("Grp1", "Grp2")
      } else {
         if (length(rows) != 2L)
            stop(mstyle$stop("Group names not of length 2."))
      }

      if (missing(cols)) {
         cols <- c("Events", "Person-Time")
      } else {
         if (length(cols) != 2L)
            stop(mstyle$stop("Outcome names not of length 2."))
      }

      dat <- array(NA_real_, dim=c(2,2,k), dimnames=list(rows, cols, slab))

      for (i in seq_len(k)) {
         tab.i <- rbind(c(x1i[i],t1i[i]), c(x2i[i],t2i[i]))
         dat[,,i] <- tab.i
      }

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
            warning(mstyle$warning(paste(sum(has.na), ifelse(sum(has.na) > 1, "tables", "table"), "with NAs omitted.")), call.=FALSE)
         }

         if (na.act == "na.fail")
            stop(mstyle$stop("Missing values in tables."))

      }

      k <- length(m1i)

      ### at least one study left?

      if (k < 1)
         stop(mstyle$stop("Processing terminated since k = 0."))

      ### row/group and column/outcome names

      if (missing(rows)) {
         rows <- c("Grp1", "Grp2")
      } else {
         if (length(rows) != 2L)
            stop(mstyle$stop("Group names not of length 2."))
      }

      if (missing(cols)) {
         cols <- c("Mean", "SD", "n")
      } else {
         if (length(cols) != 3L)
            stop(mstyle$stop("Outcome names not of length 3."))
      }

      dat <- array(NA_real_, dim=c(2,3,k), dimnames=list(rows, cols, slab))

      for (i in seq_len(k)) {
         tab.i <- rbind(c(m1i[i],sd1i[i],n1i[i]), c(m2i[i],sd2i[i],n2i[i]))
         dat[,,i] <- tab.i
      }

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
            warning(mstyle$warning(paste(sum(has.na), ifelse(sum(has.na) > 1, "tables", "table"), "with NAs omitted.")), call.=FALSE)
         }

         if (na.act == "na.fail")
            stop(mstyle$stop("Missing values in tables."))

      }

      k <- length(ri)

      ### at least one study left?

      if (k < 1)
         stop(mstyle$stop("Processing terminated since k = 0."))

      ### row/group and column/outcome names

      if (missing(rows)) {
         rows <- c("Grp")
      } else {
         if (length(rows) != 1L)
            stop(mstyle$stop("Group names not of length 1."))
      }

      if (missing(cols)) {
         cols <- c("r", "n")
      } else {
         if (length(cols) != 2L)
            stop(mstyle$stop("Outcome names not of length 2."))
      }

      dat <- array(NA_real_, dim=c(1,2,k), dimnames=list(rows, cols, slab))

      for (i in seq_len(k)) {
         tab.i <- c(ri[i],ni[i])
         dat[,,i] <- tab.i
      }

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
            warning(mstyle$warning(paste(sum(has.na), ifelse(sum(has.na) > 1, "tables", "table"), "with NAs omitted.")), call.=FALSE)
         }

         if (na.act == "na.fail")
            stop(mstyle$stop("Missing values in tables."))

      }

      k <- length(xi)

      ### at least one study left?

      if (k < 1)
         stop(mstyle$stop("Processing terminated since k = 0."))

      ### row/group and column/outcome names

      if (missing(rows)) {
         rows <- c("Grp")
      } else {
         if (length(rows) != 1L)
            stop(mstyle$stop("Group names not of length 1."))
      }

      if (missing(cols)) {
         cols <- c("Out1", "Out2")
      } else {
         if (length(cols) != 2L)
            stop(mstyle$stop("Outcome names not of length 2."))
      }

      dat <- array(NA_real_, dim=c(1,2,k), dimnames=list(rows, cols, slab))

      for (i in seq_len(k)) {
         tab.i <- c(xi[i],mi[i])
         dat[,,i] <- tab.i
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
            warning(mstyle$warning(paste(sum(has.na), ifelse(sum(has.na) > 1, "tables", "table"), "with NAs omitted.")), call.=FALSE)
         }

         if (na.act == "na.fail")
            stop(mstyle$stop("Missing values in tables."))

      }

      k <- length(xi)

      ### at least one study left?

      if (k < 1)
         stop(mstyle$stop("Processing terminated since k = 0."))

      ### row/group and column/outcome names

      if (missing(rows)) {
         rows <- c("Grp")
      } else {
         if (length(rows) != 1L)
            stop(mstyle$stop("Group names not of length 1."))
      }

      if (missing(cols)) {
         cols <- c("Events", "Person-Time")
      } else {
         if (length(cols) != 2L)
            stop(mstyle$stop("Outcome names not of length 2."))
      }

      dat <- array(NA_real_, dim=c(1,2,k), dimnames=list(rows, cols, slab))

      for (i in seq_len(k)) {
         tab.i <- c(xi[i],ti[i])
         dat[,,i] <- tab.i
      }

   }

   #########################################################################

   if (is.element(measure, c("MN","MNLN","SMN"))) {

      ### check for NAs in table data and act accordingly

      has.na <- is.na(mi) | is.na(sdi) | is.na(ni)

      if (any(has.na)) {

         not.na <- !has.na

         if (na.act == "na.omit") {
            mi   <- mi[not.na]
            sdi  <- sdi[not.na]
            ni   <- ni[not.na]
            slab <- slab[not.na]
            warning(mstyle$warning(paste(sum(has.na), ifelse(sum(has.na) > 1, "tables", "table"), "with NAs omitted.")), call.=FALSE)
         }

         if (na.act == "na.fail")
            stop(mstyle$stop("Missing values in tables."))

      }

      k <- length(ni)

      ### at least one study left?

      if (k < 1)
         stop(mstyle$stop("Processing terminated since k = 0."))

      ### row/group and column/outcome names

      if (missing(rows)) {
         rows <- c("Grp")
      } else {
         if (length(rows) != 1L)
            stop(mstyle$stop("Group names not of length 1."))
      }

      if (missing(cols)) {
         cols <- c("Mean", "SD", "n")
      } else {
         if (length(cols) != 3L)
            stop(mstyle$stop("Outcome names not of length 3."))
      }

      dat <- array(NA_real_, dim=c(1,3,k), dimnames=list(rows, cols, slab))

      for (i in seq_len(k)) {
         tab.i <- c(mi[i],sdi[i],ni[i])
         dat[,,i] <- tab.i
      }

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
            warning(mstyle$warning(paste(sum(has.na), ifelse(sum(has.na) > 1, "tables", "table"), "with NAs omitted.")), call.=FALSE)
         }

         if (na.act == "na.fail")
            stop(mstyle$stop("Missing values in tables."))

      }

      k <- length(m1i)

      ### at least one study left?

      if (k < 1)
         stop(mstyle$stop("Processing terminated since k = 0."))

      ### row/group and column/outcome names

      if (missing(rows)) {
         rows <- c("Grp")
      } else {
         if (length(rows) != 1L)
            stop(mstyle$stop("Group names not of length 1."))
      }

      if (is.element(measure, c("MC","SMCC","SMCRH","ROMC","CVRC"))) {
         if (missing(cols)) {
            cols <- c("Mean1", "Mean2", "SD1", "SD2", "n", "r")
         } else {
            if (length(cols) != 6L)
               stop(mstyle$stop("Outcome names not of length 6."))
         }
      } else {
         if (missing(cols)) {
            cols <- c("Mean1", "Mean2", "SD1", "n", "r")
         } else {
            if (length(cols) != 5L)
               stop(mstyle$stop("Outcome names not of length 5."))
         }
      }

      if (is.element(measure, c("MC","SMCC","SMCRH","ROMC","CVRC"))) {

         dat <- array(NA_real_, dim=c(1,6,k), dimnames=list(rows, cols, slab))

         for (i in seq_len(k)) {
            tab.i <- c(m1i[i],m2i[i],sd1i[i],sd2i[i],ni[i],ri[i])
            dat[,,i] <- tab.i
         }

      } else {

         dat <- array(NA_real_, dim=c(1,5,k), dimnames=list(rows, cols, slab))

         for (i in seq_len(k)) {
            tab.i <- c(m1i[i],m2i[i],sd1i[i],ni[i],ri[i])
            dat[,,i] <- tab.i
         }

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
            warning(mstyle$warning(paste(sum(has.na), ifelse(sum(has.na) > 1, "tables", "table"), "with NAs omitted.")), call.=FALSE)
         }

         if (na.act == "na.fail")
            stop(mstyle$stop("Missing values in tables."))

      }

      k <- length(ai)

      ### at least one study left?

      if (k < 1)
         stop(mstyle$stop("Processing terminated since k = 0."))

      ### row/group and column/outcome names

      if (missing(rows)) {
         rows <- c("Grp")
      } else {
         if (length(rows) != 1L)
            stop(mstyle$stop("Group names not of length 1."))
      }

      if (missing(cols)) {
         cols <- c("alpha", "m", "n")
      } else {
         if (length(cols) != 3L)
            stop(mstyle$stop("Outcome names not of length 3."))
      }

      dat <- array(NA_real_, dim=c(1,3,k), dimnames=list(rows, cols, slab))

      for (i in seq_len(k)) {
         tab.i <- c(ai[i],mi[i],ni[i])
         dat[,,i] <- tab.i
      }

   }

   #########################################################################

   return(dat)

}
