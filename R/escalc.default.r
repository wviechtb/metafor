escalc.default <- function(measure, formula, ai, bi, ci, di, n1i, n2i, x1i, x2i, t1i, t2i, m1i, m2i, sd1i, sd2i, xi, mi, ri, ti, sdi, ni, yi, vi, sei,
data, slab, subset, add=1/2, to="only0", drop00=FALSE, vtype="LS", var.names=c("yi","vi"), add.measure=FALSE, append=TRUE, replace=TRUE, digits=4, ...) {

   ### check argument specifications

   if (!is.element(measure, c("RR","OR","PETO","RD","AS","PHI","YUQ","YUY","RTET", ### 2x2 table measures
                              "PBIT","OR2D","OR2DN","OR2DL",                       ### - transformations to SMD
                              "IRR","IRD","IRSD",                                  ### two-group person-time data measures
                              "MD","SMD","SMDH","ROM",                             ### two-group mean/SD measures
                              "RPB","RBIS","D2OR","D2ORN","D2ORL",                 ### - transformations to r_PB, r_BIS, and log(OR)
                              "COR","UCOR","ZCOR",                                 ### correlations (raw and r-to-z transformed)
                              "PR","PLN","PLO","PAS","PFT",                        ### single proportions (and transformations thereof)
                              "IR","IRLN","IRS","IRFT",                            ### single-group person-time data (and transformations thereof)
                              "MN","MC","SMCC","SMCR","SMCRH","ROMC",              ### raw/standardized mean change and log(ROM) for dependent samples
                              "ARAW","AHW","ABT",                                  ### alpha (and transformations thereof)
                              "GEN")))
      stop("Unknown 'measure' specified.")

   if (!is.element(to, c("all","only0","if0all","none")))
      stop("Unknown 'to' argument specified.")

   if (any(!is.element(vtype, c("UB","LS","HO","ST","CS")), na.rm=TRUE)) ### vtype can be an entire vector, so use any() and na.rm=TRUE
      stop("Unknown 'vtype' argument specified.")

   if (add.measure) {

      if (length(var.names) == 2)
         var.names <- c(var.names, "measure")

      if (length(var.names) != 3)
         stop("Argument 'var.names' must be of length 2 or 3.")

      if (any(var.names != make.names(var.names, unique=TRUE))) {
         var.names <- make.names(var.names, unique=TRUE)
         warning(paste0("Argument 'var.names' does not contain syntactically valid variable names.\n  Variable names adjusted to: var.names = c('", var.names[1], "', '", var.names[2], "', '", var.names[3], "')."))
      }

   } else {

      if (length(var.names) == 3)
         var.names <- var.names[1:2]

      if (length(var.names) != 2)
         stop("Argument 'var.names' must be of length 2.")

      if (any(var.names != make.names(var.names, unique=TRUE))) {
         var.names <- make.names(var.names, unique=TRUE)
         warning(paste0("Argument 'var.names' does not contain syntactically valid variable names.\n  Variable names adjusted to: var.names = c('", var.names[1], "', '", var.names[2], "')."))
      }

   }

   ### check if data argument has been specified

   if (missing(data))
      data <- NULL

   ### need this at the end to check if append=TRUE can actually be done

   no.data <- is.null(data)

   ### check if data argument has been specified

   if (is.null(data)) {
      data <- sys.frame(sys.parent())
   } else {
      if (!is.data.frame(data)) {
         data <- data.frame(data)
      }
   }

   mf <- match.call()

   ### get slab and subset arguments (will be NULL when unspecified)

   mf.slab   <- mf[[match("slab",   names(mf))]]
   mf.subset <- mf[[match("subset", names(mf))]]
   slab      <- eval(mf.slab,   data, enclos=sys.frame(sys.parent()))
   subset    <- eval(mf.subset, data, enclos=sys.frame(sys.parent()))

   ### get yi (in case it has been specified)

   mf.yi <- mf[[match("yi", names(mf))]]
   yi <- eval(mf.yi, data, enclos=sys.frame(sys.parent()))

   #########################################################################
   #########################################################################
   #########################################################################

   if (is.null(yi)) {

      if (is.element(measure, c("RR","OR","RD","AS","PETO","PHI","YUQ","YUY","RTET","PBIT","OR2D","OR2DN","OR2DL"))) {

         mf.ai   <- mf[[match("ai",  names(mf))]]
         mf.bi   <- mf[[match("bi",  names(mf))]]
         mf.ci   <- mf[[match("ci",  names(mf))]]
         mf.di   <- mf[[match("di",  names(mf))]]
         mf.n1i  <- mf[[match("n1i", names(mf))]]
         mf.n2i  <- mf[[match("n2i", names(mf))]]
         ai      <- eval(mf.ai,  data, enclos=sys.frame(sys.parent()))
         bi      <- eval(mf.bi,  data, enclos=sys.frame(sys.parent()))
         ci      <- eval(mf.ci,  data, enclos=sys.frame(sys.parent()))
         di      <- eval(mf.di,  data, enclos=sys.frame(sys.parent()))
         n1i     <- eval(mf.n1i, data, enclos=sys.frame(sys.parent()))
         n2i     <- eval(mf.n2i, data, enclos=sys.frame(sys.parent()))
         if (is.null(bi)) bi <- n1i - ai
         if (is.null(di)) di <- n2i - ci

         if (!is.null(subset)) {
            ai <- ai[subset]
            bi <- bi[subset]
            ci <- ci[subset]
            di <- di[subset]
         }

         if (length(ai)==0L || length(bi)==0L || length(ci)==0L || length(di)==0L)
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")

         if (!all(length(ai) == c(length(ai),length(bi),length(ci),length(di))))
            stop("Supplied data vectors are not all of the same length.")

         if (any(c(ai, bi, ci, di) < 0, na.rm=TRUE))
            stop("One or more counts are negative.")

         ni.u <- ai + bi + ci + di ### unadjusted total sample sizes

         k <- length(ai)

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

         ### recompute group and total sample sizes (after add/to adjustment)

         n1i <- ai + bi
         n2i <- ci + di
         ni <- n1i + n2i ### ni.u computed earlier is always the 'unadjusted' total sample size

         ### compute proportions for the two groups

         p1i <- ai/n1i
         p2i <- ci/n2i

         ### (log) relative risk

         if (measure == "RR") {
            yi <- log(p1i) - log(p2i)
            vi <- 1/ai - 1/n1i + 1/ci - 1/n2i
         }

         ### (log) odds ratio

         if (is.element(measure, c("OR","OR2D","OR2DN","OR2DL"))) {
            yi <- log(p1i/(1-p1i)) - log(p2i/(1-p2i))
            vi <- 1/ai + 1/bi + 1/ci + 1/di
         }

         ### (log) odds ratio (Peto's method)

         if (measure == "PETO") {
            xt <- ai + ci ### frequency of outcome1 in both groups combined
            yt <- bi + di ### frequency of outcome2 in both groups combined
            Oi <- ai
            Ei <- xt * n1i / ni
            Vi <- xt * yt * (n1i/ni) * (n2i/ni) / (ni - 1) ### 0 when xt = 0 or yt = 0 in a table
            yi <- (Oi - Ei) / Vi                           ### then yi and vi is Inf (set to NA at end)
            vi <- 1/Vi
         }

         ### risk difference

         if (measure == "RD") {

            yi <- p1i - p2i

            if (length(vtype) == 1L)
               vtype <- rep(vtype, k)

            vi <- rep(NA_real_, k)

            mnwp1i <- sum(ai, na.rm=TRUE) / sum(n1i, na.rm=TRUE) ### sample size weighted average of proportions (same as sum(n1i*p1i)/sum(n1i))
            mnwp2i <- sum(ci, na.rm=TRUE) / sum(n2i, na.rm=TRUE) ### sample size weighted average of proportions (same as sum(n2i*p2i)/sum(n2i))

            for (i in seq_len(k)) {

               ### unbiased estimate of the sampling variance
               if (vtype[i] == "UB")
                  vi[i] <- p1i[i]*(1-p1i[i])/(n1i[i]-1) + p2i[i]*(1-p2i[i])/(n2i[i]-1)

               ### large sample approximation to the sampling variance
               if (vtype[i] == "LS")
                  vi[i] <- p1i[i]*(1-p1i[i])/n1i[i] + p2i[i]*(1-p2i[i])/n2i[i]

               ### estimator assuming homogeneity (using the average proportions)
               if (vtype[i] == "HO")
                  vi[i] <- mnwp1i*(1-mnwp1i)/n1i[i] + mnwp2i*(1-mnwp2i)/n2i[i]

            }

         }

         ### arcsine square-root risk difference

         if (measure == "AS") {
            yi <- asin(sqrt(p1i)) - asin(sqrt(p2i))
            vi <- 1/(4*n1i) + 1/(4*n2i)
         }

         ### phi coefficient (vi equation in Yule, 1912, p.603)

         if (measure == "PHI") {
            yi  <- (ai*di - bi*ci)/sqrt((ai+bi)*(ci+di)*(ai+ci)*(bi+di))
            p1. <- (ai+bi)/ni
            p2. <- (ci+di)/ni
            p.1 <- (ai+ci)/ni
            p.2 <- (bi+di)/ni
            vi  <- 1/ni * (1 - yi^2 + yi*(1+1/2*yi^2) * (p1.-p2.)*(p.1-p.2) / sqrt(p1.*p2.*p.1*p.2) - 3/4 * yi^2 * ((p1.-p2.)^2/(p1.*p2.) + (p.1-p.2)^2/(p.1*p.2)))
         }

         ### Yule's Q (vi equation in Yule, 1900, p.285, and Yule, 1912, p.593)

         if (measure == "YUQ") {
            yi <- (ai/bi)/(ci/di)
            yi <- (yi-1)/(yi+1)
            vi <- 1/4 * (1-yi^2)^2 * (1/ai + 1/bi + 1/ci + 1/di)
         }

         ### Yule's Y (vi equation in Yule, 1912, p.593)

         if (measure == "YUY") {
            yi <- (ai/bi)/(ci/di)
            yi <- (sqrt(yi)-1)/(sqrt(yi)+1)
            vi <- 1/16 * (1-yi^2)^2 * (1/ai + 1/bi + 1/ci + 1/di)
         }

         ### tetrachoric correlation

         if (measure == "RTET") {

            ### TODO: allow user to set control arguments for pmvnorm and optimizers

            ### upgrade warnings to errors (so that tables with no events or only events are skipped)
            #warn.before <- getOption("warn")
            #options(warn = 2)

            yi <- rep(NA_real_, k)
            vi <- rep(NA_real_, k)

            for (i in seq_len(k)) {

               if (is.na(ai[i]) || is.na(bi[i]) || is.na(ci[i]) || is.na(di[i]))
                  next

               res <- .rtet(ai[i], bi[i], ci[i], di[i], maxcor=.9999)

               yi[i] <- res$yi
               vi[i] <- res$vi

            }

            #options(warn = warn.before)

         }

         ### probit transformation to SMD

         if (measure == "PBIT") {
            z1i <- qnorm(p1i)
            z2i <- qnorm(p2i)
            yi <- z1i - z2i
            vi <- 2*pi*p1i*(1-p1i)*exp(z1i^2)/n1i + 2*pi*p2i*(1-p2i)*exp(z2i^2)/n2i ### from Sanchez-Meca et al. (2003) and Rosenthal (1994; Handbook chapter)
         }                                                                          ### seems to be right for stratified and cross-sectional/multinomial sampling
                                                                                    ### see code/probit directory
         ### log(OR) transformation to SMD based on logistic distribution

         if (is.element(measure, c("OR2D","OR2DL"))) {
            yi <- sqrt(3) / pi * yi
            vi <- 3 / pi^2 * vi
         }

         ### log(OR) transformation to SMD based on normal distribution (Cox & Snell method)

         if (measure == "OR2DN") {
            yi <- yi / 1.65
            vi <- vi / 1.65^2
         }

         ### sensitivity

         #if (measure == "SENS") {
         #   res <- escalc("PR", xi=ai, mi=ci, add=0, to="none", vtype=vtype)
         #   yi <- res$yi
         #   vi <- res$vi
         #}

         ### specificity

         #if (measure == "SPEC") {
         #   res <- escalc("PR", xi=di, mi=bi, add=0, to="none", vtype=vtype)
         #   yi <- res$yi
         #   vi <- res$vi
         #}

         ### [...]

         ### Note: Could in principle also compute measures commonly used in diagnostic studies.
         ### But need to take the sampling method into consideration when computing vi (so need
         ### to give this some more thought).

      }

      ######################################################################

      if (is.element(measure, c("IRR","IRD","IRSD"))) {

         mf.x1i  <- mf[[match("x1i", names(mf))]]
         mf.x2i  <- mf[[match("x2i", names(mf))]]
         mf.t1i  <- mf[[match("t1i", names(mf))]]
         mf.t2i  <- mf[[match("t2i", names(mf))]]
         x1i     <- eval(mf.x1i, data, enclos=sys.frame(sys.parent()))
         x2i     <- eval(mf.x2i, data, enclos=sys.frame(sys.parent()))
         t1i     <- eval(mf.t1i, data, enclos=sys.frame(sys.parent()))
         t2i     <- eval(mf.t2i, data, enclos=sys.frame(sys.parent()))

         if (!is.null(subset)) {
            x1i <- x1i[subset]
            x2i <- x2i[subset]
            t1i <- t1i[subset]
            t2i <- t2i[subset]
         }

         if (length(x1i)==0L || length(x2i)==0L || length(t1i)==0L || length(t2i)==0L)
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")

         if (!all(length(x1i) == c(length(x1i),length(x2i),length(t1i),length(t2i))))
            stop("Supplied data vectors are not all of the same length.")

         if (any(c(x1i, x2i) < 0, na.rm=TRUE))
            stop("One or more counts are negative.")

         if (any(c(t1i, t2i) < 0, na.rm=TRUE))
            stop("One or more person-times are negative.")

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

         ### compute rates for the two groups

         ir1i <- x1i/t1i ### t1i is the total person-time in the 1st group
         ir2i <- x2i/t2i ### t2i is the total person-time in the 2nd group

         ### (log) incidence rate ratio

         if (measure == "IRR") {
            yi <- log(ir1i) - log(ir2i)
            vi <- 1/x1i + 1/x2i
            #vi <- 1/(x1i+1/2) + 1/(x2i+1/2)
         }

         ### incidence rate difference

         if (measure == "IRD") {
            yi <- ir1i - ir2i
            vi <- ir1i/t1i + ir2i/t2i ### note: same as x1i/t1i^2 + x2i/t2i^2
         }

         ### square-root transformed incidence rate difference

         if (measure == "IRSD") {
            yi <- sqrt(ir1i) - sqrt(ir2i)
            vi <- 1/(4*t1i) + 1/(4*t2i)
         }

      }

      ######################################################################

      if (is.element(measure, c("MD","SMD","SMDH","ROM","RPB","RBIS","D2OR","D2ORN","D2ORL"))) {

         mf.m1i  <- mf[[match("m1i",  names(mf))]]
         mf.m2i  <- mf[[match("m2i",  names(mf))]]
         mf.sd1i <- mf[[match("sd1i", names(mf))]]
         mf.sd2i <- mf[[match("sd2i", names(mf))]]
         mf.n1i  <- mf[[match("n1i",  names(mf))]]
         mf.n2i  <- mf[[match("n2i",  names(mf))]]
         m1i     <- eval(mf.m1i,  data, enclos=sys.frame(sys.parent()))
         m2i     <- eval(mf.m2i,  data, enclos=sys.frame(sys.parent()))
         sd1i    <- eval(mf.sd1i, data, enclos=sys.frame(sys.parent()))
         sd2i    <- eval(mf.sd2i, data, enclos=sys.frame(sys.parent()))
         n1i     <- eval(mf.n1i,  data, enclos=sys.frame(sys.parent()))
         n2i     <- eval(mf.n2i,  data, enclos=sys.frame(sys.parent()))

         if (!is.null(subset)) {
            m1i  <- m1i[subset]
            m2i  <- m2i[subset]
            sd1i <- sd1i[subset]
            sd2i <- sd2i[subset]
            n1i  <- n1i[subset]
            n2i  <- n2i[subset]
         }

         if (length(m1i)==0L || length(m2i)==0L || length(sd1i)==0L || length(sd2i)==0L || length(n1i)==0L || length(n2i)==0L)
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")

         if (!all(length(m1i) == c(length(m1i),length(m2i),length(sd1i),length(sd2i),length(n1i),length(n2i))))
            stop("Supplied data vectors are not all of the same length.")

         if (any(c(sd1i, sd2i) < 0, na.rm=TRUE))
            stop("One or more standard deviations are negative.")

         if (any(c(n1i, n2i) < 0, na.rm=TRUE))
            stop("One or more sample sizes are negative.")

         ni.u <- n1i + n2i ### unadjusted total sample sizes

         k <- length(m1i)

         ni <- ni.u
         mi <- ni - 2
         sdpi <- sqrt(((n1i-1)*sd1i^2 + (n2i-1)*sd2i^2)/mi)
         di <- (m1i - m2i) / sdpi

         ### (raw) mean difference (with heteroscedastic variances)
         ### to use with pooled SDs, simply set sd1i = sd2i = sdpi or use vtype="HO"

         if (measure == "MD") {

            yi <- m1i - m2i

            if (length(vtype) == 1L)
               vtype <- rep(vtype, k)

            vi <- rep(NA_real_, k)

            for (i in seq_len(k)) {

               ### unbiased estimate of the sampling variance (does not assume homoscedasticity)
               if (vtype[i] == "UB" || vtype[i] == "LS")
                  vi[i] <- sd1i[i]^2/n1i[i] + sd2i[i]^2/n2i[i]

               ### estimator assuming homoscedasticity
               if (vtype[i] == "HO")
                  vi[i] <- sdpi[i]^2 * (1/n1i[i] + 1/n2i[i])

            }

         }

         ### standardized mean difference (with pooled SDs)

         if (measure == "SMD") {

            ### apply bias-correction to di values

            cmi <- .cmicalc(mi)
            #cmi <- 1
            yi <- cmi * di

            if (length(vtype) == 1L)
               vtype <- rep(vtype, k)

            vi <- rep(NA_real_, k)

            mnwyi <- sum(ni*yi, na.rm=TRUE) / sum(ni, na.rm=TRUE) ### sample size weighted average of yi's

            for (i in seq_len(k)) {

               ### unbiased estimate of the sampling variance
               if (vtype[i] == "UB")
                  vi[i] <- 1/n1i[i] + 1/n2i[i] + (1 - (mi[i]-2)/(mi[i]*cmi[i]^2)) * yi[i]^2

               ### large sample approximation to the sampling variance
               if (vtype[i] == "LS")
                  vi[i] <- 1/n1i[i] + 1/n2i[i] + yi[i]^2/(2*ni[i])

               ### estimator assuming homogeneity (using sample size weighted average of the yi's)
               if (vtype[i] == "HO")
                  vi[i] <- 1/n1i[i] + 1/n2i[i] + mnwyi^2/(2*ni[i])

            }

         }

         ### standardized mean difference (with heteroscedastic SDs)

         if (measure == "SMDH") {
            cmi <- .cmicalc(mi)
            si <- sqrt((sd1i^2 + sd2i^2)/2)
            yi <- cmi * (m1i - m2i) / si
            vi <- yi^2 * (sd1i^4 / (n1i-1) + sd2i^4 / (n2i-1)) / (2*(sd1i^2 + sd2i^2)^2) + (sd1i^2 / (n1i-1) + sd2i^2 / (n2i-1)) / ((sd1i^2 + sd2i^2)/2)
            vi <- cmi^2 * vi
            ### note: Bonett (2009) plugs in the uncorrected yi into the
            ### equation for vi; here, the corrected value is plugged in
         }

         ### ratio of means (response ratio)
         ### to use with pooled SDs, simply set sd1i = sd2i = sdpi or use vtype="HO"

         if (measure == "ROM") {

            yi <- log(m1i/m2i)

            if (length(vtype) == 1L)
               vtype <- rep(vtype, k)

            vi <- rep(NA_real_, k)

            for (i in seq_len(k)) {

               ### large sample approximation to the sampling variance (does not assume homoscedasticity)
               if (vtype[i] == "LS")
                  vi[i] <- sd1i[i]^2/(n1i[i]*m1i[i]^2) + sd2i[i]^2/(n2i[i]*m2i[i]^2)

               ### estimator assuming homoscedasticity
               if (vtype[i] == "HO")
                  vi[i] <- sdpi[i]^2/(n1i[i]*m1i[i]^2) + sdpi[i]^2/(n2i[i]*m2i[i]^2)

            }

         }

         ### point-biserial correlation obtained from the standardized mean difference
         ### this is based on Tate's model where Y|X=0 and Y|X=1 are normally distributed (with the same variance)
         ### Das Gupta (1960) describes the case where Y itself is normal, but the variance expressions therein can
         ### really only be used in some special cases (not useful in practice)

         if (is.element(measure, c("RPB","RBIS"))) {

            hi <- mi/n1i + mi/n2i
            yi <- di / sqrt(di^2 + hi) ### need this also when measure="RBIS"

            if (measure == "RPB") {    ### this only applies when measure="RPB"

               if (length(vtype) == 1L)
                  vtype <- rep(vtype, k)

               vi <- rep(NA_real_, k)

               for (i in seq_len(k)) {

                  ### estimate of the sampling variance for fixed n1i and n2i (i.e., stratified sampling)
                  if (vtype[i] == "ST" || vtype[i] == "LS")
                     vi[i] <- hi[i]^2 / (hi[i] + di[i]^2)^3 * (1/n1i[i] + 1/n2i[i] + di[i]^2/(2*ni[i]))

                  ### estimate of the sampling variance for fixed ni but random n1i and n2i (i.e., cross-sectional/multinomial sampling)
                  if (vtype[i] == "CS")
                     vi[i] <- (1-yi[i]^2)^2 * (ni[i]*yi[i]^2 / (4*n1i[i]*n2i[i]) + (2-3*yi[i]^2)/(2*ni[i])) ### from Tate (1954, 1955b)

               }

            }

         }

         ### biserial correlation obtained from the standardized mean difference (continued from above)

         if (measure == "RBIS") {
            p1i <- n1i / ni
            p2i <- n2i / ni
            zi  <- qnorm(p1i, lower.tail=FALSE)
            fzi <- dnorm(zi)
            yi  <- sqrt(p1i*p2i) / fzi * yi ### yi on the right-hand side is the point-biserial correlation from above
            #vi <- (p1i*p2i) / fzi^2 * vi   ### not correct (p1i, p2i, and fzi are random variables and vi from RBP is not correct for the bivariate normal case on which RBIS is based)
            yi.t <- ifelse(abs(yi) > 1, sign(yi), yi)
            vi  <- 1/(ni-1) * (p1i*p2i/fzi^2 - (3/2 + (1 - p1i*zi/fzi)*(1 + p2i*zi/fzi)) * yi.t^2 + yi.t^4) ### from Soper (1914)
            #vi <- 1/(ni-1) * (yi.t^4 + yi.t^2 * (p1i*p2i*zi^2/fzi^2 + (2*p1i-1)*zi/fzi - 5/2) + p1i*p2i/fzi^2) ### from Tate (1955) -- equivalent to eq. from Soper (1914)
            ### equation appears to work even if dichotomization is done based on a sample quantile value (so that p1i, p2i, and fzi are fixed by design)
         }

         ### SMD to log(OR) transformation based on logistic distribution

         if (is.element(measure, c("D2OR","D2ORL"))) {
            yi <- pi / sqrt(3) * di
            vi <- pi^2 / 3 * (1/n1i + 1/n2i + di^2/(2*ni))
         }

         ### SMD to log(OR) transformation based on normal distribution (Cox & Snell method)

         if (measure == "D2ORN") {
            yi <- 1.65 * di
            vi <- 1.65^2 * (1/n1i + 1/n2i + di^2/(2*ni))
         }

      }

      ######################################################################

      if (is.element(measure, c("COR","UCOR","ZCOR"))) {

         mf.ri   <- mf[[match("ri", names(mf))]]
         mf.ni   <- mf[[match("ni", names(mf))]]
         ri      <- eval(mf.ri, data, enclos=sys.frame(sys.parent()))
         ni      <- eval(mf.ni, data, enclos=sys.frame(sys.parent()))

         if (!is.null(subset)) {
            ri <- ri[subset]
            ni <- ni[subset]
         }

         if (length(ri)==0L || length(ni)==0L)
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")

         if (length(ri) != length(ni))
            stop("Supplied data vectors are not of the same length.")

         if (any(abs(ri) > 1, na.rm=TRUE))
            stop("One or more correlations are > 1 or < -1.")

         if (any(ni < 0, na.rm=TRUE))
            stop("One or more sample sizes are negative.")

         if (any(ni <= 4, na.rm=TRUE))
            warning("Cannot estimate sampling variance when ni <= 4.")

         ni.u <- ni ### unadjusted total sample sizes

         k <- length(ri)

         ### raw correlation

         if (measure == "COR") {
            yi <- ri
         }

         ### raw correlation with bias correction

         if (measure == "UCOR") {
            yi <- ri + ri*(1-ri^2)/(2*(ni-4))
            yi[ni <= 4] <- NA ### set corrected correlations for ni <= 4 to NA
         }

         ### sampling variances for COR or UCOR

         if (is.element(measure, c("COR","UCOR"))) {

            if (length(vtype) == 1L)
               vtype <- rep(vtype, k)

            vi <- rep(NA_real_, k)

            mnwyi <- sum(ni*yi, na.rm=TRUE) / sum(ni, na.rm=TRUE) ### sample size weighted average of yi's

            for (i in seq_len(k)) {

               ### approximation to the unbiased estimate of the sampling variance
               if (vtype[i] == "UB") {
                  vi[i] <- yi[i]^2 - 1 + (ni[i]-3)/(ni[i]-2) * ((1-ri[i]^2) + 2*(1-ri[i]^2)^2/ni[i] + 8*(1-ri[i]^2)^3/(ni[i]*(ni[i]+2)) + 48*(1-ri[i]^2)^4/(ni[i]*(ni[i]+2)*(ni[i]+4)))
               }

               ### large sample approximation to the sampling variance
               if (vtype[i] == "LS") {
                  vi[i] <- (1-yi[i]^2)^2/(ni[i]-1)
               }

               ### estimator assuming homogeneity (using sample size weighted average of the yi's)
               if (vtype[i] == "HO")
                  vi[i] <- (1-mnwyi^2)^2/(ni[i]-1)

            }

         }

         ### r-to-z transformed correlation

         if (measure == "ZCOR") {
            yi <- 1/2 * log((1+ri)/(1-ri))
            vi <- 1/(ni-3)
         }

         ### set sampling variances for ni <= 4 to NA

         vi[ni <= 4] <- NA

      }

      ######################################################################

      if (is.element(measure, c("PR","PLN","PLO","PAS","PFT"))) {

         mf.xi   <- mf[[match("xi", names(mf))]]
         mf.mi   <- mf[[match("mi", names(mf))]]
         mf.ni   <- mf[[match("ni", names(mf))]]
         xi      <- eval(mf.xi, data, enclos=sys.frame(sys.parent()))
         mi      <- eval(mf.mi, data, enclos=sys.frame(sys.parent()))
         ni      <- eval(mf.ni, data, enclos=sys.frame(sys.parent()))
         if (is.null(mi)) mi <- ni - xi

         if (!is.null(subset)) {
            xi <- xi[subset]
            mi <- mi[subset]
         }

         if (length(xi)==0L || length(mi)==0L)
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")

         if (length(xi) != length(mi))
            stop("Supplied data vectors are not all of the same length.")

         if (any(c(xi, mi) < 0, na.rm=TRUE))
            stop("One or more counts are negative.")

         ni.u <- xi + mi ### unadjusted total sample sizes

         k <- length(xi)

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

         ni  <- xi + mi
         pri <- xi/ni

         ### raw proportion

         if (measure == "PR") {

            yi <- pri

            if (length(vtype) == 1L)
               vtype <- rep(vtype, k)

            vi <- rep(NA_real_, k)

            mnwpri <- sum(xi, na.rm=TRUE) / sum(ni, na.rm=TRUE) ### sample size weighted average of proportions (same as sum(ni*pi)/sum(ni))

            for (i in seq_len(k)) {

               ### unbiased estimate of the sampling variance
               if (vtype[i] == "UB")
                  vi[i] <- pri[i]*(1-pri[i])/(ni[i]-1)

               ### large sample approximation to the sampling variance
               if (vtype[i] == "LS")
                  vi[i] <- pri[i]*(1-pri[i])/ni[i]

               ### estimator assuming homogeneity (using the average proportion)
               if (vtype[i] == "HO")
                  vi[i] <- mnwpri*(1-mnwpri)/ni[i]

            }

         }

         ### proportion with log transformation

         if (measure == "PLN") {

            yi <- log(pri)

            if (length(vtype) == 1L)
               vtype <- rep(vtype, k)

            vi <- rep(NA_real_, k)

            mnwpri <- sum(xi, na.rm=TRUE) / sum(ni, na.rm=TRUE) ### sample size weighted average of proportions (same as sum(ni*pi)/sum(ni))
            #mnwpri <- exp(sum(ni*yi)/sum(ni))                  ### alternative strategy (exp of the sample size weighted average of the log proportions)

            for (i in seq_len(k)) {

               ### large sample approximation to the sampling variance
               if (vtype[i] == "LS")
                  vi[i] <- 1/xi[i] - 1/ni[i]

               ### estimator assuming homogeneity (using the average proportion)
               if (vtype[i] == "HO")
                  vi[i] <- 1/(mnwpri*ni[i]) - 1/ni[i]

            }

         }

         ### proportion with logit (log odds) transformation

         if (measure == "PLO") {

            yi <- log(pri/(1-pri))

            if (length(vtype) == 1L)
               vtype <- rep(vtype, k)

            vi <- rep(NA_real_, k)

            mnwpri <- sum(xi, na.rm=TRUE) / sum(ni, na.rm=TRUE) ### sample size weighted average of proportions (same as sum(ni*pi)/sum(ni))
            #mnwpri <- transf.ilogit(sum(ni*yi)/sum(ni))        ### alternative strategy (inverse logit of the sample size weighted average of the logit transformed proportions)

            for (i in seq_len(k)) {

               ### large sample approximation to the sampling variance
               if (vtype[i] == "LS")
                  vi[i] <- 1/xi[i] + 1/mi[i]

               ### estimator assuming homogeneity (using the average proportion)
               if (vtype[i] == "HO")
                  vi[i] <- 1/(mnwpri*ni[i]) + 1/((1-mnwpri)*ni[i])

            }

         }

         ### proportion with arcsine square-root (angular) transformation

         if (measure == "PAS") {
            yi <- asin(sqrt(pri))
            vi <- 1/(4*ni)
         }

         ### proportion with Freeman-Tukey double arcsine transformation

         if (measure == "PFT") {
            yi <- 1/2*(asin(sqrt(xi/(ni+1))) + asin(sqrt((xi+1)/(ni+1))))
            vi <- 1/(4*ni + 2)
         }

      }

      ######################################################################

      if (is.element(measure, c("IR","IRLN","IRS","IRFT"))) {

         mf.xi   <- mf[[match("xi", names(mf))]]
         mf.ti   <- mf[[match("ti", names(mf))]]
         xi      <- eval(mf.xi, data, enclos=sys.frame(sys.parent()))
         ti      <- eval(mf.ti, data, enclos=sys.frame(sys.parent()))

         if (!is.null(subset)) {
            xi <- xi[subset]
            ti <- ti[subset]
         }

         if (length(xi)==0L || length(ti)==0L)
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")

         if (length(xi) != length(ti))
            stop("Supplied data vectors are not all of the same length.")

         if (any(xi < 0, na.rm=TRUE))
            stop("One or more counts are negative.")

         if (any(ti < 0, na.rm=TRUE))
            stop("One or more person-times are negative.")

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

         iri <- xi/ti

         ### raw incidence rate

         if (measure == "IR") {
            yi <- iri
            vi <- iri/ti ### note: same as xi/ti^2
         }

         ### log transformed incidence rate

         if (measure == "IRLN") {
            yi <- log(iri)
            vi <- 1/xi
         }

         ### square-root transformed incidence rate

         if (measure == "IRS") {
            yi <- sqrt(iri)
            vi <- 1/(4*ti)
         }

         ### incidence rate with Freeman-Tukey transformation

         if (measure == "IRFT") {
            yi <- 1/2*(sqrt(iri) + sqrt(iri+1/ti))
            vi <- 1/(4*ti)
         }

      }

      ######################################################################

      if (is.element(measure, c("MN"))) {

         mf.mi   <- mf[[match("mi",  names(mf))]]
         mf.sdi  <- mf[[match("sdi", names(mf))]]
         mf.ni   <- mf[[match("ni",  names(mf))]]
         mi      <- eval(mf.mi,  data, enclos=sys.frame(sys.parent()))
         sdi     <- eval(mf.sdi, data, enclos=sys.frame(sys.parent()))
         ni      <- eval(mf.ni,  data, enclos=sys.frame(sys.parent()))

         if (!is.null(subset)) {
            mi  <- mi[subset]
            sdi <- sdi[subset]
            ni  <- ni[subset]
         }

         if (length(mi)==0L || length(sdi)==0L || length(ni)==0L)
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")

         if (!all(length(mi) == c(length(mi),length(sdi),length(ni))))
            stop("Supplied data vectors are not all of the same length.")

         if (any(sdi < 0, na.rm=TRUE))
            stop("One or more standard deviations are negative.")

         if (any(ni < 0, na.rm=TRUE))
            stop("One or more sample sizes are negative.")

         ni.u <- ni ### unadjusted total sample sizes

         ### (raw) means

         if (measure == "MN") {
            yi <- mi
            vi <- sdi^2/ni
         }

      }

      ######################################################################

      if (is.element(measure, c("MC","SMCC","SMCR","SMCRH","ROMC"))) {

         mf.m1i  <- mf[[match("m1i",  names(mf))]]
         mf.m2i  <- mf[[match("m2i",  names(mf))]]
         mf.sd1i <- mf[[match("sd1i", names(mf))]]
         mf.sd2i <- mf[[match("sd2i", names(mf))]] ### for SMCR, do not need to supply this
         mf.ni   <- mf[[match("ni",   names(mf))]]
         mf.ri   <- mf[[match("ri",   names(mf))]]
         m1i     <- eval(mf.m1i,  data, enclos=sys.frame(sys.parent()))
         m2i     <- eval(mf.m2i,  data, enclos=sys.frame(sys.parent()))
         sd1i    <- eval(mf.sd1i, data, enclos=sys.frame(sys.parent()))
         sd2i    <- eval(mf.sd2i, data, enclos=sys.frame(sys.parent()))
         ni      <- eval(mf.ni,   data, enclos=sys.frame(sys.parent()))
         ri      <- eval(mf.ri,   data, enclos=sys.frame(sys.parent()))

         if (!is.null(subset)) {
            m1i  <- m1i[subset]
            m2i  <- m2i[subset]
            sd1i <- sd1i[subset]
            sd2i <- sd2i[subset]
            ni   <- ni[subset]
            ri   <- ri[subset]
         }

         if (is.element(measure, c("MC","SMCC","SMCRH","ROMC"))) {

            ### for these measures, need m1i, m2i, sd1i, sd2i, ni, and ri

            if (length(m1i)==0L || length(m2i)==0L || length(sd1i)==0L || length(sd2i)==0L || length(ni)==0L || length(ri)==0L)
               stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")

            if (!all(length(m1i) == c(length(m1i),length(m2i),length(sd1i),length(sd2i),length(ni),length(ri))))
               stop("Supplied data vectors are not all of the same length.")

            if (any(c(sd1i, sd2i) < 0, na.rm=TRUE))
               stop("One or more standard deviations are negative.")

         } else {

            ### for SMCR, need m1i, m2i, sd1i, ni, and ri (do not need sd2i!)

            if (length(m1i)==0L || length(m2i)==0L || length(sd1i)==0L || length(ni)==0L || length(ri)==0L)
               stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")

            if (!all(length(m1i) == c(length(m1i),length(m2i),length(sd1i),length(ni),length(ri))))
               stop("Supplied data vectors are not all of the same length.")

            if (any(sd1i < 0, na.rm=TRUE))
               stop("One or more standard deviations are negative.")

         }

         if (any(abs(ri) > 1, na.rm=TRUE))
            stop("One or more correlations are > 1 or < -1.")

         if (any(ni < 0, na.rm=TRUE))
            stop("One or more sample sizes are negative.")

         ni.u <- ni ### unadjusted total sample sizes

         ni <- ni.u
         mi <- ni - 1

         ### (raw) mean change

         if (measure == "MC") {
            yi <- m1i - m2i
            vi <- (sd1i^2 + sd2i^2 - 2*ri*sd1i*sd2i) / ni
         }

         ### standardized mean change with change score standardization (using sddi)
         ### note: does not assume homoscedasticity, since we use sddi here

         if (measure == "SMCC") {
            cmi <- .cmicalc(mi)
            sddi <- sqrt(sd1i^2 + sd2i^2 - 2*ri*sd1i*sd2i)
            yi <- cmi * (m1i - m2i) / sddi
            vi <- 1/ni + yi^2 / (2*ni)
         }

         ### standardized mean change with raw score standardization (using sd1i)
         ### note: yi does not assume homoscedasticity, but vi does

         if (measure == "SMCR") {
            cmi <- .cmicalc(mi)
            yi <- cmi * (m1i - m2i) / sd1i
            vi <- 2*(1-ri)/ni + yi^2 / (2*ni)
         }

         ### standardized mean change with raw score standardization (using sd1i)
         ### with vi computation allowing for heteroscedasticity (Bonett, 2008; and JEBS article)

         if (measure == "SMCRH") {
            cmi <- .cmicalc(mi)
            vardi <- sd1i^2 + sd2i^2 - 2*ri*sd1i*sd2i
            yi <- cmi * (m1i - m2i) / sd1i
            vi <- vardi/(sd1i^2*(ni-1)) + yi^2 / (2*(ni-1))
            vi <- cmi^2 * vi
            ### note: Bonett suggests plugging in the uncorrected yi into the
            ### equation for vi; here, the corrected value is plugged in
         }

         ### ratio of means for pre-post or matched designs (eq. 6 in Lajeunesse, 2011)
         ### to use with pooled SDs, simply set sd1i = sd2i = sdpi!

         if (measure == "ROMC") {
            yi <- log(m1i/m2i)
            vi <- sd1i^2/(ni*m1i^2) + sd2i^2/(ni*m2i^2) - 2*ri*sd1i*sd2i/(m1i*m2i*ni)
         }

      }

      ######################################################################

      if (is.element(measure, c("ARAW","AHW","ABT"))) {

         mf.ai <- mf[[match("ai", names(mf))]]
         mf.mi <- mf[[match("mi", names(mf))]]
         mf.ni <- mf[[match("ni", names(mf))]]
         ai    <- eval(mf.ai, data, enclos=sys.frame(sys.parent()))
         mi    <- eval(mf.mi, data, enclos=sys.frame(sys.parent()))
         ni    <- eval(mf.ni, data, enclos=sys.frame(sys.parent()))

         if (!is.null(subset)) {
            ai <- ai[subset]
            mi <- mi[subset]
            ni <- ni[subset]
         }

         if (length(ai)==0L || length(mi)==0L || length(ni)==0L)
            stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")

         if (!all(length(ai) == c(length(ai),length(mi),length(ni))))
            stop("Supplied data vectors are not all of the same length.")

         if (any(ai > 1, na.rm=TRUE))
            stop("One or more alpha values are > 1.")

         if (any(mi < 2, na.rm=TRUE))
            stop("One or more mi values are < 2.")

         if (any(ni < 0, na.rm=TRUE))
            stop("One or more sample sizes are negative.")

         ni.u <- ni ### unadjusted total sample sizes

         ### raw alpha values

         if (measure == "ARAW") {
            yi <- ai
            vi <- 2*mi*(1-ai)^2 / ((mi-1)*(ni-2))
         }

         ### alphas transformed with Hakstian & Whalen (1976) transformation

         if (measure == "AHW") {
            #yi <- (1-ai)^(1/3) ### technically this is the Hakstian & Whalen (1976) transformation
            yi <- 1 - (1-ai)^(1/3) ### but with this, yi remains a monotonically increasing function of ai
            vi <- 18*mi*(ni-1)*(1-ai)^(2/3) / ((mi-1)*(9*ni-11)^2)
         }

         ### alphas transformed with Bonett (2002) transformation (without bias correction)

         if (measure == "ABT") {
            #yi <- log(1-ai) - log(ni/(ni-1))
            #yi <- log(1-ai) ### technically this is the Bonett (2002) transformation
            yi <- -log(1-ai) ### but with this, yi remains a monotonically increasing function of ai
            vi <- 2*mi / ((mi-1)*(ni-2))
         }

      }

      ######################################################################

   } else {

      ### in case yi is not NULL (so user wants to convert a regular data frame to an 'escalc' object)

      ### get vi, sei, and ni

      mf.vi  <- mf[[match("vi",  names(mf))]]
      mf.sei <- mf[[match("sei", names(mf))]]
      mf.ni  <- mf[[match("ni", names(mf))]]
      vi     <- eval(mf.vi,  data, enclos=sys.frame(sys.parent()))
      sei    <- eval(mf.sei, data, enclos=sys.frame(sys.parent()))
      ni     <- eval(mf.ni,  data, enclos=sys.frame(sys.parent()))

      ### if vi is specified, this will be used (even if user has specified sei as well)
      ### otherwise, if user has specified sei, then square those values to get vi
      ### if neither is specified, then throw an error

      if (is.null(vi)) {
         if (is.null(sei)) {
            stop("Need to specify 'vi' or 'sei' argument.")
         } else {
            vi <- sei^2
         }
      }

      if (!is.null(subset)) {
         yi <- yi[subset]
         vi <- vi[subset]
         ni <- ni[subset]
      }

      if (length(yi) != length(vi))
         stop("Supplied data vectors are not of the same length.")

      if (!is.null(ni) && (length(yi) != length(ni)))
         stop("Supplied data vectors are not of the same length.")

      ni.u <- ni ### unadjusted total sample sizes

   }

   #########################################################################
   #########################################################################
   #########################################################################

   ### check for infinite values and set them to NA

   is.inf <- is.infinite(yi) | is.infinite(vi)

   if (any(is.inf)) {
      warning("Some 'yi' and/or 'vi' values equal to +-Inf. Recoded to NAs.")
      yi[is.inf] <- NA
      vi[is.inf] <- NA
   }

   ### check for NaN values and set them to NA

   is.NaN <- is.nan(yi) | is.nan(vi)

   if (any(is.NaN)) {
      yi[is.NaN] <- NA
      vi[is.NaN] <- NA
   }

   ### check for negative vi's (should not happen, but just in case)

   vi[vi < 0] <- NA

   ### add study labels if specified

   if (!is.null(slab)) {

      if (!is.null(subset))
         slab <- slab[subset]

      if (anyNA(slab))
         stop("NAs in study labels.")

      ### check if study labels are unique; if not, make them unique

      if (anyDuplicated(slab))
         slab <- make.unique(as.character(slab)) ### make.unique() only works with character vectors

      if (length(slab) != length(yi))
         stop("Study labels not of same length as data.")

      ### add slab attribute to the yi vector
      attr(yi, "slab") <- slab

   }

   ### if a subset of studies is specified (note: subsetting of other parts already done above, so yi/vi/ni.u/slab are already subsetted)

   if (!is.null(subset)) {
      if (!no.data)
         data <- data[subset,,drop=FALSE]
   }

   ### add measure attribute to the yi vector

   attr(yi, "measure") <- measure

   ### put together dataset

   if (!no.data && append) {

      ### if data argument has been specified and user wants to append

      dat <- data.frame(data)

      if (replace) {

         ### and wants to replace all values

         dat[[var.names[1]]] <- yi ### if yi variable does not exists in dat, it will be added; otherwise it will be overwritten
         dat[[var.names[2]]] <- vi ### if vi variable does not exists in dat, it will be added; otherwise it will be overwritten

         if (add.measure) {
            dat[[var.names[3]]] <- ""
            dat[[var.names[3]]][!is.na(yi)] <- measure
         }

         attr(dat[[var.names[1]]], "ni") <- ni.u

      } else {

         ### and only wants to replace any NA values

         if (is.element(var.names[1], names(dat))) { ### if yi variable is in data frame, replace NA values with newly calculated values
            is.na.yi <- is.na(dat[[var.names[1]]])
            dat[[var.names[1]]][is.na.yi] <- yi[is.na.yi]
            attributes(dat[[var.names[1]]])$ni[is.na.yi] <- ni.u[is.na.yi]
         } else {
            dat[[var.names[1]]] <- yi                ### if yi variable does not exist in dat, just add as new variable
            attr(dat[[var.names[1]]], "ni") <- ni.u
         }

         if (is.element(var.names[2], names(dat))) { ### if vi variable is in data frame, replace NA values with newly calculated values
            is.na.vi <- is.na(dat[[var.names[2]]])
            dat[[var.names[2]]][is.na.vi] <- vi[is.na.vi]
         } else {
            dat[[var.names[2]]] <- vi                ### if vi variable does not exist in dat, just add as new variable
         }

         if (add.measure) {
            if (is.element(var.names[3], names(dat))) {    ### if measure variable is in data frame, replace NA values with newly calculated values
               is.na.measure <- c(dat[[var.names[3]]] == "") & !is.na(yi)
               dat[[var.names[3]]][is.na.measure] <- measure
            } else {
               dat[[var.names[3]]] <- ""                   ### if measure variable does not exist in dat, just add as new variable
               dat[[var.names[3]]][!is.na(yi)] <- measure
            }
         }

      }

   } else {

      ### if data argument has not been specified or user does not want to append

      if (add.measure) {
         dat <- data.frame(yi, vi)
         dat$measure <- ""
         dat$measure[!is.na(yi)] <- measure
         names(dat) <- var.names
      } else {
         dat <- data.frame(yi, vi)
         names(dat) <- var.names[1:2]
      }

      attr(dat[,1], "ni") <- ni.u

   }

   attr(dat, "digits") <- digits

   ### add 'yi.names' and 'vi.names' to the first position of the corresponding attributes
   attr(dat, "yi.names") <- unique(c(var.names[1], attr(data, "yi.names"))) ### if 'yi.names' is not an attribute, attr() returns NULL, so this works fine
   attr(dat, "vi.names") <- unique(c(var.names[2], attr(data, "vi.names"))) ### if 'vi.names' is not an attribute, attr() returns NULL, so this works fine

   ### add 'out.names' back to object in case these attributes exist (if summary() has been used on the object)
   attr(dat, "sei.names")   <- attr(data, "sei.names")
   attr(dat, "zi.names")    <- attr(data, "zi.names")
   attr(dat, "ci.lb.names") <- attr(data, "ci.lb.names")
   attr(dat, "ci.ub.names") <- attr(data, "ci.ub.names")

   ### TODO: clean up attribute elements that are no longer actually part of the object

   class(dat) <- c("escalc", "data.frame")
   return(dat)

}
