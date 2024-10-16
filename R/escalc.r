escalc <- function(measure, ai, bi, ci, di, n1i, n2i, x1i, x2i, t1i, t2i, m1i, m2i, sd1i, sd2i, xi, mi, ri, ti, fi, pi, sdi, r2i, ni, yi, vi, sei,
data, slab, flip, subset, include, add=1/2, to="only0", drop00=FALSE, vtype="LS", correct=TRUE, var.names=c("yi","vi"), add.measure=FALSE, append=TRUE, replace=TRUE, digits, ...) {

   ### check argument specifications

   mstyle <- .get.mstyle()

   if (missing(measure) && missing(yi))
      stop(mstyle$stop("Must specify an effect size or outcome measure via the 'measure' argument."))

   if (!missing(yi) && missing(measure))
      measure <- "GEN"

   if (!is.character(measure))
      stop(mstyle$stop("The 'measure' argument must be a character string."))

   if (!is.element(measure, c("RR","OR","PETO","RD","AS","PHI","ZPHI","YUQ","YUY","RTET","ZTET",                # 2x2 table measures
                              "PBIT","OR2D","OR2DN","OR2DL",                                                    # 2x2 table transformations to SMDs
                              "MPRD","MPRR","MPOR","MPORC","MPPETO","MPORM",                                    # 2x2 table measures for matched pairs / pre-post data
                              "IRR","IRD","IRSD",                                                               # two-group person-time data (incidence) measures
                              "MD","SMD","SMDH","SMD1","SMD1H","ROM",                                           # two-group mean/SD measures
                              "CVR","VR",                                                                       # coefficient of variation ratio, variability ratio
                              "RPB","ZPB","RBIS","ZBIS","D2OR","D2ORN","D2ORL",                                 # two-group mean/SD transformations to r_pb, r_bis, and log(OR)
                              "COR","UCOR","ZCOR",                                                              # correlations (raw and r-to-z transformed)
                              "PCOR","ZPCOR","SPCOR","ZSPCOR",                                                  # partial and semi-partial correlations
                              "R2","ZR2",                                                                       # coefficient of determination / R^2 (raw and r-to-z transformed)
                              "PR","PLN","PLO","PRZ","PAS","PFT",                                               # single proportions (and transformations thereof)
                              "IR","IRLN","IRS","IRFT",                                                         # single-group person-time (incidence) data (and transformations thereof)
                              "MN","SMN","MNLN","CVLN","SDLN",                                                  # mean, single-group standardized mean, log(mean), log(CV), log(SD),
                              "MC","SMCC","SMCR","SMCRH","SMCRP","SMCRPH","CLESCN","AUCCN","ROMC","CVRC","VRC", # raw/standardized mean change, CLES/AUC, log(ROM), CVR, and VR for dependent samples
                              "ARAW","AHW","ABT",                                                               # alpha (and transformations thereof)
                              "REH","CLES","CLESN","AUC","AUCN",                                                # relative excess heterozygosity, common language effect size / area under the curve
                              "HR","HD",                                                                        # hazard (rate) ratios and differences
                              "GEN")))
      stop(mstyle$stop("Unknown 'measure' specified."))

   # when adding measures, remember to add measures to .setlab()

   if (!is.element(to, c("all","only0","if0all","none")))
      stop(mstyle$stop("Unknown 'to' argument specified."))

   if (any(!is.element(vtype, c("UB","LS","LS2","LS3","HO","ST","CS","AV","AV2","AVHO","H0","MAX")), na.rm=TRUE)) # vtype can be an entire vector, so use any() and na.rm=TRUE
      stop(mstyle$stop("Unknown 'vtype' argument specified."))

   if (add.measure) {

      if (length(var.names) == 2L)
         var.names <- c(var.names, "measure")

      if (length(var.names) != 3L)
         stop(mstyle$stop("Argument 'var.names' must be of length 2 or 3."))

      if (any(var.names != make.names(var.names, unique=TRUE))) {
         var.names <- make.names(var.names, unique=TRUE)
         warning(mstyle$warning(paste0("Argument 'var.names' does not contain syntactically valid variable names.\nVariable names adjusted to: var.names = c('", var.names[1], "','", var.names[2], "','", var.names[3], "').")), call.=FALSE)
      }

   } else {

      if (length(var.names) == 3L)
         var.names <- var.names[1:2]

      if (length(var.names) != 2L)
         stop(mstyle$stop("Argument 'var.names' must be of length 2."))

      if (any(var.names != make.names(var.names, unique=TRUE))) {
         var.names <- make.names(var.names, unique=TRUE)
         warning(mstyle$warning(paste0("Argument 'var.names' does not contain syntactically valid variable names.\nVariable names adjusted to: var.names = c('", var.names[1], "','", var.names[2], "').")), call.=FALSE)
      }

   }

   ### check if user is trying to use the 'formula interface' to escalc()
   ### note: if so, argument 'ai' may mistakenly be a formula, so check for that as well (further below)

   if (hasArg(formula) || hasArg(weights))
      stop(mstyle$stop("The 'formula interface' to escalc() has been deprecated."))

   ### get ... argument and check for extra/superfluous arguments

   ddd <- list(...)

   .chkdots(ddd, c("onlyo1", "addyi", "addvi"))

   ### set defaults or get onlyo1, addyi, and addvi arguments

   onlyo1  <- .chkddd(ddd$onlyo1,  FALSE, .isTRUE(ddd$onlyo1))
   addyi   <- .chkddd(ddd$addyi,   TRUE,  .isTRUE(ddd$addyi))
   addvi   <- .chkddd(ddd$addvi,   TRUE,  .isTRUE(ddd$addvi))

   ### set defaults for digits

   if (missing(digits)) {
      digits <- .set.digits(dmiss=TRUE)
   } else {
      digits <- .set.digits(digits, dmiss=FALSE)
   }

   ### check if data argument has been specified

   if (missing(data))
      data <- NULL

   ### need this at the end to check if append=TRUE can actually be done

   has.data <- !is.null(data)

   ### check if data argument has been specified

   if (is.null(data)) {
      data <- sys.frame(sys.parent())
   } else {
      if (!is.data.frame(data))
         data <- data.frame(data)
   }

   mf <- match.call()

   ### get slab, subset, and include arguments (NULL when unspecified)

   slab    <- .getx("slab",    mf=mf, data=data)
   subset  <- .getx("subset",  mf=mf, data=data)
   include <- .getx("include", mf=mf, data=data)

   ### get yi (in case it has been specified)

   yi <- .getx("yi", mf=mf, data=data)

   ### get flip (NULL if not specified)

   flip <- .getx("flip", mf=mf, data=data)

   ### for certain measures, set add=0 by default unless user explicitly set the add argument

   addval <- mf[[match("add", names(mf))]]

   if (is.element(measure, c("AS","PHI","ZPHI","RTET","ZTET","IRSD","PAS","PFT","IRS","IRFT")) && is.null(addval))
      add <- 0

   #########################################################################
   #########################################################################
   #########################################################################

   if (is.null(yi)) {

      if (is.element(measure, c("RR","OR","RD","AS","PETO","PHI","ZPHI","YUQ","YUY","RTET","ZTET","PBIT","OR2D","OR2DN","OR2DL","MPRD","MPRR","MPOR","MPORC","MPPETO","MPORM"))) {

         mf.ai <- mf[[match("ai", names(mf))]]
         if (any("~" %in% as.character(mf.ai)))
            stop(mstyle$stop("The 'formula interface' to escalc() has been deprecated."))

         ai  <- .getx("ai",  mf=mf, data=data, checknumeric=TRUE)
         bi  <- .getx("bi",  mf=mf, data=data, checknumeric=TRUE)
         ci  <- .getx("ci",  mf=mf, data=data, checknumeric=TRUE)
         di  <- .getx("di",  mf=mf, data=data, checknumeric=TRUE)
         n1i <- .getx("n1i", mf=mf, data=data, checknumeric=TRUE)
         n2i <- .getx("n2i", mf=mf, data=data, checknumeric=TRUE)
         ri  <- .getx("ri",  mf=mf, data=data, checknumeric=TRUE)
         pi  <- .getx("pi",  mf=mf, data=data, checknumeric=TRUE)

         if (!.equal.length(ai, bi, ci, di, n1i, n2i, ri, pi))
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

         if (measure == "MPORM" && !(.all.specified(ri) || .all.specified(pi)))
            stop(mstyle$stop("Need to specify also argument 'ri' (and/or 'pi') for this measure."))

         k.all <- length(ai)

         if (!is.null(subset)) {
            subset <- .chksubset(subset, k.all)
            ai <- .getsubset(ai, subset)
            bi <- .getsubset(bi, subset)
            ci <- .getsubset(ci, subset)
            di <- .getsubset(di, subset)
            ri <- .getsubset(ri, subset)
            pi <- .getsubset(pi, subset)
         }

         n1i <- ai + bi
         n2i <- ci + di

         if (any(c(ai > n1i, ci > n2i), na.rm=TRUE))
            stop(mstyle$stop("One or more event counts are larger than the corresponding group sizes."))

         if (any(c(ai, bi, ci, di) < 0, na.rm=TRUE))
            stop(mstyle$stop("One or more counts are negative."))

         if (any(c(n1i < 0, n2i < 0), na.rm=TRUE)) # note: in cross-sectional sampling, group sizes could be 0
            stop(mstyle$stop("One or more group sizes are negative."))

         if (measure == "MPORM" && !is.null(ri) && any(abs(ri) > 1, na.rm=TRUE))
            stop(mstyle$stop("One or more correlations are > 1 or < -1."))

         if (measure == "MPORM" && !is.null(pi) && any(pi < 0 | pi > 1, na.rm=TRUE))
            stop(mstyle$stop("One or more proportions are > 1 or < 0."))

         ni.u <- ai + bi + ci + di # unadjusted total sample sizes

         if (measure == "MPORM")
            ni.u <- round(ni.u / 2)

         k <- length(ai)

         ### if drop00=TRUE, set counts to NA for studies that have no events (or all events) in both arms

         if (drop00) {
            id00 <- c(ai == 0L & ci == 0L) | c(bi == 0L & di == 0L)
            id00[is.na(id00)] <- FALSE
            ai[id00] <- NA_real_
            bi[id00] <- NA_real_
            ci[id00] <- NA_real_
            di[id00] <- NA_real_
         }

         ### save unadjusted counts

         ai.u <- ai
         bi.u <- bi
         ci.u <- ci
         di.u <- di
         n1i.u <- ai + bi
         n2i.u <- ci + di

         if (to == "all") {

            ### always add to all cells in all studies

            ai <- ai + add
            ci <- ci + add
            if (!onlyo1) {
               bi <- bi + add
               di <- di + add
            }

         }

         if (to == "only0" || to == "if0all") {

            #if (onlyo1) {
            #   id0 <- c(ai == 0L | ci == 0L)
            #} else {
               id0 <- c(ai == 0L | ci == 0L | bi == 0L | di == 0L)
            #}
            id0[is.na(id0)] <- FALSE

         }

         if (to == "only0") {

            ### add to cells in studies with at least one 0 entry

            ai[id0] <- ai[id0] + add
            ci[id0] <- ci[id0] + add
            if (!onlyo1) {
               bi[id0] <- bi[id0] + add
               di[id0] <- di[id0] + add
            }

         }

         if (to == "if0all") {

            ### add to cells in all studies if there is at least one 0 entry

            if (any(id0)) {

               ai <- ai + add
               ci <- ci + add
               if (!onlyo1) {
                  bi <- bi + add
                  di <- di + add
               }

            }

         }

         ### recompute group and total sample sizes (after add/to adjustment)

         n1i <- ai + bi
         n2i <- ci + di
         ni  <- n1i + n2i # ni.u computed earlier is always the 'unadjusted' total sample size

         if (measure == "MPORM")
            ni <- round(ni / 2)

         ### compute proportions for the two groups (unadjusted and adjusted)

         p1i.u <- ai.u / n1i.u
         p2i.u <- ci.u / n2i.u
         p1i <- ai / n1i
         p2i <- ci / n2i

         ### compute sample size weighted averages of the proportions within groups (for vtype="AV")

         if (addvi) {
            mnwp1i <- .wmean(p1i, n1i, na.rm=TRUE)
            mnwp2i <- .wmean(p2i, n2i, na.rm=TRUE)
         } else {
            mnwp1i.u <- .wmean(p1i.u, n1i.u, na.rm=TRUE)
            mnwp2i.u <- .wmean(p2i.u, n2i.u, na.rm=TRUE)
         }

         ### log risk ratios

         if (measure == "RR") {

            if (addyi) {
               yi <- log(p1i) - log(p2i)
            } else {
               yi <- log(p1i.u) - log(p2i.u)
            }

            vtype <- .expand1(vtype, k)

            vi <- rep(NA_real_, k)

            if (!all(is.element(vtype, c("LS","AV"))))
               stop(mstyle$stop("For this outcome measure, 'vtype' must be either 'LS' or 'AV'."))

            for (i in seq_len(k)) {

               ### large sample approximation to the sampling variance
               if (vtype[i] == "LS") {
                  if (addvi) {
                     vi[i] <- 1/ai[i] - 1/n1i[i] + 1/ci[i] - 1/n2i[i]
                     #vi[i] <- (1-p1i[i])/(p1i[i]*n1i[i]) + (1-p2i[i])/(p2i[i]*n2i[i]) # same
                  } else {
                     vi[i] <- 1/ai.u[i] - 1/n1i.u[i] + 1/ci.u[i] - 1/n2i.u[i]
                  }
               }

               ### estimate assuming homogeneity (using the average proportions)
               if (vtype[i] == "AV") {
                  if (addvi) {
                     vi[i] <- (1-mnwp1i)/(mnwp1i*n1i[i]) + (1-mnwp2i)/(mnwp2i*n2i[i])
                  } else {
                     vi[i] <- (1-mnwp1i.u)/(mnwp1i.u*n1i.u[i]) + (1-mnwp2i.u)/(mnwp2i.u*n2i.u[i])
                  }
               }

            }

         }

         ### log odds ratio

         if (is.element(measure, c("OR","OR2D","OR2DN","OR2DL","MPORM"))) {

            if (addyi) {
               yi <- log(p1i/(1-p1i)) - log(p2i/(1-p2i))
            } else {
               yi <- log(p1i.u/(1-p1i.u)) - log(p2i.u/(1-p2i.u))
            }

            vtype <- .expand1(vtype, k)

            vi <- rep(NA_real_, k)

            if (!all(is.element(vtype, c("LS","AV"))))
               stop(mstyle$stop("For this outcome measure, 'vtype' must be either 'LS' or 'AV'."))

            for (i in seq_len(k)) {

               ### large sample approximation to the sampling variance
               if (vtype[i] == "LS") {
                  if (addvi) {
                     vi[i] <- 1/ai[i] + 1/bi[i] + 1/ci[i] + 1/di[i]
                     #vi[i] <- 1/(p1i[i]*(1-p1i[i])*n1i[i]) + 1/(p2i[i]*(1-p2i[i])*n2i[i]) # same
                  } else {
                     vi[i] <- 1/ai.u[i] + 1/bi.u[i] + 1/ci.u[i] + 1/di.u[i]
                  }
               }

               ### estimate assuming homogeneity (using the average proportions)
               if (vtype[i] == "AV") {
                  if (addvi) {
                     vi[i] <- 1/(mnwp1i*(1-mnwp1i)*n1i[i]) + 1/(mnwp2i*(1-mnwp2i)*n2i[i])
                  } else {
                     vi[i] <- 1/(mnwp1i.u*(1-mnwp1i.u)*n1i[i]) + 1/(mnwp2i.u*(1-mnwp2i.u)*n2i[i])
                  }
               }

            }

         }

         ### risk difference

         if (measure == "RD") {

            if (addyi) {
               yi <- p1i - p2i
            } else {
               yi <- p1i.u - p2i.u
            }

            vtype <- .expand1(vtype, k)

            vi <- rep(NA_real_, k)

            if (!all(is.element(vtype, c("UB","LS","AV"))))
               stop(mstyle$stop("For this outcome measure, 'vtype' must be either 'UB', 'LS', or 'AV'."))

            for (i in seq_len(k)) {

               ### unbiased estimate of the sampling variance
               if (vtype[i] == "UB") {
                  if (addvi) {
                     vi[i] <- p1i[i]*(1-p1i[i])/(n1i[i]-1) + p2i[i]*(1-p2i[i])/(n2i[i]-1)
                  } else {
                     vi[i] <- p1i.u[i]*(1-p1i.u[i])/(n1i.u[i]-1) + p2i.u[i]*(1-p2i.u[i])/(n2i.u[i]-1)
                  }
               }

               ### large sample approximation to the sampling variance
               if (vtype[i] == "LS") {
                  if (addvi) {
                     vi[i] <- p1i[i]*(1-p1i[i])/n1i[i] + p2i[i]*(1-p2i[i])/n2i[i]
                  } else {
                     vi[i] <- p1i.u[i]*(1-p1i.u[i])/n1i.u[i] + p2i.u[i]*(1-p2i.u[i])/n2i.u[i]
                  }
               }

               ### estimate assuming homogeneity (using the average proportions)
               if (vtype[i] == "AV") {
                  if (addvi) {
                     vi[i] <- mnwp1i*(1-mnwp1i)/n1i[i] + mnwp2i*(1-mnwp2i)/n2i[i]
                  } else {
                     vi[i] <- mnwp1i.u*(1-mnwp1i.u)/n1i.u[i] + mnwp2i.u*(1-mnwp2i.u)/n2i.u[i]
                  }
               }

            }

         }

         ### note: addyi and addvi only implemented for measures above

         ### log odds ratio (Peto's method)

         if (measure == "PETO") {
            xt <- ai + ci # frequency of outcome1 in both groups combined
            yt <- bi + di # frequency of outcome2 in both groups combined
            Ei <- xt * n1i / ni
            Vi <- xt * yt * (n1i/ni) * (n2i/ni) / (ni - 1) # 0 when xt = 0 or yt = 0 in a table
            yi <- (ai - Ei) / Vi                           # then yi and vi is Inf (set to NA at end)
            vi <- 1 / Vi
         }

         ### arcsine square root risk difference

         if (measure == "AS") {
            yi <- asin(sqrt(p1i)) - asin(sqrt(p2i))
            vi <- 1/(4*n1i) + 1/(4*n2i)
         }

         ### phi coefficient

         if (is.element(measure, c("PHI","ZPHI"))) {

            yi <- (ai*di - bi*ci)/sqrt((ai+bi)*(ci+di)*(ai+ci)*(bi+di))

            vtype <- .expand1(vtype, k)

            vi <- rep(NA_real_, k)

            q1i <- 1 - p1i
            q2i <- 1 - p2i
            pi1. <- (ai+bi) / ni
            pi2. <- (ci+di) / ni
            pi.1 <- (ai+ci) / ni
            pi.2 <- (bi+di) / ni

            if (!all(is.element(vtype, c("ST","LS","CS"))))
               stop(mstyle$stop("For this outcome measure, 'vtype' must be either 'ST', 'LS', or 'CS'."))

            for (i in seq_len(k)) {

               ### estimate of the sampling variance for stratified sampling
               if (vtype[i] == "ST") {
                  vi[i] <- ((n1i[i]+n2i[i])^2 * (4*n1i[i]^3*p1i[i]^2*p2i[i]*q1i[i]^2*q2i[i] +
                                                 4*n2i[i]^3*p1i[i]*p2i[i]^2*q1i[i]*q2i[i]^2 +
                                                 n1i[i]*n2i[i]^2*p2i[i]*q2i[i]*(p2i[i]*q1i[i] + p1i[i]*q2i[i])*(p2i[i]*q1i[i] + p1i[i]*(4*q1i[i] + q2i[i])) +
                                                 n1i[i]^2*n2i[i]*p1i[i]*q1i[i]*(p2i[i]*q1i[i] + p1i[i]*q2i[i])*(p1i[i]*q2i[i] + p2i[i]*(q1i[i] + 4*q2i[i])))) /
                           (4*(ai[i]+ci[i])^3*(bi[i]+di[i])^3)
               }

               ### estimate of the sampling variance for cross-sectional/multinomial sampling
               if (vtype[i] == "LS" || vtype[i] == "CS") {
                  vi[i] <- 1/ni[i] * (1 - yi[i]^2 + yi[i]*(1+1/2*yi[i]^2) * (pi1.[i]-pi2.[i])*(pi.1[i]-pi.2[i]) / sqrt(pi1.[i]*pi2.[i]*pi.1[i]*pi.2[i]) -
                                      3/4 * yi[i]^2 * ((pi1.[i]-pi2.[i])^2/(pi1.[i]*pi2.[i]) + (pi.1[i]-pi.2[i])^2/(pi.1[i]*pi.2[i]))) # Yule, 1912, p.603
               }

            }

         }

         ### Yule's Q

         if (measure == "YUQ") {
            yi <- (ai/bi) / (ci/di)
            yi <- (yi-1) / (yi+1)
            vi <- 1/4 * (1-yi^2)^2 * (1/ai + 1/bi + 1/ci + 1/di) # Yule, 1900, p.285; Yule, 1912, p.593
         }

         ### Yule's Y

         if (measure == "YUY") {
            yi <- (ai/bi) / (ci/di)
            yi <- (sqrt(yi)-1) / (sqrt(yi)+1)
            vi <- 1/16 * (1-yi^2)^2 * (1/ai + 1/bi + 1/ci + 1/di) # Yule, 1912, p.593
         }

         ### tetrachoric correlation

         if (is.element(measure, c("RTET","ZTET"))) {

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

         ### r-to-z transformation for PHI and RTET (note: NOT a variance-stabilizing transformation for these measures)

         if (is.element(measure, c("ZPHI","ZTET"))) {
            vi <- vi / (1 - ifelse(yi^2 > 1, 1, yi^2))^2
            yi <- transf.rtoz(yi)
         }

         ### probit transformation to SMD

         if (measure == "PBIT") {
            z1i <- qnorm(p1i)
            z2i <- qnorm(p2i)
            yi <- z1i - z2i
            vi <- 2*base::pi*p1i*(1-p1i)*exp(z1i^2)/n1i + 2*base::pi*p2i*(1-p2i)*exp(z2i^2)/n2i # Sanchez-Meca et al., 2003, equation 21; Rosenthal, 1994, handbook chapter
         }                                                                                      # seems to be right for stratified and cross-sectional/multinomial sampling
                                                                                                # see code/probit_transformation directory
         ### log(OR) transformation to SMD based on logistic distribution

         if (is.element(measure, c("OR2D","OR2DL"))) {
            yi <- sqrt(3) / base::pi * yi
            vi <- 3 / base::pi^2 * vi
         }

         ### log(OR) transformation to SMD based on normal distribution (Cox & Snell method)

         if (measure == "OR2DN") {
            yi <- yi / 1.65
            vi <- vi / 1.65^2
         }

         ### matched pairs / pre-post 2x2 table measures

         if (is.element(measure, c("MPRD","MPRR","MPOR"))) {
            pi12 <- bi / ni
            pi21 <- ci / ni
            pi1. <- (ai+bi) / ni
            pi.1 <- (ai+ci) / ni
         }

         if (measure == "MPRD") {
            yi <- pi1. - pi.1
            vi <- pi12*(1-pi12)/ni + 2*pi12*pi21/ni + pi21*(1-pi21)/ni
         }

         if (measure == "MPRR") {
            yi <- log(pi1.) - log(pi.1)
            vi <- (pi12 + pi21) / (ni * pi1. * pi.1)
         }

         if (measure == "MPOR") {
            yi <- log(pi1./(1-pi1.)) - log(pi.1/(1-pi.1))
            vi <- (pi12*(1-pi12) + pi21*(1-pi21) + 2*pi12*pi21) / (ni * pi1.*(1-pi1.) * pi.1*(1-pi.1))
         }

         if (measure == "MPORM") {

            ai.p <- pi * (ai.u+bi.u)
            bi.p <- ai.u - ai.p
            ci.p <- ci.u - ai.p
            di.p <- bi.u - ci.u + ai.p
            ri.p <- (ai.p*di.p - bi.p*ci.p) / sqrt((ai.p+bi.p)*(ci.p+di.p)*(ai.p+ci.p)*(bi.p+di.p))
            ri.p[ri.p < -1 | ri.p > 1] <- NA_real_
            ri <- replmiss(ri, ri.p)

            if (addvi) {
               si <- (ri * sqrt(ai * bi * ci * di) + (ai * bi)) / ni
               deltai <- ni^2 * (ni * si - ai * bi) / (ai * bi * ci * di)
               vi <- vi - 2*deltai / ni
            } else {
               si.u <- (ri * sqrt(ai.u * bi.u * ci.u * di.u) + (ai.u * bi.u)) / ni.u
               deltai.u <- ni.u^2 * (ni.u * si.u - ai.u * bi) / (ai.u * bi.u * ci.u * di.u)
               vi <- vi - 2*deltai.u / ni.u
            }
         }

         if (measure == "MPORC") {
            yi <- log(bi) - log(ci)
            vi <- 1/bi + 1/ci
         }

         if (measure == "MPPETO") {
            Ei <- (bi + ci) / 2
            Vi <- (bi + ci) / 4
            yi <- (bi - Ei) / Vi
            vi <- 1/Vi
         }

         ### Note: Could in principle also compute measures commonly used in diagnostic studies.
         ### But need to take the sampling method into consideration when computing vi (so need
         ### to give this some more thought).

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

      }

      ######################################################################

      if (is.element(measure, c("IRR","IRD","IRSD"))) {

         x1i <- .getx("x1i", mf=mf, data=data, checknumeric=TRUE)
         x2i <- .getx("x2i", mf=mf, data=data, checknumeric=TRUE)
         t1i <- .getx("t1i", mf=mf, data=data, checknumeric=TRUE)
         t2i <- .getx("t2i", mf=mf, data=data, checknumeric=TRUE)

         if (!.all.specified(x1i, x2i, t1i, t2i))
            stop(mstyle$stop("Cannot compute outcomes. Check that all of the required information is specified\n  via the appropriate arguments (i.e., x1i, x2i, t1i, t2i)."))

         if (!.equal.length(x1i, x2i, t1i, t2i))
            stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         k.all <- length(x1i)

         if (!is.null(subset)) {
            subset <- .chksubset(subset, k.all)
            x1i <- .getsubset(x1i, subset)
            x2i <- .getsubset(x2i, subset)
            t1i <- .getsubset(t1i, subset)
            t2i <- .getsubset(t2i, subset)
         }

         if (any(c(x1i, x2i) < 0, na.rm=TRUE))
            stop(mstyle$stop("One or more counts are negative."))

         if (any(c(t1i, t2i) <= 0, na.rm=TRUE))
            stop(mstyle$stop("One or more person-times are <= 0."))

         ni.u <- t1i + t2i # unadjusted total sample sizes

         k <- length(x1i)

         ### if drop00=TRUE, set counts to NA for studies that have no events in both arms

         if (drop00) {
            id00 <- c(x1i == 0L & x2i == 0L)
            id00[is.na(id00)] <- FALSE
            x1i[id00] <- NA_real_
            x2i[id00] <- NA_real_
         }

         ### save unadjusted counts

         x1i.u <- x1i
         x2i.u <- x2i

         if (to == "all") {

            ### always add to all cells in all studies

            x1i <- x1i + add
            x2i <- x2i + add

         }

         if (to == "only0" || to == "if0all") {

            id0 <- c(x1i == 0L | x2i == 0L)
            id0[is.na(id0)] <- FALSE

         }

         if (to == "only0") {

            ### add to cells in studies with at least one 0 entry

            x1i[id0] <- x1i[id0] + add
            x2i[id0] <- x2i[id0] + add

         }

         if (to == "if0all") {

            ### add to cells in all studies if there is at least one 0 entry

            if (any(id0)) {

               x1i <- x1i + add
               x2i <- x2i + add

            }

         }

         ### compute rates for the two groups (unadjusted and adjusted)
         ### t1i and t2i are the total person-times in the 1st and 2nd group

         ir1i.u <- x1i.u/t1i
         ir2i.u <- x2i.u/t2i
         ir1i <- x1i/t1i
         ir2i <- x2i/t2i

         ### log incidence rate ratio

         if (measure == "IRR") {
            if (addyi) {
               yi <- log(ir1i) - log(ir2i)
            } else {
               yi <- log(ir1i.u) - log(ir2i.u)
            }
            if (addvi) {
               vi <- 1/x1i + 1/x2i
               #vi <- 1/(x1i+1/2) + 1/(x2i+1/2)
            } else {
               vi <- 1/x1i.u + 1/x2i.u
            }
         }

         ### incidence rate difference

         if (measure == "IRD") {
            if (addyi) {
               yi <- ir1i - ir2i
            } else {
               yi <- ir1i.u - ir2i.u
            }
            if (addvi) {
               vi <- ir1i/t1i + ir2i/t2i     # same as x1i/t1i^2 + x2i/t2i^2
            } else {
               vi <- ir1i.u/t1i + ir2i.u/t2i # same as x1i.u/t1i^2 + x2i.u/t2i^2
            }
         }

         ### square root transformed incidence rate difference

         if (measure == "IRSD") {
            if (addyi) {
               yi <- sqrt(ir1i) - sqrt(ir2i)
            } else {
               yi <- sqrt(ir1i.u) - sqrt(ir2i.u)
            }
            vi <- 1/(4*t1i) + 1/(4*t2i)
         }

      }

      ######################################################################

      if (is.element(measure, c("MD","SMD","SMDH","SMD1","SMD1H","ROM","RPB","ZPB","RBIS","ZBIS","D2OR","D2ORN","D2ORL","CVR","VR"))) {

         m1i  <- .getx("m1i",  mf=mf, data=data, checknumeric=TRUE) # for VR, do not need to supply this
         m2i  <- .getx("m2i",  mf=mf, data=data, checknumeric=TRUE) # for VR, do not need to supply this
         sd1i <- .getx("sd1i", mf=mf, data=data, checknumeric=TRUE) # for SMD1, do not need to supply this
         sd2i <- .getx("sd2i", mf=mf, data=data, checknumeric=TRUE)
         n1i  <- .getx("n1i",  mf=mf, data=data, checknumeric=TRUE)
         n2i  <- .getx("n2i",  mf=mf, data=data, checknumeric=TRUE)
         di   <- .getx("di",   mf=mf, data=data, checknumeric=TRUE)
         ti   <- .getx("ti",   mf=mf, data=data, checknumeric=TRUE)
         pi   <- .getx("pi",   mf=mf, data=data, checknumeric=TRUE)

         ### for these measures, need m1i, m2i, sd1i, sd2i, n1i, and n2i (and can also specify di/ti/pi)

         if (is.element(measure, c("SMD","RPB","ZPB","RBIS","ZBIS","D2OR","D2ORN","D2ORL"))) {

            if (!.equal.length(m1i, m2i, sd1i, sd2i, n1i, n2i, di, ti, pi))
               stop(mstyle$stop("Supplied data vectors are not all of the same length."))

            ### convert pi to ti values

            ti <- replmiss(ti, .convp2t(pi, df=n1i+n2i-2))

            ### convert ti to di values

            di <- replmiss(di, ti * sqrt(1/n1i + 1/n2i))

            ### when di is available, set m1i, m2i, sd1i, and sd2i values accordingly

            m1i[!is.na(di)]  <- di[!is.na(di)]
            m2i[!is.na(di)]  <- 0
            sd1i[!is.na(di)] <- 1
            sd2i[!is.na(di)] <- 1

            if (!.all.specified(m1i, m2i, sd1i, sd2i, n1i, n2i))
               stop(mstyle$stop("Cannot compute outcomes. Check that all of the required information is specified\n  via the appropriate arguments (i.e., m1i, m2i, sd1i, sd2i, n1i, n2i (and di, ti, pi))."))

         }

         ### for these measures, need m1i, m2i, sd1i, sd2i, n1i, and n2i

         if (is.element(measure, c("MD","SMDH","SMD1H","ROM","CVR"))) {

            if (!.all.specified(m1i, m2i, sd1i, sd2i, n1i, n2i))
               stop(mstyle$stop("Cannot compute outcomes. Check that all of the required information is specified\n  via the appropriate arguments (i.e., m1i, m2i, sd1i, sd2i, n1i, n2i)."))

            if (!.equal.length(m1i, m2i, sd1i, sd2i, n1i, n2i))
               stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         }

         ### for this measure, need sd1i, sd2i, n1i, and n2i

         if (measure == "VR") {

            if (!.all.specified(sd1i, sd2i, n1i, n2i))
               stop(mstyle$stop("Cannot compute outcomes. Check that all of the required information is specified\n  via the appropriate arguments (i.e., sd1i, sd2i, n1i, n2i)."))

            if (!.equal.length(sd1i, sd2i, n1i, n2i))
               stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         }

         ### for this measure, need m1i, m2i, sd2i, n1i, and n2i

         if (measure == "SMD1") {

            if (!.all.specified(m1i, m2i, sd2i, n1i, n2i))
               stop(mstyle$stop("Cannot compute outcomes. Check that all of the required information is specified\n  via the appropriate arguments (i.e., m1i, m2i, sd2i, n1i, n2i)."))

            if (!.equal.length(m1i, m2i, sd2i, n1i, n2i))
               stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         }

         k.all <- length(n1i)

         if (!is.null(subset)) {
            subset <- .chksubset(subset, k.all)
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

         ni.u <- n1i + n2i # unadjusted total sample sizes

         k <- length(n1i)

         ni <- ni.u

         if (is.element(measure, c("SMD1","SMD1H"))) {
            mi   <- n2i - 1
            sdpi <- sd2i
            npi  <- n2i
         } else {
            mi   <- ni - 2
            sdpi <- sqrt(((n1i-1)*sd1i^2 + (n2i-1)*sd2i^2) / mi)
            npi  <- ni
         }

         di <- (m1i - m2i) / sdpi

         ### (raw) mean difference

         if (measure == "MD") {

            yi <- m1i - m2i

            vtype <- .expand1(vtype, k)

            vi <- rep(NA_real_, k)

            if (!all(is.element(vtype, c("LS","UB","HO"))))
               stop(mstyle$stop("For this outcome measure, 'vtype' must be either 'LS', 'UB', or 'HO'."))

            for (i in seq_len(k)) {

               ### unbiased estimate of the sampling variance (does not assume homoscedasticity)
               if (vtype[i] == "UB" || vtype[i] == "LS")
                  vi[i] <- sd1i[i]^2/n1i[i] + sd2i[i]^2/n2i[i]

               ### estimate assuming homoscedasticity of the variances within studies
               if (vtype[i] == "HO")
                  vi[i] <- sdpi[i]^2 * (1/n1i[i] + 1/n2i[i])

            }

         }

         ### standardized mean difference (with pooled SDs or just the SD of group 2)

         if (is.element(measure, c("SMD","SMD1"))) {

            ### apply bias-correction to di values

            cmi <- .cmicalc(mi, correct=correct)
            yi <- cmi * di

            vtype <- .expand1(vtype, k)

            vi <- rep(NA_real_, k)

            mnwyi <- .wmean(yi, ni, na.rm=TRUE) # sample size weighted average of yi's

            if (!all(is.element(vtype, c("LS","LS2","UB","AV"))))
               stop(mstyle$stop("For this outcome measure, 'vtype' must be either 'LS', 'LS2', 'UB', or 'AV'."))

            for (i in seq_len(k)) {

               ### large sample approximation to the sampling variance
               if (vtype[i] == "LS")
                  vi[i] <- 1/n1i[i] + 1/n2i[i] + yi[i]^2/(2*npi[i]) # Hedges, 1982c, equation 8; Hedges & Olkin, 1985, equation 15; see [a]

               ### alternative large sample approximation to the sampling variance
               if (vtype[i] == "LS2")
                  vi[i] <- cmi[i]^2 * (1/n1i[i] + 1/n2i[i] + di[i]^2/(2*npi[i])) # Borenstein, 2009, equation 12.17; analogous to LS2 for SMCC and SMCR; see [b]

               ### unbiased estimate of the sampling variance
               if (vtype[i] == "UB")
                  vi[i] <- 1/n1i[i] + 1/n2i[i] + (1 - (mi[i]-2)/(mi[i]*cmi[i]^2)) * yi[i]^2 # Hedges, 1983b, equation 9; see [c]

               ### estimate assuming homogeneity (using the sample size weighted average of the yi's)
               if (vtype[i] == "AV")
                  vi[i] <- 1/n1i[i] + 1/n2i[i] + mnwyi^2/(2*npi[i])

            }

         }

         ### standardized mean difference (with heteroscedastic SDs)

         if (measure == "SMDH") {

            cmi  <- .cmicalc(mi, correct=correct)
            sdpi <- sqrt((sd1i^2 + sd2i^2)/2)
            di   <- (m1i - m2i) / sdpi
            yi   <- cmi * di

            vtype <- .expand1(vtype, k)

            vi <- rep(NA_real_, k)

            if (!all(is.element(vtype, c("LS","LS2","LS3"))))
               stop(mstyle$stop("For this outcome measure, 'vtype' must be either 'LS' or 'LS2'."))

            for (i in seq_len(k)) {

               ### large sample approximation to the sampling variance
               if (vtype[i] == "LS") {
                  vi[i] <- yi[i]^2 * (sd1i[i]^4 / (n1i[i]-1) + sd2i[i]^4 / (n2i[i]-1)) / (8*sdpi[i]^4) +
                           (sd1i[i]^2 / (n1i[i]-1) + sd2i[i]^2 / (n2i[i]-1)) / sdpi[i]^2 # Bonett, 2008a, equation 8; Bonett, 2009, equation 5
                  # note: Bonett (2008a) plugs the uncorrected yi into the equation for vi; here, the corrected value is plugged in for consistency with [a]
                  #vi[i] <- cmi[i]^2 * vi[i]
               }

               ### alternative large sample approximation (replace n1i-1 and n2i-1 with n1i and n2i)
               if (vtype[i] == "LS2") {
                  #vi[i] <- sd1i[i]^2 / (n1i[i]     * sdpi[i]^2) + sd2i[i]^2 / (n2i[i]     * sdpi[i]^2) + yi[i]^2 / (8 * sdpi[i]^4) * (sd1i[i]^4 / (n1i[i]-1) + sd2i[i]^4 / (n2i[i]-1)) # based on standard application of the delta method
                  #vi[i] <- sd1i[i]^2 / ((n1i[i]-1) * sdpi[i]^2) + sd2i[i]^2 / ((n2i[i]-1) * sdpi[i]^2) + yi[i]^2 / (8 * sdpi[i]^4) * (sd1i[i]^4 / (n1i[i]-1) + sd2i[i]^4 / (n2i[i]-1)) # same as Bonett
                  vi[i] <-  sd1i[i]^2 / (n1i[i]     * sdpi[i]^2) + sd2i[i]^2 / (n2i[i]     * sdpi[i]^2) + yi[i]^2 / (8 * sdpi[i]^4) * (sd1i[i]^4 / n1i[i]     + sd2i[i]^4 / n2i[i])
               }

               ### alternative large sample approximation
               if (vtype[i] == "LS3") {
                  vi[i] <- sd1i[i]^2 / (n1i[i]     * sdpi[i]^2) + sd2i[i]^2 / (n2i[i]     * sdpi[i]^2) + yi[i]^2 / (8 * sdpi[i]^4) * (sd1i[i]^4 / (n1i[i]-1) + sd2i[i]^4 / (n2i[i]-1)) # based on standard application of the delta method
               }

            }

         }

         ### standardized mean difference standardized by SD of group 2 (with heteroscedastic SDs)

         if (measure == "SMD1H") {
            cmi <- .cmicalc(mi, correct=correct)
            yi <- cmi * di
            vi <- (sd1i^2/sd2i^2)/(n1i-1) + 1/(n2i-1) + yi^2/(2*(n2i-1)) # Bonett, 2008a, equation 12
            #vi <- cmi^2 * vi
         }

         ### ratio of means (response ratio)
         ### to use with pooled SDs, simply set sd1i = sd2i = sdpi or use vtype="HO"

         if (measure == "ROM") {

            yi <- log(m1i/m2i)

            vtype <- .expand1(vtype, k)

            vi <- rep(NA_real_, k)

            mn1wcvi <- .wmean(sd1i/m1i, n1i, na.rm=TRUE) # sample size weighted average of the coefficient of variation in group 1
            mn2wcvi <- .wmean(sd2i/m2i, n2i, na.rm=TRUE) # sample size weighted average of the coefficient of variation in group 2
            not.na  <- !(is.na(n1i) | is.na(n2i) | is.na(sd1i/m1i) | is.na(sd2i/m2i))
            mnwcvi  <- (sum(n1i[not.na]*(sd1i/m1i)[not.na]) + sum(n2i[not.na]*(sd2i/m2i)[not.na])) / sum((n1i+n2i)[not.na]) # sample size weighted average of the two CV values

            if (!all(is.element(vtype, c("LS","HO","AV","AVHO"))))
               stop(mstyle$stop("For this outcome measure, 'vtype' must be either 'LS', 'HO', 'AV', or 'AVHO'."))

            for (i in seq_len(k)) {

               ### large sample approximation to the sampling variance (does not assume homoscedasticity)
               if (vtype[i] == "LS")
                  vi[i] <- sd1i[i]^2/(n1i[i]*m1i[i]^2) + sd2i[i]^2/(n2i[i]*m2i[i]^2)

               ### estimate assuming homoscedasticity of the two variances within studies
               if (vtype[i] == "HO")
                  vi[i] <- sdpi[i]^2/(n1i[i]*m1i[i]^2) + sdpi[i]^2/(n2i[i]*m2i[i]^2)

               ### estimate using the weighted averages of the CV values
               if (vtype[i] == "AV")
                  vi[i] <- mn1wcvi^2/n1i[i] + mn2wcvi^2/n2i[i]

               ### estimate using the weighted average of two weighted averages of the CV values
               if (vtype[i] == "AVHO")
                  vi[i] <- mnwcvi^2 * (1/n1i[i] + 1/n2i[i])

            }

         }

         ### point-biserial correlation obtained from the standardized mean difference
         ### this is based on Tate's model where Y|X=0 and Y|X=1 are normally distributed (with the same variance)
         ### Das Gupta (1960) describes the case where Y itself is normal, but the variance expressions therein can
         ### really only be used in some special cases (not useful in practice)

         if (is.element(measure, c("RPB","RBIS","ZPB","ZBIS"))) {

            hi <- mi/n1i + mi/n2i
            yi <- di / sqrt(di^2 + hi) # need this also when measure="RBIS/ZBIS"

            if (is.element(measure, c("RPB","ZPB"))) { # this only applies when measure="RPB/ZPB"

               vtype <- .expand1(vtype, k)

               vi <- rep(NA_real_, k)

               if (!all(is.element(vtype, c("LS","ST","CS"))))
                  stop(mstyle$stop("For this outcome measure, 'vtype' must be either 'LS', 'ST', or 'CS'."))

               for (i in seq_len(k)) {

                  ### estimate of the sampling variance for fixed n1i and n2i (i.e., stratified sampling)
                  if (vtype[i] == "ST" || vtype[i] == "LS") {
                     vi[i] <- hi[i]^2 / (hi[i] + di[i]^2)^3 * (1/n1i[i] + 1/n2i[i] + di[i]^2/(2*ni[i]))
                     # this is consistent with escalc(measure="SMD", correct=FALSE) -> conv.delta(transf=transf.dtorpb)
                     #tmp <- escalc(measure="SMD", m1i=m1i[i], sd1i=sd1i[i], n1i=n1i[i], m2i=m2i[i], sd2i=sd2i[i], n2i=n2i[i], correct=FALSE)
                     #vi[i] <- conv.delta(yi, vi, data=tmp, transf=transf.dtorpb, replace=TRUE, n1i=n1i[i], n2i=n2i[i])$vi
                  }

                  ### estimate of the sampling variance for fixed ni but random n1i and n2i (i.e., cross-sectional/multinomial sampling)
                  if (vtype[i] == "CS")
                     vi[i] <- (1-yi[i]^2)^2 * (ni[i]*yi[i]^2 / (4*n1i[i]*n2i[i]) + (2-3*yi[i]^2)/(2*ni[i])) # Tate, 1954; Tate, 1955b

               }

            }

         }

         ### biserial correlation obtained from the standardized mean difference (continued from above)

         if (is.element(measure, c("RBIS","ZBIS"))) {
            p1i <- n1i / ni
            p2i <- n2i / ni
            zi  <- qnorm(p1i, lower.tail=FALSE)
            fzi <- dnorm(zi)
            yi  <- sqrt(p1i*p2i) / fzi * yi # yi on the right-hand side is the point-biserial correlation from above
            #vi <- (p1i*p2i) / fzi^2 * vi   # not correct (p1i, p2i, and fzi are random variables and vi from RBP is not correct for the bivariate normal case on which RBIS is based)
            yi.t <- ifelse(abs(yi) > 1, sign(yi), yi)
            vi  <- 1/(ni-1) * (p1i*p2i/fzi^2 - (3/2 + (1 - p1i*zi/fzi)*(1 + p2i*zi/fzi)) * yi.t^2 + yi.t^4) # Soper, 1914
            #vi <- 1/(ni-1) * (yi.t^4 + yi.t^2 * (p1i*p2i*zi^2/fzi^2 + (2*p1i-1)*zi/fzi - 5/2) + p1i*p2i/fzi^2) # Tate, 1955; equivalent to equation from Soper, 1914
            # equation appears to work even if dichotomization is done based on a sample quantile value (so that p1i, p2i, and fzi are fixed by design)
            # this is asymptotically consistent with escalc(measure="SMD", correct=FALSE) -> conv.delta(transf=transf.dtorbis)
            #tmp <- escalc(measure="SMD", m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i, correct=FALSE)
            #yi <- conv.delta(yi, vi, data=tmp, transf=transf.dtorbis, replace=TRUE, n1i=n1i, n2i=n2i)$yi
            #vi <- conv.delta(yi, vi, data=tmp, transf=transf.dtorbis, replace=TRUE, n1i=n1i, n2i=n2i)$vi
         }

         ### r-to-z transformation for RPB and RBIS (note: NOT a variance-stabilizing transformation for these measures)

         if (is.element(measure, c("ZPB","ZBIS"))) {
            vi <- vi / (1 - ifelse(yi^2 > 1, 1, yi^2))^2
            yi <- transf.rtoz(yi)
         }

         ### SMD to log(OR) transformation based on logistic distribution

         if (is.element(measure, c("D2OR","D2ORL"))) {
            yi <- base::pi / sqrt(3) * di
            vi <- base::pi^2 / 3 * (1/n1i + 1/n2i + di^2/(2*ni))
         }

         ### SMD to log(OR) transformation based on normal distribution (Cox & Snell method)

         if (measure == "D2ORN") {
            yi <- 1.65 * di
            vi <- 1.65^2 * (1/n1i + 1/n2i + di^2/(2*ni))
         }

         ### coefficient of variation ratio

         if (measure == "CVR") {
            if (correct) {
               yi <- log(sd1i/m1i) + 1/(2*(n1i-1)) - log(sd2i/m2i) - 1/(2*(n2i-1))
            } else {
               yi <- log(sd1i/m1i) - log(sd2i/m2i)
            }
            vi <- 1/(2*(n1i-1)) + sd1i^2/(n1i*m1i^2) + 1/(2*(n2i-1)) + sd2i^2/(n2i*m2i^2) # Nakagawa et al., 2015, equation 12, but without the '-2 rho ...' terms
         }

         ### variability ratio

         if (measure == "VR") {
            if (correct) {
               yi <- log(sd1i/sd2i) + 1/(2*(n1i-1)) - 1/(2*(n2i-1))
            } else {
               yi <- log(sd1i/sd2i)
            }
            vi <- 1/(2*(n1i-1)) + 1/(2*(n2i-1))
         }

      }

      ######################################################################

      if (is.element(measure, c("COR","UCOR","ZCOR"))) {

         ri <- .getx("ri", mf=mf, data=data, checknumeric=TRUE)
         ni <- .getx("ni", mf=mf, data=data, checknumeric=TRUE)
         ti <- .getx("ti", mf=mf, data=data, checknumeric=TRUE)
         pi <- .getx("pi", mf=mf, data=data, checknumeric=TRUE)

         if (!.equal.length(ri, ni, ti, pi))
            stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         ### convert pi to ti values

         ti <- replmiss(ti, .convp2t(pi, df=ni-2))

         ### convert ti to ri values

         ri <- replmiss(ri, ti / sqrt(ti^2 + ni-2))

         if (!.all.specified(ri, ni))
            stop(mstyle$stop("Cannot compute outcomes. Check that all of the required information is specified\n  via the appropriate arguments (i.e., ri, ni (and ti, pi))."))

         k.all <- length(ri)

         if (!is.null(subset)) {
            subset <- .chksubset(subset, k.all)
            ri <- .getsubset(ri, subset)
            ni <- .getsubset(ni, subset)
         }

         if (any(abs(ri) > 1, na.rm=TRUE))
            stop(mstyle$stop("One or more correlations are > 1 or < -1."))

         if (any(ni <= 0, na.rm=TRUE))
            stop(mstyle$stop("One or more sample sizes are <= 0."))

         if (measure != "UCOR" && any(vtype == "UB"))
            stop(mstyle$stop("Use of vtype='UB' only permitted when measure='UCOR'."))

         if (measure == "UCOR" && any(ni <= 4, na.rm=TRUE))
            warning(mstyle$warning("Cannot compute the bias-corrected correlation coefficient when ni <= 4."), call.=FALSE)

         if (measure == "ZCOR" && any(ni <= 3, na.rm=TRUE))
            warning(mstyle$warning("Cannot estimate the sampling variance when ni <= 3."), call.=FALSE)

         ni.u <- ni # unadjusted total sample sizes

         k <- length(ri)

         ### raw correlation coefficient

         if (measure == "COR")
            yi <- ri

         ### raw correlation coefficient with bias correction

         if (measure == "UCOR") {
            #yi <- ri + ri*(1-ri^2)/(2*(ni-4)) # approximation
            #yi[ni <= 4] <- NA_real_ # set corrected correlations for ni <= 4 to NA_real_
            yi <- ri * .Fcalc(1/2, 1/2, (ni-2)/2, 1-ri^2)
         }

         ### sampling variances for COR or UCOR

         if (is.element(measure, c("COR","UCOR"))) {

            vtype <- .expand1(vtype, k)

            vi <- rep(NA_real_, k)

            mnwyi <- .wmean(yi, ni, na.rm=TRUE) # sample size weighted average of yi's

            if (!all(is.element(vtype, c("LS","UB","AV"))))
               stop(mstyle$stop("For this outcome measure, 'vtype' must be either 'LS', 'UB', or 'AV'."))

            for (i in seq_len(k)) {

               ### large sample approximation to the sampling variance
               if (vtype[i] == "LS")
                  vi[i] <- (1-yi[i]^2)^2 / (ni[i]-1)

               ### unbiased estimate of the sampling variance of the bias-corrected correlation coefficient
               if (vtype[i] == "UB") {
                  #vi[i] <- yi[i]^2 - 1 + (ni[i]-3) / (ni[i]-2) * ((1-ri[i]^2) + 2*(1-ri[i]^2)^2/ni[i] + 8*(1-ri[i]^2)^3/(ni[i]*(ni[i]+2)) + 48*(1-ri[i]^2)^4/(ni[i]*(ni[i]+2)*(ni[i]+4)))
                  vi[i] <- yi[i]^2 - (1 - (ni[i]-3) / (ni[i]-2) * (1-ri[i]^2) * .Fcalc(1, 1, ni[i]/2, 1-ri[i]^2))
               }

               ### estimate assuming homogeneity (using sample size weighted average of the yi's)
               if (vtype[i] == "AV")
                  vi[i] <- (1-mnwyi^2)^2 / (ni[i]-1)

            }

         }

         ### r-to-z transformed correlation

         if (measure == "ZCOR") {
            yi <- transf.rtoz(ri)
            vi <- 1 / (ni-3)
         }

         ### set sampling variances for ni <= 3 to NA

         vi[ni <= 3] <- NA_real_

      }

      ######################################################################

      if (is.element(measure, c("PCOR","ZPCOR","SPCOR","ZSPCOR"))) {

         ri  <- .getx("ri",  mf=mf, data=data, checknumeric=TRUE)
         ti  <- .getx("ti",  mf=mf, data=data, checknumeric=TRUE)
         mi  <- .getx("mi",  mf=mf, data=data, checknumeric=TRUE)
         ni  <- .getx("ni",  mf=mf, data=data, checknumeric=TRUE)
         pi  <- .getx("pi",  mf=mf, data=data, checknumeric=TRUE)
         r2i <- .getx("r2i", mf=mf, data=data, checknumeric=TRUE)

         if (!.equal.length(ri, ti, mi, ni, pi, r2i))
            stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         ### convert pi to ti values

         ti <- replmiss(ti, .convp2t(pi, df=ni-mi-1))

         ### convert ti to ri values

         if (is.element(measure, c("PCOR","ZPCOR")))
            ri <- replmiss(ri, ti / sqrt(ti^2 + ni-mi-1))

         if (is.element(measure, c("SPCOR","ZSPCOR")))
            ri <- replmiss(ri, ti * sqrt(1-r2i) / sqrt(ni-mi-1))

         if (is.element(measure, c("PCOR","ZPCOR")) && !.all.specified(ri, mi, ni))
            stop(mstyle$stop("Cannot compute outcomes. Check that all of the required information is specified\n  via the appropriate arguments (i.e., ri, ti, mi, ni (and pi))."))

         if (is.element(measure, c("SPCOR","ZSPCOR")) && !.all.specified(ri, mi, ni, r2i))
            stop(mstyle$stop("Cannot compute outcomes. Check that all of the required information is specified\n  via the appropriate arguments (i.e., ri, ti, mi, ni, r2i (and pi))."))

         k.all <- length(ri)

         if (!is.null(subset)) {
            subset <- .chksubset(subset, k.all)
            ri  <- .getsubset(ri,  subset)
            mi  <- .getsubset(mi,  subset)
            ni  <- .getsubset(ni,  subset)
            r2i <- .getsubset(r2i, subset)
         }

         if (any(abs(ri) > 1, na.rm=TRUE))
            stop(mstyle$stop("One or more (semi-)partial correlations are > 1 or < -1."))

         if (is.element(measure, c("SPCOR","ZSPCOR")) && any(r2i > 1 | r2i < 0, na.rm=TRUE))
            stop(mstyle$stop("One or more R^2 values are > 1 or < 0."))

         if (any(ni <= 0, na.rm=TRUE))
            stop(mstyle$stop("One or more sample sizes are <= 0."))

         if (any(mi <= 0, na.rm=TRUE))
            stop(mstyle$stop("One or more mi values are <= 0."))

         if (any(ni-mi-1 <= 0, na.rm=TRUE))
            stop(mstyle$stop("One or more dfs are <= 0."))

         ni.u <- ni # unadjusted total sample sizes

         k <- length(ri)

         ### partial correlation coefficient

         if (measure == "PCOR") {

            yi <- ri

            vtype <- .expand1(vtype, k)

            vi <- rep(NA_real_, k)

            mnwyi <- .wmean(yi, ni, na.rm=TRUE) # sample size weighted average of yi's

            if (!all(is.element(vtype, c("LS","AV"))))
               stop(mstyle$stop("For this outcome measure, 'vtype' must be either 'LS' or 'AV'."))

            for (i in seq_len(k)) {

               ### large sample approximation to the sampling variance
               if (vtype[i] == "LS")
                  vi[i] <- (1 - yi[i]^2)^2 / (ni[i] - mi[i])

               ### estimate assuming homogeneity (using sample size weighted average of the yi's)
               if (vtype[i] == "AV")
                  vi[i] <- (1 - mnwyi^2)^2 / (ni[i] - mi[i])

            }

         }

         ### r-to-z transformed partial correlation

         if (measure == "ZPCOR") {
            yi <- transf.rtoz(ri)
            vi <- 1 / (ni-mi-2)
         }

         ### semi-partial correlation coefficient

         if (is.element(measure, c("SPCOR","ZSPCOR"))) {

            yi <- ri

            vtype <- .expand1(vtype, k)

            vi <- rep(NA_real_, k)

            mnwyi <- .wmean(yi, ni, na.rm=TRUE) # sample size weighted average of yi's

            if (!all(is.element(vtype, c("LS","AV"))))
               stop(mstyle$stop("For this outcome measure, 'vtype' must be either 'LS' or 'AV'."))

            for (i in seq_len(k)) {

               ### large sample approximation to the sampling variance
               if (vtype[i] == "LS")
                  vi[i] <- (r2i[i]^2 - 2*r2i[i] + (r2i[i] - yi[i]^2) + 1 - (r2i[i] - yi[i]^2)^2) / ni[i]

               ### estimate assuming homogeneity (using sample size weighted average of the yi's)
               if (vtype[i] == "AV")
                  vi[i] <- (r2i[i]^2 - 2*r2i[i] + (r2i[i] - mnwyi^2) + 1 - (r2i[i] - mnwyi^2)^2) / ni[i]

            }

         }

         ### r-to-z transformation for ZPCOR (note: NOT a variance-stabilizing transformation for this measure)

         if (measure == "ZSPCOR") {
            vi <- vi / (1 - ifelse(yi^2 > 1, 1, yi^2))^2
            yi <- transf.rtoz(yi)
         }

      }

      ######################################################################

      if (is.element(measure, c("R2","ZR2"))) {

         r2i <- .getx("r2i", mf=mf, data=data, checknumeric=TRUE)
         mi  <- .getx("mi",  mf=mf, data=data, checknumeric=TRUE)
         ni  <- .getx("ni",  mf=mf, data=data, checknumeric=TRUE)
         fi  <- .getx("fi",  mf=mf, data=data, checknumeric=TRUE)
         pi  <- .getx("pi",  mf=mf, data=data, checknumeric=TRUE)

         if (!.equal.length(r2i, mi, ni))
            stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         ### convert pi to fi values

         fi <- replmiss(fi, .convp2f(pi, df1=mi, df2=ni-mi-1))

         ### convert fi to r2i values

         r2i <- replmiss(r2i, mi*fi / (mi*fi + (ni-mi-1)))

         if (!.all.specified(r2i, mi, ni))
            stop(mstyle$stop("Cannot compute outcomes. Check that all of the required information is specified\n  via the appropriate arguments (i.e., r2i, mi, ni (and fi, pi))."))

         k.all <- length(r2i)

         if (!is.null(subset)) {
            subset <- .chksubset(subset, k.all)
            r2i <- .getsubset(r2i, subset)
            mi  <- .getsubset(mi,  subset)
            ni  <- .getsubset(ni,  subset)
         }

         if (any(r2i > 1 | r2i < 0, na.rm=TRUE))
            stop(mstyle$stop("One or more R^2 values are > 1 or < 0."))

         if (any(ni <= 0, na.rm=TRUE))
            stop(mstyle$stop("One or more sample sizes are <= 0."))

         if (any(mi <= 0, na.rm=TRUE))
            stop(mstyle$stop("One or more mi values are <= 0."))

         if (any(ni-mi- 1 <= 0, na.rm=TRUE))
            stop(mstyle$stop("One or more dfs are <= 0."))

         ni.u <- ni # unadjusted total sample sizes

         k <- length(r2i)

         ### coefficients of determination (R^2 values)

         if (measure == "R2") {

            yi <- r2i

            vtype <- .expand1(vtype, k)

            vi <- rep(NA_real_, k)

            mnwyi <- .wmean(yi, ni, na.rm=TRUE) # sample size weighted average of yi's

            if (!all(is.element(vtype, c("LS","LS2","AV"))))
               stop(mstyle$stop("For this outcome measure, 'vtype' must be either 'LS' or 'AV'."))

            for (i in seq_len(k)) {

               ### large sample approximation to the sampling variance (simplified equation)
               if (vtype[i] == "LS")
                  vi[i] <- 4 * yi[i] * (1 - yi[i])^2 / ni[i] # Kendall & Stuart, 1979, equation 27.88

               ### estimate assuming homogeneity (using sample size weighted average of the yi's)
               if (vtype[i] == "AV")
                  vi[i] <- 4 * mnwyi * (1 - mnwyi)^2 / ni[i]

               ### large sample approximation to the sampling variance (full equation)
               if (vtype[i] == "LS2")
                  vi[i] <- 4 * yi[i] * (1 - yi[i])^2 * (ni[i] - mi[i] - 1)^2 / ((ni[i]^2 - 1) * (ni[i] + 3)) # Kendall & Stuart, 1979, equation 27.87

               ### estimate assuming homogeneity (using sample size weighted average of the yi's)
               if (vtype[i] == "AV2")
                  vi[i] <- 4 * mnwyi * (1 - mnwyi)^2 * (ni[i] - mi[i] - 1)^2 / ((ni[i]^2 - 1) * (ni[i] + 3))

            }

         }

         ### r-to-z transformed coefficients of determination

         if (measure == "ZR2") {
            yi <- transf.rtoz(sqrt(r2i))
            vi <- 1 / ni # Olkin & Finn, 1995, p.162, but var(z*) is 4/n, not 16/n and here we use the 1/2 factor, so 1/n is correct
         }

      }

      ######################################################################

      if (is.element(measure, c("PR","PLN","PLO","PRZ","PAS","PFT"))) {

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

         k.all <- length(xi)

         if (!is.null(subset)) {
            subset <- .chksubset(subset, k.all)
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

         ni.u <- ni # unadjusted total sample sizes

         k <- length(xi)

         ### save unadjusted counts

         xi.u <- xi
         mi.u <- mi

         k <- length(xi)

         if (to == "all") {

            ### always add to all cells in all studies

            xi <- xi + add
            mi <- mi + add

         }

         if (to == "only0" || to == "if0all") {

            id0 <- c(xi == 0L | mi == 0L)
            id0[is.na(id0)] <- FALSE

         }

         if (to == "only0") {

            ### add to cells in studies with at least one 0 entry

            xi[id0] <- xi[id0] + add
            mi[id0] <- mi[id0] + add

         }

         if (to == "if0all") {

            ### add to cells in all studies if there is at least one 0 entry

            if (any(id0)) {

               xi <- xi + add
               mi <- mi + add

            }

         }

         ### recompute sample sizes (after add/to adjustment)

         ni <- xi + mi

         ### compute proportions (unadjusted and adjusted)

         pri.u <- xi.u/ni.u
         pri <- xi/ni

         ### raw proportion

         if (measure == "PR") {

            if (addyi) {
               yi <- pri
            } else {
               yi <- pri.u
            }

            vtype <- .expand1(vtype, k)

            vi <- rep(NA_real_, k)

            if (addvi) {
               mnwpri <- .wmean(pri, ni, na.rm=TRUE)       # sample size weighted average of proportions
            } else {
               mnwpri.u <- .wmean(pri.u, ni.u, na.rm=TRUE) # sample size weighted average of proportions
            }

            if (!all(is.element(vtype, c("LS","UB","AV"))))
               stop(mstyle$stop("For this outcome measure, 'vtype' must be either 'LS', 'UB', or 'AV'."))

            for (i in seq_len(k)) {

               ### large sample approximation to the sampling variance
               if (vtype[i] == "LS") {
                  if (addvi) {
                     vi[i] <- pri[i]*(1-pri[i])/ni[i]
                  } else {
                     vi[i] <- pri.u[i]*(1-pri.u[i])/ni.u[i]
                  }
               }

               ### unbiased estimate of the sampling variance
               if (vtype[i] == "UB") {
                  if (addvi) {
                     vi[i] <- pri[i]*(1-pri[i])/(ni[i]-1)
                  } else {
                     vi[i] <- pri.u[i]*(1-pri.u[i])/(ni.u[i]-1)
                  }
               }

               ### estimate assuming homogeneity (using the average proportion)
               if (vtype[i] == "AV") {
                  if (addvi) {
                     vi[i] <- mnwpri*(1-mnwpri)/ni[i]
                  } else {
                     vi[i] <- mnwpri.u*(1-mnwpri.u)/ni.u[i]
                  }
               }

            }

         }

         ### proportion with log transformation

         if (measure == "PLN") {

            if (addyi) {
               yi <- log(pri)
            } else {
               yi <- log(pri.u)
            }

            vtype <- .expand1(vtype, k)

            vi <- rep(NA_real_, k)

            if (addvi) {
               mnwpri <- .wmean(pri, ni, na.rm=TRUE) # sample size weighted average of proportions
               #mnwpri <- exp(.wmean(yi, ni, na.rm=TRUE)) # alternative strategy (exp of the sample size weighted average of the log proportions)
            } else {
               mnwpri.u <- .wmean(pri.u, ni.u, na.rm=TRUE) # sample size weighted average of proportions
            }

            if (!all(is.element(vtype, c("LS","AV"))))
               stop(mstyle$stop("For this outcome measure, 'vtype' must be either 'LS' or 'AV'."))

            for (i in seq_len(k)) {

               ### large sample approximation to the sampling variance
               if (vtype[i] == "LS") {
                  if (addvi) {
                     vi[i] <- 1/xi[i] - 1/ni[i]
                  } else {
                     vi[i] <- 1/xi.u[i] - 1/ni.u[i]
                  }
               }

               ### estimate assuming homogeneity (using the average proportion)
               if (vtype[i] == "AV") {
                  if (addvi) {
                     vi[i] <- 1/(mnwpri*ni[i]) - 1/ni[i]
                  } else {
                     vi[i] <- 1/(mnwpri.u*ni.u[i]) - 1/ni.u[i]
                  }
               }

            }

         }

         ### proportion with logit (log odds) transformation

         if (measure == "PLO") {

            if (addyi) {
               yi <- log(pri/(1-pri))
            } else {
               yi <- log(pri.u/(1-pri.u))
            }

            vtype <- .expand1(vtype, k)

            vi <- rep(NA_real_, k)

            if (addvi) {
               mnwpri <- .wmean(pri, ni, na.rm=TRUE) # sample size weighted average of proportions
               #mnwpri <- transf.ilogit(.wmean(yi, ni, na.rm=TRUE)) # alternative strategy (inverse logit of the sample size weighted average of the logit transformed proportions)
            } else {
               mnwpri.u <- .wmean(pri.u, ni.u, na.rm=TRUE) # sample size weighted average of proportions
            }

            if (!all(is.element(vtype, c("LS","AV"))))
               stop(mstyle$stop("For this outcome measure, 'vtype' must be either 'LS' or 'AV'."))

            for (i in seq_len(k)) {

               ### large sample approximation to the sampling variance
               if (vtype[i] == "LS") {
                  if (addvi) {
                     vi[i] <- 1/xi[i] + 1/mi[i]
                  } else {
                     vi[i] <- 1/xi.u[i] + 1/mi.u[i]
                  }
               }

               ### estimate assuming homogeneity (using the average proportion)
               if (vtype[i] == "AV") {
                  if (addvi) {
                     vi[i] <- 1/(mnwpri*ni[i]) + 1/((1-mnwpri)*ni[i])
                  } else {
                     vi[i] <- 1/(mnwpri.u*ni.u[i]) + 1/((1-mnwpri.u)*ni.u[i])
                  }
               }

            }

         }

         ### note: addyi and addvi only implemented for measures above

         ### proportion with probit transformation

         if (measure == "PRZ") {
            yi <- qnorm(pri)
            vi <- 2*base::pi*pri*(1-pri)*exp(yi^2)/ni
            # this is consistent with escalc(measure="PR") -> conv.delta(transf=qnorm)
            #tmp <- escalc(measure="PR", xi=xi, ni=ni)
            #vi <- conv.delta(yi, vi, data=tmp, transf=qnorm, replace=TRUE)$vi
         }

         ### proportion with arcsine square root (angular) transformation

         if (measure == "PAS") {
            yi <- asin(sqrt(pri))
            vi <- 1/(4*ni)
            # this is consistent with escalc(measure="PR") -> conv.delta(transf=transf.arcsin)
            #tmp <- escalc(measure="PR", xi=xi, ni=ni)
            #vi <- conv.delta(yi, vi, data=tmp, transf=transf.arcsin, replace=TRUE)$vi
         }

         ### proportion with Freeman-Tukey double arcsine transformation

         if (measure == "PFT") {
            yi <- 1/2*(asin(sqrt(xi/(ni+1))) + asin(sqrt((xi+1)/(ni+1))))
            vi <- 1/(4*ni+2)
            # this is asymptotically consistent with escalc(measure="PR") -> conv.delta(transf=transf.pft)
            #tmp <- escalc(measure="PR", xi=xi, ni=ni)
            #vi <- conv.delta(yi, vi, data=tmp, transf=transf.pft, ni=ni, replace=TRUE)$vi
         }

      }

      ######################################################################

      if (is.element(measure, c("IR","IRLN","IRS","IRFT"))) {

         xi <- .getx("xi", mf=mf, data=data, checknumeric=TRUE)
         ti <- .getx("ti", mf=mf, data=data, checknumeric=TRUE)

         if (!.all.specified(xi, ti))
            stop(mstyle$stop("Cannot compute outcomes. Check that all of the required information is specified\n  via the appropriate arguments (i.e., xi, ti)."))

         if (!.equal.length(xi, ti))
            stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         k.all <- length(xi)

         if (!is.null(subset)) {
            subset <- .chksubset(subset, k.all)
            xi <- .getsubset(xi, subset)
            ti <- .getsubset(ti, subset)
         }

         if (any(xi < 0, na.rm=TRUE))
            stop(mstyle$stop("One or more counts are negative."))

         if (any(ti <= 0, na.rm=TRUE))
            stop(mstyle$stop("One or more person-times are <= 0."))

         ni.u <- ti # unadjusted total sample sizes

         k <- length(xi)

         ### save unadjusted counts

         xi.u <- xi

         if (to == "all") {

            ### always add to all cells in all studies

            xi <- xi + add

         }

         if (to == "only0" || to == "if0all") {

            id0 <- c(xi == 0L)
            id0[is.na(id0)] <- FALSE

         }

         if (to == "only0") {

            ### add to cells in studies with at least one 0 entry

            xi[id0] <- xi[id0] + add

         }

         if (to == "if0all") {

            ### add to cells in all studies if there is at least one 0 entry

            if (any(id0)) {

               xi <- xi + add

            }

         }

         ### compute rates (unadjusted and adjusted)

         iri.u <- xi.u / ti
         iri <- xi / ti

         ### raw incidence rate

         if (measure == "IR") {
            if (addyi) {
               yi <- iri
            } else {
               yi <- iri.u
            }
            if (addvi) {
               vi <- iri / ti   # same as xi/ti^2
            } else {
               vi <- iri.u / ti # same as xi.u/ti^2
            }
         }

         ### log transformed incidence rate

         if (measure == "IRLN") {
            if (addyi) {
               yi <- log(iri)
            } else {
               yi <- log(iri.u)
            }
            if (addvi) {
               vi <- 1 / xi
            } else {
               vi <- 1 / xi.u
            }
         }

         ### square root transformed incidence rate

         if (measure == "IRS") {
            if (addyi) {
               yi <- sqrt(iri)
            } else {
               yi <- sqrt(iri.u)
            }
            vi <- 1 / (4*ti)
            # this is consistent with escalc(measure="IR") -> conv.delta(transf=sqrt)
            #tmp <- escalc(measure="IR", xi=xi, ti=ti)
            #vi <- conv.delta(yi, vi, data=tmp, transf=sqrt, replace=TRUE)$vi
         }

         ### note: addyi and addvi only implemented for measures above

         ### incidence rate with Freeman-Tukey transformation

         if (measure == "IRFT") {
            yi <- 1/2 * (sqrt(iri) + sqrt(iri+1/ti))
            vi <- 1 / (4*ti)
            # this is asymptotically consistent with escalc(measure="IR") -> conv.delta(transf=transf.irft)
            #tmp <- escalc(measure="IR", xi=xi, ti=ti)
            #vi <- conv.delta(yi, vi, data=tmp, transf=transf.irft, ti=ti, replace=TRUE)$vi
         }

      }

      ######################################################################

      if (is.element(measure, c("MN","SMN","MNLN","CVLN","SDLN"))) {

         mi  <- .getx("mi",  mf=mf, data=data, checknumeric=TRUE) # for SDLN, do not need to supply this
         sdi <- .getx("sdi", mf=mf, data=data, checknumeric=TRUE)
         ni  <- .getx("ni",  mf=mf, data=data, checknumeric=TRUE)

         ### for these measures, need mi, sdi, and ni

         if (is.element(measure, c("MN","SMN","MNLN","CVLN"))) {

            if (!.all.specified(mi, sdi, ni))
               stop(mstyle$stop("Cannot compute outcomes. Check that all of the required information is specified\n  via the appropriate arguments (i.e., mi, sdi, ni)."))

            if (!.equal.length(mi, sdi, ni))
               stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         }

         ### for this measure, need sdi and ni

         if (measure == "SDLN") {

            if (!.all.specified(sdi, ni))
               stop(mstyle$stop("Cannot compute outcomes. Check that all of the required information is specified\n  via the appropriate arguments (i.e., sdi, ni)."))

            if (!.equal.length(sdi, ni))
               stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         }

         k.all <- length(ni)

         if (!is.null(subset)) {
            subset <- .chksubset(subset, k.all)
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

         ni.u <- ni # unadjusted total sample sizes

         k <- length(ni)

         ### (raw) mean

         if (measure == "MN") {

            yi <- mi
            sdpi <- sqrt(.wmean(sdi^2, ni-1, na.rm=TRUE))

            vtype <- .expand1(vtype, k)

            vi <- rep(NA_real_, k)

            if (!all(is.element(vtype, c("LS","HO"))))
               stop(mstyle$stop("For this outcome measure, 'vtype' must be either 'LS' or 'HO'."))

            for (i in seq_len(k)) {

               ### unbiased estimate of the sampling variance
               if (vtype[i] == "LS")
                  vi[i] <- sdi[i]^2 / ni[i]

               ### estimate assuming homoscedasticity of the variances across studies
               if (vtype[i] == "HO")
                  vi[i] <- sdpi^2 / ni[i]

            }

         }

         ### single-group standardized mean

         if (measure == "SMN") {
            cmi <- .cmicalc(ni-1, correct=correct)
            yi <- cmi * mi / sdi
            vi <- 1 / ni + yi^2 / (2*ni)
         }

         ### log(mean)

         if (measure == "MNLN") {
            yi <- log(mi)
            vi <- sdi^2 / (ni*mi^2)
         }

         ### log(CV) with bias correction

         if (measure == "CVLN") {
            if (correct) {
               yi <- log(sdi/mi) + 1 / (2*(ni-1))
            } else {
               yi <- log(sdi/mi)
            }
            vi <- 1 / (2*(ni-1)) + sdi^2 / (ni*mi^2) # Nakagawa et al., 2015, but without the '-2 rho ...' term
         }

         ### log(SD) with bias correction

         if (measure == "SDLN") {
            if (correct) {
               yi <- log(sdi) + 1 / (2*(ni-1))
            } else {
               yi <- log(sdi)
            }
            vi <- 1 / (2*(ni-1))
         }

      }

      ######################################################################

      if (is.element(measure, c("MC","SMCC","SMCR","SMCRH","SMCRP","SMCRPH","CLESCN","AUCCN","ROMC","CVRC","VRC"))) {

         m1i  <- .getx("m1i",  mf=mf, data=data, checknumeric=TRUE) # for VRC, do not need to supply this
         m2i  <- .getx("m2i",  mf=mf, data=data, checknumeric=TRUE) # for VRC, do not need to supply this
         sd1i <- .getx("sd1i", mf=mf, data=data, checknumeric=TRUE)
         sd2i <- .getx("sd2i", mf=mf, data=data, checknumeric=TRUE) # for SMCR, do not need to supply this
         ri   <- .getx("ri",   mf=mf, data=data, checknumeric=TRUE)
         ni   <- .getx("ni",   mf=mf, data=data, checknumeric=TRUE)
         di   <- .getx("di",   mf=mf, data=data, checknumeric=TRUE)
         ti   <- .getx("ti",   mf=mf, data=data, checknumeric=TRUE)
         pi   <- .getx("pi",   mf=mf, data=data, checknumeric=TRUE)

         ri <- .expand1(ri, list(m1i, m2i, sd1i, sd2i, ni, di, ti, pi))

         if (is.element(measure, c("MC","SMCRH","SMCRP","SMCRPH","CLESCN","AUCCN","ROMC","CVRC"))) {

            ### for these measures, need m1i, m2i, sd1i, sd2i, ni, and ri

            if (!.all.specified(m1i, m2i, sd1i, sd2i, ri, ni))
               stop(mstyle$stop("Cannot compute outcomes. Check that all of the required information is specified\n  via the appropriate arguments (i.e., m1i, m2i, sd1i, sd2i, ni, ri)."))

            if (!.equal.length(m1i, m2i, sd1i, sd2i, ri, ni))
               stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         }

         if (measure == "SMCC") {

            ### for this measures, need m1i, m2i, sd1i, sd2i, ni, and ri (and can also specify di/ti/pi)

            if (!.equal.length(m1i, m2i, sd1i, sd2i, ri, ni, di, ti, pi))
               stop(mstyle$stop("Supplied data vectors are not all of the same length."))

            ### convert pi to ti values

            ti <- replmiss(ti, .convp2t(pi, df=ni-1))

            ### convert ti to di values

            di <- replmiss(di, ti * sqrt(1/ni))

            ### when di is available, set m1i, m2i, sd1i, sd2i, and ri values accordingly

            m1i[!is.na(di)]  <- di[!is.na(di)]
            m2i[!is.na(di)]  <- 0
            sd1i[!is.na(di)] <- 1
            sd2i[!is.na(di)] <- 1
            ri[!is.na(di)]   <- 0.5

            if (!.all.specified(m1i, m2i, sd1i, sd2i, ri, ni))
               stop(mstyle$stop("Cannot compute outcomes. Check that all of the required information is specified\n  via the appropriate arguments (i.e., m1i, m2i, sd1i, sd2i, ni, ri (and di, ti, pi))."))

         }

         if (measure == "SMCR") {

            ### for this measure, need m1i, m2i, sd1i, ni, and ri (do not need sd2i)

            if (!.all.specified(m1i, m2i, sd1i, ri, ni))
               stop(mstyle$stop("Cannot compute outcomes. Check that all of the required information is specified\n  via the appropriate arguments (i.e., m1i, m2i, sd1i, ni, ri)."))

            if (!.equal.length(m1i, m2i, sd1i, ri, ni))
               stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         }

         if (measure == "VRC") {

            ### for this measure, need sd1i, sd2i, ni, and ri

            if (!.all.specified(sd1i, sd2i, ri, ni))
               stop(mstyle$stop("Cannot compute outcomes. Check that all of the required information is specified\n  via the appropriate arguments (i.e., sd1i, sd2i, ni, ri)."))

            if (!.equal.length(sd1i, sd2i, ri, ni))
               stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         }

         k.all <- length(ni)

         if (!is.null(subset)) {
            subset <- .chksubset(subset, k.all)
            m1i  <- .getsubset(m1i,  subset)
            m2i  <- .getsubset(m2i,  subset)
            sd1i <- .getsubset(sd1i, subset)
            sd2i <- .getsubset(sd2i, subset)
            ni   <- .getsubset(ni,   subset)
            ri   <- .getsubset(ri,   subset)
         }

         if (is.element(measure, c("MC","SMCC","SMCRH","SMCRP","SMCRPH","CLESCN","AUCCN","ROMC","CVRC","VRC"))) {
            if (any(c(sd1i, sd2i) < 0, na.rm=TRUE))
               stop(mstyle$stop("One or more standard deviations are negative."))
         }

         if (measure == "SMCR") {
            if (any(sd1i < 0, na.rm=TRUE))
               stop(mstyle$stop("One or more standard deviations are negative."))
         }

         if (any(abs(ri) > 1, na.rm=TRUE))
            stop(mstyle$stop("One or more correlations are > 1 or < -1."))

         if (any(ni <= 0, na.rm=TRUE))
            stop(mstyle$stop("One or more sample sizes are <= 0."))

         ni.u <- ni # unadjusted total sample sizes

         k <- length(ni)

         ni <- ni.u
         mi <- ni - 1

         sddiffi <- sqrt(sd1i^2 + sd2i^2 - 2*ri*sd1i*sd2i) # SD of the change scores
         sdpi <- sqrt((sd1i^2+sd2i^2)/2) # pooled SD

         ### (raw) mean change

         if (measure == "MC") {
            yi <- m1i - m2i
            vi <- sddiffi^2 / ni
         }

         ### standardized mean change with change score standardization (using sddi)
         ### note: does not assume homoscedasticity, since we use sddi here

         if (measure == "SMCC") {

            cmi <- .cmicalc(mi, correct=correct)
            di <- (m1i - m2i) / sddiffi
            yi <- cmi * di

            vtype <- .expand1(vtype, k)

            vi <- rep(NA_real_, k)

            if (!all(is.element(vtype, c("LS","LS2","UB"))))
               stop(mstyle$stop("For this outcome measure, 'vtype' must be either 'LS', 'LS2', or 'UB'."))

            for (i in seq_len(k)) {

               ### large sample approximation to the sampling variance
               if (vtype[i] == "LS")
                  vi[i] <- 1/ni[i] + yi[i]^2 / (2*ni[i]) # Gibbons et al., 1993, equation 21, but using ni instead of ni-1; see [a]

               ### alternative large sample approximation to the sampling variance
               if (vtype[i] == "LS2")
                  vi[i] <- cmi[i]^2 * (1/ni[i] + di[i]^2 / (2*ni[i])) # analogous to LS2 for SMD and SMCR; see [b]

               ### unbiased estimate of the sampling variance
               if (vtype[i] == "UB")
                  vi[i] <- 1/ni[i] + (1 - (mi[i]-2)/(mi[i]*cmi[i]^2)) * yi[i]^2 # Viechtbauer, 2007d, equation 26; see [c]

            }

         }

         ### standardized mean change with raw score standardization (using sd1i)

         if (measure == "SMCR") {

            cmi <- .cmicalc(mi, correct=correct)
            di <- (m1i - m2i) / sd1i
            yi <- cmi * di

            vtype <- .expand1(vtype, k)

            vi <- rep(NA_real_, k)

            if (!all(is.element(vtype, c("LS","LS2"))))
               stop(mstyle$stop("For this outcome measure, 'vtype' must be either 'LS' or 'LS2'."))

            for (i in seq_len(k)) {

               ### large sample approximation to the sampling variance
               if (vtype[i] == "LS")
                  vi[i] <- 2*(1-ri[i])/ni[i] + yi[i]^2 / (2*ni[i]) # Becker, 1988a, equation 13

               ### alternative large sample approximation to the sampling variance
               if (vtype[i] == "LS2")
                  vi[i] <- cmi[i]^2 * (2*(1-ri[i])/ni[i] + di[i]^2 / (2*ni[i])) # corrected (!) equation from Borenstein et al., 2009; analogous to LS2 for SMD and SMCC; see [b]
                  #vi[i] <- cmi[i]^2 * 2 * (1-ri[i]) * (1/ni[i] + di[i]^2 / (2*ni[i])) # Borenstein, 2009, equation 4.28 (with J^2 multiplier) but this is incorrect

               ### unbiased estimate of the sampling variance
               if (vtype[i] == "UB") {
                  rui[i] <- ri[i] * .Fcalc(1/2, 1/2, (ni[i]-2)/2, 1-ri[i]^2) # NA when ni <= 4
                  vi[i] <- 2*(1-rui[i])/ni[i] + (1 - (mi[i]-2)/(mi[i]*cmi[i]^2)) * yi[i]^2 # Viechtbauer, 2007d, equation 37; see [c]
               }

            }

         }

         ### standardized mean change with raw score standardization (using sd1i) allowing for heteroscedasticity

         if (measure == "SMCRH") {

            cmi <- .cmicalc(mi, correct=correct)
            di <- (m1i - m2i) / sd1i
            yi <- cmi * di

            vtype <- .expand1(vtype, k)

            vi <- rep(NA_real_, k)

            if (!all(is.element(vtype, c("LS","LS2"))))
               stop(mstyle$stop("For this outcome measure, 'vtype' must be either 'LS' or 'LS2'."))

            for (i in seq_len(k)) {

               ### large sample approximation to the sampling variance
               if (vtype[i] == "LS") {
                  vi[i] <- sddiffi[i]^2/(sd1i[i]^2*(ni[i]-1)) + yi[i]^2 / (2*(ni[i]-1)) # Bonett, 2008a, equation 13
                  # note: Bonett (2008a) plugs the uncorrected yi into the equation for vi; here, the corrected value is plugged in for consistency with [a]
                  #vi <- cmi^2 * vi
               }

               ### alternative large sample approximation (replace ni-1 with ni)
               if (vtype[i] == "LS2")
                  vi[i] <- sddiffi[i]^2/(sd1i[i]^2*ni[i]) + yi[i]^2 / (2*ni[i])

            }

         }

         ### standardized mean change with raw score standardization (using (sd1i+sd2i)/2))

         if (measure == "SMCRP") {

            mi <- 2*(ni-1) / (1 + ri^2)
            cmi <- .cmicalc(mi, correct=correct)
            di <- (m1i - m2i) / sdpi
            yi <- cmi * di

            vtype <- .expand1(vtype, k)

            vi <- rep(NA_real_, k)

            if (!all(is.element(vtype, c("LS"))))
               stop(mstyle$stop("For this outcome measure, 'vtype' must be 'LS'."))

            for (i in seq_len(k)) {

               ### large sample approximation to the sampling variance
               if (vtype[i] == "LS")
                  vi[i] <- 2 * (1-ri[i]) / ni[i] + yi[i]^2 * (1 + ri[i]^2) / (4*ni[i]) # follows from Cousineau, 2020, equation 2

            }

         }

         ### standardized mean change with raw score standardization (using (sd1i+sd2i)/2)) allowing for heteroscedasticity

         if (measure == "SMCRPH") {

            mi <- 2*(ni-1) / (1 + ri^2)
            cmi <- .cmicalc(mi, correct=correct)
            di <- (m1i - m2i) / sdpi
            yi <- cmi * di

            vtype <- .expand1(vtype, k)

            vi <- rep(NA_real_, k)

            if (!all(is.element(vtype, c("LS","LS2"))))
               stop(mstyle$stop("For this outcome measure, 'vtype' must be 'LS' or 'LS2'."))

            for (i in seq_len(k)) {

               ### large sample approximation to the sampling variance
               if (vtype[i] == "LS")
                  vi[i] <- sddiffi[i]^2 / (sdpi[i]^2 * (ni[i]-1)) + yi[i]^2 * (sd1i[i]^4 + sd2i[i]^4 + 2*ri[i]^2*sd1i[i]^2*sd2i[i]^2) / (8 * sdpi[i]^4 * (ni[i]-1)) # Bonett, 2008a, equation 10

               ### alternative large sample approximation to the sampling variance (replace ni-1 with ni)
               if (vtype[i] == "LS2")
                  vi[i] <- sddiffi[i]^2 / (sdpi[i]^2 * ni[i]) + yi[i]^2 * (sd1i[i]^4 + sd2i[i]^4 + 2*ri[i]^2*sd1i[i]^2*sd2i[i]^2) / (8 * sdpi[i]^4 * ni[i])

            }

         }

         ### common language effect size / area under the curve allowing for heteroscedasticity

         if (is.element(measure, c("CLESCN","AUCCN"))) {

            di <- (m1i - m2i) / sdpi
            yi <- pnorm(di/sqrt(2))
            vi <- exp(-di^2 / 2) / (4*base::pi) * (sddiffi^2 / (sdpi^2 * (ni-1)) + di^2 * (sd1i^4 + sd2i^4 + 2*ri^2*sd1i^2*sd2i^2) / (8 * sdpi^4 * (ni-1)))

         }

         ### ratio of means for pre-post or matched designs (eq. 6 in Lajeunesse, 2011)
         ### to use with pooled SDs, simply set sd1i = sd2i = sdpi

         if (measure == "ROMC") {
            yi <- log(m1i/m2i)
            vi <- sd1i^2 / (ni*m1i^2) + sd2i^2 / (ni*m2i^2) - 2*ri*sd1i*sd2i/(m1i*m2i*ni)
         }

         ### coefficient of variation ratio for pre-post or matched designs

         if (measure == "CVRC") {
            yi <- log(sd1i/m1i) - log(sd2i/m2i)
            vi <- (1-ri^2) / (ni-1) + (m1i^2*sd2i^2 + m2i^2*sd1i^2 - 2*m1i*m2i*ri*sd1i*sd2i) / (m1i^2*m2i^2*ni)
         }

         ### variability ratio for pre-post or matched designs

         if (measure == "VRC") {
            yi <- log(sd1i/sd2i)
            vi <- (1-ri^2) / (ni-1)
         }

      }

      ######################################################################

      if (is.element(measure, c("ARAW","AHW","ABT"))) {

         ai <- .getx("ai", mf=mf, data=data, checknumeric=TRUE)
         mi <- .getx("mi", mf=mf, data=data, checknumeric=TRUE)
         ni <- .getx("ni", mf=mf, data=data, checknumeric=TRUE)

         if (!.all.specified(ai, mi, ni))
            stop(mstyle$stop("Cannot compute outcomes. Check that all of the required information is specified\n  via the appropriate arguments (i.e., ai, mi, ni)."))

         if (!.equal.length(ai, mi, ni))
            stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         k.all <- length(ai)

         if (!is.null(subset)) {
            subset <- .chksubset(subset, k.all)
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

         ni.u <- ni # unadjusted total sample sizes

         k <- length(ai)

         ### raw alpha values

         if (measure == "ARAW") {
            yi <- ai
            vi <- 2*mi*(1-ai)^2 / ((mi-1)*(ni-2))
         }

         ### alphas transformed with Hakstian & Whalen (1976) transformation

         if (measure == "AHW") {
            #yi <- (1-ai)^(1/3)    # technically this is the Hakstian & Whalen (1976) transformation
            yi <- 1 - (1-ai)^(1/3) # but with this, yi remains a monotonically increasing function of ai
            vi <- 18*mi*(ni-1)*(1-ai)^(2/3) / ((mi-1)*(9*ni-11)^2)
            #vi <- 2*mi*(1-ai)^(2/3) / (9*(mi-1)*(ni-2)) # this follows from the delta method
            # this is asymptotically consistent with escalc(measure="ARAW") -> conv.delta(transf=transf.ahw)
            #tmp <- escalc(measure="ARAW", ai=ai, mi=mi, ni=ni)
            #vi <- conv.delta(yi, vi, data=tmp, transf=transf.ahw, replace=TRUE)$vi
         }

         ### alphas transformed with Bonett (2002) transformation (without bias correction)

         if (measure == "ABT") {
            #yi <- log(1-ai) - log(ni/(ni-1))
            #yi <- log(1-ai) # technically this is the Bonett (2002) transformation
            yi <- -log(1-ai) # but with this, yi remains a monotonically increasing function of ai
            vi <- 2*mi / ((mi-1)*(ni-2))
            # this is consistent with escalc(measure="ARAW") -> conv.delta(transf=transf.abt)
            #tmp <- escalc(measure="ARAW", ai=ai, mi=mi, ni=ni)
            #vi <- conv.delta(yi, vi, data=tmp, transf=transf.abt, replace=TRUE)$vi
         }

      }

      ######################################################################

      if (measure == "REH") {

         ai <- .getx("ai", mf=mf, data=data, checknumeric=TRUE)
         bi <- .getx("bi", mf=mf, data=data, checknumeric=TRUE)
         ci <- .getx("ci", mf=mf, data=data, checknumeric=TRUE)

         if (!.all.specified(ai, bi, ci))
            stop(mstyle$stop("Cannot compute outcomes. Check that all of the required information is specified\n  via the appropriate arguments (i.e., ai, bi, ci)."))

         if (!.equal.length(ai, bi, ci))
            stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         k.all <- length(ai)

         if (!is.null(subset)) {
            subset <- .chksubset(subset, k.all)
            ai <- .getsubset(ai, subset)
            bi <- .getsubset(bi, subset)
            ci <- .getsubset(ci, subset)
         }

         if (any(ai < 0, na.rm=TRUE) || any(bi < 0, na.rm=TRUE) || any(ci < 0, na.rm=TRUE))
            stop(mstyle$stop("One or more group sizes are negative."))

         ni <- ai + bi + ci

         ni.u <- ni # unadjusted total sample sizes

         k <- length(ai)

         p0i <- ai / ni
         p1i <- bi / ni
         p2i <- ci / ni

         yi <- log(p1i) - log(2 * sqrt(p0i * p2i))
         vi <- ((1-p1i) / (4 * p0i * p2i) + 1 / p1i) / ni

      }

      ######################################################################

      if (is.element(measure, c("CLES","AUC"))) {

         ai  <- .getx("ai",  mf=mf, data=data, checknumeric=TRUE)
         n1i <- .getx("n1i", mf=mf, data=data, checknumeric=TRUE)
         n2i <- .getx("n2i", mf=mf, data=data, checknumeric=TRUE)
         mi  <- .getx("mi",  mf=mf, data=data, checknumeric=TRUE)

         if (!.all.specified(ai, n1i, n2i))
            stop(mstyle$stop("Cannot compute outcomes. Check that all of the required information is specified\n  via the appropriate arguments (i.e., ai, n1i, n2i)."))

         if (!.equal.length(ai, n1i, n2i, mi))
            stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         if (is.null(mi))
            mi <- rep(0, length(ai))

         mi[is.na(mi)] <- 0

         k.all <- length(ai)

         if (!is.null(subset)) {
            subset <- .chksubset(subset, k.all)
            ai  <- .getsubset(ai,  subset)
            n1i <- .getsubset(n1i, subset)
            n2i <- .getsubset(n2i, subset)
            mi  <- .getsubset(mi,  subset)
         }

         if (any(c(n1i, n2i) <= 0, na.rm=TRUE))
            stop(mstyle$stop("One or more group sizes are <= 0."))

         if (any(ai < 0, na.rm=TRUE) || any(ai > 1, na.rm=TRUE))
            stop(mstyle$stop("One or more AUC values are < 0 or > 1."))

         if (any(mi < 0, na.rm=TRUE) || any(mi > 1, na.rm=TRUE))
            stop(mstyle$stop("One or more 'mi' values are < 0 or > 1."))

         ni <- n1i + n2i

         ni.u <- ni # unadjusted total sample sizes

         k <- length(ai)

         yi <- ai

         vtype <- .expand1(vtype, k)

         vi <- rep(NA_real_, k)

         navgi <- (n1i+n2i)/2
         q0 <- ai*(1-ai)
         q1 <- ai/(2-ai)
         q2 <- 2*ai^2/(1+ai)

         if (!all(is.element(vtype, c("LS","LS2","H0","MAX"))))
            stop(mstyle$stop("For this outcome measure, 'vtype' must be 'LS' or 'LS2'."))

         for (i in seq_len(k)) {

            ### based on Newcombe (2006b) but using (n1i-1)*(n2i-1) in the denominator as given in Cho et al. (2019), section 2.4
            if (vtype[i] == "LS")
               vi[i] <- q0[i] / ((n1i[i]-1)*(n2i[i]-1)) * (2*navgi[i] - 1 - (3*navgi[i]-3) / ((2-ai[i])*(1+ai[i])))

            ### based on Hanley and McNeil (1982) but using (n1i-1)*(n2i-1) in the denominator and subtracting mi/4 as given in Cho et al. (2019)
            if (vtype[i] == "LS2")
               vi[i] <- (q0[i] - mi[i]/4 + (n1i[i]-1)*(q1[i]-ai[i]^2) + (n2i[i]-1)*(q2[i]-ai[i]^2)) / ((n1i[i]-1)*(n2i[i]-1))

            ### under H0: CLES=AUC=0.5 and equal variances (conservative if there are ties)
            if (vtype[i] == "H0")
               vi[i] <- (n1i[i]+n2i[i]+1)/(12*n1i[i]*n2i[i])

            ### based on sigma^2_max (eq. 7 in Bamber, 1975)
            if (vtype[i] == "MAX")
               vi[i] <- q0[i] / (min(n1i[i],n2i[i])-1)

         }

      }

      if (is.element(measure, c("CLESN","AUCN"))) {

         m1i  <- .getx("m1i",  mf=mf, data=data, checknumeric=TRUE)
         m2i  <- .getx("m2i",  mf=mf, data=data, checknumeric=TRUE)
         sd1i <- .getx("sd1i", mf=mf, data=data, checknumeric=TRUE)
         sd2i <- .getx("sd2i", mf=mf, data=data, checknumeric=TRUE)
         n1i  <- .getx("n1i",  mf=mf, data=data, checknumeric=TRUE)
         n2i  <- .getx("n2i",  mf=mf, data=data, checknumeric=TRUE)
         di   <- .getx("di",   mf=mf, data=data, checknumeric=TRUE)
         ti   <- .getx("ti",   mf=mf, data=data, checknumeric=TRUE)
         pi   <- .getx("pi",   mf=mf, data=data, checknumeric=TRUE)
         ai   <- .getx("ai",   mf=mf, data=data, checknumeric=TRUE)

         if (!.equal.length(m1i, m2i, sd1i, sd2i, n1i, n2i, di, ti, pi, ai))
            stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         if (!.all.specified(n1i, n2i))
            stop(mstyle$stop("Cannot compute outcomes. Check that all of the required information is specified\n  via the appropriate arguments."))

         k.all <- max(sapply(list(m1i, m2i, sd1i, sd2i, n1i, n2i, di, ti, pi, ai), length))

         vtype <- .expand1(vtype, k.all)

         ### if sd1i and/or sd2i have not been specified at all, set sd1i and sd2i to NA for all studies

         if (is.null(sd1i) || is.null(sd2i)) {
            sd1i <- .expand1(NA_real_, k.all)
            sd2i <- .expand1(NA_real_, k.all)
         }

         ### convert pi to ti values

         ti <- replmiss(ti, .convp2t(pi, df=n1i+n2i-2))

         ### convert ti to di values

         di <- replmiss(di, ti * sqrt(1/n1i + 1/n2i))

         ### for specified pi/ti/di values, assume homoscedasticity

         if (!is.null(di))
            vtype[!is.na(di)] <- "HO"

         ### compute di values from means and SDs (for these, do not assume homoscedasticity, unless vtype="HO")

         sdpi <- ifelse(vtype=="HO", sqrt(((n1i-1)*sd1i^2 + (n2i-1)*sd2i^2)/(n1i+n2i-2)), sqrt((sd1i^2 + sd2i^2)/2))
         di   <- replmiss(di, (m1i - m2i) / sdpi)

         ### convert di values to ai values and back (in case only ai is known, so we have di for computing vi)

         ai <- replmiss(ai, pnorm(di/sqrt(2)))
         di <- replmiss(di, qnorm(ai)*sqrt(2))

         k.all <- length(ai)

         ### if sd1i and/or sd2i is missing for a particular study, assume sd1i=sd2i=1 for that study and homoscedasticity

         sdsmiss <- is.na(sd1i) | is.na(sd2i)
         sd1i <- ifelse(sdsmiss, 1, sd1i)
         sd2i <- ifelse(sdsmiss, 1, sd2i)
         vtype[sdsmiss] <- "HO"

         if (!is.null(subset)) {
            subset <- .chksubset(subset, k.all)
            vtype  <- .getsubset(vtype,  subset)
            ai   <- .getsubset(ai,   subset)
            di   <- .getsubset(di,   subset)
            sd1i <- .getsubset(sd1i, subset)
            sd2i <- .getsubset(sd2i, subset)
            n1i  <- .getsubset(n1i,  subset)
            n2i  <- .getsubset(n2i,  subset)
         }

         if (any(c(sd1i, sd2i) < 0, na.rm=TRUE))
            stop(mstyle$stop("One or more standard deviations are negative."))

         if (any(c(n1i, n2i) <= 0, na.rm=TRUE))
            stop(mstyle$stop("One or more group sizes are <= 0."))

         if (any(ai < 0, na.rm=TRUE) || any(ai > 1, na.rm=TRUE))
            stop(mstyle$stop("One or more AUC values are < 0 or > 1."))

         ni.u <- n1i + n2i # unadjusted total sample sizes

         k <- length(ai)

         ni <- ni.u

         yi <- ai

         vi <- rep(NA_real_, k)

         if (!all(is.element(vtype, c("LS","LS2","LS3","HO"))))
            stop(mstyle$stop("For this outcome measure, 'vtype' must be 'LS' or 'HO'."))

         vri <- sd1i^2 / (sd1i^2 + sd2i^2)

         for (i in seq_len(k)) {

            ### large sample approximation to the sampling variance based on the binormal model
            if (vtype[i] == "LS") {
               vi[i] <- exp(-di[i]^2 / 2) / (8*base::pi) * (di[i]^2 * vri[i]^2 / (n1i[i]-1) + di[i]^2 * (1-vri[i])^2 / (n2i[i]-1) + 4*vri[i]/(n1i[i]-1) + 4*(1-vri[i])/(n2i[i]-1))
               # this is consistent with escalc(measure="SMDH", correct=FALSE) -> conv.delta(transf=transf.dtocles)
               #tmp <- escalc(measure="SMDH", m1i=m1i[i], sd1i=sd1i[i], n1i=n1i[i], m2i=m2i[i], sd2i=sd2i[i], n2i=n2i[i], correct=FALSE)
               #vi[i] <- conv.delta(yi, vi, data=tmp, transf=transf.dtocles, replace=TRUE)$vi
            }

            ### large sample approximation to the sampling variance based on the binormal model
            if (vtype[i] == "LS2") {
               vi[i] <- exp(-di[i]^2 / 2) / (8*base::pi) * (di[i]^2 * vri[i]^2 / (n1i[i]-0) + di[i]^2 * (1-vri[i])^2 / (n2i[i]-0) + 4*vri[i]/(n1i[i]-0) + 4*(1-vri[i])/(n2i[i]-0))
               # this is consistent with escalc(measure="SMDH", correct=FALSE, vtype="LS2") -> conv.delta(transf=transf.dtocles)
               #tmp <- escalc(measure="SMDH", m1i=m1i[i], sd1i=sd1i[i], n1i=n1i[i], m2i=m2i[i], sd2i=sd2i[i], n2i=n2i[i], correct=FALSE, vtype="LS2")
               #vi[i] <- conv.delta(yi, vi, data=tmp, transf=transf.dtocles, replace=TRUE)$vi
            }

            ### large sample approximation to the sampling variance based on the binormal model (based on standard application of the delta method)
            if (vtype[i] == "LS3") {
               vi[i] <- exp(-di[i]^2 / 2) / (8*base::pi) * (di[i]^2 * vri[i]^2 / (n1i[i]-1) + di[i]^2 * (1-vri[i])^2 / (n2i[i]-1) + 4*vri[i]/(n1i[i]-0) + 4*(1-vri[i])/(n2i[i]-0))
               # this is consistent with escalc(measure="SMDH", correct=FALSE, vtype="LS3") -> conv.delta(transf=transf.dtocles)
               #tmp <- escalc(measure="SMDH", m1i=m1i[i], sd1i=sd1i[i], n1i=n1i[i], m2i=m2i[i], sd2i=sd2i[i], n2i=n2i[i], correct=FALSE, vtype="LS3")
               #vi[i] <- conv.delta(yi, vi, data=tmp, transf=transf.dtocles, replace=TRUE)$vi
            }

            ### estimate assuming homoscedasticity of the variances within studies
            if (vtype[i] == "HO") {
               vi[i] <- exp(-di[i]^2 / 2) / (4*base::pi) * (1/n1i[i] + 1/n2i[i] + di[i]^2 / (2*(n1i[i]+n2i[i])))
               # this is consistent with escalc(measure="SMD", correct=FALSE) -> conv.delta(transf=transf.dtocles)
               #tmp <- escalc(measure="SMD", m1i=m1i[i], sd1i=sd1i[i], n1i=n1i[i], m2i=m2i[i], sd2i=sd2i[i], n2i=n2i[i], correct=FALSE)
               #vi[i] <- conv.delta(yi, vi, data=tmp, transf=transf.dtocles, replace=TRUE)$vi
            }

         }

      }

      ######################################################################

   } else {

      ### in case yi is not NULL (so user wants to convert a regular data frame to an 'escalc' object)

      ### check if yi is numeric

      if (!.is.numeric(yi))
         stop(mstyle$stop("The object/variable specified for the 'yi' argument is not numeric."))

      ### get vi, sei, and ni

      vi  <- .getx("vi",  mf=mf, data=data, checknumeric=TRUE)
      sei <- .getx("sei", mf=mf, data=data, checknumeric=TRUE)
      ni  <- .getx("ni",  mf=mf, data=data, checknumeric=TRUE)

      ### if neither vi nor sei is specified, then throw an error
      ### if only sei is specified, then square those values to get vi
      ### if vi is specified, use those values

      if (is.null(vi)) {
         if (is.null(sei)) {
            stop(mstyle$stop("Must specify 'vi' or 'sei' argument."))
         } else {
            vi <- sei^2
         }
      }

      if (!.equal.length(yi, vi, ni))
         stop(mstyle$stop("Supplied data vectors are not all of the same length."))

      k.all <- length(yi)

      ### if slab is NULL, see if we can get it from yi (subsetting is done further below; see [z])

      if (is.null(slab)) {
         slab <- attributes(yi)$slab
         if (length(slab) != k.all)
            slab <- NULL
      }

      if (!is.null(subset)) {
         subset <- .chksubset(subset, k.all)
         yi <- .getsubset(yi, subset)
         vi <- .getsubset(vi, subset)
         ni <- .getsubset(ni, subset)
      }

      ni.u <- ni # unadjusted total sample sizes

      k <- length(yi)

   }

   #########################################################################
   #########################################################################
   #########################################################################

   ### make sure yi and vi are really vectors (and not arrays)

   yi <- as.vector(yi)
   vi <- as.vector(vi)

   ### check for infinite values and set them to NA

   is.inf <- is.infinite(yi) | is.infinite(vi)

   if (any(is.inf)) {
      warning(mstyle$warning("Some 'yi' and/or 'vi' values equal to +-Inf. Recoded to NAs."), call.=FALSE)
      yi[is.inf] <- NA_real_
      vi[is.inf] <- NA_real_
   }

   ### check for NaN values and set them to NA

   is.NaN <- is.nan(yi) | is.nan(vi)

   if (any(is.NaN)) {
      yi[is.NaN] <- NA_real_
      vi[is.NaN] <- NA_real_
   }

   ### check for negative vi's (should not happen, but just in case)

   vi[vi < 0] <- NA_real_

   ### add study labels if specified

   if (!is.null(slab)) {

      if (length(slab) != k.all)
         stop(mstyle$stop(paste0("Length of the 'slab' argument (", length(slab), ") does not correspond to the size of the dataset (", k.all, ").")))

      if (is.factor(slab))
         slab <- as.character(slab)

      if (!is.null(subset))
         slab <- .getsubset(slab, subset) # [z]

      if (anyNA(slab))
         stop(mstyle$stop("NAs in study labels."))

      ### check if study labels are unique; if not, make them unique

      if (anyDuplicated(slab))
         slab <- .make.unique(slab)

   }

   ### if include/subset is NULL, set to TRUE vector

   if (is.null(include))
      include <- rep(TRUE, k.all)
   if (is.null(subset))
      subset <- rep(TRUE, k.all)

   ### turn numeric include vector into a logical vector (already done for subset)

   if (!is.null(include))
      include <- .chksubset(include, k.all, stoponk0=FALSE)

   ### apply subset to include

   include <- .getsubset(include, subset)

   ### process flip argument

   if (is.null(flip)) {
      flip <- rep(1, k.all)
   } else {
      if (is.logical(flip)) {
         flip <- .expand1(flip, k.all)
         flip <- flip %in% TRUE # so NAs are treated as FALSE
         flip <- ifelse(flip, -1, 1)
      }
   }

   if (length(flip) != k.all)
      stop(mstyle$stop(paste0("Length of the 'flip' argument (", length(flip), ") does not correspond to the size of the dataset (", k.all, ").")))

   flip <- .getsubset(flip, subset)

   yi[include] <- flip[include] * yi[include]
   vi[include] <- flip[include]^2 * vi[include]

   ### subset data frame (note: subsetting of other parts already done above, so yi/vi/ni.u/slab are already subsetted)

   if (has.data && any(!subset))
      data <- .getsubset(data, subset)

   ### put together dataset

   if (has.data && append) {

      ### if data argument has been specified and user wants to append

      dat <- data.frame(data)

      if (replace || !is.element(var.names[1], names(dat))) {
         yi.replace <- rep(TRUE, k)
      } else {
         yi.replace <- is.na(dat[[var.names[1]]])
      }

      if (replace || !is.element(var.names[2], names(dat))) {
         vi.replace <- rep(TRUE, k)
      } else {
         vi.replace <- is.na(dat[[var.names[2]]])
      }

      if (replace || !is.element(var.names[3], names(dat))) {
         measure.replace <- rep(TRUE, k)
      } else {
         measure.replace <- is.na(dat[[var.names[3]]]) | dat[[var.names[3]]] == ""
      }

      dat[[var.names[1]]][include & yi.replace] <- yi[include & yi.replace]
      dat[[var.names[2]]][include & vi.replace] <- vi[include & vi.replace]

      if (add.measure)
         dat[[var.names[3]]][!is.na(yi) & include & measure.replace] <- measure

      if (!is.null(ni.u))
         attributes(dat[[var.names[1]]])$ni[include & yi.replace] <- ni.u[include & yi.replace]

   } else {

      ### if data argument has not been specified or user does not want to append

      dat <- data.frame(yi=rep(NA_real_, k), vi=rep(NA_real_, k))
      dat$yi[include] <- yi[include]
      dat$vi[include] <- vi[include]

      if (add.measure)
         dat$measure[!is.na(yi) & include] <- measure

      attributes(dat$yi)$ni[include] <- ni.u[include]

      if (add.measure) {
         names(dat) <- var.names
      } else {
         names(dat) <- var.names[1:2]
      }

   }

   ### replace missings in measure with ""
   if (add.measure)
      dat[[var.names[3]]][is.na(dat[[var.names[3]]])] <- ""

   ### add slab attribute to the yi vector
   if (!is.null(slab))
      attr(dat[[var.names[1]]], "slab") <- slab

   ### add measure attribute to the yi vector
   attr(dat[[var.names[1]]], "measure") <- measure

   ### add digits attribute
   attr(dat, "digits") <- digits

   ### add vtype attribute
   #attr(dat, "vtype") <- vtype

   ### add 'yi.names' and 'vi.names' to the first position of the corresponding attributes (so the first is always the last one calculated/added)
   attr(dat, "yi.names") <- union(var.names[1], attr(data, "yi.names")) # if 'yi.names' is not an attribute, attr() returns NULL, so this works fine
   attr(dat, "vi.names") <- union(var.names[2], attr(data, "vi.names")) # if 'vi.names' is not an attribute, attr() returns NULL, so this works fine

   ### add 'out.names' back to object in case these attributes exist (if summary() has been used on the object)
   attr(dat, "sei.names")   <- attr(data, "sei.names")
   attr(dat, "zi.names")    <- attr(data, "zi.names")
   attr(dat, "pval.names")  <- attr(data, "pval.names")
   attr(dat, "ci.lb.names") <- attr(data, "ci.lb.names")
   attr(dat, "ci.ub.names") <- attr(data, "ci.ub.names")

   ### keep only attribute elements from yi.names and vi.names that are actually part of the object
   attr(dat, "yi.names") <- attr(dat, "yi.names")[attr(dat, "yi.names") %in% colnames(dat)]
   attr(dat, "vi.names") <- attr(dat, "vi.names")[attr(dat, "vi.names") %in% colnames(dat)]

   class(dat) <- c("escalc", "data.frame")
   return(dat)

}
