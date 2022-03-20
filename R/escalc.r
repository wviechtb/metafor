escalc <- function(measure, ai, bi, ci, di, n1i, n2i, x1i, x2i, t1i, t2i, m1i, m2i, sd1i, sd2i, xi, mi, ri, ti, sdi, r2i, ni, yi, vi, sei,
data, slab, subset, include, add=1/2, to="only0", drop00=FALSE, vtype="LS", var.names=c("yi","vi"), add.measure=FALSE, append=TRUE, replace=TRUE, digits, ...) {

   ### check argument specifications

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (missing(measure) && missing(yi))
      stop(mstyle$stop("Must specify an effect size or outcome measure via the 'measure' argument."))

   if (!missing(yi) && missing(measure))
      measure <- "GEN"

   if (!is.character(measure))
      stop(mstyle$stop("The 'measure' argument must be a character string."))

   if (!is.element(measure, c("RR","OR","PETO","RD","AS","PHI","YUQ","YUY","RTET", ### 2x2 table measures
                              "PBIT","OR2D","OR2DN","OR2DL",                       ### - transformations to SMD
                              "MPRD","MPRR","MPOR","MPORC","MPPETO","MPORM",       ### - measures for matched pairs / pre-post data
                              "IRR","IRD","IRSD",                                  ### two-group person-time data measures
                              "MD","SMD","SMDH","SMD1","ROM",                      ### two-group mean/SD measures
                              "CVR","VR",                                          ### coefficient of variation ratio, variability ratio
                              "RPB","RBIS","D2OR","D2ORN","D2ORL",                 ### - transformations to r_PB, r_BIS, and log(OR)
                              "COR","UCOR","ZCOR",                                 ### correlations (raw and r-to-z transformed)
                              "PCOR","ZPCOR","SPCOR",                              ### partial and semi-partial correlations
                              "PR","PLN","PLO","PAS","PFT",                        ### single proportions (and transformations thereof)
                              "IR","IRLN","IRS","IRFT",                            ### single-group person-time data (and transformations thereof)
                              "MN","MNLN","CVLN","SDLN","SMN",                     ### mean, log(mean), log(CV), log(SD), standardized mean
                              "MC","SMCC","SMCR","SMCRH","ROMC","CVRC","VRC",      ### raw/standardized mean change, log(ROM), CVR, and VR for dependent samples
                              "ARAW","AHW","ABT",                                  ### alpha (and transformations thereof)
                              "GEN")))
      stop(mstyle$stop("Unknown 'measure' specified."))

   if (!is.element(to, c("all","only0","if0all","none")))
      stop(mstyle$stop("Unknown 'to' argument specified."))

   if (any(!is.element(vtype, c("UB","LS","LS2","HO","ST","CS","AV","AVHO")), na.rm=TRUE)) ### vtype can be an entire vector, so use any() and na.rm=TRUE
      stop(mstyle$stop("Unknown 'vtype' argument specified."))

   if (add.measure) {

      if (length(var.names) == 2L)
         var.names <- c(var.names, "measure")

      if (length(var.names) != 3L)
         stop(mstyle$stop("Argument 'var.names' must be of length 2 or 3."))

      if (any(var.names != make.names(var.names, unique=TRUE))) {
         var.names <- make.names(var.names, unique=TRUE)
         warning(mstyle$warning(paste0("Argument 'var.names' does not contain syntactically valid variable names.\nVariable names adjusted to: var.names = c('", var.names[1], "', '", var.names[2], "', '", var.names[3], "').")), call.=FALSE)
      }

   } else {

      if (length(var.names) == 3L)
         var.names <- var.names[1:2]

      if (length(var.names) != 2L)
         stop(mstyle$stop("Argument 'var.names' must be of length 2."))

      if (any(var.names != make.names(var.names, unique=TRUE))) {
         var.names <- make.names(var.names, unique=TRUE)
         warning(mstyle$warning(paste0("Argument 'var.names' does not contain syntactically valid variable names.\nVariable names adjusted to: var.names = c('", var.names[1], "', '", var.names[2], "').")), call.=FALSE)
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

   onlyo1 <- ifelse(is.null(ddd$onlyo1), FALSE, ddd$onlyo1)
   addyi  <- ifelse(is.null(ddd$addyi),  TRUE,  ddd$addyi)
   addvi  <- ifelse(is.null(ddd$addvi),  TRUE,  ddd$addvi)

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

   ### for certain measures, set add=0 by default unless user explicitly sets the add argument

   addval <- mf[[match("add", names(mf))]]

   if (is.element(measure, c("AS","PHI","RTET","IRSD","PAS","PFT","IRS","IRFT")) && is.null(addval))
      add <- 0

   #########################################################################
   #########################################################################
   #########################################################################

   if (is.null(yi)) {

      if (is.element(measure, c("RR","OR","RD","AS","PETO","PHI","YUQ","YUY","RTET","PBIT","OR2D","OR2DN","OR2DL","MPRD","MPRR","MPOR","MPORC","MPPETO","MPORM"))) {

         mf.ai <- mf[[match("ai", names(mf))]]
         if (any("~" %in% as.character(mf.ai)))
            stop(mstyle$stop("The 'formula interface' to escalc() has been deprecated."))

         ai  <- .getx("ai",  mf=mf, data=data)
         bi  <- .getx("bi",  mf=mf, data=data)
         ci  <- .getx("ci",  mf=mf, data=data)
         di  <- .getx("di",  mf=mf, data=data)
         n1i <- .getx("n1i", mf=mf, data=data)
         n2i <- .getx("n2i", mf=mf, data=data)
         ri  <- .getx("ri",  mf=mf, data=data)

         if (!.equal.length(ai, bi, ci, di, n1i, n2i, ri))
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
            stop(mstyle$stop("Cannot compute outcomes. Check that all of the required\n  information is specified via the appropriate arguments."))

         if (measure == "MPORM" && !.all.specified(ri))
            stop(mstyle$stop("Cannot compute outcomes. Check that all of the required\n  information is specified via the appropriate arguments."))

         k.all <- length(ai)

         if (!is.null(subset)) {
            subset <- .setnafalse(subset, k=k.all)
            ai <- ai[subset]
            bi <- bi[subset]
            ci <- ci[subset]
            di <- di[subset]
         }

         n1i <- ai + bi
         n2i <- ci + di

         if (any(c(ai > n1i, ci > n2i), na.rm=TRUE))
            stop(mstyle$stop("One or more event counts are larger than the corresponding group sizes."))

         if (any(c(ai, bi, ci, di) < 0, na.rm=TRUE))
            stop(mstyle$stop("One or more counts are negative."))

         if (any(c(n1i < 0, n2i < 0), na.rm=TRUE)) # note: in cross-sectional sampling, group sizes could be 0
            stop(mstyle$stop("One or more group sizes are negative."))

         if (measure == "MPORM" && any(abs(ri) > 1, na.rm=TRUE))
               stop(mstyle$stop("One or more correlations are > 1 or < -1."))

         ni.u <- ai + bi + ci + di ### unadjusted total sample sizes

         if (measure == "MPORM")
            ni.u <- round(ni.u / 2)

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
         ni  <- n1i + n2i ### ni.u computed earlier is always the 'unadjusted' total sample size

         if (measure == "MPORM")
            ni <- round(ni / 2)

         ### compute proportions for the two groups (unadjusted and adjusted)

         p1i.u <- ai.u/n1i.u
         p2i.u <- ci.u/n2i.u
         p1i <- ai/n1i
         p2i <- ci/n2i

         ### log risk ratios

         if (measure == "RR") {
            if (addyi) {
               yi <- log(p1i) - log(p2i)
            } else {
               yi <- log(p1i.u) - log(p2i.u)
            }
            if (addvi) {
               vi <- 1/ai - 1/n1i + 1/ci - 1/n2i
            } else {
               vi <- 1/ai.u - 1/n1i.u + 1/ci.u - 1/n2i.u
            }
         }

         ### log odds ratio

         if (is.element(measure, c("OR","OR2D","OR2DN","OR2DL","MPORM"))) {
            if (addyi) {
               yi <- log(p1i/(1-p1i)) - log(p2i/(1-p2i))
            } else {
               yi <- log(p1i.u/(1-p1i.u)) - log(p2i.u/(1-p2i.u))
            }
            if (addvi) {
               vi <- 1/ai + 1/bi + 1/ci + 1/di
            } else {
               vi <- 1/ai.u + 1/bi.u + 1/ci.u + 1/di.u
            }
         }

         ### risk difference

         if (measure == "RD") {

            if (addyi) {
               yi <- p1i - p2i
            } else {
               yi <- p1i.u - p2i.u
            }
            if (length(vtype) == 1L)
               vtype <- rep(vtype, k)

            vi <- rep(NA_real_, k)

            if (addvi) {
               mnwp1i <- .wmean(p1i, n1i, na.rm=TRUE) ### sample size weighted average of proportions
               mnwp2i <- .wmean(p2i, n2i, na.rm=TRUE) ### sample size weighted average of proportions
            } else {
               mnwp1i.u <- .wmean(p1i.u, n1i.u, na.rm=TRUE) ### sample size weighted average of proportions
               mnwp2i.u <- .wmean(p2i.u, n2i.u, na.rm=TRUE) ### sample size weighted average of proportions
            }

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

               ### estimator assuming homogeneity (using the average proportions)
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
            xt <- ai + ci ### frequency of outcome1 in both groups combined
            yt <- bi + di ### frequency of outcome2 in both groups combined
            Ei <- xt * n1i / ni
            Vi <- xt * yt * (n1i/ni) * (n2i/ni) / (ni - 1) ### 0 when xt = 0 or yt = 0 in a table
            yi <- (ai - Ei) / Vi                           ### then yi and vi is Inf (set to NA at end)
            vi <- 1/Vi
         }

         ### arcsine square root risk difference

         if (measure == "AS") {
            yi <- asin(sqrt(p1i)) - asin(sqrt(p2i))
            vi <- 1/(4*n1i) + 1/(4*n2i)
         }

         ### phi coefficient

         if (measure == "PHI") {

            yi <- (ai*di - bi*ci)/sqrt((ai+bi)*(ci+di)*(ai+ci)*(bi+di))

            if (length(vtype) == 1L)
               vtype <- rep(vtype, k)

            vi <- rep(NA_real_, k)

            q1i <- 1 - p1i
            q2i <- 1 - p2i
            pi1. <- (ai+bi)/ni
            pi2. <- (ci+di)/ni
            pi.1 <- (ai+ci)/ni
            pi.2 <- (bi+di)/ni

            if (!all(is.element(vtype, c("ST","LS","CS"))))
               stop(mstyle$stop("For this outcome measure, 'vtype' must be either 'ST', 'LS', or 'CS'."))

            for (i in seq_len(k)) {

               ### estimate of the sampling variance for stratified sampling
               if (vtype[i] == "ST") {
                  vi[i] <- ((n1i[i]+n2i[i])^2*(4*n1i[i]^3*p1i[i]^2*p2i[i]*q1i[i]^2*q2i[i] + 4*n2i[i]^3*p1i[i]*p2i[i]^2*q1i[i]*q2i[i]^2 + n1i[i]*n2i[i]^2*p2i[i]*q2i[i]*(p2i[i]*q1i[i] + p1i[i]*q2i[i])*(p2i[i]*q1i[i] + p1i[i]*(4*q1i[i] + q2i[i])) + n1i[i]^2*n2i[i]*p1i[i]*q1i[i]*(p2i[i]*q1i[i] + p1i[i]*q2i[i])*(p1i[i]*q2i[i] + p2i[i]*(q1i[i] + 4*q2i[i]))))/(4*(ai[i]+ci[i])^3*(bi[i]+di[i])^3)
               }

               ### estimate of the sampling variance for cross-sectional/multinomial sampling (equation in Yule, 1912, p.603)
               if (vtype[i] == "LS" || vtype[i] == "CS") {
                  vi[i] <- 1/ni[i] * (1 - yi[i]^2 + yi[i]*(1+1/2*yi[i]^2) * (pi1.[i]-pi2.[i])*(pi.1[i]-pi.2[i]) / sqrt(pi1.[i]*pi2.[i]*pi.1[i]*pi.2[i]) - 3/4 * yi[i]^2 * ((pi1.[i]-pi2.[i])^2/(pi1.[i]*pi2.[i]) + (pi.1[i]-pi.2[i])^2/(pi.1[i]*pi.2[i])))
               }

            }

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
                                                                                    ### see code/probit_transformation directory
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

         ### matched pairs / pre-post 2x2 table measures

         if (is.element(measure, c("MPRD","MPRR","MPOR"))) {
            pi12 <- bi/ni
            pi21 <- ci/ni
            pi1. <- (ai+bi)/ni
            pi.1 <- (ai+ci)/ni
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

         x1i <- .getx("x1i", mf=mf, data=data)
         x2i <- .getx("x2i", mf=mf, data=data)
         t1i <- .getx("t1i", mf=mf, data=data)
         t2i <- .getx("t2i", mf=mf, data=data)

         if (!.all.specified(x1i, x2i, t1i, t2i))
            stop(mstyle$stop("Cannot compute outcomes. Check that all of the required\n  information is specified via the appropriate arguments."))

         if (!.equal.length(x1i, x2i, t1i, t2i))
            stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         k.all <- length(x1i)

         if (!is.null(subset)) {
            subset <- .setnafalse(subset, k=k.all)
            x1i <- x1i[subset]
            x2i <- x2i[subset]
            t1i <- t1i[subset]
            t2i <- t2i[subset]
         }

         if (any(c(x1i, x2i) < 0, na.rm=TRUE))
            stop(mstyle$stop("One or more counts are negative."))

         if (any(c(t1i, t2i) <= 0, na.rm=TRUE))
            stop(mstyle$stop("One or more person-times are <= 0."))

         ni.u <- t1i + t2i ### unadjusted total sample sizes

         k <- length(x1i)

         ### if drop00=TRUE, set counts to NA for studies that have no events in both arms

         if (drop00) {
            id00 <- c(x1i == 0L & x2i == 0L)
            id00[is.na(id00)] <- FALSE
            x1i[id00] <- NA
            x2i[id00] <- NA
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
               vi <- ir1i/t1i + ir2i/t2i ### note: same as x1i/t1i^2 + x2i/t2i^2
            } else {
               vi <- ir1i.u/t1i + ir2i.u/t2i ### note: same as x1i.u/t1i^2 + x2i.u/t2i^2
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

      if (is.element(measure, c("MD","SMD","SMDH","SMD1","ROM","RPB","RBIS","D2OR","D2ORN","D2ORL","CVR","VR"))) {

         m1i  <- .getx("m1i",  mf=mf, data=data) ### for VR, do not need to supply this
         m2i  <- .getx("m2i",  mf=mf, data=data) ### for VR, do not need to supply this
         sd1i <- .getx("sd1i", mf=mf, data=data) ### for SMD1, do not need to supply this
         sd2i <- .getx("sd2i", mf=mf, data=data)
         n1i  <- .getx("n1i",  mf=mf, data=data)
         n2i  <- .getx("n2i",  mf=mf, data=data)

         ### for these measures, need m1i, m2i, sd1i, sd2i, n1i, and n2i

         if (is.element(measure, c("MD","SMD","SMDH","ROM","RPB","RBIS","D2OR","D2ORN","D2ORL","CVR"))) {

            if (!.all.specified(m1i, m2i, sd1i, sd2i, n1i, n2i))
               stop(mstyle$stop("Cannot compute outcomes. Check that all of the required\n  information is specified via the appropriate arguments."))

            if (!.equal.length(m1i, m2i, sd1i, sd2i, n1i, n2i))
               stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         }

         ### for this measure, need sd1i, sd2i, n1i, and n2i

         if (measure == "VR") {

            if (!.all.specified(sd1i, sd2i, n1i, n2i))
               stop(mstyle$stop("Cannot compute outcomes. Check that all of the required\n  information is specified via the appropriate arguments."))

            if (!.equal.length(sd1i, sd2i, n1i, n2i))
               stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         }

         ### for this measure, need m1i, m2i, sdi, n1i, and n2i

         if (measure == "SMD1") {

            if (!.all.specified(m1i, m2i, sd2i, n1i, n2i))
               stop(mstyle$stop("Cannot compute outcomes. Check that all of the required\n  information is specified via the appropriate arguments."))

            if (!.equal.length(m1i, m2i, sd2i, n1i, n2i))
               stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         }

         k.all <- length(n1i)

         if (!is.null(subset)) {
            subset <- .setnafalse(subset, k=k.all)
            m1i  <- m1i[subset]
            m2i  <- m2i[subset]
            sd1i <- sd1i[subset]
            sd2i <- sd2i[subset]
            n1i  <- n1i[subset]
            n2i  <- n2i[subset]
         }

         if (any(c(sd1i, sd2i) < 0, na.rm=TRUE))
            stop(mstyle$stop("One or more standard deviations are negative."))

         if (any(c(n1i, n2i) <= 0, na.rm=TRUE))
            stop(mstyle$stop("One or more group sizes are <= 0."))

         ni.u <- n1i + n2i ### unadjusted total sample sizes

         k <- length(n1i)

         ni <- ni.u
         if (measure == "SMD1") {
            mi <- n2i - 1
            sdpi <- sd2i
            npi <- n2i
         } else {
            mi <- ni - 2
            sdpi <- sqrt(((n1i-1)*sd1i^2 + (n2i-1)*sd2i^2)/mi)
            npi <- ni
         }
         di <- (m1i - m2i) / sdpi

         ### (raw) mean difference (with heteroscedastic variances)
         ### to use with pooled SDs, simply set sd1i = sd2i = sdpi or use vtype="HO"

         if (measure == "MD") {

            yi <- m1i - m2i

            if (length(vtype) == 1L)
               vtype <- rep(vtype, k)

            vi <- rep(NA_real_, k)

            if (!all(is.element(vtype, c("UB","LS","HO"))))
               stop(mstyle$stop("For this outcome measure, 'vtype' must be either 'UB', 'LS', or 'HO'."))

            for (i in seq_len(k)) {

               ### unbiased estimate of the sampling variance (does not assume homoscedasticity)
               if (vtype[i] == "UB" || vtype[i] == "LS")
                  vi[i] <- sd1i[i]^2/n1i[i] + sd2i[i]^2/n2i[i]

               ### estimator assuming homoscedasticity
               if (vtype[i] == "HO")
                  vi[i] <- sdpi[i]^2 * (1/n1i[i] + 1/n2i[i])

            }

         }

         ### standardized mean difference (with pooled SDs or just the SD of group 2)

         if (is.element(measure, c("SMD","SMD1"))) {

            ### apply bias-correction to di values

            cmi <- .cmicalc(mi)
            yi <- cmi * di

            if (length(vtype) == 1L)
               vtype <- rep(vtype, k)

            vi <- rep(NA_real_, k)

            mnwyi <- .wmean(yi, ni, na.rm=TRUE) ### sample size weighted average of yi's

            if (!all(is.element(vtype, c("UB","LS","LS2","AV"))))
               stop(mstyle$stop("For this outcome measure, 'vtype' must be either 'UB', 'LS', 'LS2', or 'AV'."))

            for (i in seq_len(k)) {

               ### unbiased estimate of the sampling variance
               if (vtype[i] == "UB")
                  vi[i] <- 1/n1i[i] + 1/n2i[i] + (1 - (mi[i]-2)/(mi[i]*cmi[i]^2)) * yi[i]^2

               ### large sample approximation to the sampling variance
               if (vtype[i] == "LS")
                  vi[i] <- 1/n1i[i] + 1/n2i[i] + yi[i]^2/(2*npi[i])

               ### estimator assuming homogeneity (using sample size weighted average of the yi's)
               if (vtype[i] == "AV")
                  vi[i] <- 1/n1i[i] + 1/n2i[i] + mnwyi^2/(2*npi[i])

               ### large sample approximation to the sampling variance (using equation 4.24 from Borenstein, 2009)
               if (vtype[i] == "LS2")
                  vi[i] <- cmi[i]^2 * (1/n1i[i] + 1/n2i[i] + di[i]^2/(2*npi[i]))

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

            mn1wcvi <- .wmean(sd1i/m1i, n1i, na.rm=TRUE) ### sample size weighted average of the coefficient of variation in group 1
            mn2wcvi <- .wmean(sd2i/m2i, n2i, na.rm=TRUE) ### sample size weighted average of the coefficient of variation in group 2
            not.na  <- !(is.na(n1i) | is.na(n2i) | is.na(sd1i/m1i) | is.na(sd2i/m2i))
            mnwcvi  <- (sum(n1i[not.na]*(sd1i/m1i)[not.na]) + sum(n2i[not.na]*(sd2i/m2i)[not.na])) / sum((n1i+n2i)[not.na]) ### sample size weighted average of the two CV values

            if (!all(is.element(vtype, c("LS","HO","AV","AVHO"))))
               stop(mstyle$stop("For this outcome measure, 'vtype' must be either 'LS', 'HO', 'AV', or 'AVHO'."))

            for (i in seq_len(k)) {

               ### large sample approximation to the sampling variance (does not assume homoscedasticity)
               if (vtype[i] == "LS")
                  vi[i] <- sd1i[i]^2/(n1i[i]*m1i[i]^2) + sd2i[i]^2/(n2i[i]*m2i[i]^2)

               ### estimator assuming homoscedasticity
               if (vtype[i] == "HO")
                  vi[i] <- sdpi[i]^2/(n1i[i]*m1i[i]^2) + sdpi[i]^2/(n2i[i]*m2i[i]^2)

               ### estimator using the weighted averages of the CV values
               if (vtype[i] == "AV")
                  vi[i] <- mn1wcvi^2/n1i[i] + mn2wcvi^2/n2i[i]

               ### estimator using the weighted average of two weighted averages of the CV values
               if (vtype[i] == "AVHO")
                  vi[i] <- mnwcvi^2 * (1/n1i[i] + 1/n2i[i])

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

               if (!all(is.element(vtype, c("ST","LS","CS"))))
                  stop(mstyle$stop("For this outcome measure, 'vtype' must be either 'ST', 'LS', or 'CS'."))

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

         ### coefficient of variation ratio
         ### note: vi computed as per eq. 12 from Nakagawa et al. (2015), but without the '-2 rho ...' terms,
         ### since for normally distributed data the mean and variance (and transformations thereof) are independent

         if (measure == "CVR") {
            yi <- log(sd1i/m1i) + 1/(2*(n1i-1)) - log(sd2i/m2i) - 1/(2*(n2i-1))
            vi <- 1/(2*(n1i-1)) + sd1i^2/(n1i*m1i^2) + 1/(2*(n2i-1)) + sd2i^2/(n2i*m2i^2)
         }

         ### variability ratio

         if (measure == "VR") {
            yi <- log(sd1i/sd2i) + 1/(2*(n1i-1)) - 1/(2*(n2i-1))
            vi <- 1/(2*(n1i-1)) + 1/(2*(n2i-1))
         }

      }

      ######################################################################

      if (is.element(measure, c("COR","UCOR","ZCOR"))) {

         ri <- .getx("ri", mf=mf, data=data)
         ni <- .getx("ni", mf=mf, data=data)

         if (!.all.specified(ri, ni))
            stop(mstyle$stop("Cannot compute outcomes. Check that all of the required\n  information is specified via the appropriate arguments."))

         if (!.equal.length(ri, ni))
            stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         k.all <- length(ri)

         if (!is.null(subset)) {
            subset <- .setnafalse(subset, k=k.all)
            ri <- ri[subset]
            ni <- ni[subset]
         }

         if (any(abs(ri) > 1, na.rm=TRUE))
            stop(mstyle$stop("One or more correlations are > 1 or < -1."))

         if (any(ni <= 0, na.rm=TRUE))
            stop(mstyle$stop("One or more sample sizes are <= 0."))

         if (measure != "UCOR" && vtype == "UB")
            stop(mstyle$stop("Use of vtype='UB' only permitted when measure='UCOR'."))

         if (measure == "UCOR" && any(ni <= 4, na.rm=TRUE))
            warning(mstyle$warning("Cannot compute the bias-corrected correlation coefficient when ni <= 4."), call.=FALSE)

         if (measure == "ZCOR" && any(ni <= 3, na.rm=TRUE))
            warning(mstyle$warning("Cannot estimate the sampling variance when ni <= 3."), call.=FALSE)

         ni.u <- ni ### unadjusted total sample sizes

         k <- length(ri)

         ### raw correlation coefficient

         if (measure == "COR")
            yi <- ri

         ### raw correlation coefficient with bias correction

         if (measure == "UCOR") {
            #yi <- ri + ri*(1-ri^2)/(2*(ni-4)) ### approximation
            #yi[ni <= 4] <- NA ### set corrected correlations for ni <= 4 to NA
            yi <- ri * .Fcalc(1/2, 1/2, (ni-2)/2, 1-ri^2)
         }

         ### sampling variances for COR or UCOR

         if (is.element(measure, c("COR","UCOR"))) {

            if (length(vtype) == 1L)
               vtype <- rep(vtype, k)

            vi <- rep(NA_real_, k)

            mnwyi <- .wmean(yi, ni, na.rm=TRUE) ### sample size weighted average of yi's

            if (!all(is.element(vtype, c("UB","LS","AV"))))
               stop(mstyle$stop("For this outcome measure, 'vtype' must be either 'UB', 'LS', or 'AV'."))

            for (i in seq_len(k)) {

               ### unbiased estimate of the sampling variance of the bias-corrected correlation coefficient
               if (vtype[i] == "UB") {
                  #vi[i] <- yi[i]^2 - 1 + (ni[i]-3)/(ni[i]-2) * ((1-ri[i]^2) + 2*(1-ri[i]^2)^2/ni[i] + 8*(1-ri[i]^2)^3/(ni[i]*(ni[i]+2)) + 48*(1-ri[i]^2)^4/(ni[i]*(ni[i]+2)*(ni[i]+4)))
                  vi[i] <- yi[i]^2 - (1 - (ni[i]-3)/(ni[i]-2) * (1-ri[i]^2) * .Fcalc(1, 1, ni[i]/2, 1-ri[i]^2))
               }

               ### large sample approximation to the sampling variance
               if (vtype[i] == "LS")
                  vi[i] <- (1-yi[i]^2)^2/(ni[i]-1)

               ### estimator assuming homogeneity (using sample size weighted average of the yi's)
               if (vtype[i] == "AV")
                  vi[i] <- (1-mnwyi^2)^2/(ni[i]-1)

            }

         }

         ### r-to-z transformed correlation

         if (measure == "ZCOR") {
            yi <- 1/2 * log((1+ri)/(1-ri))
            vi <- 1/(ni-3)
         }

         ### set sampling variances for ni <= 3 to NA

         vi[ni <= 3] <- NA

      }

      ######################################################################

      if (is.element(measure, c("PCOR","ZPCOR","SPCOR"))) {

         ti  <- .getx("ti",  mf=mf, data=data)
         r2i <- .getx("r2i", mf=mf, data=data)
         mi  <- .getx("mi",  mf=mf, data=data)
         ni  <- .getx("ni",  mf=mf, data=data)

         if (measure=="PCOR" && !.all.specified(ti, ni, mi))
            stop(mstyle$stop("Cannot compute outcomes. Check that all of the required\n  information is specified via the appropriate arguments."))

         if (measure=="SPCOR" && !.all.specified(ti, ni, mi, r2i))
            stop(mstyle$stop("Cannot compute outcomes. Check that all of the required\n  information is specified via the appropriate arguments."))

         if (measure=="PCOR" && !.equal.length(ti, ni, mi))
            stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         if (measure=="SPCOR" && !.equal.length(ti, ni, mi, r2i))
            stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         k.all <- length(ti)

         if (!is.null(subset)) {
            subset <- .setnafalse(subset, k=k.all)
            ti  <- ti[subset]
            r2i <- r2i[subset]
            mi  <- mi[subset]
            ni  <- ni[subset]
         }

         if (measure=="SPCOR" && any(r2i > 1 | r2i < 0, na.rm=TRUE))
            stop(mstyle$stop("One or more R^2 values are > 1 or < 0."))

         if (any(ni <= 0, na.rm=TRUE))
            stop(mstyle$stop("One or more sample sizes are <= 0."))

         if (any(mi < 0, na.rm=TRUE))
            stop(mstyle$stop("One or more mi values are negative."))

         if (any(ni - mi - 1 <= 0, na.rm=TRUE))
            stop(mstyle$stop("One or more dfs are <= 0."))

         ni.u <- ni ### unadjusted total sample sizes

         k <- length(ti)

         ### partial correlation coefficient

         if (measure == "PCOR") {

            yi <- ti / sqrt(ti^2 + (ni - mi - 1))

            if (length(vtype) == 1L)
               vtype <- rep(vtype, k)

            vi <- rep(NA_real_, k)

            mnwyi <- .wmean(yi, ni, na.rm=TRUE) ### sample size weighted average of yi's

            if (!all(is.element(vtype, c("LS","AV"))))
               stop(mstyle$stop("For this outcome measure, 'vtype' must be either 'LS' or 'AV'."))

            for (i in seq_len(k)) {

               ### large sample approximation to the sampling variance
               if (vtype[i] == "LS")
                  vi[i] <- (1 - yi[i]^2)^2 / (ni[i] - mi[i] - 1)

               ### estimator assuming homogeneity (using sample size weighted average of the yi's)
               if (vtype[i] == "AV")
                  vi[i] <- (1 - mnwyi^2)^2 / (ni[i] - mi[i] - 1)

            }

         }

         ### r-to-z transformed partial correlation

         if (measure == "ZPCOR") {
            yi <- ti / sqrt(ti^2 + (ni - mi - 1))
            yi <- 1/2 * log((1+yi)/(1-yi))
            vi <- 1/(ni-mi-1)
         }

         ### semi-partial correlation coefficient

         if (measure == "SPCOR") {

            yi <- ti * sqrt(1 - r2i) / sqrt(ni - mi - 1)

            if (length(vtype) == 1L)
               vtype <- rep(vtype, k)

            vi <- rep(NA_real_, k)

            mnwyi <- .wmean(yi, ni, na.rm=TRUE) ### sample size weighted average of yi's

            if (!all(is.element(vtype, c("LS","AV"))))
               stop(mstyle$stop("For this outcome measure, 'vtype' must be either 'LS' or 'AV'."))

            for (i in seq_len(k)) {

               ### large sample approximation to the sampling variance
               if (vtype[i] == "LS")
                  vi[i] <- (r2i[i]^2 - 2*r2i[i] + (r2i[i] - yi[i]^2) + 1 - (r2i[i] - yi[i]^2)^2) / ni[i]

               ### estimator assuming homogeneity (using sample size weighted average of the yi's)
               if (vtype[i] == "AV")
                  vi[i] <- (r2i[i]^2 - 2*r2i[i] + (r2i[i] - mnwyi^2) + 1 - (r2i[i] - mnwyi^2)^2) / ni[i]

            }

         }

      }

      ######################################################################

      if (is.element(measure, c("PR","PLN","PLO","PAS","PFT"))) {

         xi <- .getx("xi", mf=mf, data=data)
         mi <- .getx("mi", mf=mf, data=data)
         ni <- .getx("ni", mf=mf, data=data)

         if (!.equal.length(xi, mi, ni))
            stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         ni.inc <- ni != xi + mi

         if (any(ni.inc, na.rm=TRUE))
            stop(mstyle$stop("One or more 'ni' values are not equal to 'xi + mi'."))

         mi <- replmiss(mi, ni-xi)

         if (!.all.specified(xi, mi))
            stop(mstyle$stop("Cannot compute outcomes. Check that all of the required\n  information is specified via the appropriate arguments."))

         k.all <- length(xi)

         if (!is.null(subset)) {
            subset <- .setnafalse(subset, k=k.all)
            xi <- xi[subset]
            mi <- mi[subset]
         }

         ni <- xi + mi

         if (any(xi > ni, na.rm=TRUE))
            stop(mstyle$stop("One or more event counts are larger than the corresponding group sizes."))

         if (any(c(xi, mi) < 0, na.rm=TRUE))
            stop(mstyle$stop("One or more counts are negative."))

         if (any(ni <= 0, na.rm=TRUE))
            stop(mstyle$stop("One or more group sizes are <= 0."))

         ni.u <- ni ### unadjusted total sample sizes

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

            if (length(vtype) == 1L)
               vtype <- rep(vtype, k)

            vi <- rep(NA_real_, k)

            if (addvi) {
               mnwpri <- .wmean(pri, ni, na.rm=TRUE) ### sample size weighted average of proportions
            } else {
               mnwpri.u <- .wmean(pri.u, ni.u, na.rm=TRUE) ### sample size weighted average of proportions
            }

            if (!all(is.element(vtype, c("UB","LS","AV"))))
               stop(mstyle$stop("For this outcome measure, 'vtype' must be either 'UB', 'LS', or 'AV'."))

            for (i in seq_len(k)) {

               ### unbiased estimate of the sampling variance
               if (vtype[i] == "UB") {
                  if (addvi) {
                     vi[i] <- pri[i]*(1-pri[i])/(ni[i]-1)
                  } else {
                     vi[i] <- pri.u[i]*(1-pri.u[i])/(ni.u[i]-1)
                  }
               }

               ### large sample approximation to the sampling variance
               if (vtype[i] == "LS") {
                  if (addvi) {
                     vi[i] <- pri[i]*(1-pri[i])/ni[i]
                  } else {
                     vi[i] <- pri.u[i]*(1-pri.u[i])/ni.u[i]
                  }
               }

               ### estimator assuming homogeneity (using the average proportion)
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

            if (length(vtype) == 1L)
               vtype <- rep(vtype, k)

            vi <- rep(NA_real_, k)

            if (addvi) {
               mnwpri <- .wmean(pri, ni, na.rm=TRUE) ### sample size weighted average of proportions
               #mnwpri <- exp(.wmean(yi, ni, na.rm=TRUE)) ### alternative strategy (exp of the sample size weighted average of the log proportions)
            } else {
               mnwpri.u <- .wmean(pri.u, ni.u, na.rm=TRUE) ### sample size weighted average of proportions
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

               ### estimator assuming homogeneity (using the average proportion)
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

            if (length(vtype) == 1L)
               vtype <- rep(vtype, k)

            vi <- rep(NA_real_, k)

            if (addvi) {
               mnwpri <- .wmean(pri, ni, na.rm=TRUE) ### sample size weighted average of proportions
               #mnwpri <- transf.ilogit(.wmean(yi, ni, na.rm=TRUE)) ### alternative strategy (inverse logit of the sample size weighted average of the logit transformed proportions)
            } else {
               mnwpri.u <- .wmean(pri.u, ni.u, na.rm=TRUE) ### sample size weighted average of proportions
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

               ### estimator assuming homogeneity (using the average proportion)
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

         ### proportion with arcsine square root (angular) transformation

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

         xi <- .getx("xi", mf=mf, data=data)
         ti <- .getx("ti", mf=mf, data=data)

         if (!.all.specified(xi, ti))
            stop(mstyle$stop("Cannot compute outcomes. Check that all of the required\n  information is specified via the appropriate arguments."))

         if (!.equal.length(xi, ti))
            stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         k.all <- length(xi)

         if (!is.null(subset)) {
            subset <- .setnafalse(subset, k=k.all)
            xi <- xi[subset]
            ti <- ti[subset]
         }

         if (any(xi < 0, na.rm=TRUE))
            stop(mstyle$stop("One or more counts are negative."))

         if (any(ti <= 0, na.rm=TRUE))
            stop(mstyle$stop("One or more person-times are <= 0."))

         ni.u <- ti ### unadjusted total sample sizes

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

         iri.u <- xi.u/ti
         iri <- xi/ti

         ### raw incidence rate

         if (measure == "IR") {
            if (addyi) {
               yi <- iri
            } else {
               yi <- iri.u
            }
            if (addvi) {
               vi <- iri/ti ### note: same as xi/ti^2
            } else {
               vi <- iri.u/ti ### note: same as xi.u/ti^2
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
               vi <- 1/xi
            } else {
               vi <- 1/xi.u
            }
         }

         ### square root transformed incidence rate

         if (measure == "IRS") {
            if (addyi) {
               yi <- sqrt(iri)
            } else {
               yi <- sqrt(iri.u)
            }
            vi <- 1/(4*ti)
         }

         ### note: addyi and addvi only implemented for measures above

         ### incidence rate with Freeman-Tukey transformation

         if (measure == "IRFT") {
            yi <- 1/2*(sqrt(iri) + sqrt(iri+1/ti))
            vi <- 1/(4*ti)
         }

      }

      ######################################################################

      if (is.element(measure, c("MN","MNLN","CVLN","SDLN","SMN"))) {

         mi  <- .getx("mi",  mf=mf, data=data) ### for SDLN, do not need to supply this
         sdi <- .getx("sdi", mf=mf, data=data)
         ni  <- .getx("ni",  mf=mf, data=data)

         ### for these measures, need mi, sdi, and ni

         if (is.element(measure, c("MN","MNLN","CVLN","SMN"))) {

            if (!.all.specified(mi, sdi, ni))
               stop(mstyle$stop("Cannot compute outcomes. Check that all of the required\n  information is specified via the appropriate arguments."))

            if (!.equal.length(mi, sdi, ni))
               stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         }

         ### for this measure, need sdi and ni

         if (measure == "SDLN") {

            if (!.all.specified(sdi, ni))
               stop(mstyle$stop("Cannot compute outcomes. Check that all of the required\n  information is specified via the appropriate arguments."))

            if (!.equal.length(sdi, ni))
               stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         }

         k.all <- length(ni)

         if (!is.null(subset)) {
            subset <- .setnafalse(subset, k=k.all)
            mi  <- mi[subset]
            sdi <- sdi[subset]
            ni  <- ni[subset]
         }

         if (any(sdi < 0, na.rm=TRUE))
            stop(mstyle$stop("One or more standard deviations are negative."))

         if (any(ni <= 0, na.rm=TRUE))
            stop(mstyle$stop("One or more sample sizes are <= 0."))

         if (is.element(measure, c("MNLN","CVLN")) && any(mi < 0, na.rm=TRUE))
            stop(mstyle$stop("One or more means are negative."))

         ni.u <- ni ### unadjusted total sample sizes

         k <- length(ni)

         ### (raw) mean

         if (measure == "MN") {
            yi <- mi
            vi <- sdi^2/ni
         }

         ### log(mean)

         if (measure == "MNLN") {
            yi <- log(mi)
            vi <- sdi^2/(ni*mi^2)
         }

         ### log(CV) with bias correction
         ### note: vi computed as per eq. 27 from Nakagawa et al. (2015), but without the '-2 rho ...' term,
         ### since for normally distributed data the mean and variance (and transformations thereof) are independent

         if (measure == "CVLN") {
            yi <- log(sdi/mi) + 1/(2*(ni-1))
            vi <- 1/(2*(ni-1)) + sdi^2/(ni*mi^2)
         }

         ### log(SD) with bias correction

         if (measure == "SDLN") {
            yi <- log(sdi) + 1/(2*(ni-1))
            vi <- 1/(2*(ni-1))
         }

         ### single-group standardized mean

         if (measure == "SMN") {
            cmi <- .cmicalc(ni-1)
            yi <- cmi * mi / sdi
            vi <- 1/ni + yi^2/(2*ni)
         }

      }

      ######################################################################

      if (is.element(measure, c("MC","SMCC","SMCR","SMCRH","ROMC","CVRC","VRC"))) {

         m1i  <- .getx("m1i",  mf=mf, data=data) ### for VRC, do not need to supply this
         m2i  <- .getx("m2i",  mf=mf, data=data) ### for VRC, do not need to supply this
         sd1i <- .getx("sd1i", mf=mf, data=data)
         sd2i <- .getx("sd2i", mf=mf, data=data) ### for SMCR, do not need to supply this
         ri   <- .getx("ri",   mf=mf, data=data)
         ni   <- .getx("ni",   mf=mf, data=data)

         if (is.element(measure, c("MC","SMCC","SMCRH","ROMC","CVRC"))) {

            ### for these measures, need m1i, m2i, sd1i, sd2i, ni, and ri

            if (!.all.specified(m1i, m2i, sd1i, sd2i, ni, ri))
               stop(mstyle$stop("Cannot compute outcomes. Check that all of the required\n  information is specified via the appropriate arguments."))

            if (!.equal.length(m1i, m2i, sd1i, sd2i, ni, ri))
               stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         }

         if (measure == "SMCR") {

            ### for this measure, need m1i, m2i, sd1i, ni, and ri (do not need sd2i)

            if (!.all.specified(m1i, m2i, sd1i, ni, ri))
               stop(mstyle$stop("Cannot compute outcomes. Check that all of the required\n  information is specified via the appropriate arguments."))

            if (!.equal.length(m1i, m2i, sd1i, ni, ri))
               stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         }

         if (measure == "VRC") {

            ### for this measure, need sd1i, sd2i, ni, and ri

            if (!.all.specified(sd1i, sd2i, ni, ri))
               stop(mstyle$stop("Cannot compute outcomes. Check that all of the required\n  information is specified via the appropriate arguments."))

            if (!.equal.length(sd1i, sd2i, ni, ri))
               stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         }

         k.all <- length(ni)

         if (!is.null(subset)) {
            subset <- .setnafalse(subset, k=k.all)
            m1i  <- m1i[subset]
            m2i  <- m2i[subset]
            sd1i <- sd1i[subset]
            sd2i <- sd2i[subset]
            ni   <- ni[subset]
            ri   <- ri[subset]
         }

         if (is.element(measure, c("MC","SMCC","SMCRH","ROMC","CVRC","VRC"))) {
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

         ni.u <- ni ### unadjusted total sample sizes

         k <- length(ni)

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
            di <- (m1i - m2i) / sddi
            yi <- cmi * di

            if (length(vtype) == 1L)
               vtype <- rep(vtype, k)

            vi <- rep(NA_real_, k)

            if (!all(is.element(vtype, c("LS","LS2"))))
               stop(mstyle$stop("For this outcome measure, 'vtype' must be either ''LS' or 'LS2'."))

            for (i in seq_len(k)) {

               ### large sample approximation to the sampling variance
               if (vtype[i] == "LS")
                  vi[i] <- 1/ni[i] + yi[i]^2 / (2*ni[i])

               ### large sample approximation to the sampling variance (analogous to LS2 for SMD and SMCR)
               if (vtype[i] == "LS2")
                  vi[i] <- cmi[i]^2 * (1/ni[i] + di[i]^2 / (2*ni[i]))

            }

         }

         ### standardized mean change with raw score standardization (using sd1i)
         ### note: yi does not assume homoscedasticity, but vi does

         if (measure == "SMCR") {

            cmi <- .cmicalc(mi)
            di <- (m1i - m2i) / sd1i
            yi <- cmi * di

            if (length(vtype) == 1L)
               vtype <- rep(vtype, k)

            vi <- rep(NA_real_, k)

            if (!all(is.element(vtype, c("LS","LS2"))))
               stop(mstyle$stop("For this outcome measure, 'vtype' must be either ''LS' or 'LS2'."))

            for (i in seq_len(k)) {

               ### large sample approximation to the sampling variance
               if (vtype[i] == "LS")
                  vi[i] <- 2*(1-ri[i])/ni[i] + yi[i]^2 / (2*ni[i])

               ### large sample approximation to the sampling variance (using corrected (!) equation from Borenstein et al., 2009)
               if (vtype[i] == "LS2")
                  vi[i] <- cmi[i]^2 * (2*(1-ri[i])/ni[i] + di[i]^2 / (2*ni[i]))
                  #vi[i] <- cmi[i]^2 * 2*(1-ri[i]) * (1/ni[i] + di[i]^2 / (2*ni[i])) # as in Borenstein et al. (2009) but this is incorrect

            }

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
         ### to use with pooled SDs, simply set sd1i = sd2i = sdpi

         if (measure == "ROMC") {
            yi <- log(m1i/m2i)
            vi <- sd1i^2/(ni*m1i^2) + sd2i^2/(ni*m2i^2) - 2*ri*sd1i*sd2i/(m1i*m2i*ni)
         }

         ### coefficient of variation ratio for pre-post or matched designs

         if (measure == "CVRC") {
            yi <- log(sd1i/m1i) - log(sd2i/m2i)
            vi <- (1-ri^2)/(ni-1) + (m1i^2*sd2i^2 + m2i^2*sd1i^2 - 2*m1i*m2i*ri*sd1i*sd2i) / (m1i^2*m2i^2*ni)
         }

         ### variability ratio for pre-post or matched designs

         if (measure == "VRC") {
            yi <- log(sd1i/sd2i)
            vi <- (1-ri^2)/(ni-1)
         }

      }

      ######################################################################

      if (is.element(measure, c("ARAW","AHW","ABT"))) {

         ai <- .getx("ai", mf=mf, data=data)
         mi <- .getx("mi", mf=mf, data=data)
         ni <- .getx("ni", mf=mf, data=data)

         if (!.all.specified(ai, mi, ni))
            stop(mstyle$stop("Cannot compute outcomes. Check that all of the required\n  information is specified via the appropriate arguments."))

         if (!.equal.length(ai, mi, ni))
            stop(mstyle$stop("Supplied data vectors are not all of the same length."))

         k.all <- length(ai)

         if (!is.null(subset)) {
            subset <- .setnafalse(subset, k=k.all)
            ai <- ai[subset]
            mi <- mi[subset]
            ni <- ni[subset]
         }

         if (any(ai > 1, na.rm=TRUE))
            stop(mstyle$stop("One or more alpha values are > 1."))

         if (any(mi < 2, na.rm=TRUE))
            stop(mstyle$stop("One or more mi values are < 2."))

         if (any(ni <= 0, na.rm=TRUE))
            stop(mstyle$stop("One or more sample sizes are <= 0."))

         ni.u <- ni ### unadjusted total sample sizes

         k <- length(ai)

         ### raw alpha values

         if (measure == "ARAW") {
            yi <- ai
            vi <- 2*mi*(1-ai)^2 / ((mi-1)*(ni-2))
         }

         ### alphas transformed with Hakstian & Whalen (1976) transformation

         if (measure == "AHW") {
            #yi <- (1-ai)^(1/3)    ### technically this is the Hakstian & Whalen (1976) transformation
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

      vi  <- .getx("vi",  mf=mf, data=data)
      sei <- .getx("sei", mf=mf, data=data)
      ni  <- .getx("ni",  mf=mf, data=data)

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

      if (!is.null(subset)) {
         subset <- .setnafalse(subset, k=k.all)
         yi <- yi[subset]
         vi <- vi[subset]
         ni <- ni[subset]
      }

      ni.u <- ni ### unadjusted total sample sizes

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

      if (length(slab) != k.all)
         stop(mstyle$stop("Study labels not of same length as data."))

      if (is.factor(slab))
         slab <- as.character(slab)

      if (!is.null(subset))
         slab <- slab[subset]

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

   ### turn numeric include vector into logical vector (already done for subset)

   if (!is.null(include))
      include <- .setnafalse(include, arg="include", k=k.all)

   ### apply subset to include

   include <- include[subset]

   ### subset data frame (note: subsetting of other parts already done above, so yi/vi/ni.u/slab are already subsetted)

   if (has.data && any(!subset))
      data <- data[subset,,drop=FALSE]

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
   attr(dat[[var.names[1]]], "slab") <- slab

   ### add measure attribute to the yi vector
   attr(dat[[var.names[1]]], "measure") <- measure

   ### add digits attribute
   attr(dat, "digits") <- digits

   ### add 'yi.names' and 'vi.names' to the first position of the corresponding attributes
   attr(dat, "yi.names") <- unique(c(var.names[1], attr(data, "yi.names"))) ### if 'yi.names' is not an attribute, attr() returns NULL, so this works fine
   attr(dat, "vi.names") <- unique(c(var.names[2], attr(data, "vi.names"))) ### if 'vi.names' is not an attribute, attr() returns NULL, so this works fine

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
