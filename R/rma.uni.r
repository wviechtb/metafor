rma <- rma.uni <- function(yi, vi, sei, weights, ai, bi, ci, di, n1i, n2i, x1i, x2i, t1i, t2i, m1i, m2i, sd1i, sd2i, xi, mi, ri, ti, sdi, r2i, ni, mods, #scale,
measure="GEN", intercept=TRUE,
data, slab, subset,
add=1/2, to="only0", drop00=FALSE, vtype="LS",
method="REML", weighted=TRUE, test="z", #knha=FALSE,
level=95, digits, btt, tau2, verbose=FALSE, control, ...) {

   #########################################################################

   ###### setup

   mstyle <- .get.mstyle("crayon" %in% .packages())

   ### check argument specifications
   ### (arguments "to" and "vtype" are checked inside escalc function)

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
                              "MN","MNLN","CVLN","SDLN","SMD1",                    ### mean, log(mean), log(CV), log(SD), single-group SMD
                              "MC","SMCC","SMCR","SMCRH","ROMC","CVRC","VRC",      ### raw/standardized mean change, log(ROM), CVR, and VR for dependent samples
                              "ARAW","AHW","ABT",                                  ### alpha (and transformations thereof)
                              "GEN")))
      stop(mstyle$stop("Unknown 'measure' specified."))

   if (!is.element(method, c("FE","HS","HSk","HE","DL","DLIT","GENQ","GENQM","SJ","SJIT","PM","PMM","ML","REML","EB")))
      stop(mstyle$stop("Unknown 'method' specified."))

   ### in case user specifies more than one add/to value (as one can do with rma.mh() and rma.peto())
   ### (any kind of continuity correction is directly applied to the outcomes, which are then analyzed as such)

   if (length(add) > 1L)
      add <- add[1]

   if (length(to) > 1L)
      to <- to[1]

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   if (missing(tau2))
      tau2 <- NULL

   if (missing(control))
      control <- list()

   time.start <- proc.time()

   ### get ... argument and check for extra/superfluous arguments

   ddd <- list(...)

   .chkdots(ddd, c("knha", "scale", "link", "alpha", "outlist", "onlyo1", "addyi", "addvi", "time", "skipr2", "skiphes", "att"))

   ### handle 'knha' argument from ... (note: overrides test argument)

   if (.isFALSE(ddd$knha))
      test <- "z"
   if (.isTRUE(ddd$knha))
      test <- "knha"

   if (!is.element(test, c("z", "t", "knha", "adhoc")))
      stop(mstyle$stop("Invalid option selected for 'test' argument."))

   if (!is.null(ddd$scale)) {

      scale <- ddd$scale

      #if (!inherits(scale, "formula"))
      #   stop(mstyle$stop("Must specify a formula for the 'scale' argument."))

      if (is.element(test, c("knha", "adhoc")))
         stop(mstyle$stop("Cannot use Knapp & Hartung method with location-scale models."))

      model <- "rma.ls"

   } else {

      model <- "rma.uni"

   }

   if (!is.null(ddd$link)) {
      link <- match.arg(ddd$link, c("log", "identity"))
   } else {
      link <- "log"
   }

   if (!is.null(ddd$alpha)) {
      alpha <- ddd$alpha
   } else {
      alpha <- NA
   }

   if (!is.null(ddd$att)) {
      att <- ddd$att
   } else {
      att <- NULL
   }

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

   ### set defaults for formulas

   formula.yi <- NULL
   formula.mods <- NULL
   formula.scale <- NULL

   ### set options(warn=1) if verbose > 2

   if (verbose > 2) {
      opwarn <- options(warn=1)
      on.exit(options(warn=opwarn$warn))
   }

   #########################################################################

   if (verbose && !exists(".rmspace"))
      cat("\n")

   if (verbose > 1)
      message(mstyle$message("Extracting/computing yi/vi values ..."))

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

   ### for certain measures, set add=0 by default unless user explicitly sets the add argument

   addval <- mf[[match("add", names(mf))]]

   if (is.element(measure, c("AS","PHI","RTET","IRSD","PAS","PFT","IRS","IRFT")) && is.null(addval))
      add <- 0

   ### extract yi (either NULL if not specified, a vector, a formula, or an escalc object)

   mf.yi <- mf[[match("yi", names(mf))]]
   yi <- eval(mf.yi, data, enclos=sys.frame(sys.parent())) ### NULL if user does not specify this

   ### if yi is not NULL and it is an escalc object, then use that object in place of the data argument

   if (!is.null(yi) && inherits(yi, "escalc"))
      data <- yi

   ### extract weights, slab, subset, mods, and scale values, possibly from the data frame specified via data or yi (arguments not specified are NULL)

   mf.weights <- mf[[match("weights", names(mf))]]
   mf.slab    <- mf[[match("slab",    names(mf))]]
   mf.subset  <- mf[[match("subset",  names(mf))]]
   mf.mods    <- mf[[match("mods",    names(mf))]]
   mf.scale   <- mf[[match("scale",   names(mf))]]
   weights <- eval(mf.weights, data, enclos=sys.frame(sys.parent())) ### NULL if user does not specify this
   slab    <- eval(mf.slab,    data, enclos=sys.frame(sys.parent())) ### NULL if user does not specify this
   subset  <- eval(mf.subset,  data, enclos=sys.frame(sys.parent())) ### NULL if user does not specify this
   mods    <- eval(mf.mods,    data, enclos=sys.frame(sys.parent())) ### NULL if user does not specify this
   scale   <- eval(mf.scale,   data, enclos=sys.frame(sys.parent())) ### NULL if user does not specify this

   ai <- bi <- ci <- di <- x1i <- x2i <- t1i <- t2i <- NA

   if (!is.null(subset))
      subset <- .setnafalse(subset)

   if (!is.null(yi)) {

      ### if yi is not NULL, then yi now either contains the yi values, a formula, or an escalc object

      ### if yi is a formula, extract yi and X (this overrides anything specified via the mods argument further below)

      if (inherits(yi, "formula")) {
         formula.yi <- yi
         options(na.action = "na.pass")                   ### set na.action to na.pass, so that NAs are not filtered out (we'll do that later)
         mods <- model.matrix(yi, data=data)              ### extract model matrix (now mods is no longer a formula, so [a] further below is skipped)
         attr(mods, "assign") <- NULL                     ### strip assign attribute (not needed at the moment)
         attr(mods, "contrasts") <- NULL                  ### strip contrasts attribute (not needed at the moment)
         yi <- model.response(model.frame(yi, data=data)) ### extract yi values from model frame
         options(na.action = na.act)                      ### set na.action back to na.act
         names(yi) <- NULL                                ### strip names (1:k) from yi (so res$yi is the same whether yi is a formula or not)
         intercept <- FALSE                               ### set to FALSE since formula now controls whether the intercept is included or not
      }                                                   ### note: code further below ([b]) actually checks whether intercept is included or not

      ### if yi is an escalc object, try to extract yi and vi (note that moderators must then be specified via the mods argument)

      if (inherits(yi, "escalc")) {

         if (!is.null(attr(yi, "yi.names"))) { ### if yi.names attributes is available
            yi.name <- attr(yi, "yi.names")[1] ### take the first entry to be the yi variable
         } else {                              ### if not, see if 'yi' is in the object and assume that is the yi variable
            if (!is.element("yi", names(yi)))
               stop(mstyle$stop("Cannot determine name of the 'yi' variable."))
            yi.name <- "yi"
         }
         if (!is.null(attr(yi, "vi.names"))) { ### if vi.names attributes is available
            vi.name <- attr(yi, "vi.names")[1] ### take the first entry to be the vi variable
         } else {                              ### if not, see if 'vi' is in the object and assume that is the vi variable
            if (!is.element("vi", names(yi)))
               stop(mstyle$stop("Cannot determine name of the 'vi' variable."))
            vi.name <- "vi"
         }

         ### get vi and yi variables from the escalc object (vi first, then yi)

         vi <- yi[[vi.name]]
         yi <- yi[[yi.name]]

         yi.escalc <- TRUE

      } else {

         yi.escalc <- FALSE

      }

      ### in case user passed a matrix to yi, convert it to a vector

      if (is.matrix(yi))
         yi <- as.vector(yi)

      ### number of outcomes before subsetting

      k <- length(yi)
      k.all <- k

      ### if the user has specified 'measure' to be something other than "GEN", then use that for the measure argument
      ### otherwise, if yi has a 'measure' attribute, use that to set the 'measure' argument

      if (measure == "GEN" && !is.null(attr(yi, "measure")))
         measure <- attr(yi, "measure")

      ### add measure attribute (back) to the yi vector

      attr(yi, "measure") <- measure

      ### extract vi and sei values (but only if yi wasn't an escalc object)

      if (!yi.escalc) {

         mf.vi  <- mf[[match("vi",  names(mf))]]
         mf.sei <- mf[[match("sei", names(mf))]]
         vi  <- eval(mf.vi,  data, enclos=sys.frame(sys.parent())) ### NULL if user does not specify this
         sei <- eval(mf.sei, data, enclos=sys.frame(sys.parent())) ### NULL if user does not specify this

      }

      ### extract ni argument

      mf.ni <- mf[[match("ni", names(mf))]]
      ni <- eval(mf.ni,  data, enclos=sys.frame(sys.parent())) ### NULL if user does not specify this

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

      ### in case user passes a matrix to vi, convert it to a vector
      ### note: only a row or column matrix with the right dimensions will end with the right length

      if (is.matrix(vi))
         vi <- as.vector(vi)

      ### check if user constrained vi to 0

      if (length(vi) == 1L && vi == 0) {
         vi0 <- TRUE
      } else {
         vi0 <- FALSE
      }

      ### allow easy setting of vi to a single value

      if (length(vi) == 1L)
         vi <- rep(vi, k) ### note: k is number of outcomes before subsetting

      ### check length of yi and vi

      if (length(vi) != k)
         stop(mstyle$stop("Length of 'yi' and 'vi' (or 'sei') is not the same."))

      ### if ni has not been specified but is an attribute of yi, get it

      if (is.null(ni) && !is.null(attr(yi, "ni")))
         ni <- attr(yi, "ni")

      ### check length of yi and ni (only if ni is not NULL)
      ### if there is a mismatch, then ni cannot be trusted, so set it to NULL

      if (!is.null(ni) && length(ni) != k)
         ni <- NULL

      ### if ni is now available, add it (back) as an attribute to yi

      if (!is.null(ni))
         attr(yi, "ni") <- ni

      ### note: one or more yi/vi pairs may be NA/NA (also a corresponding ni value may be NA)

      ### if slab has not been specified but is an attribute of yi, get it

      if (is.null(slab)) {

         if (!is.null(attr(yi, "slab")))
            slab <- attr(yi, "slab")

         ### check length of yi and slab (only if slab is now not NULL)
         ### if there is a mismatch, then slab cannot be trusted, so set it to NULL

         if (!is.null(slab) && length(slab) != k)
            slab <- NULL

      }

      ### subsetting of yi/vi/ni values (note: mods and slab are subsetted further below)

      if (!is.null(subset)) {

         yi <- yi[subset]
         vi <- vi[subset]
         ni <- ni[subset]

         attr(yi, "measure") <- measure ### add measure attribute back
         attr(yi, "ni")      <- ni      ### add ni attribute back

      }

   } else {

      ### if yi is NULL, try to compute yi/vi based on specified measure and supplied data

      if (is.element(measure, c("RR","OR","PETO","RD","AS","PHI","YUQ","YUY","RTET","PBIT","OR2D","OR2DN","OR2DL","MPRD","MPRR","MPOR","MPORC","MPPETO"))) {

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
         k.all <- k

         if (!is.null(subset)) {
            ai <- ai[subset]
            bi <- bi[subset]
            ci <- ci[subset]
            di <- di[subset]
         }

         dat <- escalc(measure=measure, ai=ai, bi=bi, ci=ci, di=di, add=add, to=to, drop00=drop00, vtype=vtype, onlyo1=onlyo1, addyi=addyi, addvi=addvi)

      }

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
         k.all <- k

         if (!is.null(subset)) {
            x1i <- x1i[subset]
            x2i <- x2i[subset]
            t1i <- t1i[subset]
            t2i <- t2i[subset]
         }

         dat <- escalc(measure=measure, x1i=x1i, x2i=x2i, t1i=t1i, t2i=t2i, add=add, to=to, drop00=drop00, vtype=vtype, addyi=addyi, addvi=addvi)

      }

      if (is.element(measure, c("MD","SMD","SMDH","ROM","RPB","RBIS","D2OR","D2ORN","D2ORL","CVR","VR"))) {

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
         k.all <- k

         if (!is.null(subset)) {
            m1i  <- m1i[subset]
            m2i  <- m2i[subset]
            sd1i <- sd1i[subset]
            sd2i <- sd2i[subset]
            n1i  <- n1i[subset]
            n2i  <- n2i[subset]
         }

         dat <- escalc(measure=measure, m1i=m1i, m2i=m2i, sd1i=sd1i, sd2i=sd2i, n1i=n1i, n2i=n2i, vtype=vtype)

      }

      if (is.element(measure, c("COR","UCOR","ZCOR"))) {

         mf.ri <- mf[[match("ri", names(mf))]]
         mf.ni <- mf[[match("ni", names(mf))]]
         ri <- eval(mf.ri, data, enclos=sys.frame(sys.parent()))
         ni <- eval(mf.ni, data, enclos=sys.frame(sys.parent()))

         k <- length(ri) ### number of outcomes before subsetting
         k.all <- k

         if (!is.null(subset)) {
            ri <- ri[subset]
            ni <- ni[subset]
         }

         dat <- escalc(measure=measure, ri=ri, ni=ni, vtype=vtype)

      }

      if (is.element(measure, c("PCOR","ZPCOR","SPCOR"))) {

         mf.ti  <- mf[[match("ti",  names(mf))]]
         mf.r2i <- mf[[match("r2i", names(mf))]]
         mf.mi  <- mf[[match("mi",  names(mf))]]
         mf.ni  <- mf[[match("ni",  names(mf))]]
         ti  <- eval(mf.ti,  data, enclos=sys.frame(sys.parent()))
         r2i <- eval(mf.r2i, data, enclos=sys.frame(sys.parent()))
         mi  <- eval(mf.mi,  data, enclos=sys.frame(sys.parent()))
         ni  <- eval(mf.ni,  data, enclos=sys.frame(sys.parent()))

         k <- length(ti) ### number of outcomes before subsetting
         k.all <- k

         if (!is.null(subset)) {
            ti  <- ti[subset]
            r2i <- r2i[subset]
            mi  <- mi[subset]
            ni  <- ni[subset]
         }

         dat <- escalc(measure=measure, ti=ti, r2i=r2i, mi=mi, ni=ni, vtype=vtype)

      }

      if (is.element(measure, c("PR","PLN","PLO","PAS","PFT"))) {

         mf.xi <- mf[[match("xi", names(mf))]]
         mf.mi <- mf[[match("mi", names(mf))]]
         mf.ni <- mf[[match("ni", names(mf))]]
         xi <- eval(mf.xi, data, enclos=sys.frame(sys.parent()))
         mi <- eval(mf.mi, data, enclos=sys.frame(sys.parent()))
         ni <- eval(mf.ni, data, enclos=sys.frame(sys.parent()))
         if (is.null(mi)) mi <- ni - xi

         k <- length(xi) ### number of outcomes before subsetting
         k.all <- k

         if (!is.null(subset)) {
            xi <- xi[subset]
            mi <- mi[subset]
         }

         dat <- escalc(measure=measure, xi=xi, mi=mi, add=add, to=to, vtype=vtype, addyi=addyi, addvi=addvi)

      }

      if (is.element(measure, c("IR","IRLN","IRS","IRFT"))) {

         mf.xi <- mf[[match("xi", names(mf))]]
         mf.ti <- mf[[match("ti", names(mf))]]
         xi <- eval(mf.xi, data, enclos=sys.frame(sys.parent()))
         ti <- eval(mf.ti, data, enclos=sys.frame(sys.parent()))

         k <- length(xi) ### number of outcomes before subsetting
         k.all <- k

         if (!is.null(subset)) {
            xi <- xi[subset]
            ti <- ti[subset]
         }

         dat <- escalc(measure=measure, xi=xi, ti=ti, add=add, to=to, vtype=vtype, addyi=addyi, addvi=addvi)

      }

      if (is.element(measure, c("MN","MNLN","CVLN","SDLN","SMD1"))) {

         mf.mi  <- mf[[match("mi",  names(mf))]]
         mf.sdi <- mf[[match("sdi", names(mf))]]
         mf.ni  <- mf[[match("ni",  names(mf))]]
         mi  <- eval(mf.mi,  data, enclos=sys.frame(sys.parent()))
         sdi <- eval(mf.sdi, data, enclos=sys.frame(sys.parent()))
         ni  <- eval(mf.ni,  data, enclos=sys.frame(sys.parent()))

         k <- length(ni) ### number of outcomes before subsetting
         k.all <- k

         if (!is.null(subset)) {
            mi  <- mi[subset]
            sdi <- sdi[subset]
            ni  <- ni[subset]
         }

         dat <- escalc(measure=measure, mi=mi, sdi=sdi, ni=ni, vtype=vtype)

      }

      if (is.element(measure, c("MC","SMCC","SMCR","SMCRH","ROMC","CVRC","VRC"))) {

         mf.m1i  <- mf[[match("m1i",  names(mf))]]
         mf.m2i  <- mf[[match("m2i",  names(mf))]]
         mf.sd1i <- mf[[match("sd1i", names(mf))]]
         mf.sd2i <- mf[[match("sd2i", names(mf))]]
         mf.ri   <- mf[[match("ri",   names(mf))]]
         mf.ni   <- mf[[match("ni",   names(mf))]]
         m1i  <- eval(mf.m1i,  data, enclos=sys.frame(sys.parent()))
         m2i  <- eval(mf.m2i,  data, enclos=sys.frame(sys.parent()))
         sd1i <- eval(mf.sd1i, data, enclos=sys.frame(sys.parent()))
         sd2i <- eval(mf.sd2i, data, enclos=sys.frame(sys.parent()))
         ri   <- eval(mf.ri,   data, enclos=sys.frame(sys.parent()))
         ni   <- eval(mf.ni,   data, enclos=sys.frame(sys.parent()))

         k <- length(m1i) ### number of outcomes before subsetting
         k.all <- k

         if (!is.null(subset)) {
            m1i  <- m1i[subset]
            m2i  <- m2i[subset]
            sd1i <- sd1i[subset]
            sd2i <- sd2i[subset]
            ni   <- ni[subset]
            ri   <- ri[subset]
         }

         dat <- escalc(measure=measure, m1i=m1i, m2i=m2i, sd1i=sd1i, sd2i=sd2i, ri=ri, ni=ni, vtype=vtype)

      }

      if (is.element(measure, c("ARAW","AHW","ABT"))) {

         mf.ai <- mf[[match("ai",  names(mf))]]
         mf.mi <- mf[[match("mi",  names(mf))]]
         mf.ni <- mf[[match("ni",  names(mf))]]
         ai <- eval(mf.ai,  data, enclos=sys.frame(sys.parent()))
         mi <- eval(mf.mi,  data, enclos=sys.frame(sys.parent()))
         ni <- eval(mf.ni,  data, enclos=sys.frame(sys.parent()))

         k <- length(ai) ### number of outcomes before subsetting
         k.all <- k

         if (!is.null(subset)) {
            ai <- ai[subset]
            mi <- mi[subset]
            ni <- ni[subset]
         }

         dat <- escalc(measure=measure, ai=ai, mi=mi, ni=ni, vtype=vtype)

      }

      if (is.element(measure, "GEN"))
         stop(mstyle$stop("Specify the desired outcome measure via the 'measure' argument."))

      ### note: these values are already subsetted

      yi <- dat$yi         ### one or more yi/vi pairs may be NA/NA
      vi <- dat$vi         ### one or more yi/vi pairs may be NA/NA
      ni <- attr(yi, "ni") ### unadjusted total sample sizes (ni.u in escalc)

   }

   #########################################################################

   ### allow easy setting of weights to a single value

   if (length(weights) == 1L)
      weights <- rep(weights, k) ### note: k is number of outcomes before subsetting

   ### check length of yi and weights (only if weights is not NULL)

   if (!is.null(weights) && (length(weights) != k))
      stop(mstyle$stop("Length of 'yi' and 'weights' is not the same."))

   ### subsetting of weights

   if (!is.null(subset))
      weights <- weights[subset]

   #########################################################################

   if (verbose > 1)
      message(mstyle$message("Creating model matrix ..."))

   ### convert mods formula to X matrix and set intercept equal to FALSE
   ### skipped if formula has already been specified via yi argument, since mods is then no longer a formula (see [a])

   if (inherits(mods, "formula")) {
      formula.mods <- mods
      options(na.action = "na.pass")        ### set na.action to na.pass, so that NAs are not filtered out (we'll do that later)
      mods <- model.matrix(mods, data=data) ### extract model matrix
      attr(mods, "assign") <- NULL          ### strip assign attribute (not needed at the moment)
      attr(mods, "contrasts") <- NULL       ### strip contrasts attribute (not needed at the moment)
      options(na.action = na.act)           ### set na.action back to na.act
      intercept <- FALSE                    ### set to FALSE since formula now controls whether the intercept is included or not
   }                                        ### note: code further below ([b]) actually checks whether intercept is included or not

   ### turn a vector for mods into a column vector

   if (.is.vector(mods))
      mods <- cbind(mods)

   ### turn a mods data frame into a matrix

   if (is.data.frame(mods))
      mods <- as.matrix(mods)

   ### check if model matrix contains character variables

   if (is.character(mods))
      stop(mstyle$stop("Model matrix contains character variables."))

   ### check if mods matrix has the right number of rows

   if (!is.null(mods) && nrow(mods) != k)
      stop(mstyle$stop(paste0("Number of rows in the model matrix (", nrow(mods), ") does not match length of the outcome vector (", k, ").")))

   ### for rma.ls models, get model matrix for scale part

   if (model == "rma.ls") {
      if (inherits(scale, "formula")) {
         formula.scale <- scale
         options(na.action = "na.pass")
         Z <- model.matrix(scale, data=data)
         colnames(Z)[grep("(Intercept)", colnames(Z))] <- "intrcpt"
         attr(Z, "assign") <- NULL
         attr(Z, "contrasts") <- NULL
         options(na.action = na.act)
      } else {
         Z <- scale
         if (.is.vector(Z))
            Z <- cbind(Z)
         if (is.data.frame(Z))
            Z <- as.matrix(Z)
         if (is.character(Z))
            stop(mstyle$stop("Model matrix contains character variables."))
      }
      if (nrow(Z) != k)
         stop(mstyle$stop(paste0("Number of rows in the model matrix specified via the 'scale' argument (", nrow(Z), ") does not match length of the outcome vector (", k, ").")))
   } else {
      Z <- NULL
   }

   ### generate study labels if none are specified (or none have been found in yi)

   if (verbose > 1)
      message(mstyle$message("Generating/extracting study labels ..."))

   ### study ids (1:k sequence before subsetting)

   ids <- seq_len(k)

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

      if (verbose > 1)
         message(mstyle$message("Subsetting ..."))

      mods <- mods[subset,,drop=FALSE]
      slab <- slab[subset]
      ids  <- ids[subset]
      Z    <- Z[subset,,drop=FALSE]

   }

   ### check if study labels are unique; if not, make them unique

   if (anyDuplicated(slab))
      slab <- .make.unique(slab)

   ### add slab attribute back

   attr(yi, "slab") <- slab

   ### number of outcomes after subsetting

   k <- length(yi)

   ### check for non-positive sampling variances (and set negative values to 0)

   if (any(vi <= 0, na.rm=TRUE)) {
      allvipos <- FALSE
      if (!vi0)
         warning(mstyle$warning("There are outcomes with non-positive sampling variances."), call.=FALSE)
      vi.neg <- vi < 0
      if (any(vi.neg, na.rm=TRUE)) {
         vi[vi.neg] <- 0
         warning(mstyle$warning("Negative sampling variances constrained to zero."), call.=FALSE)
      }
   } else {
      allvipos <- TRUE
   }

   ### check for (and correct?) negative/infinite weights

   if (any(weights < 0, na.rm=TRUE))
      stop(mstyle$stop("Negative weights not allowed."))

   if (any(is.infinite(weights)))
      stop(mstyle$stop("Infinite weights not allowed."))

   ### save full data (including potential NAs in yi/vi/weights/ni/mods/Z.f)

   ai.f      <- ai
   bi.f      <- bi
   ci.f      <- ci
   di.f      <- di
   x1i.f     <- x1i
   x2i.f     <- x2i
   t1i.f     <- t1i
   t2i.f     <- t2i
   yi.f      <- yi
   vi.f      <- vi
   weights.f <- weights
   ni.f      <- ni
   mods.f    <- mods
   Z.f       <- Z

   k.f <- k ### total number of observed outcomes including all NAs

   ### check for NAs and act accordingly

   has.na <- is.na(yi) | is.na(vi) | (if (is.null(mods)) FALSE else apply(is.na(mods), 1, any)) | (if (is.null(Z)) FALSE else apply(is.na(Z), 1, any)) | (if (is.null(weights)) FALSE else is.na(weights))
   not.na <- !has.na

   if (any(has.na)) {

      if (verbose > 1)
         message(mstyle$message("Handling NAs ..."))

      if (na.act == "na.omit" || na.act == "na.exclude" || na.act == "na.pass") {

         yi      <- yi[not.na]
         vi      <- vi[not.na]
         weights <- weights[not.na]
         ni      <- ni[not.na]
         mods    <- mods[not.na,,drop=FALSE]
         Z       <- Z[not.na,,drop=FALSE]
         k       <- length(yi)
         warning(mstyle$warning("Studies with NAs omitted from model fitting."), call.=FALSE)

         attr(yi, "measure") <- measure ### add measure attribute back
         attr(yi, "ni")      <- ni      ### add ni attribute back

         ### note: slab is always of the same length as the full yi vector (after subsetting), so missings are not removed and slab is not added back to yi

      }

      if (na.act == "na.fail")
         stop(mstyle$stop("Missing values in data."))

   }

   ### at least one study left?

   if (k < 1L)
      stop(mstyle$stop("Processing terminated since k = 0."))

   ### if k=1 and test != "z", set test="z" (other methods cannot be used)

   if (k == 1L && test != "z") {
      warning(mstyle$warning("Setting argument test=\"z\" since k=1."), call.=FALSE)
      test <- "z"
   }

   ### make sure that there is at least one column in X ([b])

   if (is.null(mods) && !intercept) {
      warning(mstyle$warning("Must either include an intercept and/or moderators in model.\n  Coerced intercept into the model."), call.=FALSE)
      intercept <- TRUE
   }

   ### add vector of 1s to the X matrix for the intercept (if intercept=TRUE)

   if (intercept) {
      X   <- cbind(intrcpt=rep(1,k), mods)
      X.f <- cbind(intrcpt=rep(1,k.f), mods.f)
   } else {
      X   <- mods
      X.f <- mods.f
   }

   ### drop redundant predictors
   ### note: need to save coef.na for functions that modify the data/model and then refit the model (regtest() and the
   ### various function that leave out an observation); so we can check if there are redundant/dropped predictors then

   tmp <- try(lm(yi ~ X - 1), silent=TRUE)
   if (inherits(tmp, "lm")) {
      coef.na <- is.na(coef(tmp))
   } else {
      coef.na <- rep(FALSE, NCOL(X))
   }
   if (any(coef.na)) {
      warning(mstyle$warning("Redundant predictors dropped from the model."), call.=FALSE)
      X   <- X[,!coef.na,drop=FALSE]
      X.f <- X.f[,!coef.na,drop=FALSE]
   }

   ### check whether intercept is included and if yes, move it to the first column (NAs already removed, so na.rm=TRUE for any() not necessary)

   is.int <- apply(X, 2, .is.intercept)
   if (any(is.int)) {
      int.incl <- TRUE
      int.indx <- which(is.int, arr.ind=TRUE)
      X        <- cbind(intrcpt=1,   X[,-int.indx, drop=FALSE]) ### this removes any duplicate intercepts
      X.f      <- cbind(intrcpt=1, X.f[,-int.indx, drop=FALSE]) ### this removes any duplicate intercepts
      intercept <- TRUE ### set intercept appropriately so that the predict() function works
   } else {
      int.incl <- FALSE
   }

   p <- NCOL(X) ### number of columns in X (including the intercept if it is included)

   ### check whether this is an intercept-only model

   if ((p == 1L) && .is.intercept(X)) {
      int.only <- TRUE
   } else {
      int.only <- FALSE
   }

   ### check if there are too many parameters for given k

   if (!(int.only && k == 1L)) {
      if (method == "FE") {                        ### have to estimate p parms
         if (p > k)
            stop(mstyle$stop("Number of parameters to be estimated is larger than the number of observations."))
      } else {
         if (is.numeric(tau2)) {                   ### have to estimate p parms (tau2 is fixed at value specified)
            if (p > k)
               stop(mstyle$stop("Number of parameters to be estimated is larger than the number of observations."))
         } else {
            if ((p+1) > k)                         ### have to estimate p+1 parms
               stop(mstyle$stop("Number of parameters to be estimated is larger than the number of observations."))
         }
      }
   }

   ### set/check 'btt' argument

   btt <- .set.btt(btt, p, int.incl, colnames(X))
   m <- length(btt) ### number of betas to test (m = p if all betas are tested)

   #########################################################################

   ### set default control parameters

   con <- list(verbose = FALSE,
               tau2.init = NULL,          # initial value for iterative estimators (ML, REML, EB, SJ, SJIT, DLIT)
               tau2.min = 0,              # lower bound for tau^2 value
               tau2.max = 100,            # upper bound for tau^2 value (for PM/PMM/GENQM estimators; and passed down for tau^2 CI obtained with confint())
               threshold = 10^-5,         # convergence threshold (for ML, REML, EB, SJIT, DLIT)
               tol = .Machine$double.eps^0.25, # convergence tolerance for uniroot() as used for PM, PMM, and GENQM (also used in 'll0 - ll > con$tol' check for ML/REML)
               ll0check = TRUE,           # should the 'll0 - ll > con$tol' check be conducted for ML/REML?
               maxiter = 100,             # maximum number of iterations (for ML, REML, EB, SJIT, DLIT)
               stepadj = 1,               # step size adjustment for Fisher scoring algorithm (for ML, REML, EB)
               REMLf = TRUE,              # should |X'X| term be included in the REML log likelihood?
               evtol = 1e-07,             # lower bound for eigenvalues to determine if model matrix is positive definite (also for checking if vimaxmin >= 1/con$evtol)
               alpha.init = NULL,         # initial values for scale parameters
               optimizer = "nlminb",      # optimizer to use ("optim", "nlminb", "uobyqa", "newuoa", "bobyqa", "nloptr", "nlm", "hjk", "nmk", "mads", "ucminf", "optimParallel", "constrOptim") for location-scale models
               optmethod = "BFGS",        # argument 'method' for optim() ("Nelder-Mead" and "BFGS" are sensible options)
               parallel = list(),         # parallel argument for optimParallel() (note: 'cl' argument in parallel is not passed; this is directly specified via 'cl')
               cl = NULL,                 # arguments for optimParallel()
               ncpus = 1L,                # arguments for optimParallel()
               hessianCtrl=list(r=8),     # arguments passed on to 'method.args' of hessian()
               scaleZ = TRUE)

   ### replace defaults with any user-defined values

   con.pos <- pmatch(names(control), names(con))
   con[c(na.omit(con.pos))] <- control[!is.na(con.pos)]

   if (verbose)
      con$verbose <- verbose

   verbose <- con$verbose

   ### constrain negative tau2.min values to -min(vi) (to ensure that the marginal variance is always >= 0)

   if (con$tau2.min < 0 && (-con$tau2.min > min(vi))) {
      con$tau2.min <- -min(vi)
      warning(mstyle$warning(paste0("Value of 'tau2.min' constrained to -min(vi) = ", .fcf(-min(vi), digits[["est"]]), ".")), call.=FALSE)
   }

   ### convergence indicator and change variable (for iterative estimators)

   conv <- 1
   change <- con$threshold + 1

   ### check whether model matrix is of full rank

   if (any(eigen(crossprod(X), symmetric=TRUE, only.values=TRUE)$values <= con$evtol))
      stop(mstyle$stop("Model matrix not of full rank. Cannot fit model."))

   ### check ratio of largest to smallest sampling variance
   ### note: need to exclude some special cases (0/0 = NaN, max(vi)/0 = Inf)
   ### TODO: use the condition number of diag(vi) here instead?

   vimaxmin <- max(vi) / min(vi)

   if (!is.nan(vimaxmin) && !is.infinite(vimaxmin) && vimaxmin >= 1/con$evtol)
      warning(mstyle$warning("Ratio of largest to smallest sampling variance extremely large. May not be able to obtain stable results."), call.=FALSE)

   ### iterations counter for iterative estimators (i.e., DLIT, SJIT, ML, REML, EB)
   ### (note: PM, PMM, and GENQM are also iterative, but uniroot() handles that)

   iter <- 0

   ### set some defaults

   se.tau2 <- I2 <- H2 <- QE <- QEp <- NA

   s2w <- 1

   level <- ifelse(level == 0, 1, ifelse(level >= 1, (100-level)/100, ifelse(level > .5, 1-level, level)))

   Y <- as.matrix(yi)

   #########################################################################

   ###### heterogeneity estimation for standard model (rma.uni)

   if (model == "rma.uni") {

      if (is.numeric(tau2) && method != "FE") { ### if user has fixed the tau2 value
         tau2.fix <- TRUE
         tau2.val <- tau2
      } else {
         tau2.fix <- FALSE
         tau2.val <- NA
      }

      if (verbose > 1 && !tau2.fix && method != "FE")
         message(mstyle$message("Estimating tau^2 value ...\n"))

      if (k == 1L) {
         method.sav <- method
         method = "k1" # set method to k1 so all of the stuff below is skipped
         if (!tau2.fix)
            tau2 <- 0
      }

      ### Hunter & Schmidt (HS) estimator (or k-corrected HS estimator (HSk))

      if (is.element(method, c("HS", "HSk"))) {

         if (!allvipos)
            stop(mstyle$stop(method, " estimator cannot be used when there are non-positive sampling variances in the data."))

         wi    <- 1/vi
         W     <- diag(wi, nrow=k, ncol=k)
         stXWX <- .invcalc(X=X, W=W, k=k)
         P     <- W - W %*% X %*% stXWX %*% crossprod(X,W)
         RSS   <- crossprod(Y,P) %*% Y
         if (method == "HS") {
            tau2 <- ifelse(tau2.fix, tau2.val, (RSS - k) / sum(wi))
         } else {
            tau2 <- ifelse(tau2.fix, tau2.val, (k/(k-p)*RSS - k) / sum(wi))
            ### HSk = (RSS - (k-p)) / sum(wi) * k/(k-p)
            #trP  <- sum(wi) * (k-p) / k
            #tau2 <- ifelse(tau2.fix, tau2.val, k/(k-p) * (RSS - (k-p)) / trP)
            #tau2 <- ifelse(tau2.fix, tau2.val, (RSS - (k-p)) / trP)
            #tau2 <- ifelse(tau2.fix, tau2.val, k/(k-p) * (RSS - (k-p)) / sum(wi))
         }

      }

      ### Hedges (HE) estimator (or initial value for ML, REML, EB)

      if (is.element(method, c("HE","ML","REML","EB"))) {

         stXX  <- .invcalc(X=X, W=diag(k), k=k)
         P     <- diag(k) - X %*% tcrossprod(stXX,X)
         RSS   <- crossprod(Y,P) %*% Y
         V     <- diag(vi, nrow=k, ncol=k)
         PV    <- P %*% V ### note: this is not symmetric
         trPV  <- .tr(PV) ### since PV needs to be computed anyway, can use .tr()
         tau2  <- ifelse(tau2.fix, tau2.val, (RSS - trPV) / (k-p))

      }

      ### DerSimonian-Laird (DL) estimator

      if (method == "DL") {

         if (!allvipos)
            stop(mstyle$stop("DL estimator cannot be used when there are non-positive sampling variances in the data."))

         wi    <- 1/vi
         W     <- diag(wi, nrow=k, ncol=k)
         stXWX <- .invcalc(X=X, W=W, k=k)
         P     <- W - W %*% X %*% stXWX %*% crossprod(X,W)
         RSS   <- crossprod(Y,P) %*% Y
         trP   <- .tr(P)
         tau2  <- ifelse(tau2.fix, tau2.val, (RSS - (k-p)) / trP)

      }

      ### DerSimonian-Laird (DL) estimator with iteration

      if (method == "DLIT") {

         if (is.null(con$tau2.init)) {
            tau2 <- 0
         } else {
            tau2 <- con$tau2.init
         }

         while (change > con$threshold) {

            if (verbose)
               cat(mstyle$verbose(paste("Iteration", iter, "\ttau^2 =", .fcf(tau2, digits[["var"]]), "\n")))

            iter <- iter + 1
            old2 <- tau2
            wi   <- 1/(vi + tau2)
            if (any(tau2 + vi < 0))
               stop(mstyle$stop("Some marginal variances are negative."))
            if (any(is.infinite(wi)))
               stop(mstyle$stop("Division by zero when computing the inverse variance weights."))
            W     <- diag(wi, nrow=k, ncol=k)
            stXWX <- .invcalc(X=X, W=W, k=k)
            P     <- W - W %*% X %*% stXWX %*% crossprod(X,W)
            RSS   <- crossprod(Y,P) %*% Y
            trP   <- .tr(P)
            tau2  <- ifelse(tau2.fix, tau2.val, (RSS - (k-p)) / trP)
            tau2[tau2 < con$tau2.min] <- con$tau2.min
            change <- abs(old2 - tau2)

            if (iter > con$maxiter) {
               conv <- 0
               break
            }

         }

         if (conv == 0L)
            stop(mstyle$stop("Algorithm did not converge."))

      }

      ### generalized Q-statistic estimator

      if (method == "GENQ") {

         #if (!allvipos)
         #   stop(mstyle$stop("GENQ estimator cannot be used when there are non-positive sampling variances in the data."))

         if (is.null(weights))
            stop(mstyle$stop("Must specify 'weights' when method='GENQ'."))

         A     <- diag(weights, nrow=k, ncol=k)
         stXAX <- .invcalc(X=X, W=A, k=k)
         P     <- A - A %*% X %*% stXAX %*% t(X) %*% A
         V     <- diag(vi, nrow=k, ncol=k)
         PV    <- P %*% V ### note: this is not symmetric
         trP   <- .tr(P)
         trPV  <- .tr(PV)
         RSS   <- crossprod(Y,P) %*% Y
         tau2  <- ifelse(tau2.fix, tau2.val, (RSS - trPV) / trP)

      }

      ### generalized Q-statistic estimator (median unbiased version)

      if (method == "GENQM") {

         if (is.null(weights))
            stop(mstyle$stop("Must specify 'weights' when method='GENQM'."))

         A     <- diag(weights, nrow=k, ncol=k)
         stXAX <- .invcalc(X=X, W=A, k=k)
         P     <- A - A %*% X %*% stXAX %*% t(X) %*% A
         V     <- diag(vi, nrow=k, ncol=k)
         PV    <- P %*% V ### note: this is not symmetric
         trP   <- .tr(P)

         if (!tau2.fix) {

            RSS   <- crossprod(Y,P) %*% Y

            if (.GENQ.func(con$tau2.min, P=P, vi=vi, Q=RSS, level=0, k=k, p=p, getlower=TRUE) > 0.5) {

               ### if GENQ.tau2.min is > 0.5, then estimate < tau2.min

               tau2 <- con$tau2.min

            } else {

               if (.GENQ.func(con$tau2.max, P=P, vi=vi, Q=RSS, level=0, k=k, p=p, getlower=TRUE) < 0.5) {

                  ### if GENQ.tau2.max is < 0.5, then estimate > tau2.max

                  stop(mstyle$stop("Value of 'tau2.max' too low. Try increasing 'tau2.max' or switch to another 'method'."))

               } else {

                  tau2 <- try(uniroot(.GENQ.func, interval=c(con$tau2.min, con$tau2.max), tol=con$tol, maxiter=con$maxiter, P=P, vi=vi, Q=RSS, level=0.5, k=k, p=p, getlower=FALSE, verbose=verbose, digits=digits, extendInt="no")$root, silent=TRUE)

                  if (inherits(tau2, "try-error"))
                     stop(mstyle$stop("Error in iterative search for tau2 using uniroot()."))

               }

            }

         } else {

            tau2 <- tau2.val

         }

         wi <- 1/(vi + tau2)

      }

      ### Sidik-Jonkman (SJ) estimator

      if (method == "SJ") {

         if (is.null(con$tau2.init)) {
            tau2.0 <- c(var(yi) * (k-1)/k)
         } else {
            tau2.0 <- con$tau2.init
         }

         wi    <- 1/(vi + tau2.0)
         W     <- diag(wi, nrow=k, ncol=k)
         stXWX <- .invcalc(X=X, W=W, k=k)
         P     <- W - W %*% X %*% stXWX %*% crossprod(X,W)
         RSS   <- crossprod(Y,P) %*% Y
         V     <- diag(vi, nrow=k, ncol=k)
         PV    <- P %*% V ### note: this is not symmetric
         tau2  <- ifelse(tau2.fix, tau2.val, tau2.0 * RSS / (k-p))

      }

      ### Sidik-Jonkman (SJ) estimator with iteration

      if (method == "SJIT") {

         if (is.null(con$tau2.init)) {
            tau2 <- var(yi) * (k-1)/k
         } else {
            tau2 <- con$tau2.init
         }

         tau2.0 <- tau2

         while (change > con$threshold) {

            if (verbose)
               cat(mstyle$verbose(paste("Iteration", iter, "\ttau^2 =", .fcf(tau2, digits[["var"]]), "\n")))

            iter <- iter + 1
            old2 <- tau2

            wi     <- 1/(vi + tau2)
            W      <- diag(wi, nrow=k, ncol=k)
            stXWX  <- .invcalc(X=X, W=W, k=k)
            P      <- W - W %*% X %*% stXWX %*% crossprod(X,W)
            RSS    <- crossprod(Y,P) %*% Y
            V      <- diag(vi, nrow=k, ncol=k)
            PV     <- P %*% V ### note: this is not symmetric
            tau2   <- ifelse(tau2.fix, tau2.val, tau2 * RSS / (k-p))
            change <- abs(old2 - tau2)

            if (iter > con$maxiter) {
               conv <- 0
               break
            }

         }

         if (conv == 0L)
            stop(mstyle$stop("Algorithm did not converge."))

      }

      ### Paule-Mandel (PM) estimator

      if (method == "PM") {

         if (!allvipos)
            stop(mstyle$stop("PM estimator cannot be used when there are non-positive sampling variances in the data."))

         if (!tau2.fix) {

            if (.QE.func(con$tau2.min, Y=Y, vi=vi, X=X, k=k, objective=0) < k-p) {

               tau2 <- con$tau2.min

            } else {

               if (.QE.func(con$tau2.max, Y=Y, vi=vi, X=X, k=k, objective=0) > k-p) {

                     stop(mstyle$stop("Value of 'tau2.max' too low. Try increasing 'tau2.max' or switch to another 'method'."))

               } else {

                  tau2 <- try(uniroot(.QE.func, interval=c(con$tau2.min, con$tau2.max), tol=con$tol, maxiter=con$maxiter, Y=Y, vi=vi, X=X, k=k, objective=k-p, verbose=verbose, digits=digits, extendInt="no")$root, silent=TRUE)

                  if (inherits(tau2, "try-error"))
                     stop(mstyle$stop("Error in iterative search for tau2 using uniroot()."))

               }

            }

            #W <- diag(wi, nrow=k, ncol=k)
            #stXWX <- .invcalc(X=X, W=W, k=k)
            #P <- W - W %*% X %*% stXWX %*% crossprod(X,W) ### needed for se.tau2 computation below (not when using the simpler equation)

         } else {

            tau2 <- tau2.val

         }

         wi <- 1/(vi + tau2)

      }

      ### Paule-Mandel (PM) estimator (median unbiased version)

      if (method == "PMM") {

         if (!allvipos)
            stop(mstyle$stop("PMM estimator cannot be used when there are non-positive sampling variances in the data."))

         if (!tau2.fix) {

            if (.QE.func(con$tau2.min, Y=Y, vi=vi, X=X, k=k, objective=0) < qchisq(0.5, df=k-p)) {

               tau2 <- con$tau2.min

            } else {

               if (.QE.func(con$tau2.max, Y=Y, vi=vi, X=X, k=k, objective=0) > qchisq(0.5, df=k-p)) {

                  stop(mstyle$stop("Value of 'tau2.max' too low. Try increasing 'tau2.max' or switch to another 'method'."))

               } else {

                  tau2 <- try(uniroot(.QE.func, interval=c(con$tau2.min, con$tau2.max), tol=con$tol, maxiter=con$maxiter, Y=Y, vi=vi, X=X, k=k, objective=qchisq(0.5, df=k-p), verbose=verbose, digits=digits, extendInt="no")$root, silent=TRUE)

                  if (inherits(tau2, "try-error"))
                     stop(mstyle$stop("Error in iterative search for tau2. Try increasing 'tau2.max' or switch to another 'method'."))

               }

            }

         } else {

            tau2 <- tau2.val

         }

         wi <- 1/(vi + tau2)

      }

      ### maximum-likelihood (ML), restricted maximum-likelihood (REML), and empirical Bayes (EB) estimators

      if (is.element(method, c("ML","REML","EB"))) {

         if (is.null(con$tau2.init)) {         ### check if user specified initial value for tau2
            tau2 <- max(0, tau2, con$tau2.min) ### if not, use HE estimate (or con$tau2.min) as initial estimate for tau2
         } else {
            tau2 <- con$tau2.init              ### if yes, use value specified by user
         }

         while (change > con$threshold) {

            if (verbose)
               cat(mstyle$verbose(paste(mstyle$verbose(paste("Iteration", iter, "\ttau^2 =", .fcf(tau2, digits[["var"]]), "\n")))))

            iter <- iter + 1
            old2 <- tau2
            wi   <- 1/(vi + tau2)
            if (any(tau2 + vi < 0))
               stop(mstyle$stop("Some marginal variances are negative."))
            if (any(is.infinite(wi)))
               stop(mstyle$stop("Division by zero when computing the inverse variance weights."))
            W     <- diag(wi, nrow=k, ncol=k)
            stXWX <- .invcalc(X=X, W=W, k=k)
            P     <- W - W %*% X %*% stXWX %*% crossprod(X,W)

            if (method == "ML") {
               PP  <- P %*% P
               adj <- (crossprod(Y,PP) %*% Y - sum(wi)) / sum(wi^2)
            }
            if (method == "REML") {
               PP  <- P %*% P
               adj <- (crossprod(Y,PP) %*% Y - .tr(P)) / .tr(PP)
            }
            if (method == "EB") {
               adj <- (crossprod(Y,P) %*% Y * k/(k-p) - k) / sum(wi)
            }

            adj <- adj * con$stepadj ### apply (user-defined) step adjustment

            while (tau2 + adj < con$tau2.min) ### use step-halving if necessary
               adj <- adj / 2

            tau2   <- ifelse(tau2.fix, tau2.val, tau2 + adj)
            change <- abs(old2 - tau2)

            if (iter > con$maxiter) {
               conv <- 0
               break
            }

         }

         if (conv == 0L)
            stop(mstyle$stop("Fisher scoring algorithm did not converge. See 'help(rma)' for possible remedies."))

         ### check if ll is larger when tau^2 = 0 (only if ll0check=TRUE and only possible/sensible if allvipos and !tau2.fix)
         ### note: this doesn't catch the case where tau^2 = 0 is a local maximum

         if (is.element(method, c("ML","REML")) && con$ll0check && allvipos && !tau2.fix) {

            wi    <- 1/vi
            W     <- diag(wi, nrow=k, ncol=k)
            stXWX <- .invcalc(X=X, W=W, k=k)
            beta  <- stXWX %*% crossprod(X,W) %*% Y
            RSS   <- sum(wi*(yi - X %*% beta)^2)
            if (method == "ML")
               ll0 <- -1/2 * (k)   * log(2*base::pi) - 1/2 * sum(log(vi)) - 1/2 * RSS
            if (method == "REML")
               ll0 <- -1/2 * (k-p) * log(2*base::pi) - 1/2 * sum(log(vi)) - 1/2 * determinant(crossprod(X,W) %*% X, logarithm=TRUE)$modulus - 1/2 * RSS

            wi    <- 1/(vi + tau2)
            if (any(tau2 + vi < 0))
               stop(mstyle$stop("Some marginal variances are negative."))
            if (any(is.infinite(wi)))
               stop(mstyle$stop("Division by zero when computing the inverse variance weights."))
            W     <- diag(wi, nrow=k, ncol=k)
            stXWX <- .invcalc(X=X, W=W, k=k)
            beta  <- stXWX %*% crossprod(X,W) %*% Y
            RSS   <- sum(wi*(yi - X %*% beta)^2)
            if (method == "ML")
               ll <- -1/2 * (k)   * log(2*base::pi) - 1/2 * sum(log(vi + tau2)) - 1/2 * RSS
            if (method == "REML")
               ll <- -1/2 * (k-p) * log(2*base::pi) - 1/2 * sum(log(vi + tau2)) - 1/2 * determinant(crossprod(X,W) %*% X, logarithm=TRUE)$modulus - 1/2 * RSS

            if (ll0 - ll > con$tol && tau2 > con$threshold) {
               warning(mstyle$warning("Fisher scoring algorithm may have gotten stuck at a local maximum.\n  Setting tau^2 = 0. Check the profile likelihood plot with profile()."), call.=FALSE)
               tau2 <- 0
            }

         }

         ### need to run this so that wi and P are based on the final tau^2 value

         wi <- 1/(vi + tau2)
         if (any(tau2 + vi < 0))
            stop(mstyle$stop("Some marginal variances are negative."))
         if (any(is.infinite(wi)))
            stop(mstyle$stop("Division by zero when computing the inverse variance weights."))
         W <- diag(wi, nrow=k, ncol=k)
         stXWX <- .invcalc(X=X, W=W, k=k)
         P <- W - W %*% X %*% stXWX %*% crossprod(X,W)

      }

      ### make sure that tau2 is >= con$tau2.min

      tau2 <- max(con$tau2.min, c(tau2))

      ### check if any marginal variances are negative (only possible if user has changed tau2.min)

      if (!is.na(tau2) && any(tau2 + vi < 0))
         stop(mstyle$stop("Some marginal variances are negative."))

      ### verbose output upon convergence for ML/REML/EB estimators

      if (verbose && is.element(method, c("ML","REML","EB"))) {
         cat(mstyle$verbose(paste("Iteration", iter, "\ttau^2 =", .fcf(tau2, digits[["var"]]), "\n")))
         cat(mstyle$verbose(paste("Fisher scoring algorithm converged after", iter, "iterations.\n")))
      }

      ### standard error of the tau^2 estimators (also when the user has fixed/specified a tau^2 value)
      ### see notes.pdf and note: .tr(P%*%P) = sum(P*t(P)) = sum(P*P) (since P is symmetric)

      if (method == "HS")
         se.tau2 <- sqrt(1/sum(wi)^2 * (2*(k-p) + 4*max(tau2,0)*.tr(P) + 2*max(tau2,0)^2*sum(P*P))) ### note: wi = 1/vi
      if (method == "HSk")
         se.tau2 <- k/(k-p) * sqrt(1/sum(wi)^2 * (2*(k-p) + 4*max(tau2,0)*.tr(P) + 2*max(tau2,0)^2*sum(P*P)))
         #se.tau2 <- sqrt(1/trP^2 * (2*(k-p) + 4*max(tau2,0)*.tr(P) + 2*max(tau2,0)^2*sum(P*P))) # trP <- sum(wi) * (k-p) / k
      if (method == "HE")
         se.tau2 <- sqrt(1/(k-p)^2 * (2*sum(PV*t(PV)) + 4*max(tau2,0)*trPV + 2*max(tau2,0)^2*(k-p)))
      if (method == "DL" || method == "DLIT")
         se.tau2 <- sqrt(1/trP^2 * (2*(k-p) + 4*max(tau2,0)*trP + 2*max(tau2,0)^2*sum(P*P)))
      if (method == "GENQ" || method == "GENQM")
         se.tau2 <- sqrt(1/trP^2 * (2*sum(PV*t(PV)) + 4*max(tau2,0)*sum(PV*P) + 2*max(tau2,0)^2*sum(P*P)))
      if (method == "SJ")
         se.tau2 <- sqrt(tau2.0^2/(k-p)^2 * (2*sum(PV*t(PV)) + 4*max(tau2,0)*sum(PV*P) + 2*max(tau2,0)^2*sum(P*P)))
      if (method == "ML")
         se.tau2 <- sqrt(2/sum(wi^2)) ### note: wi = 1/(vi + tau2) for ML, REML, EB, PM, and SJIT
      if (method == "REML")
         se.tau2 <- sqrt(2/sum(P*P))
      if (method == "EB" || method == "PM" || method == "PMM" || method == "SJIT") {
         #V  <- diag(vi, nrow=k, ncol=k)
         #PV <- P %*% V ### note: this is not symmetric
         #se.tau2 <- sqrt((k/(k-p))^2 / sum(wi)^2 * (2*sum(PV*t(PV)) + 4*max(tau2,0)*sum(PV*P) + 2*max(tau2,0)^2*sum(P*P)))
         se.tau2 <- sqrt(2*k^2/(k-p) / sum(wi)^2) ### these two equations are actually identical, but this one is much simpler
      }

      if (k == 1L)
         method <- method.sav

   }

   #########################################################################

   ###### parameter estimation for location-scale model (rma.ls)

   if (model == "rma.ls") {

      if (!is.element(method, c("ML", "REML")))
         stop(mstyle$stop("Location-scale model can only be fitted with ML or REML estimation."))

      tau2.fix <- FALSE

      if (is.numeric(tau2))
         warning(mstyle$warning("Argument 'tau2' ignored for location-scale models."), call.=FALSE)

      ### get optimizer arguments from control argument

      optimizer  <- match.arg(con$optimizer, c("optim","nlminb","uobyqa","newuoa","bobyqa","nloptr","nlm","hjk","nmk","mads","ucminf","optimParallel","constrOptim"))
      optmethod  <- match.arg(con$optmethod, c("Nelder-Mead","BFGS","CG","L-BFGS-B","SANN","Brent"))
      parallel   <- con$parallel
      cl         <- con$cl
      ncpus      <- con$ncpus
      optcontrol <- control[is.na(con.pos)] ### get arguments that are control arguments for optimizer

      if (length(optcontrol) == 0L)
         optcontrol <- list()

      ### if control argument 'ncpus' is larger than 1, automatically switch to optimParallel optimizer

      if (ncpus > 1L)
         optimizer <- "optimParallel"

      ### when using an identify link, force optimizer to constrOptim

      if (link == "identity")
         optimizer <- "constrOptim"

      if (link == "log" && optimizer == "constrOptim")
         stop(mstyle$stop("Cannot use 'constrOptim' optimizer when using a log link."))

      reml <- ifelse(method=="REML", TRUE, FALSE)

      ### set NLOPT_LN_BOBYQA as the default algorithm for nloptr optimizer
      ### and by default use a relative convergence criterion of 1e-8 on the function value

      if (optimizer=="nloptr" && !is.element("algorithm", names(optcontrol)))
         optcontrol$algorithm <- "NLOPT_LN_BOBYQA"

      if (optimizer=="nloptr" && !is.element("ftol_rel", names(optcontrol)))
         optcontrol$ftol_rel <- 1e-8

      ### for mads, set trace=FALSE and tol=1e-6 by default

      if (optimizer=="mads" && !is.element("trace", names(optcontrol)))
         optcontrol$trace <- FALSE

      if (optimizer=="mads" && !is.element("tol", names(optcontrol)))
         optcontrol$tol <- 1e-6

      ### check that the required packages are installed

      if (is.element(optimizer, c("uobyqa","newuoa","bobyqa"))) {
         if (!requireNamespace("minqa", quietly=TRUE))
            stop(mstyle$stop("Please install the 'minqa' package to use this optimizer."))
      }

      if (optimizer == "nloptr") {
         if (!requireNamespace("nloptr", quietly=TRUE))
            stop(mstyle$stop("Please install the 'nloptr' package to use this optimizer."))
      }

      if (is.element(optimizer, c("hjk","nmk","mads"))) {
         if (!requireNamespace("dfoptim", quietly=TRUE))
            stop(mstyle$stop("Please install the 'dfoptim' package to use this optimizer."))
      }

      if (optimizer == "ucminf") {
         if (!requireNamespace("ucminf", quietly=TRUE))
            stop(mstyle$stop("Please install the 'ucminf' package to use this optimizer."))
      }

      if (optimizer == "optimParallel") {
         if (!requireNamespace("optimParallel", quietly=TRUE))
            stop(mstyle$stop("Please install the 'optimParallel' package to use this optimizer."))
      }

      if (!isTRUE(ddd$skiphes) && !requireNamespace("numDeriv", quietly=TRUE))
         stop(mstyle$stop("Please install the 'numDeriv' package to compute the Hessian."))

      ### drop redundant predictors

      tmp <- try(lm(yi ~ Z - 1), silent=TRUE)
      if (inherits(tmp, "lm")) {
         coef.na.Z <- is.na(coef(tmp))
      } else {
         coef.na.Z <- rep(FALSE, NCOL(Z))
      }
      if (any(coef.na.Z)) {
         warning(mstyle$warning("Redundant predictors dropped from the scale model."), call.=FALSE)
         Z   <- Z[,!coef.na.Z,drop=FALSE]
         Z.f <- Z.f[,!coef.na.Z,drop=FALSE]
      }

      ### check whether intercept is included and if yes, move it to the first column (NAs already removed, so na.rm=TRUE for any() not necessary)

      is.int <- apply(Z, 2, .is.intercept)
      if (any(is.int)) {
         Z.int.incl <- TRUE
         int.indx   <- which(is.int, arr.ind=TRUE)
         Z          <- cbind(intrcpt=1,   Z[,-int.indx, drop=FALSE]) ### this removes any duplicate intercepts
         Z.f        <- cbind(intrcpt=1, Z.f[,-int.indx, drop=FALSE]) ### this removes any duplicate intercepts
         Z.intercept <- TRUE ### set intercept appropriately so that the predict() function works
      } else {
         Z.int.incl <- FALSE
      }

      q <- NCOL(Z) ### number of columns in Z (including the intercept if it is included)

      ### check whether model matrix is of full rank

      if (any(eigen(crossprod(Z), symmetric=TRUE, only.values=TRUE)$values <= con$evtol))
         stop(mstyle$stop("Model matrix for scale part of the model not of full rank. Cannot fit model."))

      ### check whether this is an intercept-only model

      is.int <- apply(Z, 2, .is.intercept)
      if ((q == 1L) && is.int) {
         Z.int.only <- TRUE
      } else {
         Z.int.only <- FALSE
      }

      ### checks on alpha argument

      if (missing(alpha) || is.null(alpha) || all(is.na(alpha))) {
         alpha <- rep(NA, q)
      } else {
         if (length(alpha) == 1L)
            alpha <- rep(alpha, q)
         if (length(alpha) != q)
            stop(mstyle$stop(paste0("Length of 'alpha' argument (", length(alpha), ") does not match actual number of parameters (", q, ").")))
      }

      ### rescale Z matrix (only for models with moderators and models including an intercept term and when alpha[1] is not fixed)

      if (!Z.int.only && Z.int.incl && con$scaleZ && is.na(alpha[1])) {
         Zsave <- Z
         meanZ <- colMeans(Z[, 2:q, drop=FALSE])
         sdZ   <- apply(Z[, 2:q, drop=FALSE], 2, sd) ### consider using colSds() from matrixStats package
         is.d  <- apply(Z, 2, .is.dummy) ### is each column a dummy variable (i.e., only 0s and 1s)?
         mZ    <- rbind(c(intrcpt=1, -1*ifelse(is.d[-1], 0, meanZ/sdZ)), cbind(0, diag(ifelse(is.d[-1], 1, 1/sdZ), nrow=length(is.d)-1, ncol=length(is.d)-1)))
         imZ   <- try(suppressWarnings(solve(mZ)), silent=TRUE)
         Z[,!is.d] <- apply(Z[, !is.d, drop=FALSE], 2, scale) ### rescale the non-dummy variables
         if (any(!is.na(alpha))) {
            if (inherits(imZ, "try-error"))
               stop(mstyle$stop("Unable to rescale starting values for the scale parameters."))
            alpha <- diag(imZ) * alpha
         }
      } else {
         mZ <- NULL
      }

      if (k == 1L && Z.int.only) {
         if (link == "log")
            con$alpha.init <- -10000
         if (link == "identity")
            con$alpha.init <- 0.00001
      }

      ### set/transform/check alpha.init

      if (verbose > 1)
         message(mstyle$message("Extracting/computing initial values ..."))

      if (is.null(con$alpha.init)) {

         fit <- suppressWarnings(rma(yi, vi, mods=X, intercept=FALSE, method="HE"))
         tmp <- rstandard(fit)

         if (link == "log") {

            tmp <- rma(log(tmp$resid^2), 4/tmp$resid^2*tmp$se^2, mods=Z, intercept=FALSE, method="FE")
            #tmp <- rma(log(tmp$resid^2), tmp$se^2, mods=Z, intercept=FALSE, method="FE")
            #tmp <- rma(log(tmp$resid^2), 1, mods=Z, intercept=FALSE, method="FE")
            alpha.init <- coef(tmp)

         }

         if (link == "identity") {

            #tmp <- rma(tmp$resid^2, 4*tmp$resid^2*tmp$se^2, mods=Z, intercept=FALSE, method="FE")
            tmp <- rma(tmp$resid^2, tmp$se^2, mods=Z, intercept=FALSE, method="FE")
            #tmp <- rma(tmp$resid^2, 1, mods=Z, intercept=FALSE, method="FE")
            alpha.init <- coef(tmp)
            if (any(Z %*% alpha.init < 0))
               alpha.init <- ifelse(is.int, fit$tau2+.01, 0)
            if (any(Z %*% alpha.init < 0))
               stop(mstyle$stop("Unable to find suitable starting values for the scale parameters."))

         }

      } else {

         alpha.init <- con$alpha.init

         if (!Z.int.only && Z.int.incl && con$scaleZ && is.na(alpha[1])) {
            if (inherits(imZ, "try-error"))
               stop(mstyle$stop("Unable to rescale starting values for the scale parameters."))
            alpha.init <- c(imZ %*% cbind(alpha.init))
         }

         if (link == "identity" && any(Z %*% alpha.init < 0))
            stop(mstyle$stop("Starting values for the scale parameters lead to one or more negative tau^2 values."))

      }

      if (length(alpha.init) != q)
         stop(mstyle$stop(paste0("Length of 'alpha.init' argument (", length(alpha.init), ") does not match actual number of parameters (", q, ").")))

      if (anyNA(alpha.init))
         stop(mstyle$stop("No missing values allowed in 'alpha.init'."))

      if (verbose > 1)
         message(mstyle$message("Estimating scale parameters ...\n"))

      ### obtain initial values for beta (only need this when optimizing over beta and alpha jointly)

      #wi <- 1/vi
      #W  <- diag(wi, nrow=k, ncol=k)
      #stXWX <- .invcalc(X=X, W=W, k=k)
      #beta.init <- stXWX %*% crossprod(X,W) %*% Y

      if (is.element(optimizer, c("optim","constrOptim"))) {
         par.arg <- "par"
         ctrl.arg <- ", control=optcontrol"
      }
      if (optimizer=="nlminb") {
         par.arg <- "start"
         ctrl.arg <- ", control=optcontrol"
      }
      if (is.element(optimizer, c("uobyqa","newuoa","bobyqa"))) {
         par.arg <- "par"
         optimizer <- paste0("minqa::", optimizer) ### need to use this since loading nloptr masks bobyqa() and newuoa() functions
         ctrl.arg <- ", control=optcontrol"
      }
      if (optimizer=="nloptr") {
         par.arg <- "x0"
         optimizer <- paste0("nloptr::nloptr") ### need to use this due to requireNamespace()
         ctrl.arg <- ", opts=optcontrol"
      }
      if (optimizer=="nlm") {
         par.arg <- "p" ### because of this, must use argument name pX for p (number of columns in X matrix)
         ctrl.arg <- paste(names(optcontrol), unlist(optcontrol), sep="=", collapse=", ")
         if (nchar(ctrl.arg) != 0L)
            ctrl.arg <- paste0(", ", ctrl.arg)
      }
      if (is.element(optimizer, c("hjk","nmk","mads"))) {
         par.arg <- "par"
         optimizer <- paste0("dfoptim::", optimizer) ### need to use this so that the optimizers can be found
         ctrl.arg <- ", control=optcontrol"
      }
      if (optimizer=="ucminf") {
         par.arg <- "par"
         optimizer <- paste0("ucminf::ucminf") ### need to use this due to requireNamespace()
         ctrl.arg <- ", control=optcontrol"
      }
      if (optimizer=="optimParallel") {

         par.arg <- "par"
         optimizer <- paste0("optimParallel::optimParallel") ### need to use this due to requireNamespace()
         ctrl.arg <- ", control=optcontrol, parallel=parallel"

         parallel$cl <- NULL

         if (is.null(cl)) {

            ncpus <- as.integer(ncpus)

            if (ncpus < 1L)
               stop(mstyle$stop("Control argument 'ncpus' must be >= 1."))

            cl <- parallel::makePSOCKcluster(ncpus)
            on.exit(parallel::stopCluster(cl), add=TRUE)

         } else {

            if (!inherits(cl, "SOCKcluster"))
               stop(mstyle$stop("Specified cluster is not of class 'SOCKcluster'."))

         }

         parallel$cl <- cl

         if (is.null(parallel$forward))
            parallel$forward <- FALSE

         if (is.null(parallel$loginfo)) {
            if (verbose) {
               parallel$loginfo <- TRUE
            } else {
               parallel$loginfo <- FALSE
            }
         }

      }

      #return(list(con=con, optimizer=optimizer, optmethod=optmethod, optcontrol=optcontrol, ctrl.arg=ctrl.arg))

      if (link == "log") {

         optcall <- paste(optimizer, "(", par.arg, "=alpha.init, .ll.rma.ls, ", ifelse(optimizer=="optim", "method=optmethod, ", ""),
                                                       "yi=yi, vi=vi, X=X, Z=Z, reml=reml, k=k, pX=p, alpha.val=alpha, verbose=verbose, digits=digits,
                                                       #hessian=TRUE,
                                                       REMLf=con$REMLf, link=link, mZ=mZ", ctrl.arg, ")\n", sep="")
         #return(optcall)
         if (verbose) {
            opt.res <- try(eval(parse(text=optcall)), silent=!verbose)
         } else {
            opt.res <- try(suppressWarnings(eval(parse(text=optcall))), silent=!verbose)
         }

      }

      if (link == "identity") {

         opt.res <- try(constrOptim(theta=alpha.init, f=.ll.rma.ls, grad=NULL, ui=Z, ci=rep(0,k),
                                    yi=yi, vi=vi, X=X, Z=Z, reml=reml, k=k, pX=p, alpha.val=alpha, verbose=verbose, digits=digits,
                                    REMLf=con$REMLf, link=link, mZ=mZ), silent=!verbose)

      }

      #return(opt.res)

      if (optimizer == "optimParallel::optimParallel" && verbose) {
         tmp <- capture.output(print(opt.res$loginfo))
         .print.output(tmp, mstyle$verbose)
      }

      if (inherits(opt.res, "try-error"))
         stop(mstyle$stop("Error during the optimization. Use verbose=TRUE and see help(rma) for more details on the optimization routines."))

      ### convergence checks

      if (is.element(optimizer, c("optim","constrOptim","nlminb","dfoptim::hjk","dfoptim::nmk","optimParallel::optimParallel")) && opt.res$convergence != 0)
         stop(mstyle$stop(paste0("Optimizer (", optimizer, ") did not achieve convergence (convergence = ", opt.res$convergence, ").")))

      if (is.element(optimizer, c("dfoptim::mads")) && opt.res$convergence > optcontrol$tol)
         stop(mstyle$stop(paste0("Optimizer (", optimizer, ") did not achieve convergence (convergence = ", opt.res$convergence, ").")))

      if (is.element(optimizer, c("minqa::uobyqa","minqa::newuoa","minqa::bobyqa")) && opt.res$ierr != 0)
         stop(mstyle$stop(paste0("Optimizer (", optimizer, ") did not achieve convergence (ierr = ", opt.res$ierr, ").")))

      if (optimizer=="nloptr::nloptr" && !(opt.res$status >= 1 && opt.res$status <= 4))
         stop(mstyle$stop(paste0("Optimizer (", optimizer, ") did not achieve convergence (status = ", opt.res$status, ").")))

      if (optimizer=="ucminf::ucminf" && !(opt.res$convergence == 1 || opt.res$convergence == 2))
         stop(mstyle$stop(paste0("Optimizer (", optimizer, ") did not achieve convergence (convergence = ", opt.res$convergence, ").")))

      if (verbose > 2) {
         cat("\n")
         tmp <- capture.output(print(opt.res))
         .print.output(tmp, mstyle$verbose)
      }

      ### copy estimated values to 'par' so code below works

      if (optimizer=="nloptr::nloptr")
         opt.res$par <- opt.res$solution
      if (optimizer=="nlm")
         opt.res$par <- opt.res$estimate

      ### replace fixed alpha values in opt.res$par

      opt.res$par <- ifelse(is.na(alpha), opt.res$par, alpha)

      ### try to compute vcov matrix for scale parameter estimates

      H <- NA
      vb.alpha <- matrix(NA_real_, nrow=q, ncol=q)
      hest <- is.na(alpha)

      if (any(hest) && !isTRUE(ddd$skiphes)) {

         if (verbose > 1)
            message(mstyle$message("\nComputing Hessian ..."))

         H <- try(numDeriv::hessian(func=.ll.rma.ls, x=opt.res$par, method.args=con$hessianCtrl, yi=yi, vi=vi, X=X,
                                    Z=Z, reml=reml, k=k, pX=p, alpha.val=alpha, verbose=FALSE, digits=digits, REMLf=con$REMLf, link=link, mZ=mZ), silent=TRUE)

         if (inherits(H, "try-error")) {

            warning(mstyle$warning("Error when trying to compute Hessian."), call.=FALSE)

         } else {

            H.hest <- H[hest, hest, drop=FALSE]
            iH.hest <- try(suppressWarnings(chol2inv(chol(H.hest))), silent=TRUE)

            if (inherits(iH.hest, "try-error") || any(is.na(iH.hest)) || any(is.infinite(iH.hest))) {

               warning(mstyle$warning("Error when trying to invert Hessian."), call.=FALSE)

            } else {

               vb.alpha[hest, hest] <- iH.hest

            }

         }

      }

      ### get scale parameter estimates

      alpha.val <- alpha
      alpha <- cbind(opt.res$par)

      ### scale back alpha and vb.alpha

      if (!Z.int.only && Z.int.incl && con$scaleZ && is.na(alpha.val[1])) {
         alpha <- mZ %*% alpha
         vb.alpha[!hest,] <- 0
         vb.alpha[,!hest] <- 0
         vb.alpha <- mZ %*% vb.alpha %*% t(mZ)
         vb.alpha[!hest,] <- NA
         vb.alpha[,!hest] <- NA
         Z <- Zsave
      }

      ### set/check 'att' argument

      att <- .set.btt(att, q, Z.int.incl, colnames(Z))
      m.alpha <- length(att) ### number of alphas to test (m = z if all betas are tested)

      ### QM calculation

      QM.alpha <- try(as.vector(t(alpha)[att] %*% chol2inv(chol(vb.alpha[att,att])) %*% alpha[att]), silent=TRUE)

      if (inherits(QM.alpha, "try-error"))
         QM.alpha <- NA

      se.alpha <- sqrt(diag(vb.alpha))

      rownames(alpha) <- rownames(vb.alpha) <- colnames(vb.alpha) <- colnames(Z)

      names(se.alpha) <- NULL
      zval.alpha  <- c(alpha/se.alpha)

      if (is.element(test, c("t"))) {
         dfs.alpha <- k-q
         QM.alpha <- QM.alpha / m.alpha
         QMp.alpha  <- pf(QM.alpha, df1=m.alpha, df2=dfs.alpha, lower.tail=FALSE)
         pval.alpha <- 2*pt(abs(zval.alpha), df=dfs.alpha, lower.tail=FALSE)
         crit <- qt(level/2, df=dfs.alpha, lower.tail=FALSE)
      } else {
         dfs.alpha <- NA
         QMp.alpha  <- pchisq(QM.alpha, df=m.alpha, lower.tail=FALSE)
         pval.alpha  <- 2*pnorm(abs(zval.alpha), lower.tail=FALSE)
         crit <- qnorm(level/2, lower.tail=FALSE)
      }

      ci.lb.alpha <- c(alpha - crit * se.alpha)
      ci.ub.alpha <- c(alpha + crit * se.alpha)

      if (link == "log")
         tau2 <- exp(as.vector(Z %*% alpha))
      if (link == "identity")
         tau2 <- as.vector(Z %*% alpha)

   }

   ### fixed-effects model (note: set tau2 to zero even when tau2 is specified)

   if (method == "FE")
      tau2 <- 0

   #########################################################################

   ###### model fitting, test statistics, and confidence intervals

   if (verbose > 1)
      message(mstyle$message("\nModel fitting ..."))

   wi <- 1/(vi + tau2)
   W  <- diag(wi, nrow=k, ncol=k)
   M  <- diag(vi + tau2, nrow=k, ncol=k)

   if (weighted) {

      #########################
      ### weighted analysis ###
      #########################

      ### fit model with weighted estimation

      if (is.null(weights) || is.element(test, c("knha","adhoc"))) {

         ### if no weights are specified, use default inverse variance weights, that is, 1/vi or 1/(vi + tau2)
         ### also, even with weights, if test="knha" or "adhoc", need to run this to get RSS.knha

         ### if any vi = 0 and tau^2 is estimated to be 0 (or is set to 0 for a FE model), then get Inf for wi

         if (any(is.infinite(wi)))
            stop(mstyle$stop("Division by zero when computing the inverse variance weights."))

         stXWX <- .invcalc(X=X, W=W, k=k)
         beta  <- stXWX %*% crossprod(X,W) %*% Y
         vb    <- stXWX
         RSS.f <- sum(wi*(yi - X %*% beta)^2)
         #P     <- W - W %*% X %*% stXWX %*% crossprod(X,W)
         #RSS.f <- crossprod(Y,P) %*% Y
         RSS.knha <- RSS.f

      }

      if (!is.null(weights)) {

         ### if weights are specified, use them (note: RSS.f is recomputed if test="knha" or "adhoc")

         A     <- diag(weights, nrow=k, ncol=k)
         stXAX <- .invcalc(X=X, W=A, k=k)
         beta  <- stXAX %*% crossprod(X,A) %*% Y
         vb    <- stXAX %*% t(X) %*% A %*% M %*% A %*% X %*% stXAX
         RSS.f <- sum(wi*(yi - X %*% beta)^2)
         #P     <- W - W %*% X %*% stXAX %*% t(X) %*% A - A %*% X %*% stXAX %*% t(X) %*% W + A %*% X %*% stXAX %*% t(X) %*% W %*% X %*% stXAX %*% t(X) %*% A
         #RSS.f <- crossprod(Y,P) %*% Y

      }

      #return(list(beta=beta, vb=vb, se=sqrt(diag(vb)), RSS.f=RSS.f))

      ### calculate scaling factor for Knapp & Hartung method
      ### note: catch cases where RSS.knha is extremely small, which is probably due to all yi being equal
      ### then set s2w to 0 (to avoid the strange looking output we would obtain if we don't do this)

      if (is.element(test, c("knha","adhoc"))) {

         if (RSS.knha <= .Machine$double.eps) {
            s2w <- 0
         } else {
            s2w <- RSS.knha / (k-p)
         }

      }

   } else {

      ###########################
      ### unweighted analysis ###
      ###########################

      ### fit model with unweighted estimation
      ### note: 1) if user has specified weights, they are ignored
      ###       2) but if method="GENQ/GENQM", they were used to estimate tau^2

      stXX  <- .invcalc(X=X, W=diag(k), k=k)
      beta  <- stXX %*% crossprod(X,Y)
      vb    <- tcrossprod(stXX,X) %*% M %*% X %*% stXX
      RSS.f <- sum(wi*(yi - X %*% beta)^2)
      #P     <- W - W %*% X %*% tcrossprod(stXX,X) - X %*% stXX %*% crossprod(X,W) + X %*% stXX %*% crossprod(X,W) %*% X %*% tcrossprod(stXX,X)
      #RSS.f <- crossprod(Y,P) %*% Y

      ### calculate scaling factor for Knapp & Hartung method

      if (is.element(test, c("knha","adhoc"))) {

         if (any(is.infinite(wi)))
            stop(mstyle$stop("Division by zero when computing the inverse variance weights."))

         stXWX     <- .invcalc(X=X, W=W, k=k)
         beta.knha <- stXWX %*% crossprod(X,W) %*% Y
         RSS.knha  <- sum(wi*(yi - X %*% beta.knha)^2)
         #P         <- W - W %*% X %*% stXWX %*% crossprod(X,W)
         #RSS.knha  <- c(crossprod(Y,P) %*% Y)

         if (RSS.knha <= .Machine$double.eps) {
            s2w <- 0
         } else {
            s2w <- RSS.knha / (k-p)
         }

      }

   }

   if (verbose > 1)
      message(mstyle$message("Conducting tests of the fixed effects ..."))

   ### the Knapp & Hartung method as described in the literature is for random/mixed-effects models

   if (method == "FE" && is.element(test, c("knha","adhoc")))
      warning(mstyle$warning("Knapp & Hartung method is not meant to be used in the context of FE models."), call.=FALSE)

   ### Knapp & Hartung method with ad-hoc correction so that the scale factor is always >= 1

   if (test == "adhoc")
      s2w[s2w < 1] <- 1

   ### for Knapp & Hartung method, apply scaling to vb

   vb <- s2w * vb

   ### QM calculation

   QM <- try(as.vector(t(beta)[btt] %*% chol2inv(chol(vb[btt,btt])) %*% beta[btt]), silent=TRUE)

   if (inherits(QM, "try-error"))
      QM <- NA

   rownames(beta) <- rownames(vb) <- colnames(vb) <- colnames(X) ### may not be needed, but what's the harm?

   se <- sqrt(diag(vb))
   names(se) <- NULL
   zval <- c(beta/se)

   if (is.element(test, c("knha","adhoc","t"))) {
      dfs <- k-p
      QM  <- QM / m
      if (dfs > 0) {
         QMp  <- pf(QM, df1=m, df2=dfs, lower.tail=FALSE)
         pval <- 2*pt(abs(zval), df=dfs, lower.tail=FALSE)
         crit <- qt(level/2, df=dfs, lower.tail=FALSE)
      } else {
         QMp  <- NaN
         pval <- NaN
         crit <- NaN
      }
   } else {
      dfs  <- NA
      QMp  <- pchisq(QM, df=m, lower.tail=FALSE)
      pval <- 2*pnorm(abs(zval), lower.tail=FALSE)
      crit <- qnorm(level/2, lower.tail=FALSE)
   }

   ci.lb <- c(beta - crit * se)
   ci.ub <- c(beta + crit * se)

   #########################################################################

   ### heterogeneity test (Wald-type test of the extra coefficients in the saturated model)

   if (verbose > 1)
      message(mstyle$message("Conducting heterogeneity test ..."))

   if (allvipos) {

      ### heterogeneity test (always uses inverse variance method)
      ### note: this is unaffected by the 'weighted' argument, since under H0, the same parameters are estimated
      ### and weighted estimation provides the most efficient estimates; therefore, also any arbitrary weights
      ### specified by the user are not relevant here (different from what the metan command in Stata does!)
      ### see also: Chen, Z., Ng, H. K. T., & Nadarajah, S. (2014). A note on Cochran test for homogeneity in
      ### one-way ANOVA and meta-analysis. Statistical Papers, 55(2), 301-310. This shows that the weights used
      ### are not relevant.

      if (k > p) {

         wi    <- 1/vi
         W.FE  <- diag(wi, nrow=k, ncol=k) ### note: ll.REML below involves W, so cannot overwrite W
         stXWX <- .invcalc(X=X, W=W.FE, k=k)
         P     <- W.FE - W.FE %*% X %*% stXWX %*% crossprod(X,W.FE) ### need P below for calculation of I^2
         QE    <- max(0, c(crossprod(Y,P) %*% Y))
         #beta.FE <- stXWX %*% crossprod(X,W.FE) %*% Y
         #QE    <- max(0, sum(wi*(yi - X %*% beta.FE)^2))
         QEp   <- pchisq(QE, df=k-p, lower.tail=FALSE)

         ### calculation of 'typical' sampling variance

         #vt <- (k-1) / (sum(wi) - sum(wi^2)/sum(wi)) ### this only applies to the RE model
         #vt <- 1/mean(wi) ### harmonic mean of vi's (see Takkouche et al., 1999)
         vt <- (k-p) / .tr(P)

         ### calculation of I^2 and H^2

         if (method == "FE") {
            I2 <- max(0, 100 * (QE - (k-p)) / QE)
            H2 <- QE / (k-p)
         } else {
            I2 <- 100 * mean(tau2) / (vt + mean(tau2)) ### must use mean(tau2) in case tau2 is vector from location-scale model
            H2 <- mean(tau2) / vt + 1                  ### must use mean(tau2) in case tau2 is vector from location-scale model
         }

      } else {

         QE  <- 0
         QEp <- 1
         I2  <- 0
         H2  <- 1
         vt  <- 0

      }

   } else {

      if (!vi0)
         warning(mstyle$warning(paste0("Cannot compute ", ifelse(int.only, "Q", "QE"), "-test, I^2, or H^2 when there are non-positive sampling variances in the data.")), call.=FALSE)

      vt <- NA

   }

   #########################################################################

   ### compute pseudo R^2 statistic for mixed-effects models with an intercept (only for rma.uni models)

   if (!int.only && int.incl && method != "FE" && model == "rma.uni" && !isTRUE(ddd$skipr2)) {

      if (verbose > 1) {
         message(mstyle$message("Fitting RE model for R^2 computation ..."))
         res.RE <- try(rma.uni(yi, vi, weights=weights, method=method, weighted=weighted, test=test, verbose=ifelse(verbose, TRUE, FALSE), control=con, digits=digits), silent=FALSE)
      } else {
         res.RE <- try(suppressWarnings(rma.uni(yi, vi, weights=weights, method=method, weighted=weighted, test=test, verbose=ifelse(verbose, TRUE, FALSE), control=con, digits=digits)), silent=TRUE)
      }

      if (!inherits(res.RE, "try-error")) {

         tau2.RE <- res.RE$tau2

         if (identical(tau2.RE,0)) {
            R2 <- 0
         } else {
            R2 <- max(0, 100 * (tau2.RE - tau2) / tau2.RE)
         }

      } else {

         R2 <- NA

      }

   } else {

      R2 <- NULL

   }

   #########################################################################

   ###### fit statistics

   if (verbose > 1)
      message(mstyle$message("Computing fit statistics and log likelihood ..."))

   ### note: tau2 is not counted as a parameter when it was fixed by the user (same for fixed alpha values)
   parms <- p + ifelse(model == "rma.uni", ifelse(method == "FE" || tau2.fix, 0, 1), sum(is.na(alpha.val)))

   ll.ML    <- -1/2 * (k) * log(2*base::pi) - 1/2 * sum(log(vi + tau2)) - 1/2 * RSS.f
   ll.REML  <- -1/2 * (k-p) * log(2*base::pi) + ifelse(con$REMLf, 1/2 * determinant(crossprod(X), logarithm=TRUE)$modulus, 0) +
               -1/2 * sum(log(vi + tau2)) - 1/2 * determinant(crossprod(X,W) %*% X, logarithm=TRUE)$modulus - 1/2 * RSS.f

   if (k > p) {
      dev.ML <- -2 * (ll.ML - sum(dnorm(yi, mean=yi, sd=sqrt(vi), log=TRUE)))
   } else {
      dev.ML <- 0
   }
   AIC.ML    <- -2 * ll.ML   + 2*parms
   BIC.ML    <- -2 * ll.ML   +   parms * log(k)
   AICc.ML   <- -2 * ll.ML   + 2*parms * max(k, parms+2) / (max(k, parms+2) - parms - 1)
   dev.REML  <- -2 * (ll.REML - 0) ### saturated model has ll = 0 when using the full REML likelihood
   AIC.REML  <- -2 * ll.REML + 2*parms
   BIC.REML  <- -2 * ll.REML +   parms * log(k-p)
   AICc.REML <- -2 * ll.REML + 2*parms * max(k-p, parms+2) / (max(k-p, parms+2) - parms - 1)

   fit.stats <- matrix(c(ll.ML, dev.ML, AIC.ML, BIC.ML, AICc.ML, ll.REML, dev.REML, AIC.REML, BIC.REML, AICc.REML), ncol=2, byrow=FALSE)
   dimnames(fit.stats) <- list(c("ll","dev","AIC","BIC","AICc"), c("ML","REML"))
   fit.stats <- data.frame(fit.stats)

   #########################################################################

   ###### prepare output

   if (verbose > 1)
      message(mstyle$message("Preparing output ..."))

   p.eff <- p
   k.eff <- k

   if (is.null(ddd$outlist)) {

      res <- list(b=beta, beta=beta, se=se, zval=zval, pval=pval, ci.lb=ci.lb, ci.ub=ci.ub, vb=vb,
                  tau2=tau2, se.tau2=se.tau2, tau2.fix=tau2.fix, tau2.f=tau2,
                  I2=I2, H2=H2, R2=R2, vt=vt,
                  QE=QE, QEp=QEp, QM=QM, QMp=QMp,
                  k=k, k.f=k.f, k.eff=k.eff, k.all=k.all, p=p, p.eff=p.eff, parms=parms,
                  int.only=int.only, int.incl=int.incl, intercept=intercept, allvipos=allvipos, coef.na=coef.na,
                  yi=yi, vi=vi, X=X, weights=weights, yi.f=yi.f, vi.f=vi.f, X.f=X.f, weights.f=weights.f, M=M,
                  ai.f=ai.f, bi.f=bi.f, ci.f=ci.f, di.f=di.f,
                  x1i.f=x1i.f, x2i.f=x2i.f, t1i.f=t1i.f, t2i.f=t2i.f, ni=ni, ni.f=ni.f,
                  ids=ids, not.na=not.na, subset=subset, slab=slab, slab.null=slab.null,
                  measure=measure, method=method, model=model, weighted=weighted,
                  test=test, dfs=dfs, s2w=s2w, btt=btt, m=m,
                  digits=digits, level=level, control=control, verbose=verbose,
                  add=add, to=to, drop00=drop00,
                  fit.stats=fit.stats,
                  formula.yi=formula.yi, formula.mods=formula.mods,
                  version=packageVersion("metafor"), call=mf)

   }

   if (!is.null(ddd$outlist)) {

      if (ddd$outlist == "minimal") {
         res <- list(b=beta, beta=beta, se=se, zval=zval, pval=pval, ci.lb=ci.lb, ci.ub=ci.ub, vb=vb,
                     tau2=tau2, se.tau2=se.tau2, tau2.fix=tau2.fix,
                     I2=I2, H2=H2, R2=R2,
                     QE=QE, QEp=QEp, QM=QM, QMp=QMp,
                     k=k, k.eff=k.eff, p=p, p.eff=p.eff, parms=parms,
                     int.only=int.only,
                     measure=measure, method=method, model=model,
                     test=test, dfs=dfs, btt=btt, m=m,
                     digits=digits,
                     fit.stats=fit.stats)
      } else {
         res <- eval(parse(text=paste0("list(", ddd$outlist, ")")))
      }

   }

   if (model == "rma.ls") {

      res$alpha          <- alpha
      res$vb.alpha       <- vb.alpha
      res$se.alpha       <- se.alpha
      res$zval.alpha     <- zval.alpha
      res$pval.alpha     <- pval.alpha
      res$ci.lb.alpha    <- ci.lb.alpha
      res$ci.ub.alpha    <- ci.ub.alpha
      res$alpha.fix      <- !is.na(alpha.val)
      res$q              <- q
      res$alphas         <- q
      res$link           <- link
      res$Z              <- Z
      res$Z.f            <- Z.f
      res$tau2.f         <- rep(NA, k.f)
      res$tau2.f[not.na] <- tau2
      res$att            <- att
      res$m.alpha        <- m.alpha
      res$dfs.alpha      <- dfs.alpha
      res$QM.alpha       <- QM.alpha
      res$QMp.alpha      <- QMp.alpha
      res$formula.scale  <- formula.scale
      res$Z.int.incl     <- Z.int.incl
      res$Z.int.only     <- Z.int.only
      res$H              <- H

   }

   time.end <- proc.time()
   res$time <- unname(time.end - time.start)[3]

   if (.isTRUE(ddd$time))
      .print.time(res$time)

   if (verbose || .isTRUE(ddd$time))
      cat("\n")

   if (model == "rma.ls") {
      class(res) <- c("rma.ls", "rma.uni", "rma")
   } else {
      class(res) <- c("rma.uni", "rma")
   }
   return(res)

}
