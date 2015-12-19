rma <- rma.uni <- function(yi, vi, sei, weights, ai, bi, ci, di, n1i, n2i, x1i, x2i, t1i, t2i, m1i, m2i, sd1i, sd2i, xi, mi, ri, ti, sdi, ni, mods, #scale,
measure="GEN", intercept=TRUE,
data, slab, subset,
add=1/2, to="only0", drop00=FALSE, vtype="LS",
method="REML", weighted=TRUE, knha=FALSE,
level=95, digits=4, btt, tau2, verbose=FALSE, control) {

   #########################################################################

   ###### setup

   ### check for some incorrect argument specifications
   ### (arguments "to" and "vtype" are checked inside escalc function)

   if (!is.element(measure, c("GEN",
                              "RR","OR","PETO","RD","AS","PHI","YUQ","YUY","RTET", ### 2x2 table measures
                              "PBIT","OR2D","OR2DN","OR2DL",                       ### - transformations to SMD
                              "IRR","IRD","IRSD",                                  ### two-group person-time data measures
                              "MD","SMD","SMDH","ROM",                             ### two-group mean/SD measures
                              "RPB","RBIS","D2OR","D2ORN","D2ORL",                 ### - transformations to r_PB, r_BIS, and log(OR)
                              "COR","UCOR","ZCOR",                                 ### correlations (raw and r-to-z transformed)
                              "PR","PLN","PLO","PAS","PFT",                        ### single proportions (and transformations thereof)
                              "IR","IRLN","IRS","IRFT",                            ### single-group person-time data (and transformations thereof)
                              "MN","MC","SMCC","SMCR","SMCRH","ROMC",              ### raw/standardized mean change and log(ROM) for dependent samples
                              "ARAW","AHW","ABT")))                                ### alpha (and transformations thereof)
      stop("Unknown 'measure' specified.")

   if (!is.element(method, c("FE","HS","HE","DL","GENQ","SJ","ML","REML","EB","DLIT","SJIT","PM")))
      stop("Unknown 'method' specified.")

   ### in case user specifies more than one add/to value (as one can do with rma.mh() and rma.peto())
   ### (any kind of continuity correction is directly applied to the outcomes, which are then analyzed as such)

   if (length(add) > 1)
      add <- add[1]

   if (length(to) > 1)
      to <- to[1]

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop("Unknown 'na.action' specified under options().")

   if (missing(tau2))
      tau2 <- NULL

   if (missing(control))
      control <- list()

   if (!(is.logical(knha) || (is.character(knha) && is.element(knha, c("adhoc", "tdist")))))
      stop("Invalid option selected for 'knha' argument.")

   very.verbose <- ifelse(!is.logical(verbose) && verbose > 1, TRUE, FALSE)

   #########################################################################

   if (very.verbose)
      message("Extracting/computing yi/vi values ...")

   ### check if data argument has been specified

   if (missing(data))
      data <- NULL

   if (is.null(data)) {
      data <- sys.frame(sys.parent())
   } else {
      if (!is.data.frame(data)) {
         data <- data.frame(data)
      }
   }

   ### extract yi (either NULL if not specified, a vector, a formula, or an escalc object)

   mf <- match.call()
   mf.yi <- mf[[match("yi", names(mf))]]
   yi <- eval(mf.yi, data, enclos=sys.frame(sys.parent())) ### NULL if user does not specify this

   ### if yi is not NULL and it is an escalc object, then use that object in place of the data argument

   if (!is.null(yi) && is.element("escalc", class(yi)))
      data <- yi

   ### extract weights, slab, subset, mods, and tau2 values, possibly from the data frame specified via data or yi (arguments not specified are NULL)

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

   is.formula <- FALSE

   #if (measure == "GEN") { ### this way, one *must* use measure="GEN" when specifying yi/vi directly
   if (!is.null(yi)) {      ### this way, one can specify yi/vi directly, but still set 'measure' to something else

      ### if yi is not NULL, then yi now either contains the yi values, a formula, or an escalc object

      ### if yi is a formula, extract yi and X (this overrides anything specified via the mods argument further below)

      if (is.element("formula", class(yi))) {
         options(na.action = "na.pass")                   ### set na.action to na.pass, so that NAs are not filtered out (we'll do that later)
         mods <- model.matrix(yi, data=data)              ### extract model matrix (now mods is no longer a formula, so part further below is skipped)
         attr(mods, "assign") <- NULL                     ### strip assign attribute (not needed at the moment)
         yi <- model.response(model.frame(yi, data=data)) ### extract dependent variable from model frame (the yi values)
         options(na.action = na.act)                      ### set na.action back to na.act
         names(yi) <- NULL                                ### strip names (1:k) from yi (so res$yi is the same whether yi is a formula or not)
         intercept <- FALSE                               ### set to FALSE since formula now controls whether the intercept is included or not
         is.formula <- TRUE                               ### note: code further below actually checks whether intercept is included or not
      }

      ### if yi is an escalc object, try to extract yi and vi (note that moderators must then be specified via the mods argument)

      if (is.element("escalc", class(yi))) {

         if (!is.null(attr(yi, "yi.names"))) { ### if yi.names attributes is available
            yi.name <- attr(yi, "yi.names")[1] ### take the first entry to be the yi variable
         } else {                              ### if not, see if 'yi' is in the object and assume that is the yi variable
            if (!is.element("yi", names(yi)))
               stop("Cannot determine name of the 'yi' variable.")
            yi.name <- "yi"
         }
         if (!is.null(attr(yi, "vi.names"))) { ### if vi.names attributes is available
            vi.name <- attr(yi, "vi.names")[1] ### take the first entry to be the vi variable
         } else {                              ### if not, see if 'vi' is in the object and assume that is the vi variable
            if (!is.element("vi", names(yi)))
               stop("Cannot determine name of the 'vi' variable.")
            vi.name <- "vi"
         }

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

      mf.ni <- mf[[match("ni",  names(mf))]]
      ni <- eval(mf.ni,  data, enclos=sys.frame(sys.parent())) ### NULL if user does not specify this

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

      ### in case user passes a matrix to vi, convert it to a vector

      if (is.matrix(vi))
         vi <- as.vector(vi)

      ### allow easy setting of vi to a single value

      if (length(vi) == 1)
         vi <- rep(vi, k) ### note: k is number of outcomes before subsetting

      ### check length of yi and vi

      if (length(vi) != k)
         stop("Length of 'yi' and 'vi' (or 'sei') vectors are not the same.")

      ### if ni has not been specified but is an attribute of yi, get it

      if (is.null(ni) && !is.null(attr(yi, "ni")))
         ni <- attr(yi, "ni")

      ### check length of yi and ni (only if ni is not NULL)
      ### if there is a mismatch, then ni cannot be trusted, so set it to NULL

      if (!is.null(ni) && length(ni) != k)
         ni <- NULL
         #stop("Length of 'yi' and 'ni' vectors are not the same.")

      ### if ni is now available, add it (back) as an attribute to yi

      if (!is.null(ni))
         attr(yi, "ni") <- ni

      ### note: one or more yi/vi pairs may be NA/NA (also a corresponding ni value may be NA)

      ### if slab has not been specified but is an attribute of yi, get it

      if (is.null(slab)) {

         if (!is.null(attr(yi, "slab")))
            slab <- attr(yi, "slab")

         ### check length of yi and slab (only if slab is not NULL)
         ### if there is a mismatch, then slab cannot be trusted, so set it to NULL

         if (is.null(slab) && length(slab) != k)
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

      if (is.element(measure, c("RR","OR","PETO","RD","AS","PHI","YUQ","YUY","RTET","PBIT","OR2D","OR2DN","OR2DL"))) {

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

         k <- length(ai) ### number of outcomes before subsetting

         if (!is.null(subset)) {
            ai <- ai[subset]
            bi <- bi[subset]
            ci <- ci[subset]
            di <- di[subset]
         }

         dat <- escalc(measure=measure, ai=ai, bi=bi, ci=ci, di=di, add=add, to=to, drop00=drop00, vtype=vtype)

      }

      if (is.element(measure, c("IRR","IRD","IRSD"))) {

         mf.x1i  <- mf[[match("x1i", names(mf))]]
         mf.x2i  <- mf[[match("x2i", names(mf))]]
         mf.t1i  <- mf[[match("t1i", names(mf))]]
         mf.t2i  <- mf[[match("t2i", names(mf))]]
         x1i     <- eval(mf.x1i, data, enclos=sys.frame(sys.parent()))
         x2i     <- eval(mf.x2i, data, enclos=sys.frame(sys.parent()))
         t1i     <- eval(mf.t1i, data, enclos=sys.frame(sys.parent()))
         t2i     <- eval(mf.t2i, data, enclos=sys.frame(sys.parent()))

         k <- length(x1i) ### number of outcomes before subsetting

         if (!is.null(subset)) {
            x1i <- x1i[subset]
            x2i <- x2i[subset]
            t1i <- t1i[subset]
            t2i <- t2i[subset]
         }

         dat <- escalc(measure=measure, x1i=x1i, x2i=x2i, t1i=t1i, t2i=t2i, add=add, to=to, drop00=drop00, vtype=vtype)

      }

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

         k <- length(m1i) ### number of outcomes before subsetting

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

         mf.ri   <- mf[[match("ri", names(mf))]]
         mf.ni   <- mf[[match("ni", names(mf))]]
         ri      <- eval(mf.ri, data, enclos=sys.frame(sys.parent()))
         ni      <- eval(mf.ni, data, enclos=sys.frame(sys.parent()))

         k <- length(ri) ### number of outcomes before subsetting

         if (!is.null(subset)) {
            ri <- ri[subset]
            ni <- ni[subset]
         }

         dat <- escalc(measure=measure, ri=ri, ni=ni, vtype=vtype)

      }

      if (is.element(measure, c("PR","PLN","PLO","PAS","PFT"))) {

         mf.xi   <- mf[[match("xi", names(mf))]]
         mf.mi   <- mf[[match("mi", names(mf))]]
         mf.ni   <- mf[[match("ni", names(mf))]]
         xi      <- eval(mf.xi, data, enclos=sys.frame(sys.parent()))
         mi      <- eval(mf.mi, data, enclos=sys.frame(sys.parent()))
         ni      <- eval(mf.ni, data, enclos=sys.frame(sys.parent()))
         if (is.null(mi)) mi <- ni - xi

         k <- length(xi) ### number of outcomes before subsetting

         if (!is.null(subset)) {
            xi <- xi[subset]
            mi <- mi[subset]
         }

         dat <- escalc(measure=measure, xi=xi, mi=mi, add=add, to=to, vtype=vtype)

      }

      if (is.element(measure, c("IR","IRLN","IRS","IRFT"))) {

         mf.xi   <- mf[[match("xi", names(mf))]]
         mf.ti   <- mf[[match("ti", names(mf))]]
         xi      <- eval(mf.xi, data, enclos=sys.frame(sys.parent()))
         ti      <- eval(mf.ti, data, enclos=sys.frame(sys.parent()))

         k <- length(xi) ### number of outcomes before subsetting

         if (!is.null(subset)) {
            xi <- xi[subset]
            ti <- ti[subset]
         }

         dat <- escalc(measure=measure, xi=xi, ti=ti, add=add, to=to, vtype=vtype)

      }

      if (is.element(measure, c("MN"))) {

         mf.mi   <- mf[[match("mi",  names(mf))]]
         mf.sdi  <- mf[[match("sdi", names(mf))]]
         mf.ni   <- mf[[match("ni",  names(mf))]]
         mi      <- eval(mf.mi,  data, enclos=sys.frame(sys.parent()))
         sdi     <- eval(mf.sdi, data, enclos=sys.frame(sys.parent()))
         ni      <- eval(mf.ni,  data, enclos=sys.frame(sys.parent()))

         k <- length(mi) ### number of outcomes before subsetting

         if (!is.null(subset)) {
            mi  <- mi[subset]
            sdi <- sdi[subset]
            ni  <- ni[subset]
         }

         dat <- escalc(measure=measure, mi=mi, sdi=sdi, ni=ni, vtype=vtype)

      }

      if (is.element(measure, c("MC","SMCC","SMCR","SMCRH","ROMC"))) {

         mf.m1i  <- mf[[match("m1i",  names(mf))]]
         mf.m2i  <- mf[[match("m2i",  names(mf))]]
         mf.sd1i <- mf[[match("sd1i", names(mf))]]
         mf.sd2i <- mf[[match("sd2i", names(mf))]]
         mf.ri   <- mf[[match("ri",   names(mf))]]
         mf.ni   <- mf[[match("ni",   names(mf))]]
         m1i     <- eval(mf.m1i,  data, enclos=sys.frame(sys.parent()))
         m2i     <- eval(mf.m2i,  data, enclos=sys.frame(sys.parent()))
         sd1i    <- eval(mf.sd1i, data, enclos=sys.frame(sys.parent()))
         sd2i    <- eval(mf.sd2i, data, enclos=sys.frame(sys.parent()))
         ri      <- eval(mf.ri,   data, enclos=sys.frame(sys.parent()))
         ni      <- eval(mf.ni,   data, enclos=sys.frame(sys.parent()))

         k <- length(m1i) ### number of outcomes before subsetting

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

         mf.ai   <- mf[[match("ai",  names(mf))]]
         mf.mi   <- mf[[match("mi",  names(mf))]]
         mf.ni   <- mf[[match("ni",  names(mf))]]
         ai      <- eval(mf.ai,  data, enclos=sys.frame(sys.parent()))
         mi      <- eval(mf.mi,  data, enclos=sys.frame(sys.parent()))
         ni      <- eval(mf.ni,  data, enclos=sys.frame(sys.parent()))

         k <- length(ai) ### number of outcomes before subsetting

         if (!is.null(subset)) {
            ai <- ai[subset]
            mi <- mi[subset]
            ni <- ni[subset]
         }

         dat <- escalc(measure=measure, ai=ai, mi=mi, ni=ni, vtype=vtype)

      }

      if (is.element(measure, "GEN"))
         stop("Specify the desired outcome measure via the 'measure' argument.")

      ### note: these values are already subsetted

      yi <- dat$yi         ### one or more yi/vi pairs may be NA/NA
      vi <- dat$vi         ### one or more yi/vi pairs may be NA/NA
      ni <- attr(yi, "ni") ### unadjusted total sample sizes (ni.u in escalc)

   }

   #########################################################################

   ### allow easy setting of weights to a single value

   if (length(weights) == 1)
      weights <- rep(weights, k) ### note: k is number of outcomes before subsetting

   ### check length of yi and weights (only if weights is not NULL)

   if (!is.null(weights) && (length(weights) != k))
      stop("Length of 'yi' and 'weights' vectors are not the same.")

   ### subsetting of weights

   if (!is.null(subset))
      weights <- weights[subset]

   #########################################################################

   if (very.verbose)
      message("Creating model matrix ...")

   ### convert mods formula to X matrix and set intercept equal to FALSE
   ### skipped if formula has already been specified via yi argument, since mods is then no longer a formula

   if (class(mods) == "formula") {
      options(na.action = "na.pass")        ### set na.action to na.pass, so that NAs are not filtered out (we'll do that later)
      mods <- model.matrix(mods, data=data) ### extract model matrix
      attr(mods, "assign") <- NULL          ### strip assign attribute (not needed at the moment)
      options(na.action = na.act)           ### set na.action back to na.act
      intercept <- FALSE                    ### set to FALSE since formula now controls whether the intercept is included or not
      is.formula <- TRUE                    ### note: code further below actually checks whether intercept is included or not
   }

   ### turn a row vector for mods into a column vector

   if (is.vector(mods))
      mods <- cbind(mods)

   ### turn a mods data frame into a matrix

   if (is.data.frame(mods))
      mods <- as.matrix(mods)

   ### check if model matrix contains character variables

   if (is.character(mods))
      stop("Model matrix contains character variables.")

   ### check if mods matrix has the right number of rows

   if (!is.null(mods) && (nrow(mods) != k))
      stop("Number of rows of the model matrix does not match length of 'yi' argument.")

   ### in case scale is a formula, get model matrix for it

   if (class(scale) == "formula") {
      options(na.action = "na.pass")
      Z <- model.matrix(scale, data=data)
      attr(Z, "assign") <- NULL
      options(na.action = na.act)
      model <- "rma.tau2"
      if (nrow(Z) != k)
         stop("Number of rows of the model matrix for tau2 does not match length of 'yi' argument.")
   } else {
      Z <- NULL
      model <- "rma.uni"
   }

   ### study ids (1:k sequence before subsetting)

   ids <- seq_len(k)

   ### generate study labels if none are specified (or none have been found in yi)

   if (very.verbose)
      message("Generating/extracting study labels ...")

   if (is.null(slab)) {

      slab.null <- TRUE
      slab      <- ids

   } else {

      if (anyNA(slab))
         stop("NAs in study labels.")

      if (length(slab) != k)
         stop("Study labels not of same length as data.")

      slab.null <- FALSE

   }

   ### if a subset of studies is specified

   if (!is.null(subset)) {

      if (very.verbose)
         message("Subsetting ...")

      mods <- mods[subset,,drop=FALSE]
      slab <- slab[subset]
      ids  <- ids[subset]
      Z    <- Z[subset,,drop=FALSE]

   }

   ### check if study labels are unique; if not, make them unique

   if (anyDuplicated(slab))
      slab <- make.unique(as.character(slab)) ### make.unique() only works with character vectors

   ### add slab attribute back

   attr(yi, "slab") <- slab

   ### number of outcomes after subsetting

   k <- length(yi)

   ### check for non-positive sampling variances (and set negative values to 0)

   if (any(vi <= 0, na.rm=TRUE)) {
      allvipos <- FALSE
      warning("There are outcomes with non-positive sampling variances.")
      vi.neg <- vi < 0
      if (any(vi.neg, na.rm=TRUE)) {
         vi[vi.neg] <- 0
         warning("Negative sampling variances constrained to zero.")
      }
   } else {
      allvipos <- TRUE
   }

   ### check for (and correct?) negative/infinite weights

   if (any(weights < 0, na.rm=TRUE))
      stop("Negative weights not allowed.")

   if (any(is.infinite(weights)))
      stop("Infinite weights not allowed.")

   ### save full data (including potential NAs in yi/vi and/or mods)

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
   #Z.f      <- Z ### don't need this at the moment

   k.f <- k ### total number of observed outcomes including all NAs (on yi/vi and/or mods)

   ### check for NAs and act accordingly

   YVXZW.na <- is.na(yi) | is.na(vi) | if (is.null(mods)) FALSE else apply(is.na(mods), 1, any) | if (is.null(Z)) FALSE else apply(is.na(Z), 1, any) | if (is.null(weights)) FALSE else is.na(weights)

   if (any(YVXZW.na)) {

      if (very.verbose)
         message("Handling NAs ...")

      not.na <- !YVXZW.na

      if (na.act == "na.omit" || na.act == "na.exclude" || na.act == "na.pass") {

         yi      <- yi[not.na]
         vi      <- vi[not.na]
         weights <- weights[not.na]
         ni      <- ni[not.na]
         mods    <- mods[not.na,,drop=FALSE]
         Z       <- Z[not.na,,drop=FALSE]
         k       <- length(yi)
         warning("Studies with NAs omitted from model fitting.")

         attr(yi, "measure") <- measure ### add measure attribute back
         attr(yi, "ni")      <- ni      ### add ni attribute back

         ### note: slab is always of the same length as the full yi vector (after subsetting), so missings are not removed and slab is not added back to yi

      }

      if (na.act == "na.fail")
         stop("Missing values in data.")

   } else {
      not.na <- rep(TRUE, k)
   }

   ### at least one study left?

   if (k < 1)
      stop("Processing terminated since k = 0.")

   ### check if k is equal to 1

   if (k == 1) {
      method <- "FE"
      knha   <- FALSE
   }

   ### make sure that there is at least one column in X

   if (is.null(mods) && !intercept) {
      warning("Must either include an intercept and/or moderators in model.\n  Coerced intercept into the model.")
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

   ### check whether intercept is included and if yes, move it to the first column (NAs already removed, so na.rm=TRUE for any() not necessary)
   ### note that a column may turn into an 'intercept' due to the removal of NAs (e.g., X=cbind(c(1,1,1,1), c(1,2,1,1), c(5,NA,2,8))), so does
   ### this have consequences for how X.f gets used in other functions?

   is.int <- apply(X, 2, .is.int.func)
   if (any(is.int)) {
      int.incl <- TRUE
      int.indx <- which(is.int, arr.ind=TRUE)
      X        <- cbind(intrcpt=1,   X[,-int.indx, drop=FALSE]) ### this removes any duplicate intercepts
      X.f      <- cbind(intrcpt=1, X.f[,-int.indx, drop=FALSE]) ### this removes any duplicate intercepts
      if (is.formula)
         intercept <- TRUE ### set intercept appropriately so that the predict() function works
   } else {
      int.incl <- FALSE
   }

   ### drop redundant predictors
   ### note: need to save coef.na for functions that modify the data/model and then refit the model (regtest() and the
   ### various function that leave out an observation); so we can check if there are redundant/dropped predictors then

   tmp <- lm(yi ~ X - 1)
   coef.na <- is.na(coef(tmp))
   if (any(coef.na)) {
      warning("Redundant predictors dropped from the model.")
      X   <- X[,!coef.na,drop=FALSE]
      X.f <- X.f[,!coef.na,drop=FALSE]
   }

   p <- NCOL(X) ### number of columns in X (including the intercept if it is included)

   ### check whether this is an intercept-only model

   if ((p == 1L) && (all(sapply(X, identical, 1)))) {
      int.only <- TRUE
   } else {
      int.only <- FALSE
   }

   ### check if there are too many parameters for given k

   if (method == "FE") {                        ### have to estimate p parms
      if (p > k)
         stop("Number of parameters to be estimated is larger than the number of observations.")
   } else {
      if (is.numeric(tau2)) {                   ### have to estimate p parms (tau2 is fixed at value specified)
         if (p > k)
            stop("Number of parameters to be estimated is larger than the number of observations.")
      } else {
         if ((p+1) > k)                         ### have to estimate p+1 parms
            stop("Number of parameters to be estimated is larger than the number of observations.")
      }
   }

   ### set/check 'btt' argument

   btt <- .set.btt(btt, p, int.incl)
   m <- length(btt) ### number of betas to test (m = p if all betas are tested)

   #########################################################################

   ### set default control parameters

   con <- list(verbose = FALSE,
               tau2.init = NULL,     # initial value for iterative estimators (ML, REML, EB, SJ, SJIT, DLIT) or scale parameters
               tau2.min = 0,         # lower bound for tau^2 value
               tau2.max = 100,       # upper bound for tau^2 value (only needed for PM estimator; and passed down for tau^2 CI obtained with confint())
               threshold = 10^-5,    # convergence threshold (for ML, REML, EB, SJIT, DLIT; also used for PM)
               maxiter = 100,        # maximum number of iterations (for ML, REML, EB, SJIT, DLIT)
               stepadj = 1,          # step size adjustment for Fisher scoring algorithm (for ML, REML, EB)
               REMLf = TRUE,         # should |X'X| term be included in the REML log likelihood?
               tol = 1e-07,          # tolerance for checking whether X'X is of full rank
               optimizer = "nlminb", # optimizer to use ("optim", "nlminb", "uobyqa", "newuoa", "bobyqa", "nloptr") for location-scale model
               optmethod = "BFGS")   # argument 'method' for optim() ("Nelder-Mead" and "BFGS" are sensible options)

   ### replace defaults with any user-defined values

   con.pos <- pmatch(names(control), names(con))
   con[c(na.omit(con.pos))] <- control[!is.na(con.pos)]

   if (verbose)
      con$verbose <- verbose

   verbose <- con$verbose

   ### convergence indicator and change variable (for iterative estimators)

   conv <- 1
   change <- con$threshold + 1

   #########################################################################

   ### check whether model matrix is of full rank

   if (any(eigen(crossprod(X), symmetric=TRUE, only.values=TRUE)$values <= con$tol))
      stop("Model matrix not of full rank. Cannot fit model.")

   iter <- 0 ### iterations counter

   ### set some defaults

   se.tau2 <- I2 <- H2 <- QE <- QEp <- NA

   s2w <- 1

   alpha <- ifelse(level > 1, (100-level)/100, 1-level)

   Y <- as.matrix(yi)

   #########################################################################

   ###### heterogeneity estimation

   if (model == "rma.tau2") {

      ### code for location-scale model

      if (method != "ML" && method != "REML")
         stop("Location-scale model can only be fitted with ML or REML estimation.")

      if (!weighted)
         stop("Cannot use weighted estimation to fit location-scale model.")

      if (!is.null(weights))
         stop("Cannot use user-defined weights for location-scale model.")

      if (any(eigen(crossprod(Z), symmetric=TRUE, only.values=TRUE)$values <= con$tol))
         stop("Model matrix for scale part of the model not of full rank. Cannot fit model.")

      optimizer  <- match.arg(con$optimizer, c("optim","nlminb","uobyqa","newuoa","bobyqa","nloptr"))
      optmethod  <- con$optmethod
      optcontrol <- control[is.na(con.pos)] ### get arguments that are control arguments for optimizer

      if (length(optcontrol) == 0)
         optcontrol <- list()

      ### set NLOPT_LN_BOBYQA as the default algorithm for nloptr optimizer
      ### and by default use a relative convergence criterion of 1e-8 on the function value

      if (optimizer=="nloptr" && !is.element("algorithm", names(optcontrol)))
         optcontrol$algorithm <- "NLOPT_LN_BOBYQA"

      if (optimizer=="nloptr" && !is.element("ftol_rel", names(optcontrol)))
         optcontrol$ftol_rel <- 1e-8

      reml <- ifelse(method=="REML", TRUE, FALSE)

      if (is.element(optimizer, c("uobyqa","newuoa","bobyqa"))) {
         if (!requireNamespace("minqa", quietly=TRUE))
            stop("Please install the 'minqa' package to use this optimizer.")
      }

      if (optimizer == "nloptr") {
         if (!requireNamespace("nloptr", quietly=TRUE))
            stop("Please install the 'nloptr' package to use this optimizer.")
      }

      if (!requireNamespace("numDeriv", quietly=TRUE))
         stop("Please install the 'numDeriv' package to fit a location-scale model.")

      p.tau2 <- NCOL(Z)

      ### check if length of tau2.init matches number of parameters

      if (!is.null(con$tau2.init)) {
         if (length(con$tau2.init) != p.tau2)
            stop(paste("Length of 'tau2.init' argument (", length(con$tau2.init), ") does not match actual number of parameters (", p.tau2, ").", sep=""))
      } else {
         con$tau2.init <- rep(0.0001, p.tau2) ### need a better way to set default initial value(s)
      }

      ### obtain initial values for b (only need this when optimizing over b.fe and b.tau2 jointly)

      #wi <- 1/vi
      #W  <- diag(wi, nrow=k, ncol=k)
      #stXWX <- .invcalc(X=X, W=W, k=k)
      #b     <- stXWX %*% crossprod(X,W) %*% Y

      if (optimizer=="optim")
         par.arg <- "par"
      if (optimizer=="nlminb")
         par.arg <- "start"
      if (is.element(optimizer, c("uobyqa","newuoa","bobyqa"))) {
         par.arg <- "par"
         optimizer <- paste0("minqa::", optimizer) ### need to use this since loading nloptr masks bobyqa() and newuoa() functions
      }
      if (optimizer=="nloptr") {
         par.arg <- "x0"
         optimizer <- paste0("nloptr::nloptr") ### need to use this due to requireNamespace()
      }

      optcall <- paste(optimizer, "(", par.arg, "=con$tau2.init, .ll.rma.tau2, ",
                                                    ifelse(optimizer=="optim", "method=optmethod, ", ""),
                                                    "yi=yi, vi=vi, X=X, Z=Z, reml=reml,
                                                    k=k, p=p,
                                                    verbose=verbose, digits=digits, REMLf=con$REMLf, ",
                                                    ifelse(optimizer=="nloptr::nloptr", "opts=optcontrol)", "control=optcontrol)"), sep="")
      #return(optcall)
      opt.res <- try(eval(parse(text=optcall)), silent=!verbose)
      #return(opt.res)

      if (inherits(opt.res, "try-error"))
         stop("Error during optimization for location-scale model.")

      if (is.element(optimizer, c("optim","nlminb")) && opt.res$convergence != 0)
         stop(paste0("Optimizer (", optimizer, ") did not achieve convergence (convergence = ", opt.res$convergence, ")."))

      if (is.element(optimizer, c("minqa::uobyqa","minqa::newuoa","minqa::bobyqa")) && opt.res$ierr != 0)
         stop(paste0("Optimizer (", optimizer, ") did not achieve convergence (ierr = ", opt.res$ierr, ")."))

      if (optimizer=="nloptr::nloptr" && !(opt.res$status >= 1 && opt.res$status <= 4))
         stop(paste0("Optimizer (", optimizer, ") did not achieve convergence (status = ", opt.res$status, ")."))

      if (verbose) {
         cat("\n")
         print(opt.res)
      }

      ### copy 'solution' to 'par' so code below works (all other optimizers use 'par' for the solution)

      if (optimizer=="nloptr::nloptr")
         opt.res$par <- opt.res$solution

      ### try to compute vcov matrix for scale parameters

      vb.tau2 <- matrix(NA, nrow=p.tau2, ncol=p.tau2)
      se.tau2 <- rep(NA, p.tau2)

      h <- try(numDeriv::hessian(.ll.rma.tau2, x=opt.res$par, yi=yi, vi=vi, X=X, Z=Z, reml=reml, k=k, p=p, verbose=FALSE, digits=digits, REMLf=con$REMLf), silent=TRUE)

      if (!inherits(h, "try-error")) {
         chol.h <- try(chol(h), silent=TRUE)
         if (!inherits(chol.h, "try-error")) {
            vb.tau2 <- chol2inv(chol.h)
            se.tau2 <- sqrt(diag(vb.tau2))
         }
      }

      ### scale parameter estimates

      b.tau2 <- cbind(opt.res$par)

      colnames(Z)[grep("(Intercept)", colnames(Z))] <- "intrcpt"
      rownames(b.tau2) <- rownames(vb.tau2) <- colnames(vb.tau2) <- colnames(Z)

      crit <- qnorm(alpha/2, lower.tail=FALSE)

      names(se.tau2) <- NULL
      zval.tau2  <- c(b.tau2/se.tau2)
      pval.tau2  <- 2*pnorm(abs(zval.tau2), lower.tail=FALSE)
      ci.lb.tau2 <- c(b.tau2 - crit * se.tau2)
      ci.ub.tau2 <- c(b.tau2 + crit * se.tau2)

      tau2 <- exp(as.vector(Z %*% b.tau2))
      tau2.fix <- FALSE

   } else {

      ### if user has NOT specified a formula for tau2

      if (is.numeric(tau2)) { ### if user has fixed the tau2 value
         tau2.fix <- TRUE
         tau2.val <- tau2
      } else {
         tau2.fix <- FALSE
         tau2.val <- NA
      }

      if (very.verbose && !tau2.fix)
         message("Estimating tau^2 value ...")

      ### Hunter & Schmidt (HS) estimator

      if (method == "HS") {

         if (!allvipos)
            stop("HS estimator cannot be used with non-positive sampling variances.")

         wi    <- 1/vi
         W     <- diag(wi, nrow=k, ncol=k)
         stXWX <- .invcalc(X=X, W=W, k=k)
         P     <- W - W %*% X %*% stXWX %*% crossprod(X,W)
         RSS   <- crossprod(Y,P) %*% Y
         tau2  <- ifelse(tau2.fix, tau2.val, (RSS-k)/sum(wi))

      }

      ### Hedges (HE) estimator (or initial value for ML, REML, EB)

      if (is.element(method, c("HE","ML","REML","EB"))) {

         stXX  <- .invcalc(X=X, W=diag(k), k=k)
         P     <- diag(k) - X %*% tcrossprod(stXX,X)
         RSS   <- crossprod(Y,P) %*% Y
         V     <- diag(vi, nrow=k, ncol=k)
         PV    <- P %*% V ### careful: is not symmetric
         trPV  <- .tr(PV) ### since PV needs to be computed anyway, can use .tr()
         tau2  <- ifelse(tau2.fix, tau2.val, (RSS - trPV) / (k-p))

      }

      ### DerSimonian-Laird (DL) estimator

      if (method == "DL") {

         if (!allvipos)
            stop("DL estimator cannot be used with non-positive sampling variances.")

         wi    <- 1/vi
         W     <- diag(wi, nrow=k, ncol=k)
         stXWX <- .invcalc(X=X, W=W, k=k)
         P     <- W - W %*% X %*% stXWX %*% crossprod(X,W)
         RSS   <- crossprod(Y,P) %*% Y
         trP   <- .tr(P)
         tau2  <- ifelse(tau2.fix, tau2.val, (RSS - (k-p)) / trP)

      }

      if (method == "GENQ") {

         #if (!allvipos)
         #   stop("GENQ estimator cannot be used with non-positive sampling variances.")

         if (is.null(weights))
            stop("Must specify 'weights' when method='GENQ'.")

         A     <- diag(weights, nrow=k, ncol=k)
         stXAX <- .invcalc(X=X, W=A, k=k)
         P     <- A - A %*% X %*% stXAX %*% t(X) %*% A
         RSS   <- crossprod(Y,P) %*% Y
         V     <- diag(vi, nrow=k, ncol=k)
         PV    <- P %*% V ### careful: is not symmetric
         trP   <- .tr(P)
         trPV  <- .tr(PV)
         tau2  <- ifelse(tau2.fix, tau2.val, (RSS - trPV) / trP)

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
         PV    <- P %*% V ### careful: is not symmetric
         tau2  <- ifelse(tau2.fix, tau2.val, tau2.0 * RSS / (k-p))

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
               cat("Iteration", iter, "\ttau^2 =", formatC(tau2, format="f", digits=digits), "\n")

            iter     <- iter + 1
            tau2.old <- tau2
            wi       <- 1/(vi + tau2)
            if (any(tau2 + vi < 0))
               stop("Some marginal variances are negative.")
            if (any(is.infinite(wi)))
               stop("Division by zero when computing the inverse variance weights.")
            W        <- diag(wi, nrow=k, ncol=k)
            stXWX    <- .invcalc(X=X, W=W, k=k)
            P        <- W - W %*% X %*% stXWX %*% crossprod(X,W)
            RSS      <- crossprod(Y,P) %*% Y
            trP      <- .tr(P)
            tau2     <- ifelse(tau2.fix, tau2.val, (RSS - (k-p)) / trP)
            tau2[tau2 < con$tau2.min] <- con$tau2.min
            change   <- abs(tau2.old - tau2)

            if (iter > con$maxiter) {
               conv <- 0
               break
            }

         }

         if (conv == 0L)
            stop("Algorithm did not converge.")

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
               cat("Iteration", iter, "\ttau^2 =", formatC(tau2, format="f", digits=digits), "\n")

            iter     <- iter + 1
            tau2.old <- tau2

            wi     <- 1/(vi + tau2)
            W      <- diag(wi, nrow=k, ncol=k)
            stXWX  <- .invcalc(X=X, W=W, k=k)
            P      <- W - W %*% X %*% stXWX %*% crossprod(X,W)
            RSS    <- crossprod(Y,P) %*% Y
            V      <- diag(vi, nrow=k, ncol=k)
            PV     <- P %*% V ### careful: is not symmetric
            tau2   <- ifelse(tau2.fix, tau2.val, tau2 * RSS / (k-p))
            change <- abs(tau2.old - tau2)

            if (iter > con$maxiter) {
               conv <- 0
               break
            }

         }

         if (conv == 0L)
            stop("Algorithm did not converge.")

      }

      ### Paule-Mandel (PM) estimator

      if (method == "PM") {

         if (!allvipos)
            stop("PM estimator cannot be used with non-positive sampling variances.")

         if (.QE.func(con$tau2.min, Y=Y, vi=vi, X=X, k=k, objective=k-p) < con$tau2.min) {

            tau2 <- con$tau2.min

         } else {

            tau2 <- ifelse(tau2.fix, tau2.val, try(uniroot(.QE.func, interval=c(con$tau2.min, con$tau2.max), tol=con$threshold, maxiter=con$maxiter, Y=Y, vi=vi, X=X, k=k, objective=k-p, verbose=verbose, digits=digits, extendInt="upX")$root, silent=TRUE))

            if (!is.numeric(tau2))
               stop("Error in iterative search for tau2. Try increasing 'tau2.max' or switch to another 'method'.")

         }

         wi    <- 1/(vi + tau2)
         W     <- diag(wi, nrow=k, ncol=k)
         stXWX <- .invcalc(X=X, W=W, k=k)
         P     <- W - W %*% X %*% stXWX %*% crossprod(X,W) ### needed for se.tau2 computation below

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
               cat("Iteration", iter, "\ttau^2 =", formatC(tau2, format="f", digits=digits), "\n")

            iter     <- iter + 1
            tau2.old <- tau2
            wi       <- 1/(vi + tau2)
            if (any(tau2 + vi < 0))
               stop("Some marginal variances are negative.")
            if (any(is.infinite(wi)))
               stop("Division by zero when computing the inverse variance weights.")
            W        <- diag(wi, nrow=k, ncol=k)
            stXWX    <- .invcalc(X=X, W=W, k=k)
            P        <- W - W %*% X %*% stXWX %*% crossprod(X,W)

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

            while (tau2 + adj < con$tau2.min) { ### use step-halving if necessary
               #print(adj)
               adj <- adj / 2
            }

            tau2   <- ifelse(tau2.fix, tau2.val, tau2 + adj)
            change <- abs(tau2.old - tau2)

            if (iter > con$maxiter) {
               conv <- 0
               break
            }

         }

         if (conv == 0L)
            stop("Fisher scoring algorithm did not converge. See 'help(rma)' for possible remedies.")

      }

      ### make sure that tau2 is >= con$tau2.min

      tau2 <- max(con$tau2.min, c(tau2))

      ### check if any marginal variances negative (only possible if user has changed tau2.min)

      if (any(tau2 + vi < 0))
         stop("Some marginal variances are negative.")

      ### verbose output upon convergence for ML/REML/EB estimators

      if (verbose && is.element(method, c("ML","REML","EB"))) {
         cat("Iteration", iter, "\ttau^2 =", formatC(tau2, format="f", digits=digits), "\n")
         cat("Fisher scoring algorithm converged after", iter, "iterations.\n")
      }

      ### standard error of the tau^2 estimators (also when the user has fixed/specified a tau^2 value)
      ### see notes.pdf and note: .tr(P%*%P) = sum(P*t(P)) = sum(P*P) (since P is symmetric)

      if (method == "HS") {
         se.tau2 <- sqrt(1/sum(wi)^2 * (2*(k-p) + 4*max(tau2,0)*.tr(P) + 2*max(tau2,0)^2*sum(P*P)))
      }
      if (method == "HE") {
         se.tau2 <- sqrt(1/(k-p)^2 * (2*sum(PV*t(PV)) + 4*max(tau2,0)*trPV + 2*max(tau2,0)^2*(k-p)))
      }
      if (method == "DL" || method == "DLIT") {
         se.tau2 <- sqrt(1/trP^2 * (2*(k-p) + 4*max(tau2,0)*trP + 2*max(tau2,0)^2*sum(P*P)))
      }
      if (method == "GENQ") {
         se.tau2 <- sqrt(1/trP^2 * (2*sum(PV*t(PV)) + 4*max(tau2,0)*sum(PV*P) + 2*max(tau2,0)^2*sum(P*P)))
      }
      if (method == "SJ") {
         se.tau2 <- sqrt(tau2.0^2/(k-p)^2 * (2*sum(PV*t(PV)) + 4*max(tau2,0)*sum(PV*P) + 2*max(tau2,0)^2*sum(P*P)))
      }
      if (method == "ML") {
         se.tau2 <- sqrt(2/sum(wi^2))
      }
      if (method == "REML") {
         se.tau2 <- sqrt(2/sum(P*P))
      }
      if (method == "EB" || method == "PM" || method == "SJIT") {
         V  <- diag(vi, nrow=k, ncol=k)
         PV <- P %*% V ### careful: is not symmetric
         #se.tau2 <- sqrt((k/(k-p))^2 / sum(wi)^2 * (2*sum(PV*t(PV)) + 4*max(tau2,0)*sum(PV*P) + 2*max(tau2,0)^2*sum(P*P)))
         se.tau2 <- sqrt(2*k^2/(k-p) / sum(wi)^2) ### these two equations are actually identical, but this one is much simpler
      }

   }

   ### fixed-effects model (note: set tau2 to zero even when tau2 is specified)

   if (method == "FE")
      tau2 <- 0

   #########################################################################

   ###### model fitting, test statistics, and confidence intervals

   if (very.verbose)
      message("Model fitting ...")

   wi <- 1/(vi + tau2)
   W  <- diag(wi, nrow=k, ncol=k)
   M  <- diag(vi + tau2, nrow=k, ncol=k)

   if (weighted) {

      #########################
      ### weighted analysis ###
      #########################

      ### fit model with weighted estimation

      if (is.null(weights)) {

         ### if no weights are specified, use default inverse variance weights, that is, 1/vi or 1/(vi + tau2)

         ### if any vi = 0 and tau^2 is estimated to be 0 (or is set to 0 for a FE model), then get Inf for wi

         if (any(is.infinite(wi)))
            stop("Division by zero when computing the inverse variance weights.")

         stXWX <- .invcalc(X=X, W=W, k=k)
         b     <- stXWX %*% crossprod(X,W) %*% Y
         vb    <- stXWX
         RSS.f <- sum(wi*(yi - X %*% b)^2)
         #P     <- W - W %*% X %*% stXWX %*% crossprod(X,W)
         #RSS.f <- crossprod(Y,P) %*% Y

      } else {

         ### if weights are specified, use them

         A     <- diag(weights, nrow=k, ncol=k)
         stXAX <- .invcalc(X=X, W=A, k=k)
         b     <- stXAX %*% crossprod(X,A) %*% Y
         vb    <- stXAX %*% t(X) %*% A %*% M %*% A %*% X %*% stXAX
         RSS.f <- sum(wi*(yi - X %*% b)^2)
         #P     <- W - W %*% X %*% stXAX %*% t(X) %*% A - A %*% X %*% stXAX %*% t(X) %*% W + A %*% X %*% stXAX %*% t(X) %*% W %*% X %*% stXAX %*% t(X) %*% A
         #RSS.f <- crossprod(Y,P) %*% Y

      }

      #return(list(b=b, vb=vb, se=sqrt(diag(vb)), RSS.f=RSS.f))

      ### calculate scaling factor for Knapp & Hartung method
      ### note: catch cases where RSS.f is extremely small, which is probably due to all yi being equal
      ### then set s2w to 0 (to avoid the strange looking output we would obtain if we don't do this)

      if ((is.logical(knha) && knha) || is.character(knha)) {

         if (RSS.f <= .Machine$double.eps) {
            s2w <- 0
         } else {
            s2w <- RSS.f / (k-p)
         }

      }

   } else {

      ###########################
      ### unweighted analysis ###
      ###########################

      ### fit model with unweighted estimation
      ### note: 1) if user has specified weights, they are ignored
      ###       2) but if method="GENQ", they were used to estimate tau^2

      stXX <- .invcalc(X=X, W=diag(k), k=k)
      b    <- stXX %*% crossprod(X,Y)
      vb   <- tcrossprod(stXX,X) %*% M %*% X %*% stXX
      #P    <- W - W %*% X %*% tcrossprod(stXX,X) - X %*% stXX %*% crossprod(X,W) + X %*% stXX %*% crossprod(X,W) %*% X %*% tcrossprod(stXX,X)
      #RSS.f <- crossprod(Y,P) %*% Y
      RSS.f <- sum(wi*(yi - X %*% b)^2)

      ### calculate scaling factor for Knapp & Hartung method

      if ((is.logical(knha) && knha) || is.character(knha)) {

         stXWX     <- .invcalc(X=X, W=W, k=k)
         b.knha    <- stXWX %*% crossprod(X,W) %*% Y
         #P        <- W - W %*% X %*% stXWX %*% crossprod(X,W)
         #RSS.knha <- c(crossprod(Y,P) %*% Y)
         RSS.knha  <- sum(wi*(yi - X %*% b.knha)^2)

         if (RSS.knha <= .Machine$double.eps) {
            s2w <- 0
         } else {
            s2w <- RSS.knha / (k-p)
         }

      }

   }

   ### the Knapp & Hartung method as described in the literature is for random/mixed-effects models

   if (method == "FE" && ((is.logical(knha) && knha) || is.character(knha)))
      warning("Knapp & Hartung method is not meant to be used in the context of FE models.")

   ### Knapp & Hartung method with ad-hoc correction so that the scale factor is always >= 1

   if (knha == "adhoc")
      s2w[s2w < 1] <- 1

   ### to force use of a t-distribution with no adjustment to the SEs (for testing purposes)

   if (knha == "tdist")
      s2w <- 1

   ### for Knapp & Hartung method, apply scaling to vb

   vb <- s2w * vb

   ### QM calculation

   QM <- try(c(t(b)[btt] %*% chol2inv(chol(vb[btt,btt])) %*% b[btt]), silent=TRUE)

   if (inherits(QM, "try-error"))
      QM <- NA

   rownames(b) <- rownames(vb) <- colnames(vb) <- colnames(X) ### may not be needed, but what's the harm?

   se <- sqrt(diag(vb))
   names(se) <- NULL
   zval <- c(b/se)

   if ((is.logical(knha) && knha) || is.character(knha)) {
      dfs  <- k-p
      QM   <- QM / m
      if (dfs > 0) {
         QMp  <- pf(QM, df1=m, df2=dfs, lower.tail=FALSE)
         pval <- 2*pt(abs(zval), df=dfs, lower.tail=FALSE)
         crit <- qt(alpha/2, df=dfs, lower.tail=FALSE)
      } else {
         QMp  <- NaN
         pval <- NaN
         crit <- NaN
      }
   } else {
      dfs  <- NA
      QMp  <- pchisq(QM, df=m, lower.tail=FALSE)
      pval <- 2*pnorm(abs(zval), lower.tail=FALSE)
      crit <- qnorm(alpha/2, lower.tail=FALSE)
   }

   ci.lb <- c(b - crit * se)
   ci.ub <- c(b + crit * se)

   #########################################################################

   ### heterogeneity test (Wald-type test of the extra coefficients in the saturated model)

   if (very.verbose)
      message("Heterogeneity testing ...")

   if (allvipos) {

      ### heterogeneity test (always uses inverse variance method)
      ### note: this is unaffected by the 'weighted' argument, since under H0, the same parameters are estimated
      ### and weighted estimation provides the most efficient estimates; therefore, also any arbitrary weights
      ### specified by the user are not relevant here (different from what the metan command in Stata does!)
      ### see also: Chen, Z., Ng, H. K. T., & Nadarajah, S. (2014). A note on Cochran test for homogeneity in
      ### one-way ANOVA and meta-analysis. Statistical Papers, 55(2), 301-310. This shows that the weights used
      ### are not relevant.

      wi    <- 1/vi
      W.FE  <- diag(wi, nrow=k, ncol=k) ### care: ll.REML below involves W, so cannot overwrite W
      stXWX <- .invcalc(X=X, W=W.FE, k=k)
      P     <- W.FE - W.FE %*% X %*% stXWX %*% crossprod(X,W.FE) ### need P below for calculation of I^2
      QE    <- max(0, c(crossprod(Y,P) %*% Y))
      #b.FE  <- stXWX %*% crossprod(X,W.FE) %*% Y
      #QE    <- max(0, sum(wi*(yi - X %*% b.FE)^2))
      QEp   <- ifelse(k-p >= 1, pchisq(QE, df=k-p, lower.tail=FALSE), 1)

      ### calculation of I^2 and H^2 in the RE/ME model

      if (method != "FE") {
         #vi.avg <- (k-1) / (sum(wi) - sum(wi^2)/sum(wi)) ### this only applies to the RE model
         #vi.avg <- 1/mean(wi) ### harmonic mean of vi's (see Takkouche et al., 1999)
         vi.avg  <- (k-p) / .tr(P)
         I2      <- 100 * mean(tau2) / (vi.avg + mean(tau2)) ### must use mean(tau2) in case tau2 is vector from location-scale model
         H2      <- mean(tau2) / vi.avg + 1                  ### must use mean(tau2) in case tau2 is vector from location-scale model
      }

   } else {

      warning(paste0("Cannot compute ", ifelse(int.only, "Q", "QE"), "-test, I^2, or H^2 with non-positive sampling variances."))

   }

   #########################################################################

   ### compute pseudo R^2 statistic for mixed-effects models with an intercept

   if (!int.only && int.incl && method != "FE" && model == "rma.uni") {

      if (very.verbose) {
         message("Fitting RE model for R^2 computation ...")
         res.RE <- try(rma(yi, vi, weights=weights, method=method, weighted=weighted, knha=knha, verbose=ifelse(verbose, TRUE, FALSE), control=con, digits=digits), silent=FALSE)
      } else {
         res.RE <- try(suppressWarnings(rma(yi, vi, weights=weights, method=method, weighted=weighted, knha=knha, verbose=ifelse(verbose, TRUE, FALSE), control=con, digits=digits)), silent=FALSE)
      }

      if (!inherits(res.RE, "try-error")) {

         tau2.RE <- res.RE$tau2

         if (identical(tau2.RE,0)) {
            R2 <- NA
         } else {
            R2 <- round(max(0, 100 * (tau2.RE - tau2) / tau2.RE), 2)
         }

      } else {

         R2 <- NA

      }

   } else {

      R2 <- NULL

   }

   #########################################################################

   ###### fit statistics

   if (very.verbose)
      message("Computing fit statistics and log likelihood ...")

   ### note: tau2 is not counted as a parameter when it was fixed by the user
   parms <- p + ifelse(model == "rma.uni", ifelse(method == "FE" || tau2.fix, 0, 1), p.tau2)

   ll.ML    <- -1/2 * (k)   * log(2*base::pi)                                                                                 - 1/2 * sum(log(vi + tau2))                                                                   - 1/2 * RSS.f
   ll.REML  <- -1/2 * (k-p) * log(2*base::pi) + ifelse(con$REMLf, 1/2 * determinant(crossprod(X), logarithm=TRUE)$modulus, 0) - 1/2 * sum(log(vi + tau2)) - 1/2 * determinant(crossprod(X,W) %*% X, logarithm=TRUE)$modulus - 1/2 * RSS.f
   #ll.REML <- -1/2 * (k-p) * log(2*base::pi) +                   1/2 * log(det(crossprod(X)))                                - 1/2 * sum(log(vi + tau2)) - 1/2 * log(det(crossprod(X,W) %*% X))                            - 1/2 * RSS.f
   #ll.REML <- -1/2 * (k-p) * log(2*base::pi)                                                                                 - 1/2 * sum(log(vi + tau2)) - 1/2 * log(det(crossprod(X,W) %*% X))                            - 1/2 * RSS.f

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

   if (very.verbose)
      message("Preparing output ...")

   p.eff <- p
   k.eff <- k

   res <- list(b=b, se=se, zval=zval, pval=pval, ci.lb=ci.lb, ci.ub=ci.ub, vb=vb,
               tau2=tau2, se.tau2=se.tau2, tau2.fix=tau2.fix,
               k=k, k.f=k.f, k.eff=k.eff, p=p, p.eff=p.eff, parms=parms, m=m,
               QE=QE, QEp=QEp, QM=QM, QMp=QMp, I2=I2, H2=H2, R2=R2,
               int.only=int.only, int.incl=int.incl, allvipos=allvipos, coef.na=coef.na,
               yi=yi, vi=vi, X=X, weights=weights, yi.f=yi.f, vi.f=vi.f, X.f=X.f, weights.f=weights.f, M=M,
               ai.f=ai.f, bi.f=bi.f, ci.f=ci.f, di.f=di.f,
               x1i.f=x1i.f, x2i.f=x2i.f, t1i.f=t1i.f, t2i.f=t2i.f, ni=ni, ni.f=ni.f,
               ids=ids, not.na=not.na, subset=subset, slab=slab, slab.null=slab.null,
               measure=measure, method=method, weighted=weighted, knha=knha, dfs=dfs, s2w=s2w, btt=btt, intercept=intercept, digits=digits, level=level, control=control, verbose=verbose,
               add=add, to=to, drop00=drop00,
               fit.stats=fit.stats, version=packageVersion("metafor"), model=model, call=mf)

   if (model == "rma.tau2") {

      res$b.tau2     <- b.tau2
      res$vb.tau2    <- vb.tau2
      res$se.tau2    <- se.tau2
      res$zval.tau2  <- zval.tau2
      res$pval.tau2  <- pval.tau2
      res$ci.lb.tau2 <- ci.lb.tau2
      res$ci.ub.tau2 <- ci.ub.tau2
      res$Z <- Z

   }

   class(res) <- c("rma.uni", "rma")
   return(res)

}
