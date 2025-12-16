rma.glmm <- function(ai, bi, ci, di, n1i, n2i, x1i, x2i, t1i, t2i, xi, mi, ti, ni, mods,
measure, data, slab, subset,
add=1/2, to="only0", drop00=TRUE, intercept=TRUE,
model="UM.FS", method="ML", coding=1/2, cor=FALSE,
test="z", level=95, btt, nAGQ=7, verbose=FALSE, digits, control, ...) {

   #########################################################################

   ###### setup

   mstyle <- .get.mstyle()

   ### check argument specifications
   ### (arguments "to" and "vtype" are checked inside escalc function)

   if (missing(measure))
      stop(mstyle$stop("Must specify the 'measure' argument."))

   if (!is.element(measure, c("OR","IRR","PLO","IRLN", "PR","RR","RD","PLN")))
      stop(mstyle$stop("Unknown 'measure' specified."))

   if (!is.element(method, c("FE","EE","CE","ML")))
      stop(mstyle$stop("Unknown 'method' specified."))

   if (!is.element(coding, c(1/2, 1, 0)))
      stop(mstyle$stop("Unknown 'coding' option specified."))

   ### in case the user specified more than one add/to value (as one can do with rma.mh() and rma.peto())
   ### (never apply any kind of continuity correction to the data used in the actual model fitting for models implemented in this function)

   if (length(add) > 1L)
      add <- add[1]

   if (length(to) > 1L)
      to <- to[1]

   ### model argument only relevant for 2x2 table data (measure="OR") and for 2-group rate data (measure="IRR")
   ### UM.FS/UM.RS = unconditional GLMM with fixed/random study effects (logistic or poisson mixed-effects model with fixed/random intercepts)
   ### CM.EL/CM.AL = conditional GLMM (exact/approximate) (hypergeometric or conditional logistic model)
   ### BV/MV       = bi/multivariate model (logistic or poisson mixed-effects model with unstructured covariance matrix) -- not implemented

   if (!is.element(model, c("UM.FS","UM.RS","CM.EL","CM.AL")))
      stop(mstyle$stop("Unknown 'model' specified."))

   ### no need for CM.AL for IRR -- use CM.EL

   if (model == "CM.AL" && measure == "IRR")
      model <- "CM.EL"

   ### check if user changed model for measures where this is not relevant; if so, issue a warning

   if (is.element(measure, c("PLO","PR","PLN","IRLN")) && !is.null(match.call()$model))
      warning(mstyle$warning("Argument 'model' not relevant for this outcome measure."), call.=FALSE)

   ### warning about experimental measures

   if (!is.element(measure, c("OR","IRR","PLO","IRLN")))
      warning(mstyle$warning("The use of this 'measure' is experimental - treat results with caution."), call.=FALSE)

   if (is.element(model, c("CM.EL","CM.AL")) && is.element(measure, c("RR","RD")))
      stop(mstyle$stop("Cannot use this measure with model='CM.EL' or model='CM.AL'."))

   na.act <- getOption("na.action")
   on.exit(options(na.action=na.act), add=TRUE)

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   if (missing(control))
      control <- list()

   time.start <- proc.time()

   ### get ... argument and check for extra/superfluous arguments

   ddd <- list(...)

   .chkdots(ddd, c("vtype", "tdist", "outlist", "onlyo1", "addyi", "addvi", "time", "retdat", "family", "retfit", "skiphet", "i2def", "link"))

   if (is.null(ddd$vtype)) {
      vtype <- "LS"
   } else {
      vtype <- ddd$vtype
   }

   ### handle 'tdist' argument from ... (note: overrides test argument)

   if (isFALSE(ddd$tdist))
      test <- "z"
   if (isTRUE(ddd$tdist))
      test <- "t"

   if (!is.element(test, c("z", "t")))
      stop(mstyle$stop("Unknown option specified for the 'test' argument."))

   ### set defaults or get 'onlyo1', 'addyi', and 'addvi' arguments

   onlyo1 <- .chkddd(ddd$onlyo1, FALSE)
   addyi  <- .chkddd(ddd$addyi,  TRUE)
   addvi  <- .chkddd(ddd$addvi,  TRUE)

   ### set default for 'i2def'

   i2def <- .chkddd(ddd$i2def, "1")

   ### set defaults for 'digits'

   if (missing(digits)) {
      digits <- .set.digits(dmiss=TRUE)
   } else {
      digits <- .set.digits(digits, dmiss=FALSE)
   }

   ### set default for formula.mods

   formula.mods <- NULL

   ### set options(warn=1) if verbose > 2

   if (verbose > 2) {
      opwarn <- options(warn=1)
      on.exit(options(warn=opwarn$warn), add=TRUE)
   }

   if (is.null(ddd$link)) {
      if (measure=="OR" || measure=="PLO")
         link <- "logit"
      if (measure=="RR" || measure=="PLN")
         link <- "log"
      if (measure=="RD" || measure=="PR")
         link <- "identity"
      if (measure=="IRR" || measure=="IRLN")
         link <- "log"
   } else {
      link <- ddd$link
   }

   #########################################################################

   if (verbose) .space()

   if (verbose > 1)
      message(mstyle$message("Extracting the data and computing yi/vi values ..."))

   ### check if the 'data' argument was specified

   if (missing(data))
      data <- NULL

   if (is.null(data)) {
      data <- sys.frame(sys.parent())
   } else {
      if (!is.data.frame(data))
         data <- data.frame(data)
   }

   mf <- match.call()

   ### extract slab, subset, and mods values, possibly from the data frame specified via data (arguments not specified are NULL)

   slab   <- .getx("slab",   mf=mf, data=data)
   subset <- .getx("subset", mf=mf, data=data)
   mods   <- .getx("mods",   mf=mf, data=data)

   ai <- bi <- ci <- di <- x1i <- x2i <- t1i <- t2i <- xi <- mi <- ti <- NA_real_

   ### calculate yi and vi values

   if (is.element(measure, c("OR","RR","RD"))) {

      ai  <- .getx("ai",  mf=mf, data=data, checknumeric=TRUE)
      bi  <- .getx("bi",  mf=mf, data=data, checknumeric=TRUE)
      ci  <- .getx("ci",  mf=mf, data=data, checknumeric=TRUE)
      di  <- .getx("di",  mf=mf, data=data, checknumeric=TRUE)
      n1i <- .getx("n1i", mf=mf, data=data, checknumeric=TRUE)
      n2i <- .getx("n2i", mf=mf, data=data, checknumeric=TRUE)

      if (is.null(bi)) bi <- n1i - ai
      if (is.null(di)) di <- n2i - ci

      k <- length(ai) # number of outcomes before subsetting
      k.all <- k

      if (!is.null(subset)) {
         subset <- .chksubset(subset, k)
         ai <- .getsubset(ai, subset)
         bi <- .getsubset(bi, subset)
         ci <- .getsubset(ci, subset)
         di <- .getsubset(di, subset)
      }

      args <- list(ai=ai, bi=bi, ci=ci, di=di, add=add, to=to, drop00=drop00, onlyo1=onlyo1, addyi=addyi, addvi=addvi)

   }

   if (is.element(measure, c("IRR"))) {

      x1i <- .getx("x1i", mf=mf, data=data, checknumeric=TRUE)
      x2i <- .getx("x2i", mf=mf, data=data, checknumeric=TRUE)
      t1i <- .getx("t1i", mf=mf, data=data, checknumeric=TRUE)
      t2i <- .getx("t2i", mf=mf, data=data, checknumeric=TRUE)

      k <- length(x1i) # number of outcomes before subsetting
      k.all <- k

      if (!is.null(subset)) {
         subset <- .chksubset(subset, k)
         x1i <- .getsubset(x1i, subset)
         x2i <- .getsubset(x2i, subset)
         t1i <- .getsubset(t1i, subset)
         t2i <- .getsubset(t2i, subset)
      }

      args <- list(x1i=x1i, x2i=x2i, t1i=t1i, t2i=t2i, add=add, to=to, drop00=drop00, onlyo1=onlyo1, addyi=addyi, addvi=addvi)

   }

   if (is.element(measure, c("PLO","PR","PLN"))) {

      xi <- .getx("xi", mf=mf, data=data, checknumeric=TRUE)
      mi <- .getx("mi", mf=mf, data=data, checknumeric=TRUE)
      ni <- .getx("ni", mf=mf, data=data, checknumeric=TRUE)

      if (is.null(mi)) mi <- ni - xi

      k <- length(xi) # number of outcomes before subsetting
      k.all <- k

      if (!is.null(subset)) {
         subset <- .chksubset(subset, k)
         xi <- .getsubset(xi, subset)
         mi <- .getsubset(mi, subset)
      }

      args <- list(xi=xi, mi=mi, add=add, to=to, onlyo1=onlyo1, addyi=addyi, addvi=addvi)

   }

   if (is.element(measure, c("IRLN"))) {

      xi <- .getx("xi", mf=mf, data=data, checknumeric=TRUE)
      ti <- .getx("ti", mf=mf, data=data, checknumeric=TRUE)

      k <- length(xi) # number of outcomes before subsetting
      k.all <- k

      if (!is.null(subset)) {
         subset <- .chksubset(subset, k)
         xi <- .getsubset(xi, subset)
         ti <- .getsubset(ti, subset)
      }

      args <- list(xi=xi, ti=ti, add=add, to=to, onlyo1=onlyo1, addyi=addyi, addvi=addvi)

   }

   args <- c(args, list(measure=measure, vtype=vtype))

   dat <- .do.call(escalc, args)

   yi <- dat$yi         # one or more yi/vi pairs may be NA/NA (note: yi/vi pairs that are NA/NA may still have 'valid' table data)
   vi <- dat$vi         # one or more yi/vi pairs may be NA/NA (note: yi/vi pairs that are NA/NA may still have 'valid' table data)
   ni <- attr(yi, "ni") # unadjusted total sample sizes (ni.u in escalc)

   ### study ids (1:k sequence before subsetting)

   ids <- seq_len(k)

   #########################################################################

   if (verbose > 1)
      message(mstyle$message("Creating the model matrix ..."))

   ### convert mods formula to X matrix and set intercept equal to FALSE

   if (inherits(mods, "formula")) {
      formula.mods <- mods
      if (.is.tilde1(formula.mods)) { # needed so 'mods = ~ 1' without 'data' specified works
         mods <- matrix(1, nrow=k, ncol=1)
         intercept <- FALSE
      } else {
         options(na.action = "na.pass")        # set na.action to na.pass, so that NAs are not filtered out (we'll do that later)
         mods <- model.matrix(mods, data=data) # extract the model matrix
         attr(mods, "assign") <- NULL          # strip the 'assign' attribute (not used at the moment)
         options(na.action = na.act)           # set na.action back to na.act
         intercept <- FALSE                    # set 'intercept' to FALSE since the formula now controls whether the intercept is included
      }
   }

   ### turn a vector for mods into a column vector

   if (.is.vector(mods))
      mods <- cbind(mods)

   ### turn a mods data frame into a matrix

   if (is.data.frame(mods))
      mods <- as.matrix(mods)

   ### check if the model matrix contains character variables

   if (is.character(mods))
      stop(mstyle$stop("The model matrix contains character variables."))

   ### check if the 'mods' matrix has the right number of rows

   if (!is.null(mods) && nrow(mods) != k)
      stop(mstyle$stop(paste0("Number of rows in the model matrix (", nrow(mods), ") does not match the length of the the outcome vector (", k, ").")))

   ### generate study labels if none are specified

   if (verbose > 1)
      message(mstyle$message("Generating/extracting the study labels ..."))

   if (is.null(slab)) {

      slab.null <- TRUE
      slab      <- ids

   } else {

      if (anyNA(slab))
         stop(mstyle$stop("NAs in study labels."))

      if (length(slab) != k)
         stop(mstyle$stop(paste0("Length of the 'slab' argument (", length(slab), ") does not correspond to the size of the dataset (", k, ").")))

      if (is.factor(slab))
         slab <- as.character(slab)

      slab.null <- FALSE

   }

   ### if a subset of studies is specified (note: tables, yi/vi, and ni are already subsetted above)

   if (!is.null(subset)) {

      if (verbose > 1)
         message(mstyle$message("Subsetting ..."))

      mods <- .getsubset(mods, subset)
      slab <- .getsubset(slab, subset)
      ids  <- .getsubset(ids,  subset)

   }

   ### check if the study labels are unique; if not, make them unique

   if (anyDuplicated(slab))
      slab <- .make.unique(slab)

   ### add the 'slab' attribute back to 'yi'

   attr(yi, "slab") <- slab

   k <- length(yi) # number of tables/outcomes after subsetting (can all still include NAs)

   ### if drop00=TRUE, set counts to NA for studies that have no events (or all events) in both arms (corresponding yi/vi will also be NA/NA then)

   if (is.element(measure, c("OR","RR","RD"))) {
      if (drop00) {
         id00 <- c(ai == 0L & ci == 0L) | c(bi == 0L & di == 0L)
         id00[is.na(id00)] <- FALSE
         ai[id00] <- NA_real_
         bi[id00] <- NA_real_
         ci[id00] <- NA_real_
         di[id00] <- NA_real_
      }
   }

   if (is.element(measure, c("IRR"))) {
      if (drop00) {
         id00 <- c(x1i == 0L & x2i == 0L)
         id00[is.na(id00)] <- FALSE
         x1i[id00] <- NA_real_
         x2i[id00] <- NA_real_
      }
   }

   ### save the full data (including potential NAs in table data, yi/vi/ni/mods) (after subsetting)

   outdat.f <- list(ai=ai, bi=bi, ci=ci, di=di, x1i=x1i, x2i=x2i, t1i=t1i, t2i=t2i, xi=xi, mi=mi, ni=ni, ti=ti)

   yi.f <- yi
   vi.f <- vi
   ni.f <- ni
   mods.f <- mods

   k.f <- k # total number of tables/outcomes and rows in the model matrix (including all NAs)

   ### check for NAs in tables (and corresponding mods) and act accordingly

   if (is.element(measure, c("OR","RR","RD"))) {

      has.na <- is.na(ai) | is.na(bi) | is.na(ci) | is.na(di) | (if (is.null(mods)) FALSE else apply(is.na(mods), 1, any))
      not.na <- !has.na

      if (any(has.na)) {

         if (verbose > 1)
            message(mstyle$message("Handling NAs in the table data ..."))

         if (na.act == "na.omit" || na.act == "na.exclude" || na.act == "na.pass") {
            ai   <- ai[not.na]
            bi   <- bi[not.na]
            ci   <- ci[not.na]
            di   <- di[not.na]
            mods <- mods[not.na,,drop=FALSE]
            k    <- length(ai)
            warning(mstyle$warning(paste(sum(has.na), ifelse(sum(has.na) > 1, "studies", "study"), "with NAs omitted from model fitting.")), call.=FALSE)
         }

         if (na.act == "na.fail")
            stop(mstyle$stop("Missing values in studies."))

      }

   }

   if (is.element(measure, "IRR")) {

      has.na <- is.na(x1i) | is.na(x2i) | is.na(t1i) | is.na(t2i) | (if (is.null(mods)) FALSE else apply(is.na(mods), 1, any))
      not.na <- !has.na

      if (any(has.na)) {

         if (verbose > 1)
            message(mstyle$message("Handling NAs in the table data ..."))

         if (na.act == "na.omit" || na.act == "na.exclude" || na.act == "na.pass") {
            x1i  <- x1i[not.na]
            x2i  <- x2i[not.na]
            t1i  <- t1i[not.na]
            t2i  <- t2i[not.na]
            mods <- mods[not.na,,drop=FALSE]
            k    <- length(x1i)
            warning(mstyle$warning(paste(sum(has.na), ifelse(sum(has.na) > 1, "studies", "study"), "with NAs omitted from model fitting.")), call.=FALSE)
         }

         if (na.act == "na.fail")
            stop(mstyle$stop("Missing values in studies."))

      }

   }

   if (is.element(measure, c("PLO","PR","PLN"))) {

      has.na <- is.na(xi) | is.na(mi) | (if (is.null(mods)) FALSE else apply(is.na(mods), 1, any))
      not.na <- !has.na

      if (any(has.na)) {

         if (verbose > 1)
            message(mstyle$message("Handling NAs in the table data ..."))

         if (na.act == "na.omit" || na.act == "na.exclude" || na.act == "na.pass") {
            xi   <- xi[not.na]
            mi   <- mi[not.na]
            mods <- mods[not.na,,drop=FALSE]
            k    <- length(xi)
            warning(mstyle$warning(paste(sum(has.na), ifelse(sum(has.na) > 1, "studies", "study"), "with NAs omitted from model fitting.")), call.=FALSE)
         }

         if (na.act == "na.fail")
            stop(mstyle$stop("Missing values in studies."))

      }

   }

   if (is.element(measure, "IRLN")) {

      has.na <- is.na(xi) | is.na(ti) | (if (is.null(mods)) FALSE else apply(is.na(mods), 1, any))
      not.na <- !has.na

      if (any(has.na)) {

         if (verbose > 1)
            message(mstyle$message("Handling NAs in the table data ..."))

         if (na.act == "na.omit" || na.act == "na.exclude" || na.act == "na.pass") {
            xi   <- xi[not.na]
            ti   <- ti[not.na]
            mods <- mods[not.na,,drop=FALSE]
            k    <- length(xi)
            warning(mstyle$warning(paste(sum(has.na), ifelse(sum(has.na) > 1, "studies", "study"), "with NAs omitted from model fitting.")), call.=FALSE)
         }

         if (na.act == "na.fail")
            stop(mstyle$stop("Missing values in studies."))

      }

   }

   ### note: k   = number of tables (and corresponding rows of 'mods') after removing NAs
   ###       k.f = total number of tables/outcomes and rows in the model matrix (including all NAs) stored in .f elements

   ### at least one study left?

   if (k < 1L)
      stop(mstyle$stop("Processing terminated since k = 0."))

   ### check for NAs in yi/vi and act accordingly (yi/vi pair can be NA/NA if add=0 is used)
   ### note: if a table was removed because of NAs in mods, must also remove the corresponding yi/vi pair;
   ###       also, must use mods.f here, since NAs in mods were already removed above (and need a separate
   ###       mods.yi element, so that dimensions of the model matrix and vi are guaranteed to match up)

   mods.yi <- mods.f
   yivi.na <- is.na(yi) | is.na(vi) | (if (is.null(mods.yi)) FALSE else apply(is.na(mods.yi), 1, any))
   not.na.yivi <- !yivi.na

   if (any(yivi.na)) {

      if (verbose > 1)
         message(mstyle$message("Handling NAs in yi/vi ..."))

      if (na.act == "na.omit" || na.act == "na.exclude" || na.act == "na.pass") {

         yi <- yi[not.na.yivi]
         ni <- ni[not.na.yivi]
         vi <- vi[not.na.yivi]
         mods.yi <- mods.f[not.na.yivi,,drop=FALSE]
         warning(mstyle$warning("Some yi/vi values are NA."), call.=FALSE)

         attr(yi, "measure") <- measure # add measure attribute back
         attr(yi, "ni")      <- ni      # add ni attribute back

      }

      if (na.act == "na.fail")
         stop(mstyle$stop("Missing yi/vi values."))

   }

   k.yi <- length(yi) # number of yi/vi pairs that are not NA

   ### make sure that there is at least one column in X

   if (is.null(mods) && !intercept) {
      warning(mstyle$warning("Must either include an intercept and/or moderators in the model.\nCoerced an intercept into the model."), call.=FALSE)
      intercept <- TRUE
   }

   if (!is.null(mods) && ncol(mods) == 0L) {
      warning(mstyle$warning("Cannot fit model with an empty model matrix. Coerced an intercept into the model."), call.=FALSE)
      intercept <- TRUE
   }

   ### add vector of 1s to the X matrix for the intercept (if intercept=TRUE)

   if (intercept) {
      X    <- cbind(intrcpt=rep(1,k), mods)
      X.f  <- cbind(intrcpt=rep(1,k.f), mods.f)
      X.yi <- cbind(intrcpt=rep(1,k.yi), mods.yi)
   } else {
      X    <- mods
      X.f  <- mods.f
      X.yi <- mods.yi
   }

   ### drop redundant predictors
   ### note: yi may have become shorter than X due to the omission of NAs, so just use a fake yi vector here

   tmp <- lm(rep(0,k) ~ X - 1)
   coef.na <- is.na(coef(tmp))
   if (any(coef.na)) {
      warning(mstyle$warning("Redundant predictors dropped from the model."), call.=FALSE)
      X    <- X[,!coef.na,drop=FALSE]
      X.f  <- X.f[,!coef.na,drop=FALSE]
   }

   ### need to do this separately for X.yi, since model matrix may have fewer rows due to removal of NA/NA pairs for yi/vi

   tmp <- lm(yi ~ X.yi - 1)
   coef.na <- is.na(coef(tmp))
   if (any(coef.na))
      X.yi <- X.yi[,!coef.na,drop=FALSE]

   ### check whether the intercept is included and if yes, move it to the first column (NAs already removed, so na.rm=TRUE for any() not necessary)

   is.int <- apply(X, 2, .is.intercept)
   if (any(is.int)) {
      int.incl <- TRUE
      int.indx <- which(is.int, arr.ind=TRUE)
      X        <- cbind(intrcpt=1,   X[,-int.indx, drop=FALSE]) # note: this removes any duplicate intercepts
      X.f      <- cbind(intrcpt=1, X.f[,-int.indx, drop=FALSE]) # note: this removes any duplicate intercepts
      intercept <- TRUE # set intercept appropriately so that the predict() function works
   } else {
      int.incl <- FALSE
   }

   ### need to do this separately for X.yi, since model matrix may have fewer rows due to removal of NA/NA pairs for yi/vi

   is.int <- apply(X.yi, 2, .is.intercept)
   if (any(is.int)) {
      int.indx <- which(is.int, arr.ind=TRUE)
      X.yi     <- cbind(intrcpt=1, X.yi[,-int.indx, drop=FALSE]) # note: this removes any duplicate intercepts
   }

   p <- NCOL(X) # number of columns in X (including the intercept if it is included)

   ### note: number of columns in X.yi may be lower than p; but computation of I^2 below is based on p

   ### make sure variable names in X are unique

   colnames(X) <- colnames(X.f) <- .make.unique(colnames(X))

   ### check whether this is an intercept-only model

   if ((p == 1L) && .is.intercept(X)) {
      int.only <- TRUE
   } else {
      int.only <- FALSE
   }

   ### check if there are too many parameters for given k

   if (is.element(method, c("FE","EE","CE")) && p > k)
      stop(mstyle$stop("Number of parameters to be estimated is larger than the number of observations."))
   if (!is.element(method, c("FE","EE","CE")) && (p+1) > k)
      stop(mstyle$stop("Number of parameters to be estimated is larger than the number of observations."))

   ### set/check 'btt' argument

   btt <- .set.btt(btt, p, int.incl, colnames(X))
   m <- length(btt) # number of betas to test (m = p if all betas are tested)

   #########################################################################

   ### set defaults for control parameters

   con <- list(verbose = FALSE,            # also passed on to glm/glmer/optim/nlminb/minqa (uobyqa/newuoa/bobyqa)
               package="lme4",             # package for fitting logistic mixed-effects models ("lme4", "GLMMadaptive", "glmmTMB")
               optimizer = "nlminb",       # optimizer to use for CM.EL+OR ("optim","nlminb","uobyqa","newuoa","bobyqa","nloptr","nlm","hjk","nmk","mads","ucminf","lbfgsb3c","subplex","BBoptim","optimParallel","clogit","clogistic","Rcgmin","Rvmmin")
               optmethod = "BFGS",         # argument 'method' for optim() ("Nelder-Mead" and "BFGS" are sensible options)
               parallel = list(),          # parallel argument for optimParallel() (note: 'cl' argument in parallel is not passed; this is directly specified via 'cl')
               cl = NULL,                  # arguments for optimParallel()
               ncpus = 1L,                 # arguments for optimParallel()
               scaleX = TRUE,              # whether non-dummy variables in the X matrix should be rescaled before model fitting
               evtol = 1e-07,              # lower bound for eigenvalues to determine if the model matrix is positive definite
               dnchgcalc = "dFNCHypergeo", # method for calculating dnchg ("dFNCHypergeo" from BiasedUrn package or "dnoncenhypergeom")
               dnchgprec = 1e-10,          # precision for dFNCHypergeo()
               hesspack = "numDeriv",      # package for computing the Hessian (numDeriv or pracma)
               tau2tol = 1e-04)            # for "CM.EL" + "ML", threshold for treating tau^2 values as effectively equal to 0

   ### replace defaults with any user-defined values

   con.pos <- pmatch(names(control), names(con))
   con[c(na.omit(con.pos))] <- control[!is.na(con.pos)]

   if (verbose)
      con$verbose <- verbose

   verbose <- con$verbose

   optimizer <- match.arg(con$optimizer, c("optim","nlminb","uobyqa","newuoa","bobyqa","nloptr","nlm","hjk","nmk","mads","ucminf","lbfgsb3c","subplex","BBoptim","optimParallel","clogit","clogistic","Nelder-Mead","BFGS","CG","L-BFGS-B","SANN","Brent","Rcgmin","Rvmmin"))
   optmethod <- match.arg(con$optmethod, c("Nelder-Mead","BFGS","CG","L-BFGS-B","SANN","Brent"))
   if (optimizer %in% c("Nelder-Mead","BFGS","CG","L-BFGS-B","SANN","Brent")) {
      optmethod <- optimizer
      optimizer <- "optim"
   }
   package  <- match.arg(con$package, c("lme4","GLMMadaptive","glmmTMB"))
   parallel <- con$parallel
   cl       <- con$cl
   ncpus    <- con$ncpus

   if (con$dnchgcalc != "dnoncenhypergeom" && con$dnchgcalc != "dFNCHypergeo")
      stop(mstyle$stop("Unknown dnchgcalc method specified."))

   if (is.element(optimizer, c("clogit","clogistic")) && method == "ML")
      stop(mstyle$stop("Cannot use 'clogit' or 'clogistic' with method='ML'."))

   if (package == "lme4" && is.element(measure, c("OR","RR","RD","IRR")) && model == "UM.RS" && method == "ML" && nAGQ > 1) {
      warning(mstyle$warning("Not possible to fit RE/ME model='UM.RS' with nAGQ > 1 with glmer(). nAGQ automatically set to 1."), call.=FALSE)
      nAGQ <- 1
   }

   ### if control argument 'ncpus' is larger than 1, automatically switch to the 'optimParallel' optimizer

   if (ncpus > 1L)
      optimizer <- "optimParallel"

   pos.optCtrl <- pmatch(names(control), "optCtrl", nomatch=0)
   if (sum(pos.optCtrl) > 0) {
      optCtrl <- control[[which(pos.optCtrl == 1)]]
   } else {
      optCtrl <- list()
   }

   ### set NLOPT_LN_BOBYQA as the default algorithm for nloptr optimizer
   ### and by default use a relative convergence criterion of 1e-8 on the function value

   if (optimizer == "nloptr" && !is.element("algorithm", names(optCtrl)))
      optCtrl$algorithm <- "NLOPT_LN_BOBYQA"

   if (optimizer == "nloptr" && !is.element("ftol_rel", names(optCtrl)))
      optCtrl$ftol_rel <- 1e-8

   ### for mads, set trace=FALSE and tol=1e-6 by default

   if (optimizer == "mads" && !is.element("trace", names(optCtrl)))
      optCtrl$trace <- FALSE

   if (optimizer == "mads" && !is.element("tol", names(optCtrl)))
      optCtrl$tol <- 1e-6

   ### for subplex, set reltol=1e-8 by default (the default in subplex() is .Machine$double.eps)

   if (optimizer == "subplex" && !is.element("reltol", names(optCtrl)))
      optCtrl$reltol <- 1e-8

   ### for BBoptim, set trace=FALSE by default

   if (optimizer == "BBoptim" && !is.element("trace", names(optCtrl)))
      optCtrl$trace <- FALSE

   if (optimizer == "optim") {
      con.pos <- pmatch(names(optCtrl), "REPORT", nomatch=0) # set REPORT to 1 if it is not already set by the user
      if (sum(con.pos) > 0) {
         names(optCtrl)[which(con.pos == 1)] <- "REPORT"
      } else {
         optCtrl$REPORT <- 1
      }
      optCtrl$trace <- con$verbose # trace for optim is a non-negative integer
   }

   if (optimizer == "nlminb")
      optCtrl$trace <- ifelse(con$verbose > 0, 1, 0) # set trace to 1, so information is printed every iteration

   if (is.element(optimizer, c("uobyqa", "newuoa", "bobyqa")))
      optCtrl$iprint <- ifelse(con$verbose > 0, 3, 0) # set iprint to 3 for maximum information

   pos.clogitCtrl <- pmatch(names(control), "clogitCtrl", nomatch=0)
   if (sum(pos.clogitCtrl) > 0) {
      clogitCtrl <- control[[which(pos.clogitCtrl == 1)]]
   } else {
      clogitCtrl <- list()
   }

   pos.clogisticCtrl <- pmatch(names(control), "clogisticCtrl", nomatch=0)
   if (sum(pos.clogisticCtrl) > 0) {
      clogisticCtrl <- control[[which(pos.clogisticCtrl == 1)]]
   } else {
      clogisticCtrl <- list()
   }

   pos.glmCtrl <- pmatch(names(control), "glmCtrl", nomatch=0)
   if (sum(pos.glmCtrl) > 0) {
      glmCtrl <- control[[which(pos.glmCtrl == 1)]]
   } else {
      glmCtrl <- list()
   }
   glmCtrl$trace <- ifelse(con$verbose > 0, TRUE, FALSE) # trace for glmCtrl is logical

   pos.glmerCtrl <- pmatch(names(control), "glmerCtrl", nomatch=0)
   if (sum(pos.glmerCtrl) > 0) {
      glmerCtrl <- control[[which(pos.glmerCtrl == 1)]]
   } else {
      glmerCtrl <- list()
   }

   pos.intCtrl <- pmatch(names(control), "intCtrl", nomatch=0)
   if (sum(pos.intCtrl) > 0) {
      intCtrl <- control[[which(pos.intCtrl == 1)]]
   } else {
      intCtrl <- list()
   }
   con.pos <- pmatch(names(intCtrl), "lower", nomatch=0)
   if (sum(con.pos) > 0) {
      names(intCtrl)[which(con.pos == 1)] <- "lower"
   } else {
      intCtrl$lower <- -Inf
   }
   con.pos <- pmatch(names(intCtrl), "upper", nomatch=0)
   if (sum(con.pos) > 0) {
      names(intCtrl)[which(con.pos == 1)] <- "upper"
   } else {
      intCtrl$upper <- Inf
   }
   con.pos <- pmatch(names(intCtrl), "subdivisions", nomatch=0)
   if (sum(con.pos) > 0) {
      names(intCtrl)[which(con.pos == 1)] <- "subdivisions"
   } else {
      intCtrl$subdivisions <- 100L
   }
   con.pos <- pmatch(names(intCtrl), "rel.tol", nomatch=0)
   if (sum(con.pos) > 0) {
      names(intCtrl)[which(con.pos == 1)] <- "rel.tol"
   } else {
      intCtrl$rel.tol <- .Machine$double.eps^0.25
   }

   pos.hessianCtrl <- pmatch(names(control), "hessianCtrl", nomatch=0)
   if (sum(pos.hessianCtrl) > 0) {
      hessianCtrl <- control[[which(pos.hessianCtrl == 1)]]
   } else {
      hessianCtrl <- list(r=16)
   }

   #return(list(verbose=verbose, optimizer=optimizer, dnchgcalc=con$dnchgcalc, dnchgprec=con$dnchgprec, optCtrl=optCtrl, glmCtrl=glmCtrl, glmerCtrl=glmerCtrl, intCtrl=intCtrl, hessianCtrl=hessianCtrl))

   #########################################################################

   ### check that the required packages are installed

   if (is.element(measure, c("OR","RR","RD","IRR"))) {
      if ((model == "UM.FS" && method == "ML") || (model == "UM.RS") || (model == "CM.AL" && method == "ML") || (model == "CM.EL" && method == "ML")) {
         if (!requireNamespace(package, quietly=TRUE))
            stop(mstyle$stop(paste0("Please install the '", package, "' package to fit this model.")))
      }
   }

   if (is.element(measure, c("PLO","PR","PLN","IRLN")) && method == "ML") {
      if (!requireNamespace(package, quietly=TRUE))
         stop(mstyle$stop(paste0("Please install the '", package, "' package to fit this model.")))
   }

   if (measure == "OR" && model == "CM.EL") {

      if (is.element(optimizer, c("uobyqa","newuoa","bobyqa"))) {
         if (!requireNamespace("minqa", quietly=TRUE))
            stop(mstyle$stop("Please install the 'minqa' package to fit this model."))
      }

      if (is.element(optimizer, c("nloptr","ucminf","lbfgsb3c","subplex","optimParallel"))) {
         if (!requireNamespace(optimizer, quietly=TRUE))
            stop(mstyle$stop(paste0("Please install the '", optimizer, "' package to use this optimizer.")))
      }

      if (is.element(optimizer, c("hjk","nmk","mads"))) {
         if (!requireNamespace("dfoptim", quietly=TRUE))
            stop(mstyle$stop("Please install the 'dfoptim' package to use this optimizer."))
      }

      if (optimizer == "BBoptim") {
         if (!requireNamespace("BB", quietly=TRUE))
            stop(mstyle$stop("Please install the 'BB' package to use this optimizer."))
      }

      if (is.element(optimizer, c("Rcgmin","Rvmmin"))) {
         if (!requireNamespace("optimx", quietly=TRUE))
            stop(mstyle$stop(paste0("Please install the 'optimx' package to use this optimizer.")))
      }

      if (is.element(optimizer, c("optim","nlminb","uobyqa","newuoa","bobyqa","nloptr","nlm","hjk","nmk","mads","ucminf","lbfgsb3c","subplex","BBoptim","optimParallel","Rcgmin","Rvmmin"))) {
         con$hesspack <- match.arg(con$hesspack, c("numDeriv","pracma","calculus"))
         if (!requireNamespace(con$hesspack, quietly=TRUE))
            stop(mstyle$stop(paste0("Please install the '", con$hesspack, "' package to fit this model.")))
         if (con$dnchgcalc == "dFNCHypergeo") {
            if (!requireNamespace("BiasedUrn", quietly=TRUE))
               stop(mstyle$stop("Please install the 'BiasedUrn' package to fit this model."))
         }
      }

      if (optimizer == "clogit") {
         if (!requireNamespace("survival", quietly=TRUE))
            stop(mstyle$stop("Please install the 'survival' package to fit this model."))
         coxph   <- survival::coxph
         Surv    <- survival::Surv
         clogit  <- survival::clogit
         strata  <- survival::strata
      }

      if (optimizer == "clogistic") {
         if (!requireNamespace("Epi", quietly=TRUE))
            stop(mstyle$stop("Please install the 'Epi' package to fit this model."))
      }

   }

   ### check whether the model matrix is of full rank

   if (!.chkpd(crossprod(X), tol=con$evtol))
      stop(mstyle$stop("Model matrix not of full rank. Cannot fit model."))

   #########################################################################
   #########################################################################
   #########################################################################

   se.tau2 <- ci.lb.tau2 <- ci.ub.tau2 <- I2 <- H2 <- QE <- QEp <- NA_real_
   se.warn <- FALSE

   rho <- NA_real_

   level <- .level(level)

   ###### model fitting, test statistics, and confidence intervals

   ### upgrade warnings to errors (for some testing)
   #o.warn <- getOption("warn")
   #on.exit(options(warn = o.warn), add=TRUE)
   #options(warn = 2)

   ### rescale X matrix (only for models with moderators and models including an intercept term)

   if (!int.only && int.incl && con$scaleX) {
      Xsave <- X
      meanX <- colMeans(X[, 2:p, drop=FALSE])
      sdX   <- apply(X[, 2:p, drop=FALSE], 2, sd)          # consider using colSds() from matrixStats package
      is.d  <- apply(X, 2, .is.dummy)                      # is each column a dummy variable (i.e., only 0s and 1s)?
      X[,!is.d] <- apply(X[, !is.d, drop=FALSE], 2, scale) # rescale the non-dummy variables
   }

   #########################################################################
   #########################################################################
   #########################################################################

   ### two group outcomes (odds ratios and incidence rate ratios)

   if (is.element(measure, c("OR","RR","RD","IRR"))) {

      ######################################################################

      if (is.element(model, c("UM.FS","UM.RS"))) {

         ### prepare data for the unconditional models

         if (is.element(measure, c("OR","RR","RD"))) {               #                           xi   mi   study   group1  group2  group12  offset  intrcpt  mod1
            dat.grp <- cbind(xi=c(rbind(ai,ci)), mi=c(rbind(bi,di))) # grp-level outcome data    ai   bi   i       1       0       +1/2     NULL    1        x1i
                                                                     #                           ci   di   i       0       1       -1/2     NULL    0        0
            if (is.null(ddd$family)) {
               if (measure == "OR")
                  dat.fam <- binomial(link=link)
               if (measure == "RR")
                  dat.fam <- binomial(link=link)
               if (measure == "RD")
                  #dat.fam <- eval(parse(text="binomial(link=\"identity\")"))
                  dat.fam <- binomial(link=link)
            } else {
               dat.fam <- ddd$family
            }
            dat.off <- NULL
         }

         if (is.element(measure, c("IRR"))) {                        #                           xi   ti   study   group1  group2  group12  offset  intrcpt  mod1
            dat.grp <- c(rbind(x1i,x2i))                             # grp-level outcome data    x1i  t1i  i       1       0       +1/2     t1i     1        x1i
                                                                     # log(ti) for offset        x2i  t2i  i       0       1       -1/2     t2i     0        0
            if (is.null(ddd$family)) {
               dat.fam <- poisson(link=link)
            } else {
               dat.fam <- ddd$family
            }
            dat.off <- log(c(rbind(t1i,t2i)))
         }

         group1  <- rep(c(1,0), times=k)                             # group dummy for 1st group (ai,bi for group 1)
         group2  <- rep(c(0,1), times=k)                             # group dummy for 2nd group (ci,di for group 2) (not really needed)
         group12 <- rep(c(1/2,-1/2), times=k)                        # group dummy with +- 1/2 coding
         study   <- factor(rep(seq_len(k), each=2L))                 # study factor
         const   <- cbind(rep(1,2*k))                                # intercept for random study effects model

         X.fit   <- X[rep(seq(k), each=2L),,drop=FALSE]              # duplicate each row in X (drop=FALSE, so column names are preserved)
         X.fit   <- cbind(group1*X.fit[,,drop=FALSE])                # then multiply by group1 dummy (intercept, if included, becomes the group1 dummy)

         if (coding == 1/2)
            group <- group12
         if (coding == 1)
            group <- group1
         if (coding == 0)
            group <- group2

         rownames(X.fit) <- seq_len(2*k)

         if (isTRUE(ddd$retdat))
            return(list(dat.grp=dat.grp, X.fit=X.fit, study=study, dat.off = if (!is.null(dat.off)) dat.off else NULL, const=const, group1=group1, group2=group2, group12=group12, group=group, dat.fam=dat.fam))

         ###################################################################

         ####################################################
         ### unconditional model with fixed study effects ###
         ####################################################

         if (model == "UM.FS") {

            ### fit FE model

            if (verbose)
               message(mstyle$message("Fitting the FE model ..."))

            if (k > 1) {
               res.FE <- try(glm(dat.grp ~ -1 + X.fit + study, offset=dat.off, family=dat.fam, control=glmCtrl), silent=!verbose)
            } else {
               res.FE <- try(glm(dat.grp ~ -1 + X.fit + const, offset=dat.off, family=dat.fam, control=glmCtrl), silent=!verbose)
            }

            if (inherits(res.FE, "try-error"))
               stop(mstyle$stop(paste0("Cannot fit FE model", ifelse(verbose, ".", " (set 'verbose=TRUE' to obtain further details)."))))

            ### log-likelihood

            #ll.FE <- with(data.frame(dat.grp), sum(dbinom(xi, xi+mi, predict(res.FE, type="response"), log=TRUE))) # model has a NULL offset
            #ll.FE <- with(data.frame(dat.grp), sum(dpois(xi, predict(res.FE, type="response"), log=TRUE)))         # offset already incorporated into predict()
            ll.FE <- c(logLik(res.FE)) # same as above

            ### fit saturated FE model (= QE model)

            QEconv <- FALSE
            ll.QE <- NA_real_

            if (!isTRUE(ddd$skiphet)) {

               if (k > 1 && verbose)
                  message(mstyle$message("Fitting the saturated model ..."))

               if (k > 1) {
                  X.QE   <- model.matrix(~ -1 + X.fit + study + study:group1)
                  res.QE <- try(glm(dat.grp ~ -1 + X.QE, offset=dat.off, family=dat.fam, control=glmCtrl), silent=!verbose)
               } else {
                  res.QE <- res.FE
               }

               if (inherits(res.QE, "try-error")) {

                  warning(mstyle$warning(paste0("Cannot fit saturated model", ifelse(verbose, ".", " (set 'verbose=TRUE' to obtain further details)."))), call.=FALSE)

               } else {

                  QEconv <- TRUE

                  ### log-likelihood

                  #ll.QE <- with(data.frame(dat.grp), sum(dbinom(xi, xi+mi, xi/(xi+mi), log=TRUE))) # model has a NULL offset
                  #ll.QE <- with(data.frame(dat.grp), sum(dpois(xi, xi, log=TRUE)))                 # offset not relevant for saturated model
                  ll.QE  <- c(logLik(res.QE)) # same as above

                  ### extract coefficients and variance-covariance matrix for Wald-type test for heterogeneity
                  #b2.QE <- cbind(na.omit(coef(res.QE)[-seq_len(k+p)]))                          # coef() still includes aliased coefficients as NAs, so filter them out
                  b2.QE  <- cbind(coef(res.QE, complete=FALSE)[-seq_len(k+p)])                   # aliased coefficients are removed by coef() when complete=FALSE
                  vb2.QE <- vcov(res.QE, complete=FALSE)[-seq_len(k+p),-seq_len(k+p),drop=FALSE] # aliased coefficients are removed by vcov() when complete=FALSE

               }

            }

            if (method == "ML") {

               ### fit ML model

               if (verbose)
                  message(mstyle$message("Fitting the ML model ..."))

               if (package == "lme4") {
                  if (verbose) {
                     res.ML <- try(lme4::glmer(dat.grp ~ -1 + X.fit + study + (group - 1 | study), offset=dat.off, family=dat.fam, nAGQ=nAGQ, verbose=verbose, control=do.call(lme4::glmerControl, glmerCtrl)), silent=!verbose)
                  } else {
                     res.ML <- suppressMessages(try(lme4::glmer(dat.grp ~ -1 + X.fit + study + (group - 1 | study), offset=dat.off, family=dat.fam, nAGQ=nAGQ, verbose=verbose, control=do.call(lme4::glmerControl, glmerCtrl)), silent=!verbose))
                  }
               }

               if (package == "GLMMadaptive") {
                  if (is.element(measure, c("OR","RR","RD"))) {
                     dat.mm <- data.frame(xi=dat.grp[,"xi"], mi=dat.grp[,"mi"], study=study, group=group)
                     res.ML <- try(GLMMadaptive::mixed_model(cbind(xi,mi) ~ -1 + X.fit + study, random = ~ group - 1 | study, data=dat.mm, family=dat.fam, nAGQ=nAGQ, verbose=verbose, control=glmerCtrl), silent=!verbose)
                  } else {
                     dat.mm <- data.frame(xi=dat.grp, study=study, group=group)
                     res.ML <- try(GLMMadaptive::mixed_model(xi ~ -1 + X.fit + study + offset(dat.off), random = ~ group - 1 | study, data=dat.mm, family=dat.fam, nAGQ=nAGQ, verbose=verbose, control=glmerCtrl), silent=!verbose)
                  }
               }

               if (package == "glmmTMB") {
                  if (verbose) {
                     res.ML <- try(glmmTMB::glmmTMB(dat.grp ~ -1 + X.fit + study + (group - 1 | study), offset=dat.off, family=dat.fam, verbose=verbose, data=NULL, control=do.call(glmmTMB::glmmTMBControl, glmerCtrl)), silent=!verbose)
                  } else {
                     res.ML <- suppressMessages(try(glmmTMB::glmmTMB(dat.grp ~ -1 + X.fit + study + (group - 1 | study), offset=dat.off, family=dat.fam, verbose=verbose, data=NULL, control=do.call(glmmTMB::glmmTMBControl, glmerCtrl)), silent=!verbose))
                  }
               }

               if (inherits(res.ML, "try-error"))
                  stop(mstyle$stop(paste0("Cannot fit ML model", ifelse(verbose, ".", " (set 'verbose=TRUE' to obtain further details)."))))

               ### log-likelihood

               #ll.ML <- with(data.frame(dat.grp), sum(dbinom(xi, xi+mi, fitted(res.ML), log=TRUE))) # not correct (since it does not incorporate the random effects; same as ll.FE if tau^2=0)
               #ll.ML <- with(data.frame(dat.grp), sum(dbinom(xi, xi+mi, plogis(qlogis(fitted(res.ML)) + group12*unlist(ranef(res.ML))), log=TRUE))) # not correct (since one really has to integrate; same as ll.FE if tau^2=0)
               #ll.ML <- c(logLik(res.ML)) # this is not the same as ll.FE when tau^2 = 0 (not sure why)
               if (package == "lme4") {
                  if (is.na(ll.QE)) {
                     ll.ML <- c(logLik(res.ML))
                  } else {
                     ll.ML <- ll.QE - 1/2 * deviance(res.ML) # this makes ll.ML comparable to ll.FE (same as ll.FE when tau^2=0)

                  }
               } else {
                  ll.ML <- c(logLik(res.ML)) # not 100% sure how comparable this is to ll.FE when tau^2 = 0 (seems correct for glmmTMB)
               }
            }

            #return(list(res.FE, res.QE, res.ML, ll.FE=ll.FE, ll.QE=ll.QE, ll.ML=ll.ML))
            #res.FE <- res[[1]]; res.QE <- res[[2]]; res.ML <- res[[3]]

            if (is.element(method, c("FE","EE","CE"))) {
               beta   <- cbind(coef(res.FE)[seq_len(p)])
               vb     <- vcov(res.FE)[seq_len(p),seq_len(p),drop=FALSE]
               tau2   <- 0
               sigma2 <- NA_real_
               parms  <- p + k
               p.eff  <- p + k
               k.eff  <- 2*k
            }

            if (method == "ML") {

               if (package  == "lme4") {
                  beta   <- cbind(lme4::fixef(res.ML)[seq_len(p)])
                  vb     <- as.matrix(vcov(res.ML))[seq_len(p),seq_len(p),drop=FALSE]
                  tau2   <- lme4::VarCorr(res.ML)[[1]][1]
               }

               if (package == "GLMMadaptive") {
                  beta   <- cbind(GLMMadaptive::fixef(res.ML)[seq_len(p)])
                  vb     <- as.matrix(vcov(res.ML))[seq_len(p),seq_len(p),drop=FALSE]
                  tau2   <- res.ML$D[1,1]
               }

               if (package  == "glmmTMB") {
                  beta   <- cbind(glmmTMB::fixef(res.ML)$cond[seq_len(p)])
                  vb     <- as.matrix(vcov(res.ML)$cond)[seq_len(p),seq_len(p),drop=FALSE]
                  tau2   <- glmmTMB::VarCorr(res.ML)[[1]][[1]][[1]]
               }

               sigma2 <- NA_real_
               parms  <- p + k + 1
               p.eff  <- p + k
               k.eff  <- 2*k

            }

            #return(list(beta=beta, vb=vb, tau2=tau2, sigma2=sigma2, parms=parms, p.eff=p.eff, k.eff=k.eff, b2.QE=b2.QE, vb2.QE=vb2.QE))

         }

         ###################################################################

         #####################################################
         ### unconditional model with random study effects ###
         #####################################################

         if (model == "UM.RS") {

            ### fit FE model

            if (verbose)
               message(mstyle$message("Fitting the FE model ..."))

            if (package == "lme4") {
               if (verbose) {
                  res.FE <- try(lme4::glmer(dat.grp ~ -1 + X.fit + const + (1 | study), offset=dat.off, family=dat.fam, nAGQ=nAGQ, verbose=verbose, control=do.call(lme4::glmerControl, glmerCtrl)), silent=!verbose)
               } else {
                  res.FE <- suppressMessages(try(lme4::glmer(dat.grp ~ -1 + X.fit + const + (1 | study), offset=dat.off, family=dat.fam, nAGQ=nAGQ, verbose=verbose, control=do.call(lme4::glmerControl, glmerCtrl)), silent=!verbose))
               }
            }

            if (package == "GLMMadaptive") {
               if (is.element(measure, c("OR","RR","RD"))) {
                  dat.mm <- data.frame(xi=dat.grp[,"xi"], mi=dat.grp[,"mi"], study=study, const=const)
                  res.FE <- try(GLMMadaptive::mixed_model(cbind(xi,mi) ~ -1 + X.fit + const, random = ~ 1 | study, data=dat.mm, family=dat.fam, nAGQ=nAGQ, verbose=verbose, control=glmerCtrl), silent=!verbose)
               } else {
                  dat.mm <- data.frame(xi=dat.grp, study=study, const=const)
                  res.FE <- try(GLMMadaptive::mixed_model(xi ~ -1 + X.fit + const + offset(dat.off), random = ~ 1 | study, data=dat.mm, family=dat.fam, nAGQ=nAGQ, verbose=verbose, control=glmerCtrl), silent=!verbose)
               }
            }

            if (package == "glmmTMB") {
               if (verbose) {
                  res.FE <- try(glmmTMB::glmmTMB(dat.grp ~ -1 + X.fit + const + (1 | study), offset=dat.off, family=dat.fam, verbose=verbose, data=NULL, control=do.call(glmmTMB::glmmTMBControl, glmerCtrl)), silent=!verbose)
               } else {
                  res.FE <- suppressMessages(try(glmmTMB::glmmTMB(dat.grp ~ -1 + X.fit + const + (1 | study), offset=dat.off, family=dat.fam, verbose=verbose, data=NULL, control=do.call(glmmTMB::glmmTMBControl, glmerCtrl)), silent=!verbose))
               }
            }

            if (inherits(res.FE, "try-error"))
               stop(mstyle$stop(paste0("Cannot fit FE model", ifelse(verbose, ".", " (set 'verbose=TRUE' to obtain further details)."))))

            ### log-likelihood

            ll.FE <- c(logLik(res.FE))

            ### fit saturated FE model (= QE model)
            ### notes: 1) must remove aliased terms before fitting (for GLMMadaptive to work)
            ###        2) use the sigma^2 value from the FE model as the starting value for the study-level random effect

            QEconv <- FALSE
            ll.QE <- NA_real_

            if (!isTRUE(ddd$skiphet)) {

               if (k > 1 && verbose)
                  message(mstyle$message("Fitting the saturated model ..."))

               if (k > 1) {

                  X.QE   <- model.matrix(~ -1 + X.fit + const + study:group1)
                  res.QE <- try(glm(dat.grp ~ -1 + X.QE, offset=dat.off, family=dat.fam, control=glmCtrl), silent=TRUE)
                  X.QE   <- X.QE[,!is.na(coef(res.QE)),drop=FALSE]

                  if (package == "lme4") {
                     if (verbose) {
                        res.QE <- try(lme4::glmer(dat.grp ~ -1 + X.QE + (1 | study), offset=dat.off, family=dat.fam, start=c(sqrt(lme4::VarCorr(res.FE)[[1]][1])), nAGQ=nAGQ, verbose=verbose, control=do.call(lme4::glmerControl, glmerCtrl)), silent=!verbose)
                     } else {
                        res.QE <- suppressMessages(try(lme4::glmer(dat.grp ~ -1 + X.QE + (1 | study), offset=dat.off, family=dat.fam, start=c(sqrt(lme4::VarCorr(res.FE)[[1]][1])), nAGQ=nAGQ, verbose=verbose, control=do.call(lme4::glmerControl, glmerCtrl)), silent=!verbose))
                     }
                  }

                  if (package == "GLMMadaptive") {
                     glmerCtrl$max_coef_value <- 50
                     if (is.element(measure, c("OR","RR","RD"))) {
                        dat.mm <- data.frame(xi=dat.grp[,"xi"], mi=dat.grp[,"mi"], study=study)
                        res.QE <- try(GLMMadaptive::mixed_model(cbind(xi,mi) ~ -1 + X.QE, random = ~ 1 | study, data=dat.mm, family=dat.fam, nAGQ=nAGQ, verbose=verbose, control=glmerCtrl, initial_values=list(D=matrix(res.FE$D[1,1]))), silent=!verbose)
                     } else {
                        dat.mm <- data.frame(xi=dat.grp, study=study)
                        res.QE <- try(GLMMadaptive::mixed_model(xi ~ -1 + X.QE + offset(dat.off), random = ~ 1 | study, data=dat.mm, family=dat.fam, nAGQ=nAGQ, verbose=verbose, control=glmerCtrl), silent=!verbose)
                     }
                  }

                  if (package == "glmmTMB") {
                     if (verbose) {
                        res.QE <- try(glmmTMB::glmmTMB(dat.grp ~ -1 + X.QE + (1 | study), offset=dat.off, family=dat.fam, start=list(theta=sqrt(glmmTMB::VarCorr(res.FE)[[1]][[1]][[1]])), verbose=verbose, data=NULL, control=do.call(glmmTMB::glmmTMBControl, glmerCtrl)), silent=!verbose)
                     } else {
                        res.QE <- suppressMessages(try(glmmTMB::glmmTMB(dat.grp ~ -1 + X.QE + (1 | study), offset=dat.off, family=dat.fam, start=list(theta=sqrt(glmmTMB::VarCorr(res.FE)[[1]][[1]][[1]])), verbose=verbose, data=NULL, control=do.call(glmmTMB::glmmTMBControl, glmerCtrl)), silent=!verbose))
                     }
                  }

               } else {
                  res.QE <- res.FE
               }

               if (inherits(res.QE, "try-error")) {

                  warning(mstyle$warning(paste0("Cannot fit saturated model", ifelse(verbose, ".", " (set 'verbose=TRUE' to obtain further details)."))), call.=FALSE)

               } else {

                  QEconv <- TRUE

                  ### log-likelihood

                  ll.QE <- c(logLik(res.QE))

                  ### extract coefficients and variance-covariance matrix for Wald-type test for heterogeneity (aliased coefficients are already removed)

                  if (package == "lme4") {
                     b2.QE  <- cbind(lme4::fixef(res.QE)[-seq_len(p+1)])
                     vb2.QE <- as.matrix(vcov(res.QE))[-seq_len(p+1),-seq_len(p+1),drop=FALSE]
                  }

                  if (package == "GLMMadaptive") {
                     b2.QE  <- cbind(GLMMadaptive::fixef(res.QE)[-seq_len(p+1)])
                     vb2.QE <- as.matrix(vcov(res.QE))[-seq_len(p+1),-seq_len(p+1),drop=FALSE]
                     vb2.QE <- vb2.QE[-nrow(vb2.QE), -ncol(vb2.QE)]
                  }

                  if (package == "glmmTMB") {
                     b2.QE  <- cbind(glmmTMB::fixef(res.QE)$cond[-seq_len(p+1)])
                     vb2.QE <- as.matrix(vcov(res.QE)$cond)[-seq_len(p+1),-seq_len(p+1),drop=FALSE]
                  }

               }

            }

            if (method == "ML") {

               ### fit ML model
               ### notes: 1) not recommended alternative: using group1 instead of group12 for the random effect (since that forces the variance in group 2 to be lower)
               ###        2) this approach is okay if we also allow group1 random effect and intercepts to correlate (in fact, this is identical to the bivariate model)
               ###        3) start=c(sqrt(lme4::VarCorr(res.FE)[[1]][1])) has no effect, since the start value for tau^2 is not specified (and using 0 is probably not ideal for that)

               if (verbose)
                  message(mstyle$message("Fitting the ML model ..."))

               if (package == "lme4") {
                  if (verbose) {
                     if (cor) {
                        res.ML <- try(lme4::glmer(dat.grp ~ -1 + X.fit + const + (group | study),  offset=dat.off, family=dat.fam, nAGQ=nAGQ, verbose=verbose, control=do.call(lme4::glmerControl, glmerCtrl)), silent=!verbose)
                     } else {
                        res.ML <- try(lme4::glmer(dat.grp ~ -1 + X.fit + const + (group || study), offset=dat.off, family=dat.fam, nAGQ=nAGQ, verbose=verbose, control=do.call(lme4::glmerControl, glmerCtrl)), silent=!verbose)
                     }
                  } else {
                     if (cor) {
                        res.ML <- suppressMessages(try(lme4::glmer(dat.grp ~ -1 + X.fit + const + (group | study),  offset=dat.off, family=dat.fam, nAGQ=nAGQ, verbose=verbose, control=do.call(lme4::glmerControl, glmerCtrl)), silent=!verbose))
                     } else {
                        res.ML <- suppressMessages(try(lme4::glmer(dat.grp ~ -1 + X.fit + const + (group || study), offset=dat.off, family=dat.fam, nAGQ=nAGQ, verbose=verbose, control=do.call(lme4::glmerControl, glmerCtrl)), silent=!verbose))
                     }
                  }
               }

               if (package == "GLMMadaptive") {
                  if (is.element(measure, c("OR","RR","RD"))) {
                     dat.mm <- data.frame(xi=dat.grp[,"xi"], mi=dat.grp[,"mi"], study=study, const=const, group=group)
                     if (cor) {
                        res.ML <- try(GLMMadaptive::mixed_model(cbind(xi,mi) ~ -1 + X.fit + const, random = ~ group | study,  data=dat.mm, family=dat.fam, nAGQ=nAGQ, verbose=verbose, control=glmerCtrl), silent=!verbose)
                     } else {
                        res.ML <- try(GLMMadaptive::mixed_model(cbind(xi,mi) ~ -1 + X.fit + const, random = ~ group || study, data=dat.mm, family=dat.fam, nAGQ=nAGQ, verbose=verbose, control=glmerCtrl), silent=!verbose)
                     }
                  } else {
                     dat.mm <- data.frame(xi=dat.grp, study=study, const=const, group=group)
                     if (cor) {
                        res.ML <- try(GLMMadaptive::mixed_model(xi ~ -1 + X.fit + const + offset(dat.off), random = ~ group | study,  data=dat.mm, family=dat.fam, nAGQ=nAGQ, verbose=verbose, control=glmerCtrl), silent=!verbose)
                     } else {
                        res.ML <- try(GLMMadaptive::mixed_model(xi ~ -1 + X.fit + const + offset(dat.off), random = ~ group || study, data=dat.mm, family=dat.fam, nAGQ=nAGQ, verbose=verbose, control=glmerCtrl), silent=!verbose)
                     }
                  }
               }

               if (package == "glmmTMB") {
                  if (verbose) {
                     if (cor) {
                        res.ML <- try(glmmTMB::glmmTMB(dat.grp ~ -1 + X.fit + const + (group | study),                   offset=dat.off, family=dat.fam, verbose=verbose, data=NULL, control=do.call(glmmTMB::glmmTMBControl, glmerCtrl)), silent=!verbose)
                     } else {
                        res.ML <- try(glmmTMB::glmmTMB(dat.grp ~ -1 + X.fit + const + (1 | study) + (group - 1 | study), offset=dat.off, family=dat.fam, verbose=verbose, data=NULL, control=do.call(glmmTMB::glmmTMBControl, glmerCtrl)), silent=!verbose)
                     }
                  } else {
                     if (cor) {
                        res.ML <- suppressMessages(try(glmmTMB::glmmTMB(dat.grp ~ -1 + X.fit + const + (group | study),                   offset=dat.off, family=dat.fam, verbose=verbose, data=NULL, control=do.call(glmmTMB::glmmTMBControl, glmerCtrl)), silent=!verbose))
                     } else {
                        res.ML <- suppressMessages(try(glmmTMB::glmmTMB(dat.grp ~ -1 + X.fit + const + (1 | study) + (group - 1 | study), offset=dat.off, family=dat.fam, verbose=verbose, data=NULL, control=do.call(glmmTMB::glmmTMBControl, glmerCtrl)), silent=!verbose))
                     }
                  }
               }

               if (inherits(res.ML, "try-error"))
                  stop(mstyle$stop(paste0("Cannot fit ML model", ifelse(verbose, ".", " (set 'verbose=TRUE' to obtain further details)."))))

               ### log-likelihood

               ll.ML <- c(logLik(res.ML))

            }

            #return(list(res.FE, res.QE, res.ML, ll.FE=ll.FE, ll.QE=ll.QE, ll.ML=ll.ML))
            #res.FE <- res[[1]]; res.QE <- res[[2]]; res.ML <- res[[3]]

            if (is.element(method, c("FE","EE","CE"))) {

               tau2 <- 0

               if (package == "lme4") {
                  beta   <- cbind(lme4::fixef(res.FE)[seq_len(p)])
                  vb     <- as.matrix(vcov(res.FE))[seq_len(p),seq_len(p),drop=FALSE]
                  sigma2 <- lme4::VarCorr(res.FE)[[1]][1]
               }

               if (package == "GLMMadaptive") {
                  beta   <- cbind(GLMMadaptive::fixef(res.FE)[seq_len(p)])
                  vb     <- as.matrix(vcov(res.FE))[seq_len(p),seq_len(p),drop=FALSE]
                  sigma2 <- res.FE$D[1,1]
               }

               if (package == "glmmTMB") {
                  beta   <- cbind(glmmTMB::fixef(res.FE)$cond[seq_len(p)])
                  vb     <- as.matrix(vcov(res.FE)$cond)[seq_len(p),seq_len(p),drop=FALSE]
                  sigma2 <- glmmTMB::VarCorr(res.FE)[[1]][[1]][[1]]
               }

               parms  <- p + 1 + 1
               p.eff  <- p + 1
               k.eff  <- 2*k
            }

            if (method == "ML") {

               if (package == "lme4") {
                  beta   <- cbind(lme4::fixef(res.ML)[seq_len(p)])
                  vb     <- as.matrix(vcov(res.ML))[seq_len(p),seq_len(p),drop=FALSE]
                  if (cor) {
                     tau2   <- lme4::VarCorr(res.ML)[[1]][2,2]
                     sigma2 <- lme4::VarCorr(res.ML)[[1]][1,1]
                     rho    <- lme4::VarCorr(res.ML)[[1]][1,2] / sqrt(tau2 * sigma2)
                  } else {
                     tau2   <- lme4::VarCorr(res.ML)[[2]][1]
                     sigma2 <- lme4::VarCorr(res.ML)[[1]][1]
                  }
               }

               if (package == "GLMMadaptive") {
                  beta   <- cbind(GLMMadaptive::fixef(res.ML)[seq_len(p)])
                  vb     <- as.matrix(vcov(res.ML))[seq_len(p),seq_len(p),drop=FALSE]
                  tau2   <- res.ML$D[2,2]
                  sigma2 <- res.ML$D[1,1]
                  if (cor)
                     rho <- res.ML$D[1,2] / sqrt(tau2 * sigma2)
               }

               if (package == "glmmTMB") {
                  beta   <- cbind(glmmTMB::fixef(res.ML)$cond[seq_len(p)])
                  vb     <- as.matrix(vcov(res.ML)$cond)[seq_len(p),seq_len(p),drop=FALSE]
                  if (cor) {
                     tau2   <- glmmTMB::VarCorr(res.ML)[[1]][[1]][2,2]
                     sigma2 <- glmmTMB::VarCorr(res.ML)[[1]][[1]][1,1]
                     rho    <- glmmTMB::VarCorr(res.ML)[[1]][[1]][1,2] / sqrt(tau2 * sigma2)
                  } else {
                     tau2   <- glmmTMB::VarCorr(res.ML)[[1]][[2]][[1]]
                     sigma2 <- glmmTMB::VarCorr(res.ML)[[1]][[1]][[1]]
                  }
               }

               parms  <- p + 1 + 2
               p.eff  <- p + 1
               k.eff  <- 2*k

            }

            #return(list(beta=beta, vb=vb, tau2=tau2, sigma2=sigma2, parms=parms, p.eff=p.eff, k.eff=k.eff, b2.QE=b2.QE, vb2.QE=vb2.QE))

         }

         ###################################################################

      }

      ######################################################################

      if ((measure=="IRR" && model == "CM.EL") || (measure=="OR" && model=="CM.AL") || (measure=="OR" && model=="CM.EL")) {

         ### prepare data for the conditional models

         if (measure == "OR") {
            dat.grp <- cbind(xi=ai, mi=ci)   # conditional outcome data (number of cases in group 1 conditional on total number of cases)
            dat.off <- log((ai+bi)/(ci+di))  # log(n1i/n2i) for offset
         }

         if (measure == "IRR") {
            dat.grp <- cbind(xi=x1i, mi=x2i) # conditional outcome data (number of events in group 1 conditional on total number of events)
            dat.off <- log(t1i/t2i)          # log(t1i/t1i) for offset
         }

         study <- factor(seq_len(k))         # study factor
         X.fit <- X

         if (isTRUE(ddd$retdat))
            return(list(dat.grp=dat.grp, X.fit=X.fit, study=study, dat.off = if (!is.null(dat.off)) dat.off else NULL))

         ###################################################################

         ###############################################################
         ### conditional model (approx. ll for ORs / exact for IRRs) ###
         ###############################################################

         ### fit FE model

         if (verbose)
            message(mstyle$message("Fitting the FE model ..."))

         res.FE <- try(glm(dat.grp ~ -1 + X.fit, offset=dat.off, family=binomial, control=glmCtrl), silent=!verbose)

         if (inherits(res.FE, "try-error"))
            stop(mstyle$stop(paste0("Cannot fit FE model", ifelse(verbose, ".", " (set 'verbose=TRUE' to obtain further details)."))))

         ### log-likelihood

         #ll.FE <- with(data.frame(dat.grp), sum(dbinom(xi, xi+mi, predict(res.FE, type="response"), log=TRUE))) # offset already incorporated into predict()
         #ll.FE <- with(data.frame(dat.grp), sum(dpois(xi, predict(res.FE, type="response"), log=TRUE)))         # offset already incorporated into predict()
         ll.FE  <- c(logLik(res.FE)) # same as above

         ### fit saturated FE model (= QE model)

         QEconv <- FALSE
         ll.QE <- NA_real_

         if (!isTRUE(ddd$skiphet)) {

            if (k > 1 && verbose)
               message(mstyle$message("Fitting the saturated model ..."))

            if (k > 1) {
               X.QE   <- model.matrix(~ -1 + X.fit + study)
               res.QE <- try(glm(dat.grp ~ -1 + X.QE, offset=dat.off, family=binomial, control=glmCtrl), silent=!verbose)
            } else {
               res.QE <- res.FE
            }

            if (inherits(res.QE, "try-error")) {

               warning(mstyle$warning(paste0("Cannot fit saturated model", ifelse(verbose, ".", " (set 'verbose=TRUE' to obtain further details)."))), call.=FALSE)

            } else {

               QEconv <- TRUE

               ### log-likelihood

               #ll.QE <- with(data.frame(dat.grp), sum(dbinom(xi, xi+mi, xi/(xi+mi), log=TRUE))) # offset not relevant for saturated model
               #ll.QE <- with(data.frame(dat.grp), sum(dpois(xi, xi, log=TRUE)))                 # offset not relevant for saturated model
               ll.QE  <- c(logLik(res.QE)) # same as above

               ### extract coefficients and variance-covariance matrix for Wald-type test for heterogeneity

               #b2.QE <- cbind(na.omit(coef(res.QE)[-seq_len(p)]))                        # coef() still includes aliased coefficients as NAs, so filter them out
               b2.QE  <- cbind(coef(res.QE, complete=FALSE)[-seq_len(p)])                 # aliased coefficients are removed by coef() when complete=FALSE
               vb2.QE <- vcov(res.QE, complete=FALSE)[-seq_len(p),-seq_len(p),drop=FALSE] # aliased coefficients are removed by vcov() when complete=FALSE

            }

            #return(list(res.FE, res.QE, ll.FE, ll.QE))
            #res.FE <- res[[1]]; res.QE <- res[[2]]

         }

         if (method == "ML") {

            ### fit ML model
            ### notes: 1) suppressMessages to suppress the 'one random effect per observation' warning

            if (verbose)
               message(mstyle$message("Fitting the ML model ..."))

            if (package == "lme4") {
               if (verbose) {
                  res.ML <- try(lme4::glmer(dat.grp ~ -1 + X.fit + (1 | study), offset=dat.off, family=binomial, nAGQ=nAGQ, verbose=verbose, control=do.call(lme4::glmerControl, glmerCtrl)), silent=!verbose)
               } else {
                  res.ML <- suppressMessages(try(lme4::glmer(dat.grp ~ -1 + X.fit + (1 | study), offset=dat.off, family=binomial, nAGQ=nAGQ, verbose=verbose, control=do.call(lme4::glmerControl, glmerCtrl)), silent=!verbose))
               }
            }

            if (package == "GLMMadaptive") {
               dat.mm <- data.frame(xi=dat.grp[,"xi"], mi=dat.grp[,"mi"], study=study)
               res.ML <- try(GLMMadaptive::mixed_model(cbind(xi,mi) ~ -1 + X.fit + offset(dat.off), random = ~ 1 | study, data=dat.mm, family=binomial, nAGQ=nAGQ, verbose=verbose, control=glmerCtrl), silent=!verbose)
            }

            if (package == "glmmTMB") {
               if (verbose) {
                  res.ML <- try(glmmTMB::glmmTMB(dat.grp ~ -1 + X.fit + (1 | study), offset=dat.off, family=binomial, verbose=verbose, data=NULL, control=do.call(glmmTMB::glmmTMBControl, glmerCtrl)), silent=!verbose)
               } else {
                  res.ML <- suppressMessages(try(glmmTMB::glmmTMB(dat.grp ~ -1 + X.fit + (1 | study), offset=dat.off, family=binomial, verbose=verbose, data=NULL, control=do.call(glmmTMB::glmmTMBControl, glmerCtrl)), silent=!verbose))
               }
            }

            if (inherits(res.ML, "try-error"))
               stop(mstyle$stop(paste0("Cannot fit ML model", ifelse(verbose, ".", " (set 'verbose=TRUE' to obtain further details)."))))

            ### log-likelihood

            if (package == "lme4") {
               if (is.na(ll.QE)) {
                  ll.ML <- c(logLik(res.ML))
               } else {
                  if (verbose) {
                     ll.ML <- ll.QE - 1/2 * deviance(res.ML) # this makes ll.ML comparable to ll.FE (same as ll.FE when tau^2=0)
                  } else {
                     ll.ML <- ll.QE - 1/2 * suppressWarnings(deviance(res.ML)) # suppressWarnings() to suppress 'Warning in sqrt(object$devResid()) : NaNs produced'
                  }
               }
            } else {
               ll.ML <- c(logLik(res.ML)) # not 100% sure how comparable this is to ll.FE when tau^2 = 0 (seems correct for glmmTMB)
            }

         }

         #return(list(res.FE, res.QE, res.ML, ll.FE=ll.FE, ll.QE=ll.QE, ll.ML=ll.ML))
         #res.FE <- res[[1]]; res.QE <- res[[2]]; res.ML <- res[[3]]

         if (is.element(method, c("FE","EE","CE"))) {
            beta   <- cbind(coef(res.FE)[seq_len(p)])
            vb     <- vcov(res.FE)[seq_len(p),seq_len(p),drop=FALSE]
            tau2   <- 0
            sigma2 <- NA_real_
            parms  <- p
            p.eff  <- p
            k.eff  <- k
         }

         if (method == "ML") {

            if (package == "lme4") {
               beta   <- cbind(lme4::fixef(res.ML)[seq_len(p)])
               vb     <- as.matrix(vcov(res.ML))[seq_len(p),seq_len(p),drop=FALSE]
               tau2   <- lme4::VarCorr(res.ML)[[1]][1]
            }

            if (package == "GLMMadaptive") {
               beta   <- cbind(GLMMadaptive::fixef(res.ML)[seq_len(p)])
               vb     <- as.matrix(vcov(res.ML))[seq_len(p),seq_len(p),drop=FALSE]
               tau2   <- res.ML$D[1,1]
            }

            if (package == "glmmTMB") {
               beta   <- cbind(glmmTMB::fixef(res.ML)$cond[seq_len(p)])
               vb     <- as.matrix(vcov(res.ML)$cond)[seq_len(p),seq_len(p),drop=FALSE]
               tau2   <- glmmTMB::VarCorr(res.ML)[[1]][[1]][[1]]
            }

            sigma2 <- NA_real_
            parms  <- p + 1
            p.eff  <- p
            k.eff  <- k

         }

         #return(list(beta=beta, vb=vb, tau2=tau2, sigma2=sigma2, parms=parms, p.eff=p.eff, k.eff=k.eff, b2.QE=b2.QE, vb2.QE=vb2.QE))

         ###################################################################

      }

      if (measure=="OR" && model=="CM.EL") {

         ####################################################
         ### conditional model (exact likelihood for ORs) ###
         ####################################################

         if (verbose)
            message(mstyle$message("Fitting the FE model ..."))

         if (is.element(optimizer, c("optim","nlminb","uobyqa","newuoa","bobyqa","nloptr","nlm","hjk","nmk","mads","ucminf","lbfgsb3c","subplex","BBoptim","optimParallel","Rcgmin","Rvmmin"))) {

            if (optimizer == "optim") {
               par.arg <- "par"
               ctrl.arg <- ", control=optCtrl"
            }

            if (optimizer == "nlminb") {
               par.arg <- "start"
               ctrl.arg <- ", control=optCtrl"
            }

            if (is.element(optimizer, c("uobyqa","newuoa","bobyqa"))) {
               par.arg <- "par"
               optimizer <- paste0("minqa::", optimizer)
               ctrl.arg <- ", control=optCtrl"
            }

            if (optimizer == "nloptr") {
               par.arg <- "x0"
               optimizer <- paste0("nloptr::nloptr")
               ctrl.arg <- ", opts=optCtrl"
            }

            if (optimizer == "nlm") {
               par.arg <- "p"
               ctrl.arg <- paste(names(optCtrl), unlist(optCtrl), sep="=", collapse=", ")
               if (nchar(ctrl.arg) != 0L)
                  ctrl.arg <- paste0(", ", ctrl.arg)
            }

            if (is.element(optimizer, c("hjk","nmk","mads"))) {
               par.arg <- "par"
               optimizer <- paste0("dfoptim::", optimizer)
               ctrl.arg <- ", control=optCtrl"
            }

            if (is.element(optimizer, c("ucminf","lbfgsb3c","subplex"))) {
               par.arg <- "par"
               optimizer <- paste0(optimizer, "::", optimizer)
               ctrl.arg <- ", control=optCtrl"
            }

            if (optimizer == "BBoptim") {
               par.arg <- "par"
               optimizer <- "BB::BBoptim"
               ctrl.arg <- ", quiet=TRUE, control=optCtrl"
            }

            if (optimizer == "Rcgmin") {
               par.arg <- "par"
               optimizer <- "optimx::Rcgmin"
               ctrl.arg <- ", gr='grnd', control=optCtrl"
               #ctrl.arg <- ", control=optCtrl"
            }

            if (optimizer == "Rvmmin") {
               par.arg <- "par"
               optimizer <- "optimx::Rvmmin"
               ctrl.arg <- ", gr='grnd', control=optCtrl"
               #ctrl.arg <- ", control=optCtrl"
            }

            if (optimizer == "optimParallel") {

               par.arg <- "par"
               optimizer <- paste0("optimParallel::optimParallel")
               ctrl.arg <- ", control=optCtrl, parallel=parallel"

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

            ### fit FE model
            ### notes: 1) this routine uses direct optimization over the non-central hypergeometric distribution
            ###        2) start values from CM.AL model (res.FE) and tau^2=0 (random=FALSE)
            ###        3) no integration needed for FE model, so intCtrl is not actually relevant
            ###        4) results can be sensitive to the scaling of moderators

            optcall <- paste0(optimizer, "(", par.arg, "=c(coef(res.FE)[seq_len(p)], 0),
               .dnchg, ", ifelse(optimizer=="optim", "method=optmethod, ", ""),
               "ai=ai, bi=bi, ci=ci, di=di, X.fit=X.fit, random=FALSE,
               verbose=verbose, digits=digits,
               dnchgcalc=con$dnchgcalc, dnchgprec=con$dnchgprec, intCtrl=intCtrl", ctrl.arg, ")\n")

            #return(optcall)

            if (verbose) {
               res.FE <- try(eval(str2lang(optcall)), silent=!verbose)
            } else {
               res.FE <- try(suppressWarnings(eval(str2lang(optcall))), silent=!verbose)
            }

            #return(res.FE)

            if (optimizer == "optimParallel::optimParallel" && verbose) {
               tmp <- capture.output(print(res.FE$loginfo))
               .print.output(tmp, mstyle$verbose)
            }

            if (inherits(res.FE, "try-error"))
               stop(mstyle$stop(paste0("Cannot fit FE model", ifelse(verbose, ".", " (set 'verbose=TRUE' to obtain further details)."))))

            ### convergence checks

            if (is.element(optimizer, c("optim","nlminb","dfoptim::hjk","dfoptim::nmk","lbfgsb3c::lbfgsb3c","subplex::subplex","BB::BBoptim","optimx::Rcgmin","optimx::Rvmmin","optimParallel::optimParallel")) && res.FE$convergence != 0)
               stop(mstyle$stop(paste0("Cannot fit FE model. Optimizer (", optimizer, ") did not achieve convergence (convergence = ", res.FE$convergence, ").")))

            if (is.element(optimizer, c("dfoptim::mads")) && res.FE$convergence > optCtrl$tol)
               stop(mstyle$stop(paste0("Cannot fit FE model. Optimizer (", optimizer, ") did not achieve convergence (convergence = ", res.FE$convergence, ").")))

            if (is.element(optimizer, c("minqa::uobyqa","minqa::newuoa","minqa::bobyqa")) && res.FE$ierr != 0)
               stop(mstyle$stop(paste0("Cannot fit FE model. Optimizer (", optimizer, ") did not achieve convergence (ierr = ", res.FE$ierr, ").")))

            if (optimizer=="nloptr::nloptr" && !(res.FE$status >= 1 && res.FE$status <= 4))
               stop(mstyle$stop(paste0("Cannot fit FE model. Optimizer (", optimizer, ") did not achieve convergence (status = ", res.FE$status, ").")))

            if (optimizer=="ucminf::ucminf" && !(res.FE$convergence == 1 || res.FE$convergence == 2))
               stop(mstyle$stop(paste0("Cannot fit FE model. Optimizer (", optimizer, ") did not achieve convergence (convergence = ", res.FE$convergence, ").")))

            if (verbose > 2) {
               cat("\n")
               tmp <- capture.output(print(res.FE))
               .print.output(tmp, mstyle$verbose)
            }

            ### copy estimated values to 'par'

            if (optimizer=="nloptr::nloptr")
               res.FE$par <- res.FE$solution
            if (optimizer=="nlm")
               res.FE$par <- res.FE$estimate

            res.FE$par <- unname(res.FE$par)

            if (verbose > 1)
               message(mstyle$message("Computing the Hessian ..."))

            if (con$hesspack == "numDeriv")
               h.FE <- numDeriv::hessian(.dnchg, x=res.FE$par, method.args=hessianCtrl, ai=ai, bi=bi, ci=ci, di=di, X.fit=X.fit, random=FALSE, verbose=verbose, digits=digits, dnchgcalc=con$dnchgcalc, dnchgprec=con$dnchgprec)
            if (con$hesspack == "pracma")
               h.FE <- pracma::hessian(.dnchg, x0=res.FE$par, ai=ai, bi=bi, ci=ci, di=di, X.fit=X.fit, random=FALSE, verbose=verbose, digits=digits, dnchgcalc=con$dnchgcalc, dnchgprec=con$dnchgprec)
            if (con$hesspack == "calculus")
               h.FE <- calculus::hessian(.dnchg, var=res.FE$par, params=list(ai=ai, bi=bi, ci=ci, di=di, X.fit=X.fit, random=FALSE, verbose=verbose, digits=digits, dnchgcalc=con$dnchgcalc, dnchgprec=con$dnchgprec))
            #return(list(res.FE=res.FE, h.FE=h.FE))

            ### log-likelihood

            if (is.element(optimizer, c("optim","dfoptim::hjk","dfoptim::nmk","dfoptim::mads","ucminf::ucminf","lbfgsb3c::lbfgsb3c","subplex::subplex","BB::BBoptim","optimx::Rcgmin","optimx::Rvmmin","optimParallel::optimParallel")))
               ll.FE <- -1 * res.FE$value
            if (is.element(optimizer, c("nlminb","nloptr::nloptr")))
               ll.FE <- -1 * res.FE$objective
            if (is.element(optimizer, c("minqa::uobyqa","minqa::newuoa","minqa::bobyqa")))
               ll.FE <- -1 * res.FE$fval
            if (optimizer == "nlm")
               ll.FE <- -1 * res.FE$minimum

            ### fit saturated FE model (= QE model)
            ### notes: 1) must figure out which terms are aliased in saturated model and remove those terms before fitting
            ###        2) start values from CM.AL model (res.QE) and tau^2=0 (random=FALSE)
            ###        3) so only try to fit saturated model if this was possible with CM.AL
            ###        4) no integration needed for FE model, so intCtrl is not relevant

            if (QEconv) { # QEconv is FALSE when skiphet=TRUE so this then also gets skipped automatically

               if (k > 1 && verbose)
                  message(mstyle$message("Fitting the saturated model ..."))

               if (k > 1) {

                  b.QE <- coef(res.QE, complete=TRUE) # res.QE is from CM.AL model
                  is.aliased <- is.na(b.QE)
                  b.QE <- b.QE[!is.aliased]
                  X.QE <- X.QE[,!is.aliased,drop=FALSE]

                  optcall <- paste0(optimizer, "(", par.arg, "=c(b.QE, 0),
                     .dnchg, ", ifelse(optimizer=="optim", "method=optmethod, ", ""),
                     "ai=ai, bi=bi, ci=ci, di=di, X.fit=X.QE, random=FALSE,
                     verbose=verbose, digits=digits,
                     dnchgcalc=con$dnchgcalc, dnchgprec=con$dnchgprec, intCtrl=intCtrl", ctrl.arg, ")\n")

                  #return(optcall)

                  if (verbose) {
                     res.QE <- try(eval(str2lang(optcall)), silent=!verbose)
                  } else {
                     res.QE <- try(suppressWarnings(eval(str2lang(optcall))), silent=!verbose)
                  }

                  #return(res.QE)

                  if (optimizer == "optimParallel::optimParallel" && verbose) {
                     tmp <- capture.output(print(res.QE$loginfo))
                     .print.output(tmp, mstyle$verbose)
                  }

                  if (inherits(res.QE, "try-error")) {
                     warning(mstyle$warning(paste0("Cannot fit saturated model", ifelse(verbose, ".", " (set 'verbose=TRUE' to obtain further details)."))), call.=FALSE)
                     QEconv <- FALSE
                     ll.QE <- NA_real_
                  }

                  ### convergence checks

                  if (QEconv && is.element(optimizer, c("optim","nlminb","dfoptim::hjk","dfoptim::nmk","lbfgsb3c::lbfgsb3c","subplex::subplex","BB::BBoptim","optimx::Rcgmin","optimx:Rvmmin","optimParallel::optimParallel")) && res.QE$convergence != 0) {
                     warning(mstyle$warning(paste0("Cannot fit saturated model. Optimizer (", optimizer, ") did not achieve convergence (convergence = ", res.QE$convergence, ").")), call.=FALSE)
                     QEconv <- FALSE
                     ll.QE <- NA_real_
                  }

                  if (QEconv && is.element(optimizer, c("dfoptim::mads")) && res.QE$convergence > optCtrl$tol) {
                     warning(mstyle$warning(paste0("Cannot fit saturated model. Optimizer (", optimizer, ") did not achieve convergence (convergence = ", res.QE$convergence, ").")), call.=FALSE)
                     QEconv <- FALSE
                     ll.QE <- NA_real_
                  }

                  if (QEconv && is.element(optimizer, c("minqa::uobyqa","minqa::newuoa","minqa::bobyqa")) && res.QE$ierr != 0) {
                     warning(mstyle$warning(paste0("Cannot fit saturated model. Optimizer (", optimizer, ") did not achieve convergence (ierr = ", res.QE$ierr, ").")), call.=FALSE)
                     QEconv <- FALSE
                     ll.QE <- NA_real_
                  }

                  if (QEconv && optimizer=="nloptr::nloptr" && !(res.QE$status >= 1 && res.QE$status <= 4)) {
                     warning(mstyle$warning(paste0("Cannot fit saturated model. Optimizer (", optimizer, ") did not achieve convergence (status = ", res.QE$status, ").")), call.=FALSE)
                     QEconv <- FALSE
                     ll.QE <- NA_real_
                  }

                  if (QEconv && optimizer=="ucminf::ucminf" && !(res.QE$convergence == 1 || res.QE$convergence == 2)) {
                     warning(mstyle$warning(paste0("Cannot fit saturated model. Optimizer (", optimizer, ") did not achieve convergence (convergence = ", res.QE$convergence, ").")), call.=FALSE)
                     QEconv <- FALSE
                     ll.QE <- NA_real_
                  }

                  if (verbose > 2) {
                     cat("\n")
                     tmp <- capture.output(print(res.QE))
                     .print.output(tmp, mstyle$verbose)
                  }

                  ### copy estimated values to 'par'

                  if (QEconv && optimizer=="nloptr::nloptr")
                     res.QE$par <- res.QE$solution
                  if (QEconv && optimizer=="nlm")
                     res.QE$par <- res.QE$estimate

                  res.QE$par <- unname(res.QE$par)

                  if (QEconv) {
                     if (verbose > 1)
                        message(mstyle$message("Computing the Hessian ..."))
                     if (con$hesspack == "numDeriv")
                        h.QE <- numDeriv::hessian(.dnchg, x=res.QE$par, method.args=hessianCtrl, ai=ai, bi=bi, ci=ci, di=di, X.fit=X.QE, random=FALSE, verbose=verbose, digits=digits, dnchgcalc=con$dnchgcalc, dnchgprec=con$dnchgprec)
                     if (con$hesspack == "pracma")
                        h.QE <- pracma::hessian(.dnchg, x0=res.QE$par, ai=ai, bi=bi, ci=ci, di=di, X.fit=X.QE, random=FALSE, verbose=verbose, digits=digits, dnchgcalc=con$dnchgcalc, dnchgprec=con$dnchgprec)
                     if (con$hesspack == "calculus")
                        h.QE <- calculus::hessian(.dnchg, var=res.QE$par, params=list(ai=ai, bi=bi, ci=ci, di=di, X.fit=X.QE, random=FALSE, verbose=verbose, digits=digits, dnchgcalc=con$dnchgcalc, dnchgprec=con$dnchgprec))
                  }

               } else {
                  res.QE <- res.FE
                  h.QE   <- h.FE
               }

               #return(list(res.QE, h.QE))

            }

            if (k > 1 && QEconv) {

               ### log-likelihood

               if (is.element(optimizer, c("optim","dfoptim::hjk","dfoptim::nmk","dfoptim::mads","ucminf::ucminf","lbfgsb3c::lbfgsb3c","subplex::subplex","BB::BBoptim","optimx::Rcgmin","optimx::Rvmmin","optimParallel::optimParallel")))
                  ll.QE <- -1 * res.QE$value
               if (is.element(optimizer, c("nlminb","nloptr::nloptr")))
                  ll.QE <- -1 * res.QE$objective
               if (is.element(optimizer, c("minqa::uobyqa","minqa::newuoa","minqa::bobyqa")))
                  ll.QE <- -1 * res.QE$fval
               if (optimizer == "nlm")
                  ll.QE <- -1 * res.QE$minimum

               ### extract coefficients and variance-covariance matrix for Wald-type test for heterogeneity

               #return(res.QE)
               b2.QE    <- res.QE$par                                  # recall: aliased coefficients are already removed
               hessian  <- h.QE                                        # take hessian from hessian() (again, aliased coefs are already removed)
               #hessian <- res.QE$hessian                              # take hessian from optim() (again, aliased coefs are already removed)
               p.QE     <- length(b2.QE)                               # how many parameters are left in saturated model?
               b2.QE    <- b2.QE[-p.QE]                                # remove last element (for tau^2, constrained to 0)
               hessian  <- hessian[-p.QE,-p.QE,drop=FALSE]             # remove last row/column (for tau^2, constrained to 0)
               p.QE     <- length(b2.QE)                               # how many parameters are now left?
               is.0     <- colSums(hessian == 0L) == p.QE              # any columns in hessian entirely composed of 0s?
               b2.QE    <- b2.QE[!is.0]                                # keep coefficients where this is not the case
               hessian  <- hessian[!is.0,!is.0,drop=FALSE]             # keep parts of hessian where this is not the case
               b2.QE    <- cbind(b2.QE[-seq_len(p)])                   # remove first p coefficients
               h.A      <- hessian[seq_len(p),seq_len(p),drop=FALSE]   # upper left part of hessian
               h.B      <- hessian[seq_len(p),-seq_len(p),drop=FALSE]  # upper right part of hessian
               h.C      <- hessian[-seq_len(p),seq_len(p),drop=FALSE]  # lower left part of hessian
               h.D      <- hessian[-seq_len(p),-seq_len(p),drop=FALSE] # lower right part of hessian (of which we need the inverse)
               chol.h.A <- try(chol(h.A), silent=!verbose)             # see if h.A can be inverted with chol()
               if (inherits(chol.h.A, "try-error") || anyNA(chol.h.A)) {
                  warning(mstyle$warning("Cannot invert the Hessian for the saturated model."), call.=FALSE)
                  QE.Wld <- NA_real_
               } else {
                  Ivb2.QE  <- h.D-h.C%*%chol2inv(chol.h.A)%*%h.B       # inverse of the inverse of the lower right part
                  QE.Wld   <- c(t(b2.QE) %*% Ivb2.QE %*% b2.QE)        # Wald statistic (note: this approach only requires taking the inverse of h.A)
               }                                                       # see: https://en.wikipedia.org/wiki/Invertible_matrix#Blockwise_inversion

               #vb2.QE <- chol2inv(chol(hessian))[-seq_len(p),-seq_len(p),drop=FALSE] # take inverse, then take part relevant for QE test
               #QE.Wld <- c(t(b2.QE) %*% chol2inv(chol(vb2.QE)) %*% b2.QE)

            }

         }

         if (is.element(optimizer, c("clogit","clogistic"))) {

            ### fit FE model
            ### notes: 1) this routine uses either clogit() from the survival package or clogistic() from the Epi package
            ###        2) the dataset must be in group-level and IPD format (i.e., not in the conditional format)
            ###        3) if the studies are large, the IPD dataset may also be very large, and R may run out of memory
            ###        4) for larger datasets, run time is often excessive (and may essentially freeze R)
            ###        5) suppressMessages for clogit() to suppress the 'beta may be infinite' warning

            ### prepare IPD dataset                                                                                                      #                   study  event  group1  intrcpt  moderator
                                                                                                                                         #                   i      1      1       1        x1i       (repeated ai times)
            event   <- unlist(lapply(seq_len(k), function(i) c(rep.int(1,ai[i]), rep.int(0,bi[i]), rep.int(1,ci[i]), rep.int(0,di[i])))) # event dummy       i      0      1       1        x1i       (repeated bi times)
            group1  <- unlist(lapply(seq_len(k), function(i) c(rep.int(1,ai[i]), rep.int(1,bi[i]), rep.int(0,ci[i]), rep.int(0,di[i])))) # group1 dummy      i      1      0       0        0         (repeated ci times)
            study.l <- factor(rep(seq_len(k), times=ni))        # study factor                                                                               i      0      0       0        0         (repeated di times)
            X.fit.l <- X[rep(seq_len(k), times=ni),,drop=FALSE] # repeat each row in X ni times each
            X.fit.l <- cbind(group1*X.fit.l)                    # multiply by group1 dummy (including intercept, which becomes the group1 dummy)
            const   <- rep(1,length(event))

            if (isTRUE(ddd$retdat))
               return(data.frame(event, group1, study.l, X.fit.l, const))

            ### fit FE model

            if (k > 1) {

               if (optimizer == "clogit") {
                  args.clogit <- clogitCtrl
                  args.clogit$formula <- event ~ X.fit.l + strata(study.l)
                  res.FE <- try(do.call(clogit, args.clogit), silent=!verbose)
               }
               if (optimizer == "clogistic") {
                  args.clogistic <- clogisticCtrl
                  args.clogistic$formula <- event ~ X.fit.l
                  args.clogistic$strata <- study.l
                  res.FE <- try(do.call(Epi::clogistic, args.clogistic), silent=!verbose)
               }

            } else {
               stop(mstyle$stop(paste0("Cannot use '", optimizer, "' optimizer when k=1.")))
            }

            if (inherits(res.FE, "try-error"))
               stop(mstyle$stop(paste0("Cannot fit FE model", ifelse(verbose, ".", " (set 'verbose=TRUE' to obtain further details)."))))

            ### fit saturated FE model (= QE model)
            ### notes: 1) must figure out which terms are aliased in saturated model and remove those terms before fitting
            ###        2) fixed effects part does not include 'study' factor, since this is incorporated into the strata
            ###        3) however, for calculating the log-likelihood, we need to go back to the conditional data, so we need to reconstruct X.QE (the study.l:group1 coefficients are the study coefficients)

            if (QEconv) { # QEconv is FALSE when skiphet=TRUE so this then also gets skipped automatically

               if (verbose)
                  message(mstyle$message("Fitting the saturated model ..."))

               b.QE <- coef(res.QE, complete=TRUE) # res.QE is from CM.AL model
               is.aliased <- is.na(b.QE)

               X.QE.l <- model.matrix(~ -1 + X.fit.l + study.l:group1)
               X.QE.l <- X.QE.l[,!is.aliased,drop=FALSE]
               X.QE   <- X.QE[,!is.aliased,drop=FALSE]

               if (optimizer == "clogit") {
                  args.clogit <- clogitCtrl
                  args.clogit$formula <- event ~ X.QE.l + strata(study.l)
                  #args.clogit$method <- "efron" # c("exact", "approximate", "efron", "breslow")
                  if (verbose) {
                     res.QE <- try(do.call(clogit, args.clogit), silent=!verbose)
                  } else {
                     res.QE <- try(suppressWarnings(do.call(clogit, args.clogit)), silent=!verbose)
                  }
               }

               if (optimizer == "clogistic") {
                  args.clogistic <- clogisticCtrl
                  args.clogistic$formula <- event ~ X.QE.l
                  args.clogistic$strata <- study.l
                  res.QE <- try(do.call(Epi::clogistic, args.clogistic), silent=!verbose)
               }

               if (inherits(res.QE, "try-error"))
                  stop(mstyle$stop(paste0("Cannot fit saturated model", ifelse(verbose, ".", " (set 'verbose=TRUE' to obtain further details)."))))

               ### log-likelihood

               ll.FE <- -1 * .dnchg(c(cbind(coef(res.FE)),0), ai=ai, bi=bi, ci=ci, di=di, X.fit=X.fit, random=FALSE, verbose=verbose, digits=digits, dnchgcalc=con$dnchgcalc, dnchgprec=con$dnchgprec)
               ll.QE <- -1 * .dnchg(c(cbind(coef(res.QE)),0), ai=ai, bi=bi, ci=ci, di=di, X.fit=X.QE,  random=FALSE, verbose=verbose, digits=digits, dnchgcalc=con$dnchgcalc, dnchgprec=con$dnchgprec)

               ### extract coefficients and variance-covariance matrix for Wald-type test for heterogeneity

               b2.QE  <- cbind(coef(res.QE)[-seq_len(p)])                 # aliased coefficients are already removed
               vb2.QE <- vcov(res.QE)[-seq_len(p),-seq_len(p),drop=FALSE] # aliased coefficients are already removed

            }

         }

         #return(list(res.FE, res.QE, ll.FE=ll.FE, ll.QE=ll.QE))
         #res.FE <- res[[1]]; res.QE <- res[[2]]

         if (method == "ML") {

            ### fit ML model
            ### notes: 1) cannot use clogit() or clogistic() for this (do not allow for the addition of random effects)
            ###        2) mclogit() from mclogit package may be an alternative (but it only provides a PQL method)
            ###        3) start values from CM.AL model (add 0.01 to tau^2 estimate, in case estimate of tau^2 is 0)
            ###        4) optimization involves integration, so intCtrl is relevant
            ###        5) results can be sensitive to the scaling of moderators

            if (verbose)
               message(mstyle$message("Fitting the ML model ..."))

            optcall <- paste0(optimizer, "(", par.arg, "=c(beta, log(tau2+0.01)),
               .dnchg, ", ifelse(optimizer=="optim", "method=optmethod, ", ""),
               "ai=ai, bi=bi, ci=ci, di=di, X.fit=X.fit, random=TRUE,
               verbose=verbose, digits=digits,
               dnchgcalc=con$dnchgcalc, dnchgprec=con$dnchgprec, intCtrl=intCtrl", ctrl.arg, ")\n")

            #return(optcall)

            if (verbose) {
               res.ML <- try(eval(str2lang(optcall)), silent=!verbose)
            } else {
               res.ML <- try(suppressWarnings(eval(str2lang(optcall))), silent=!verbose)
            }

            #return(res.ML)

            if (optimizer == "optimParallel::optimParallel" && verbose) {
               tmp <- capture.output(print(res.ML$loginfo))
               .print.output(tmp, mstyle$verbose)
            }

            if (inherits(res.ML, "try-error"))
               stop(mstyle$stop(paste0("Cannot fit ML model", ifelse(verbose, ".", " (set 'verbose=TRUE' to obtain further details)."))))

            ### convergence checks

            if (is.element(optimizer, c("optim","nlminb","dfoptim::hjk","dfoptim::nmk","lbfgsb3c::lbfgsb3c","subplex::subplex","BB::BBoptim","optimx::Rcgmin","optimx::Rvmmin","optimParallel::optimParallel")) && res.ML$convergence != 0)
               stop(mstyle$stop(paste0("Cannot fit ML model. Optimizer (", optimizer, ") did not achieve convergence (convergence = ", res.ML$convergence, ").")))

            if (is.element(optimizer, c("dfoptim::mads")) && res.ML$convergence > optCtrl$tol)
               stop(mstyle$stop(paste0("Cannot fit ML model. Optimizer (", optimizer, ") did not achieve convergence (convergence = ", res.ML$convergence, ").")))

            if (is.element(optimizer, c("minqa::uobyqa","minqa::newuoa","minqa::bobyqa")) && res.ML$ierr != 0)
               stop(mstyle$stop(paste0("Cannot fit ML model. Optimizer (", optimizer, ") did not achieve convergence (ierr = ", res.ML$ierr, ").")))

            if (optimizer=="nloptr::nloptr" && !(res.ML$status >= 1 && res.ML$status <= 4))
               stop(mstyle$stop(paste0("Cannot fit ML model. Optimizer (", optimizer, ") did not achieve convergence (status = ", res.ML$status, ").")))

            if (optimizer=="ucminf::ucminf" && !(res.ML$convergence == 1 || res.ML$convergence == 2))
               stop(mstyle$stop(paste0("Cannot fit ML model. Optimizer (", optimizer, ") did not achieve convergence (convergence = ", res.ML$convergence, ").")))

            if (verbose > 2) {
               cat("\n")
               tmp <- capture.output(print(res.ML))
               .print.output(tmp, mstyle$verbose)
            }

            ### copy estimated values to 'par'

            if (optimizer=="nloptr::nloptr")
               res.ML$par <- res.ML$solution
            if (optimizer=="nlm")
               res.ML$par <- res.ML$estimate

            res.ML$par <- unname(res.ML$par)

            if (verbose > 1)
               message(mstyle$message("Computing the Hessian ..."))

            tau2eff0 <- exp(res.ML$par[p+1]) < con$tau2tol

            if (tau2eff0)
               method <- "T0"

            if (con$hesspack == "numDeriv")
               h.ML <- numDeriv::hessian(.dnchg, x=res.ML$par, method.args=hessianCtrl, ai=ai, bi=bi, ci=ci, di=di, X.fit=X.fit, random=!tau2eff0, verbose=verbose, digits=digits, dnchgcalc=con$dnchgcalc, dnchgprec=con$dnchgprec, intCtrl=intCtrl)
            if (con$hesspack == "pracma")
               h.ML <- pracma::hessian(.dnchg, x0=res.ML$par, ai=ai, bi=bi, ci=ci, di=di, X.fit=X.fit, random=!tau2eff0, verbose=verbose, digits=digits, dnchgcalc=con$dnchgcalc, dnchgprec=con$dnchgprec, intCtrl=intCtrl)
            if (con$hesspack == "calculus")
               h.ML <- calculus::hessian(.dnchg, var=res.ML$par, params=list(ai=ai, bi=bi, ci=ci, di=di, X.fit=X.fit, random=!tau2eff0, verbose=verbose, digits=digits, dnchgcalc=con$dnchgcalc, dnchgprec=con$dnchgprec, intCtrl=intCtrl))
            #return(list(res.ML, h.ML))

            ### log-likelihood

            if (is.element(optimizer, c("optim","dfoptim::hjk","dfoptim::nmk","dfoptim::mads","ucminf::ucminf","lbfgsb3c::lbfgsb3c","subplex::subplex","BB::BBoptim","optimx::Rcgmin","optimx:Rvmmin","optimParallel::optimParallel")))
               ll.ML <- -1 * res.ML$value
            if (is.element(optimizer, c("nlminb","nloptr::nloptr")))
               ll.ML <- -1 * res.ML$objective
            if (is.element(optimizer, c("minqa::uobyqa","minqa::newuoa","minqa::bobyqa")))
               ll.ML <- -1 * res.ML$fval
            if (optimizer == "nlm")
               ll.ML <- -1 * res.ML$minimum

         }

         #return(list(res.FE, res.QE, res.ML, ll.FE=ll.FE, ll.QE=ll.QE, ll.ML=ll.ML))
         #res.FE <- res[[1]]; res.QE <- res[[2]]; res.ML <- res[[3]]

         if (is.element(method, c("FE","EE","CE","T0"))) {

            if (!is.element(optimizer, c("clogit","clogistic"))) {
               beta <- cbind(res.FE$par[seq_len(p)])
               chol.h <- try(chol(h.FE[seq_len(p),seq_len(p)]), silent=!verbose)    # see if Hessian can be inverted with chol()
               if (inherits(chol.h, "try-error") || anyNA(chol.h)) {
                  if (anyNA(chol.h))
                     stop(mstyle$stop(paste0("Cannot invert the Hessian for the ", ifelse(method == "T0", "ML", method), " model.")))
                  warning(mstyle$warning("Choleski factorization of Hessian failed. Trying inversion via QR decomposition."), call.=FALSE)
                  vb <- try(qr.solve(h.FE[seq_len(p),seq_len(p)]), silent=!verbose) # see if Hessian can be inverted with qr.solve()
                  if (inherits(vb, "try-error"))
                     stop(mstyle$stop(paste0("Cannot invert the Hessian for the ", ifelse(method == "T0", "ML", method), " model.")))
               } else {
                  vb <- chol2inv(chol.h)
               }
            }
            if (is.element(optimizer, c("clogit","clogistic"))) {
               beta <- cbind(coef(res.FE)[seq_len(p)])
               vb <- vcov(res.FE)[seq_len(p),seq_len(p),drop=FALSE]
            }
            tau2   <- 0
            sigma2 <- NA_real_
            parms  <- p
            p.eff  <- p
            k.eff  <- k

         }

         if (method == "ML") {

            beta <- cbind(res.ML$par[seq_len(p)])
            chol.h <- try(chol(h.ML), silent=!verbose)      # see if Hessian can be inverted with chol()
            if (inherits(chol.h, "try-error") || anyNA(chol.h)) {
               if (anyNA(chol.h))
                  stop(mstyle$stop("Cannot invert the Hessian for the ML model."))
               warning(mstyle$warning("Choleski factorization of Hessian failed. Trying inversion via QR decomposition."), call.=FALSE)
               vb.f <- try(qr.solve(h.ML), silent=!verbose) # see if Hessian can be inverted with qr.solve()
               if (inherits(vb.f, "try-error"))
                  stop(mstyle$stop("Cannot invert the Hessian for the ML model."))
            } else {
               vb.f <- chol2inv(chol.h)
            }
            vb <- vb.f[seq_len(p),seq_len(p),drop=FALSE]
            if (any(diag(vb) <= 0))
               stop(mstyle$stop("Cannot compute var-cov matrix of the fixed effects."))
            tau2   <- exp(res.ML$par[p+1])
            sigma2 <- NA_real_
            parms  <- p + 1
            p.eff  <- p
            k.eff  <- k
            if (vb.f[p+1,p+1] >= 0) {
               se.tau2 <- sqrt(vb.f[p+1,p+1]) * tau2 # delta rule: vb[p+1,p+1] is the variance of log(tau2), so vb[p+1,p+1] * tau2^2 is the variance of exp(log(tau2))
               crit <- qnorm(level/2, lower.tail=FALSE)
               ci.lb.tau2 <- exp(res.ML$par[p+1] - crit * sqrt(vb.f[p+1,p+1]))
               ci.ub.tau2 <- exp(res.ML$par[p+1] + crit * sqrt(vb.f[p+1,p+1]))
            }

         }

         if (is.element(method, c("ML","T0"))) {

            tmp <- try(rma.uni(measure="PETO", ai=ai, bi=bi, ci=ci, di=di, add=0, mods=X.fit, intercept=FALSE, skipr2=TRUE), silent=TRUE)

            if (!inherits(tmp, "try-error")) {
               gvar1 <- det(vcov(tmp))
               gvar2 <- det(vb)
               ratio <- (gvar1 / gvar2)^(1/(2*m))
               if (!is.na(ratio) && ratio >= 100) {
                  warning(mstyle$warning("Standard errors of fixed effects appear to be unusually small. Treat results with caution."), call.=FALSE)
                  se.warn <- TRUE
               }
               if (!is.na(ratio) && ratio <= 1/100) {
                  warning(mstyle$warning("Standard errors of fixed effects appear to be unusually large. Treat results with caution."), call.=FALSE)
                  se.warn <- TRUE
               }
            }

         }

         if (method == "T0") {

            tau2    <- exp(res.ML$par[p+1])
            parms   <- p + 1
            se.tau2 <- 0
            method  <- "ML"

         }

         #return(list(beta=beta, vb=vb, tau2=tau2, sigma2=sigma2, parms=parms, p.eff=p.eff, k.eff=k.eff, b2.QE=b2.QE, vb2.QE=vb2.QE))

      }

   }

   #########################################################################
   #########################################################################
   #########################################################################

   ### one group outcomes (log odds and log transformed rates)

   if (is.element(measure, c("PLO","PR","PLN","IRLN"))) {

      ### prepare data

      if (is.element(measure, c("PLO","PR","PLN"))) {
         dat.grp <- cbind(xi=xi,mi=mi)
         if (is.null(ddd$family)) {
            if (measure == "PLO")
               dat.fam <- binomial(link=link)
               #dat.fam <- binomial(link="probit")
            if (measure == "PR")
               #dat.fam <- eval(parse(text="binomial(link=\"identity\")"))
               dat.fam <- binomial(link=link)
            if (measure == "PLN")
               dat.fam <- binomial(link=link)
         } else {
            dat.fam <- ddd$family
         }
         dat.off <- NULL
      }

      if (is.element(measure, c("IRLN"))) {
         dat.grp <- xi
         if (is.null(ddd$family)) {
            dat.fam <- poisson(link=link)
         } else {
            dat.fam <- ddd$family
         }
         dat.off <- log(ti)
      }

      study <- factor(seq_len(k)) # study factor
      X.fit <- X

      if (isTRUE(ddd$retdat))
         return(list(dat.grp=dat.grp, X.fit=X.fit, study=study, dat.off = if (!is.null(dat.off)) dat.off else NULL, dat.fam=dat.fam))

      ### fit FE model

      if (verbose)
         message(mstyle$message("Fitting the FE model ..."))

      res.FE <- try(glm(dat.grp ~ -1 + X.fit, offset=dat.off, family=dat.fam, control=glmCtrl), silent=!verbose)

      if (inherits(res.FE, "try-error"))
         stop(mstyle$stop(paste0("Cannot fit FE model", ifelse(verbose, ".", " (set 'verbose=TRUE' to obtain further details)."))))

      ### log-likelihood

      #ll.FE <- with(data.frame(dat.grp), sum(dbinom(xi, xi+mi, predict(res.FE, type="response"), log=TRUE))) # model has a NULL offset
      #ll.FE <- with(data.frame(dat.grp), sum(dpois(xi, predict(res.FE, type="response"), log=TRUE)))         # offset already incorporated into predict()
      ll.FE <- c(logLik(res.FE)) # same as above

      ### fit saturated FE model (= QE model)
      ### notes: 1) suppressWarnings() to suppress warning "glm.fit: fitted probabilities numerically 0 or 1 occurred"

      QEconv <- FALSE
      ll.QE <- NA_real_

      if (!isTRUE(ddd$skiphet)) {

         if (k > 1 && verbose)
            message(mstyle$message("Fitting the saturated model ..."))

         if (k > 1) {
            X.QE <- model.matrix(~ -1 + X.fit + study)
            if (verbose) {
               res.QE <- try(glm(dat.grp ~ -1 + X.QE, offset=dat.off, family=dat.fam, control=glmCtrl), silent=!verbose)
            } else {
               res.QE <- try(suppressWarnings(glm(dat.grp ~ -1 + X.QE, offset=dat.off, family=dat.fam, control=glmCtrl)), silent=!verbose)
            }
         } else {
            res.QE <- res.FE
         }

         if (inherits(res.QE, "try-error")) {

            warning(mstyle$warning(paste0("Cannot fit saturated model", ifelse(verbose, ".", " (set 'verbose=TRUE' to obtain further details)."))), call.=FALSE)

         } else {

            QEconv <- TRUE

            ### log-likelihood

            #ll.QE <- with(data.frame(dat.grp), sum(dbinom(xi, xi+mi, xi/(xi+mi), log=TRUE))) # model has a NULL offset
            #ll.QE <- with(data.frame(dat.grp), sum(dpois(xi, xi, log=TRUE)))                 # offset not relevant for saturated model
            ll.QE <- c(logLik(res.QE)) # same as above

            ### extract coefficients and variance-covariance matrix for Wald-type test for heterogeneity

            #b2.QE <- cbind(na.omit(coef(res.QE)[-seq_len(p)]))                        # coef() still includes aliased coefficients as NAs, so filter them out
            b2.QE  <- cbind(coef(res.QE, complete=FALSE)[-seq_len(p)])                 # aliased coefficients are removed by coef() when complete=FALSE
            vb2.QE <- vcov(res.QE, complete=FALSE)[-seq_len(p),-seq_len(p),drop=FALSE] # aliased coefficients are removed by vcov() when complete=FALSE

         }

      }

      if (method == "ML") {

         ### fit ML model
         ### notes: 1) suppressMessages to suppress the 'one random effect per observation' warning

         if (verbose)
            message(mstyle$message("Fitting the ML model ..."))

         if (package == "lme4") {
            if (verbose) {
               res.ML <- try(lme4::glmer(dat.grp ~ -1 + X.fit + (1 | study), offset=dat.off, family=dat.fam, nAGQ=nAGQ, verbose=verbose, control=do.call(lme4::glmerControl, glmerCtrl)), silent=!verbose)
            } else {
               res.ML <- suppressMessages(try(lme4::glmer(dat.grp ~ -1 + X.fit + (1 | study), offset=dat.off, family=dat.fam, nAGQ=nAGQ, verbose=verbose, control=do.call(lme4::glmerControl, glmerCtrl)), silent=!verbose))
            }
         }

         if (package == "GLMMadaptive") {
            if (is.element(measure, c("PLO","PR","PLN"))) {
               dat.mm <- data.frame(xi=dat.grp[,"xi"], mi=dat.grp[,"mi"], study=study)
               res.ML <- try(GLMMadaptive::mixed_model(cbind(xi,mi) ~ -1 + X.fit, random = ~ 1 | study, data=dat.mm, family=dat.fam, nAGQ=nAGQ, verbose=verbose, control=glmerCtrl), silent=!verbose)
            } else {
               dat.mm <- data.frame(xi=dat.grp, study=study)
               res.ML <- try(GLMMadaptive::mixed_model(xi ~ -1 + X.fit + offset(dat.off), random = ~ 1 | study, data=dat.mm, family=dat.fam, nAGQ=nAGQ, verbose=verbose, control=glmerCtrl), silent=!verbose)
            }
         }

         if (package == "glmmTMB") {
            if (verbose) {
               res.ML <- try(glmmTMB::glmmTMB(dat.grp ~ -1 + X.fit + (1 | study), offset=dat.off, family=dat.fam, verbose=verbose, data=NULL, control=do.call(glmmTMB::glmmTMBControl, glmerCtrl)), silent=!verbose)
            } else {
               res.ML <- suppressMessages(try(glmmTMB::glmmTMB(dat.grp ~ -1 + X.fit + (1 | study), offset=dat.off, family=dat.fam, verbose=verbose, data=NULL, control=do.call(glmmTMB::glmmTMBControl, glmerCtrl)), silent=!verbose))
            }
         }

         if (inherits(res.ML, "try-error"))
            stop(mstyle$stop(paste0("Cannot fit ML model", ifelse(verbose, ".", " (set 'verbose=TRUE' to obtain further details)."))))

         #return(res.ML)

         ### log-likelihood

         #ll.ML <- with(data.frame(dat.grp), sum(dbinom(xi, xi+mi, fitted(res.ML), log=TRUE))) # not correct (since it does not incorporate the random effects; same as ll.FE if tau^2=0)
         #ll.ML <- with(data.frame(dat.grp), sum(dbinom(xi, xi+mi, plogis(qlogis(fitted(res.ML)) + group12*unlist(ranef(res.ML))), log=TRUE))) # not correct (since one really has to integrate; same as ll.FE if tau^2=0)
         #ll.ML <- with(data.frame(dat.grp), sum(dbinom(xi, xi+mi, plogis(predict(res.ML))))) # not correct (since one really has to integrate; same as ll.FE if tau^2=0)
         #ll.ML <- c(logLik(res.ML)) # this is not the same as ll.FE when tau^2 = 0 (not sure why)
         if (package == "lme4") {
            ll.ML <- ll.QE - 1/2 * deviance(res.ML) # this makes ll.ML comparable to ll.FE (same as ll.FE when tau^2=0)
         } else {
            ### FIXME: When using GLMMadaptive, ll is not comparable for FE model when tau^2 = 0
            ll.ML <- c(logLik(res.ML))
         }

      }

      #return(list(res.FE, res.QE, res.ML, ll.FE=ll.FE, ll.QE=ll.QE, ll.ML=ll.ML))
      #res.FE <- res[[1]]; res.QE <- res[[2]]; res.ML <- res[[3]]

      if (is.element(method, c("FE","EE","CE"))) {
         beta   <- cbind(coef(res.FE)[seq_len(p)])
         vb     <- vcov(res.FE)[seq_len(p),seq_len(p),drop=FALSE]
         tau2   <- 0
         sigma2 <- NA_real_
         parms  <- p
         p.eff  <- p
         k.eff  <- k
      }

      if (method == "ML") {

         if (package == "lme4") {
            beta   <- cbind(lme4::fixef(res.ML)[seq_len(p)])
            vb     <- as.matrix(vcov(res.ML))[seq_len(p),seq_len(p),drop=FALSE]
            tau2   <- lme4::VarCorr(res.ML)[[1]][1]
         }

         if (package == "GLMMadaptive") {
            beta   <- cbind(GLMMadaptive::fixef(res.ML)[seq_len(p)])
            vb     <- as.matrix(vcov(res.ML))[seq_len(p),seq_len(p),drop=FALSE]
            tau2   <- res.ML$D[1,1]
         }

         if (package == "glmmTMB") {
            beta   <- cbind(glmmTMB::fixef(res.ML)$cond[seq_len(p)])
            vb     <- as.matrix(vcov(res.ML)$cond)[seq_len(p),seq_len(p),drop=FALSE]
            tau2   <- glmmTMB::VarCorr(res.ML)[[1]][[1]][[1]]
         }

         sigma2 <- NA_real_
         parms  <- p + 1
         p.eff  <- p
         k.eff  <- k

      }

      #return(list(beta=beta, vb=vb, tau2=tau2, sigma2=sigma2, parms=parms, p.eff=p.eff, k.eff=k.eff, b2.QE=b2.QE, vb2.QE=vb2.QE))

   }

   #########################################################################
   #########################################################################
   #########################################################################

   ### heterogeneity tests (Wald-type and likelihood ratio tests of the extra coefficients in the saturated model)

   if (verbose > 1)
      message(mstyle$message("Conducting the heterogeneity tests ..."))

   if (k > 1 && QEconv) {

      ### for OR + CM.EL + NOT clogit/clogistic, QE.Wld is already calculated, so skip this part then

      if (!(measure == "OR" && model == "CM.EL" && !is.element(optimizer, c("clogit","clogistic")))) {

         if (nrow(vb2.QE) > 0) {

            chol.h <- try(chol(vb2.QE), silent=!verbose) # see if Hessian can be inverted with chol()

            if (inherits(chol.h, "try-error") || anyNA(chol.h)) {
               warning(mstyle$warning("Cannot invert the Hessian for the saturated model."), call.=FALSE)
               QE.Wld <- NA_real_
            } else {
               QE.Wld <- try(c(t(b2.QE) %*% chol2inv(chol.h) %*% b2.QE), silent=!verbose)
               if (inherits(QE.Wld, "try-error")) {
                  warning(mstyle$warning("Cannot invert the Hessian for the saturated model."), call.=FALSE)
                  QE.Wld <- NA_real_
               }
            }

         } else {

            QE.Wld <- 0 # if vb2.QE has 0x0 dims, then fitted model is the saturated model and QE.Wld must be 0

         }

      }

      QE.LRT <- -2 * (ll.FE - ll.QE)

      QE.Wld[QE.Wld <= 0] <- 0
      QE.LRT[QE.LRT <= 0] <- 0

      #QE.df <- length(b2.QE) # removed coefficients are not counted if dfs are determined like this
      QE.df <- k-p            # this yields always the same dfs regardless of how many coefficients are removed
      if (QE.df > 0L) {
         QEp.Wld <- pchisq(QE.Wld, df=QE.df, lower.tail=FALSE)
         QEp.LRT <- pchisq(QE.LRT, df=QE.df, lower.tail=FALSE)
      } else {
         QEp.Wld <- 1
         QEp.LRT <- 1
      }

   } else {

      QE.Wld  <- NA_real_
      QE.LRT  <- NA_real_
      QEp.Wld <- NA_real_
      QEp.LRT <- NA_real_
      QE.df   <- NA_integer_

   }

   ### calculation of I^2 and H^2

   wi    <- 1/vi
   W     <- .diag(wi)
   stXWX <- .invcalc(X=X.yi, W=W, k=k.yi)
   P     <- W - W %*% X.yi %*% stXWX %*% crossprod(X.yi,W)
   if (i2def == "1")
      vt <- (k.yi-p) / .tr(P)
   if (i2def == "2")
      vt <- 1/mean(wi) # harmonic mean of vi's (see Takkouche et al., 1999)
   #vt <- (k-1) / (sum(wi) - sum(wi^2)/sum(wi)) # this only applies to the RE model
   I2  <- 100 * tau2 / (vt + tau2)
   H2  <- tau2 / vt + 1

   ### testing of the fixed effects in the model

   if (verbose > 1)
      message(mstyle$message("Conducting the tests of the fixed effects ..."))

   chol.h <- try(chol(vb[btt,btt]), silent=!verbose) # see if Hessian can be inverted with chol()

   if (inherits(chol.h, "try-error") || anyNA(chol.h)) {
      warning(mstyle$warning("Cannot invert the Hessian for the QM-test."), call.=FALSE)
      QM <- NA_real_
   } else {
      QM <- as.vector(t(beta)[btt] %*% chol2inv(chol.h) %*% beta[btt])
   }

   ### scale back beta and vb

   if (!int.only && int.incl && con$scaleX) {
      mX <- rbind(c(intrcpt=1, -1*ifelse(is.d[-1], 0, meanX/sdX)), cbind(0, diag(ifelse(is.d[-1], 1, 1/sdX), nrow=length(is.d)-1, ncol=length(is.d)-1)))
      beta <- mX %*% beta
      vb <- mX %*% vb %*% t(mX)
      X <- Xsave
   }

   ### ddf calculation

   if (test == "t") {
      ddf <- k-p
   } else {
      ddf <- NA_integer_
   }

   ### abbreviate certain coefficient names

   if (isTRUE(ddd$abbrev)) {
      tmp <- colnames(X)
      tmp <- gsub("relevel(factor(", "", tmp, fixed=TRUE)
      tmp <- gsub("\\), ref = \"[[:alnum:]]*\")", "", tmp)
      tmp <- gsub("poly(", "", tmp, fixed=TRUE)
      tmp <- gsub(", degree = [[:digit:]], raw = TRUE)", "^", tmp)
      tmp <- gsub(", degree = [[:digit:]], raw = T)", "^", tmp)
      tmp <- gsub(", degree = [[:digit:]])", "^", tmp)
      tmp <- gsub("rcs\\([[:alnum:]]*, [[:digit:]]\\)", "", tmp)
      tmp <- gsub("factor(", "", tmp, fixed=TRUE)
      tmp <- gsub("I(", "", tmp, fixed=TRUE)
      tmp <- gsub(")", "", tmp, fixed=TRUE)
      colnames(X) <- tmp
   }

   rownames(beta) <- rownames(vb) <- colnames(vb) <- colnames(X.f) <- colnames(X)

   ve <- diag(vb)
   se <- ifelse(ve >= 0, sqrt(ve), NA_real_)
   names(se) <- NULL
   zval <- c(beta/se)

   if (test == "t") {
      QM   <- QM / m
      QMdf <- c(m, k-p)
      QMp  <- if (QMdf[2] > 0) pf(QM, df1=QMdf[1], df2=QMdf[2], lower.tail=FALSE) else NA_real_
      pval <- if (ddf > 0) 2*pt(abs(zval), df=ddf, lower.tail=FALSE) else rep(NA_real_, p)
      crit <- if (ddf > 0) qt(level/2, df=ddf, lower.tail=FALSE) else rep(NA_real_, p)
   } else {
      QMdf <- c(m, NA_integer_)
      QMp  <- pchisq(QM, df=QMdf[1], lower.tail=FALSE)
      pval <- 2*pnorm(abs(zval), lower.tail=FALSE)
      crit <- qnorm(level/2, lower.tail=FALSE)
   }

   ci.lb <- c(beta - crit * se)
   ci.ub <- c(beta + crit * se)

   #return(list(beta=beta, se=se, zval=zval, ci.lb=ci.lb, ci.ub=ci.ub, vb=vb, tau2=tau2, QM=QM, QMp=QMp))

   #########################################################################

   ###### fit statistics

   if (verbose > 1)
      message(mstyle$message("Computing fit statistics and log-likelihood ..."))

   ll.ML     <- ifelse(is.element(method, c("FE","EE","CE")), ll.FE, ll.ML)
   ll.REML   <- NA_real_
   dev.ML    <- -2 * (ll.ML - ll.QE)
   AIC.ML    <- -2 * ll.ML + 2*parms
   BIC.ML    <- -2 * ll.ML +   parms * log(k.eff)
   AICc.ML   <- -2 * ll.ML + 2*parms * max(k.eff, parms+2) / (max(k.eff, parms+2) - parms - 1)
   dev.REML  <- NA_real_
   AIC.REML  <- NA_real_
   BIC.REML  <- NA_real_
   AICc.REML <- NA_real_

   fit.stats <- matrix(c(ll.ML, dev.ML, AIC.ML, BIC.ML, AICc.ML, ll.REML, dev.REML, AIC.REML, BIC.REML, AICc.REML), ncol=2, byrow=FALSE)
   dimnames(fit.stats) <- list(c("ll","dev","AIC","BIC","AICc"), c("ML","REML"))
   fit.stats <- data.frame(fit.stats)

   #########################################################################

   ###### prepare output

   if (verbose > 1)
      message(mstyle$message("Preparing the output ..."))

   weighted  <- TRUE

   if (is.null(ddd$outlist) || ddd$outlist == "nodata") {

      outdat <- list(ai=ai, bi=bi, ci=ci, di=di, x1i=x1i, x2i=x2i, t1i=t1i, t2i=t2i, xi=xi, mi=mi, ti=ti)

      res <- list(b=beta, beta=beta, se=se, zval=zval, pval=pval, ci.lb=ci.lb, ci.ub=ci.ub, vb=vb,
                  tau2=tau2, se.tau2=se.tau2, sigma2=sigma2, rho=rho, ci.lb.tau2=ci.lb.tau2, ci.ub.tau2=ci.ub.tau2,
                  I2=I2, H2=H2, vt=vt,
                  QE.Wld=QE.Wld, QEp.Wld=QEp.Wld, QE.LRT=QE.LRT, QEp.LRT=QEp.LRT, QE.df=QE.df, QM=QM, QMdf=QMdf, QMp=QMp,
                  k=k, k.f=k.f, k.yi=k.yi, k.eff=k.eff, k.all=k.all, p=p, p.eff=p.eff, parms=parms,
                  int.only=int.only, int.incl=int.incl, intercept=intercept,
                  yi=yi, vi=vi, X=X, yi.f=yi.f, vi.f=vi.f, X.f=X.f,
                  chksumyi=digest::digest(as.vector(yi)), chksumvi=digest::digest(as.vector(vi)), chksumX=digest::digest(X),
                  outdat.f=outdat.f, outdat=outdat, ni=ni, ni.f=ni.f,
                  ids=ids, not.na=not.na, subset=subset, not.na.yivi=not.na.yivi, slab=slab, slab.null=slab.null,
                  measure=measure, method=method, model=model, weighted=weighted,
                  test=test, dfs=ddf, ddf=ddf, btt=btt, m=m,
                  digits=digits, level=level, control=control, verbose=verbose,
                  add=add, to=to, drop00=drop00,
                  fit.stats=fit.stats, se.warn=se.warn,
                  formula.yi=NULL, formula.mods=formula.mods, version=packageVersion("metafor"), call=mf)

      if (is.null(ddd$outlist))
         res <- append(res, list(data=data), which(names(res) == "fit.stats"))

   } else {

      if (ddd$outlist == "minimal") {
         res <- list(b=beta, beta=beta, se=se, zval=zval, pval=pval, ci.lb=ci.lb, ci.ub=ci.ub, vb=vb,
                     tau2=tau2, se.tau2=se.tau2, sigma2=sigma2,
                     I2=I2, H2=H2,
                     QE.Wld=QE.Wld, QEp.Wld=QEp.Wld, QE.LRT=QE.LRT, QEp.LRT=QEp.LRT, QE.df=QE.df, QEp=QEp, QM=QM, QMdf=QMdf, QMp=QMp,
                     k=k, k.eff=k.eff, p=p, p.eff=p.eff, parms=parms,
                     int.only=int.only,
                     chksumyi=digest::digest(as.vector(yi)), chksumvi=digest::digest(as.vector(vi)), chksumX=digest::digest(X),
                     measure=measure, method=method, model=model,
                     test=test, dfs=ddf, ddf=ddf, btt=btt, m=m,
                     digits=digits, level=level,
                     fit.stats=fit.stats)
      } else {
         res <- eval(str2lang(paste0("list(", ddd$outlist, ")")))
      }

   }

   if (isTRUE(ddd$retfit)) {
      res$res.FE <- res.FE
      if (!isTRUE(ddd$skiphet))
         res$res.QE <- res.QE
      if (method == "ML")
         res$res.ML <- res.ML
   }

   time.end <- proc.time()
   res$time <- unname(time.end - time.start)[3]

   if (isTRUE(ddd$time))
      .print.time(res$time)

   if (verbose || isTRUE(ddd$time))
      cat("\n")

   class(res) <- c("rma.glmm", "rma")
   return(res)

}
