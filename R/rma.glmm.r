rma.glmm <- function(ai, bi, ci, di, n1i, n2i, x1i, x2i, t1i, t2i, xi, mi, ti, ni, mods,
measure, intercept=TRUE,
data, slab, subset,
add=1/2, to="only0", drop00=TRUE, vtype="LS",
model="UM.FS", method="ML", tdist=FALSE, # weighted=TRUE, ### change tdist to test="t" (or "z")?
level=95, digits=4, btt, nAGQ=7, verbose=FALSE, control) { # tau2,

   #########################################################################

   ###### setup

   ### check argument specifications
   ### (arguments "to" and "vtype" are checked inside escalc function)

   if (missing(measure))
      stop("Need to specify 'measure' argument.")

   if (!is.element(measure, c("OR","IRR","PLO","IRLN")))
      stop("Unknown 'measure' specified.")

   if (!is.element(method, c("FE","ML")))
      stop("Unknown 'method' specified.")

   ### in case user specifies more than one add/to value (as one can do with rma.mh() and rma.peto())
   ### (never apply any kind of continuity correction to the data used in the actual model fitting for models implemented in this function)

   if (length(add) > 1)
      add <- add[1]

   if (length(to) > 1)
      to <- to[1]

   ### model argument only relevant for 2x2 table data (measure="OR") and for 2-group rate data (measure="IRR")
   ### UM.FS/UM.RS = unconditional GLMM with fixed/random study effects (logistic or poisson mixed-effects model with fixed/random intercepts)
   ### CM.EL/CM.AL = conditional GLMM (exact/approximate) (hypergeometric or conditional logistic model)
   ### BV/MV       = bi/multivariate model (logistic or poisson mixed-effects model with unstructured covariance matrix) -- not implemented

   if (!is.element(model, c("UM.FS","UM.RS","CM.EL","CM.AL")))
      stop("Unknown 'model' specified.")

   ### no need for CM.AL for IRR -- use CM.EL

   if (model == "CM.AL" && measure == "IRR")
      model <- "CM.EL"

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop("Unknown 'na.action' specified under options().")

   if (missing(control))
      control <- list()

   if (is.element(measure, c("OR","IRR")) && model == "UM.RS" && method == "ML" && nAGQ > 1) {
      warning("Currently not possible to fit RE/ME model='UM.RS' with nAGQ > 1. nAGQ automatically set to 1.")
      nAGQ <- 1
   }

   knha <- tdist

   #########################################################################

   if (verbose > 1)
      message("Extracting data and computing yi/vi values ...")

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

   ### extract slab, subset, and mods values, possibly from the data frame specified via data (arguments not specified are NULL)

   mf <- match.call()
   mf.slab   <- mf[[match("slab",   names(mf))]]
   mf.subset <- mf[[match("subset", names(mf))]]
   mf.mods   <- mf[[match("mods",   names(mf))]]
   slab   <- eval(mf.slab,   data, enclos=sys.frame(sys.parent()))
   subset <- eval(mf.subset, data, enclos=sys.frame(sys.parent()))
   mods   <- eval(mf.mods,   data, enclos=sys.frame(sys.parent()))

   ai <- bi <- ci <- di <- x1i <- x2i <- t1i <- t2i <- xi <- mi <- ti <- ni <- NA

   ### calculate yi and vi values

   if (is.element(measure, "OR")) {

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

   if (is.element(measure, "IRR")) {

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

   if (is.element(measure, "PLO")) {

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

   if (is.element(measure, "IRLN")) {

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

   yi <- dat$yi         ### one or more yi/vi pairs may be NA/NA (note: yi/vi pairs that are NA/NA may still have 'valid' table data)
   vi <- dat$vi         ### one or more yi/vi pairs may be NA/NA (note: yi/vi pairs that are NA/NA may still have 'valid' table data)
   ni <- attr(yi, "ni") ### unadjusted total sample sizes (ni.u in escalc)

   ### study ids (1:k sequence before subsetting)

   ids <- seq_len(k)

   #########################################################################

   if (verbose > 1)
      message("Creating model matrix ...")

   ### convert mods formula to X matrix and set intercept equal to FALSE

   is.formula <- inherits(mods, "formula")
   if (is.formula) {
      options(na.action = "na.pass")        ### set na.action to na.pass, so that NAs are not filtered out (we'll do that later)
      mods <- model.matrix(mods, data=data) ### extract model matrix
      attr(mods, "assign") <- NULL          ### strip assign attribute (not needed at the moment)
      options(na.action = na.act)           ### set na.action back to na.act
      intercept <- FALSE                    ### set to FALSE since formula now controls whether the intercept is included or not
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
      stop("Number of rows of the model matrix does not match length of the outcome data.")

   ### generate study labels if none are specified

   if (verbose > 1)
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

   ### if a subset of studies is specified (note: tables, yi/vi, and ni are already subsetted above)

   if (!is.null(subset)) {

      if (verbose > 1)
         message("Subsetting ...")

      mods <- mods[subset,,drop=FALSE]
      slab <- slab[subset]
      ids  <- ids[subset]

   }

   ### check if study labels are unique; if not, make them unique

   if (anyDuplicated(slab))
      slab <- make.unique(as.character(slab)) ### make.unique() only works with character vectors

   ### add slab attribute back

   attr(yi, "slab") <- slab

   k <- length(yi) ### number of tables/outcomes after subsetting (can all still include NAs)

   ### if drop00=TRUE, set counts to NA for studies that have no events (or all events) in both arms (corresponding yi/vi will also be NA/NA then)

   if (measure=="OR") {
      if (drop00) {
         id00 <- c(ai == 0L & ci == 0L) | c(bi == 0L & di == 0L)
         id00[is.na(id00)] <- FALSE
         ai[id00] <- NA
         bi[id00] <- NA
         ci[id00] <- NA
         di[id00] <- NA
      }
   }

   if (measure=="IRR") {
      if (drop00) {
         id00 <- c(x1i == 0L & x2i == 0L)
         id00[is.na(id00)] <- FALSE
         x1i[id00] <- NA
         x2i[id00] <- NA
      }
   }

   ### save full data (including potential NAs in table data, yi/vi, and/or mods) (after subsetting)

   ai.f   <- ai
   bi.f   <- bi
   ci.f   <- ci
   di.f   <- di
   x1i.f  <- x1i
   x2i.f  <- x2i
   t1i.f  <- t1i
   t2i.f  <- t2i
   xi.f   <- xi
   mi.f   <- mi
   ti.f   <- ti
   yi.f   <- yi
   vi.f   <- vi
   ni.f   <- ni
   mods.f <- mods

   k.f <- k ### total number of tables/outcomes and rows in the model matrix (including all NAs)

   ### check for NAs in tables (and corresponding mods) and act accordingly

   if (is.element(measure, "OR")) {

      aibicidimods.na <- is.na(ai) | is.na(bi) | is.na(ci) | is.na(di) | if (is.null(mods)) FALSE else apply(is.na(mods), 1, any)

      if (any(aibicidimods.na)) {

         if (verbose > 1)
            message("Handling NAs in table data ...")

         not.na <- !aibicidimods.na

         if (na.act == "na.omit" || na.act == "na.exclude" || na.act == "na.pass") {
            ai   <- ai[not.na]
            bi   <- bi[not.na]
            ci   <- ci[not.na]
            di   <- di[not.na]
            mods <- mods[not.na,,drop=FALSE]
            k    <- length(ai)
            warning("Studies with NAs omitted from model fitting.")
         }

         if (na.act == "na.fail")
            stop("Missing values in studies.")

      } else {
         not.na <- rep(TRUE, k)
      }

   }

   if (is.element(measure, "IRR")) {

      x1ix2it1it2imods.na <- is.na(x1i) | is.na(x2i) | is.na(t1i) | is.na(t2i) | if (is.null(mods)) FALSE else apply(is.na(mods), 1, any)

      if (any(x1ix2it1it2imods.na)) {

         if (verbose > 1)
            message("Handling NAs in table data ...")

         not.na <- !x1ix2it1it2imods.na

         if (na.act == "na.omit" || na.act == "na.exclude" || na.act == "na.pass") {
            x1i  <- x1i[not.na]
            x2i  <- x2i[not.na]
            t1i  <- t1i[not.na]
            t2i  <- t2i[not.na]
            mods <- mods[not.na,,drop=FALSE]
            k    <- length(x1i)
            warning("Studies with NAs omitted from model fitting.")
         }

         if (na.act == "na.fail")
            stop("Missing values in studies.")

      } else {
         not.na <- rep(TRUE, k)
      }

   }

   if (is.element(measure, "PLO")) {

      ximimods.na <- is.na(xi) | is.na(mi) | if (is.null(mods)) FALSE else apply(is.na(mods), 1, any)

      if (any(ximimods.na)) {

         if (verbose > 1)
            message("Handling NAs in table data ...")

         not.na <- !ximimods.na

         if (na.act == "na.omit" || na.act == "na.exclude" || na.act == "na.pass") {
            xi   <- xi[not.na]
            mi   <- mi[not.na]
            mods <- mods[not.na,,drop=FALSE]
            k    <- length(xi)
            warning("Studies with NAs omitted from model fitting.")
         }

         if (na.act == "na.fail")
            stop("Missing values in studies.")

      } else {
         not.na <- rep(TRUE, k)
      }

   }

   if (is.element(measure, "IRLN")) {

      xitimods.na <- is.na(xi) | is.na(ti) | if (is.null(mods)) FALSE else apply(is.na(mods), 1, any)

      if (any(xitimods.na)) {

         if (verbose > 1)
            message("Handling NAs in table data ...")

         not.na <- !xitimods.na

         if (na.act == "na.omit" || na.act == "na.exclude" || na.act == "na.pass") {
            xi   <- xi[not.na]
            ti   <- ti[not.na]
            mods <- mods[not.na,,drop=FALSE]
            k    <- length(xi)
            warning("Studies with NAs omitted from model fitting.")
         }

         if (na.act == "na.fail")
            stop("Missing values in studies.")

      } else {
         not.na <- rep(TRUE, k)
      }

   }

   ### note: k   = number of tables (and corresponding rows of 'mods') after removing NAs
   ###       k.f = total number of tables/outcomes and rows in the model matrix (including all NAs) stored in .f elements

   ### at least one study left?

   if (k < 1)
      stop("Processing terminated since k = 0.")

   ### check for NAs in yi/vi and act accordingly (yi/vi pair can be NA/NA is add=0 is used)
   ### note: if a table was removed because of NAs in mods, must also remove the corresponding yi/vi pair;
   ###       also, must use mods.f here, since NAs in mods were already removed above (and need a separate
   ###       mods.yi element, so that dimensions of the model matrix and vi are guaranteed to match up)

   mods.yi <- mods.f
   yivi.na <- is.na(yi) | is.na(vi) | if (is.null(mods.yi)) FALSE else apply(is.na(mods.yi), 1, any)

   if (any(yivi.na)) {

      if (verbose > 1)
         message("Handling NAs in yi/vi ...")

      not.na.yivi <- !yivi.na

      if (na.act == "na.omit" || na.act == "na.exclude" || na.act == "na.pass") {

         yi <- yi[not.na.yivi]
         ni <- ni[not.na.yivi]
         vi <- vi[not.na.yivi]
         mods.yi <- mods.f[not.na.yivi,,drop=FALSE]
         warning("Some yi/vi values are NA.")

         attr(yi, "measure") <- measure ### add measure attribute back
         attr(yi, "ni")      <- ni      ### add ni attribute back

      }

      if (na.act == "na.fail")
         stop("Missing yi/vi values.")

   } else {
      not.na.yivi <- rep(TRUE, k)
   }

   k.yi <- length(yi) ### number of yi/vi pairs that are not NA

   ### make sure that there is at least one column in X

   if (is.null(mods) && !intercept) {
      warning("Must either include an intercept and/or moderators in model.\n  Coerced intercept into the model.")
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
   ### careful: yi may have become shorter than X due to the omission of NAs, so just use a fake yi vector here

   tmp <- lm(rep(0,k) ~ X - 1)
   coef.na <- is.na(coef(tmp))
   if (any(coef.na)) {
      warning("Redundant predictors dropped from the model.")
      X    <- X[,!coef.na,drop=FALSE]
      X.f  <- X.f[,!coef.na,drop=FALSE]
   }

   ### need to do this separately for X.yi, since model matrix may have fewer rows due to removal of NA/NA pairs for yi/vi

   tmp <- lm(yi ~ X.yi - 1)
   coef.na <- is.na(coef(tmp))
   if (any(coef.na))
      X.yi <- X.yi[,!coef.na,drop=FALSE]

   ### check whether intercept is included and if yes, move it to the first column (NAs already removed, so na.rm=TRUE for any() not necessary)

   is.int <- apply(X, 2, .is.int.func)
   if (any(is.int)) {
      int.incl <- TRUE
      int.indx <- which(is.int, arr.ind=TRUE)
      X        <- cbind(intrcpt=1,   X[,-int.indx, drop=FALSE]) ### note: this removes any duplicate intercepts
      X.f      <- cbind(intrcpt=1, X.f[,-int.indx, drop=FALSE]) ### note: this removes any duplicate intercepts
      if (is.formula)
         intercept <- TRUE ### set intercept appropriately so that the predict() function works
   } else {
      int.incl <- FALSE
   }

   ### need to do this separately for X.yi, since model matrix may have fewer rows due to removal of NA/NA pairs for yi/vi

   is.int <- apply(X.yi, 2, .is.int.func)
   if (any(is.int)) {
      int.indx <- which(is.int, arr.ind=TRUE)
      X.yi     <- cbind(intrcpt=1, X.yi[,-int.indx, drop=FALSE]) ### note: this removes any duplicate intercepts
   }

   p <- NCOL(X) ### number of columns in X (including the intercept if it is included)

   ### note: number of columns in X.yi may be lower than p; but computation of I^2 below is based on p

   ### check whether this is an intercept-only model

   if ((p == 1L) && (all(sapply(X, identical, 1)))) {
      int.only <- TRUE
   } else {
      int.only <- FALSE
   }

   ### check if there are too many parameters for given k

   if (method == "FE" && p > k)
      stop("Number of parameters to be estimated is larger than the number of observations.")
   if (method != "FE" && (p+1) > k)
      stop("Number of parameters to be estimated is larger than the number of observations.")

   ### set/check 'btt' argument

   btt <- .set.btt(btt, p, int.incl)
   m <- length(btt) ### number of betas to test (m = p if all betas are tested)

   #########################################################################

   ### set default control parameters

   con <- list(verbose = FALSE,            # also passed on to glm/glmer/optim/nlminb/minqa (uobyqa/newuoa/bobyqa)
               optimizer = "optim",        # optimizer to use ("optim", "nlminb", "uobyqa", "newuoa", "bobyqa", "clogit", "clogistic")
               optmethod = "BFGS",         # argument 'method' for optim() ("Nelder-Mead" and "BFGS" are sensible options)
               scale = TRUE,               # should non-dummy variables in the X matrix be rescaled before model fitting?
               tol = 1e-07,                # lower bound for eigenvalues to determine if var-cov matrix is positive definite
               dnchgcalc = "dFNCHypergeo", # method for calculating dnchg ("dFNCHypergeo" from BiasedUrn package or "dnoncenhypergeom")
               dnchgprec = 1e-10)          # precision for dFNCHypergeo()

   ### replace defaults with any user-defined values

   con.pos <- pmatch(names(control), names(con))
   con[c(na.omit(con.pos))] <- control[!is.na(con.pos)]

   if (verbose)
      con$verbose <- verbose

   verbose <- con$verbose

   pos.optCtrl <- pmatch(names(control), "optCtrl", nomatch=0)
   if (sum(pos.optCtrl) > 0) {
      optCtrl <- control[[which(pos.optCtrl == 1)]]
   } else {
      optCtrl <- list()
   }

   if (con$optimizer == "optim") {
      con.pos <- pmatch(names(optCtrl), "REPORT", nomatch=0) ### set REPORT to 1 if it is not already set by the user
      if (sum(con.pos) > 0) {
         optCtrl[which(con.pos == 1)] <- 1
         names(optCtrl)[which(con.pos == 1)] <- "REPORT"
      } else {
         optCtrl$REPORT <- 1
      }
      optCtrl$trace <- con$verbose ### trace for optim is a non-negative integer
   }

   if (con$optimizer == "nlminb")
      optCtrl$trace <- ifelse(con$verbose > 0, 1, 0) ### set trace to 1, so information is printed every iteration

   if (is.element(con$optimizer, c("uobyqa", "newuoa", "bobyqa")))
      optCtrl$iprint <- ifelse(con$verbose > 0, 3, 0) ### set iprint to 3 for maximum information

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
   glmCtrl$trace <- ifelse(con$verbose > 0, TRUE, FALSE) ### trace for glmCtrl is logical

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
      intCtrl[which(con.pos == 1)] <- -Inf
      names(intCtrl)[which(con.pos == 1)] <- "lower"
   } else {
      intCtrl$lower <- -Inf
   }
   con.pos <- pmatch(names(intCtrl), "upper", nomatch=0)
   if (sum(con.pos) > 0) {
      intCtrl[which(con.pos == 1)] <- Inf
      names(intCtrl)[which(con.pos == 1)] <- "upper"
   } else {
      intCtrl$upper <- Inf
   }
   con.pos <- pmatch(names(intCtrl), "subdivisions", nomatch=0)
   if (sum(con.pos) > 0) {
      intCtrl[which(con.pos == 1)] <- 100L
      names(intCtrl)[which(con.pos == 1)] <- "subdivisions"
   } else {
      intCtrl$subdivisions <- 100L
   }
   con.pos <- pmatch(names(intCtrl), "rel.tol", nomatch=0)
   if (sum(con.pos) > 0) {
      intCtrl[which(con.pos == 1)] <- .Machine$double.eps^0.25
      names(intCtrl)[which(con.pos == 1)] <- "rel.tol"
   } else {
      intCtrl$rel.tol <- .Machine$double.eps^0.25
   }

   pos.hessianCtrl <- pmatch(names(control), "hessianCtrl", nomatch=0)
   if (sum(pos.hessianCtrl) > 0) {
      hessianCtrl <- control[[which(pos.hessianCtrl == 1)]]
   } else {
      hessianCtrl <- list(r=8)
   }

   #return(list(verbose=verbose, optimizer=con$optimizer, dnchgcalc=con$dnchgcalc, dnchgprec=con$dnchgprec, optCtrl=optCtrl, glmCtrl=glmCtrl, glmerCtrl=glmerCtrl, intCtrl=intCtrl, hessianCtrl=hessianCtrl))

   if (!is.element(con$optimizer, c("optim","nlminb","uobyqa","newuoa","bobyqa","clogit","clogistic")))
      stop("Unknown optimizer specified.")

   if (con$dnchgcalc != "dnoncenhypergeom" && con$dnchgcalc != "dFNCHypergeo")
      stop("Unknown dnchgcalc method specified.")

   if (is.element(con$optimizer, c("clogit", "clogistic")) && method == "ML")
      stop("Cannot use 'clogit' or 'clogistic' with method='ML'.")

   #########################################################################

   ### check that required packages are available

   if (is.element(measure, c("OR","IRR"))) {
      if ((model == "UM.FS" && method == "ML") || (model == "UM.RS") || (model == "CM.AL" && method == "ML") || (model == "CM.EL" && method == "ML")) {
         if (!requireNamespace("lme4", quietly=TRUE))
            stop("Please install the 'lme4' package to fit this model.")
      }
   }

   if (is.element(measure, c("PLO","IRLN")) && method == "ML") {
      if (!requireNamespace("lme4", quietly=TRUE))
         stop("Please install the 'lme4' package to fit this model.")
   }

   if (measure == "OR" && model == "CM.EL") {
      if (is.element(con$optimizer, c("uobyqa", "newuoa", "bobyqa"))) {
         if (!requireNamespace("minqa", quietly=TRUE))
            stop("Please install the 'minqa' package to fit this model.")
         minqa <- get(con$optimizer, envir=loadNamespace("minqa"))
         con$optimizer <- "minqa"
      }
      if (con$optimizer == "optim" || con$optimizer == "nlminb" || con$optimizer == "minqa") {
         if (!requireNamespace("numDeriv", quietly=TRUE))
            stop("Please install the 'numDeriv' package to fit this model.")
         if (con$dnchgcalc == "dFNCHypergeo") {
            if (!requireNamespace("BiasedUrn", quietly=TRUE))
               stop("Please install the 'BiasedUrn' package to fit this model.")
         }
      }
      if (con$optimizer == "clogit") {
         if (!requireNamespace("survival", quietly=TRUE))
            stop("Please install the 'survival' package to fit this model.")
         coxph <- survival::coxph
         Surv  <- survival::Surv
      }
      if (con$optimizer == "clogistic") {
         if (!requireNamespace("Epi", quietly=TRUE))
            stop("Please install the 'Epi' package to fit this model.")
      }
   }

   ### check whether model matrix is of full rank

   if (any(eigen(crossprod(X), symmetric=TRUE, only.values=TRUE)$values <= con$tol))
      stop("Model matrix not of full rank. Cannot fit model.")

   #########################################################################
   #########################################################################
   #########################################################################

   se.tau2 <- I2 <- H2 <- QE <- QEp <- NA

   alpha <- ifelse(level > 1, (100-level)/100, 1-level)

   ###### model fitting, test statistics, and confidence intervals

   ### upgrade warnings to errors (for some testing)
   #o.warn <- getOption("warn")
   #on.exit(options(warn = o.warn))
   #options(warn = 2)

   ### rescale X matrix (only for models with moderators and models including an intercept term)

   if (!int.only && int.incl && con$scale) {
      Xsave <- X
      meanX <- colMeans(X[, 2:p, drop=FALSE])
      sdX   <- apply(X[, 2:p, drop=FALSE], 2, sd) ### consider using colSds() from matrixStats package
      is.d  <- apply(X, 2, function(x) all(sapply(x, identical, 0) | sapply(x, identical, 1))) ### is each column a dummy variable (i.e., only 0s and 1s)?
      X[,!is.d] <- apply(X[, !is.d, drop=FALSE], 2, scale) ### rescale the non-dummy variables
   }

   #########################################################################
   #########################################################################
   #########################################################################

   ### two group outcomes (odds ratios and incidence rate ratios)

   if (is.element(measure, c("OR","IRR"))) {

      ######################################################################

      if (is.element(model, c("UM.FS","UM.RS"))) {

         ### prepare grp-level data for the unconditional models

         if (measure == "OR") {                                      ###                           xi   mi   study   group1  group2  group12  offset  intrcpt  mod1
            dat.grp <- cbind(xi=c(rbind(ai,ci)), mi=c(rbind(bi,di))) ### grp-level outcome data    ai   bi   i       1       0       +1/2     NULL    1        x1i
            dat.fam <- binomial                                      ###                           ci   di   i       0       1       -1/2     NULL    0        0
            dat.off <- NULL
         }

         if (measure == "IRR") {                                     ###                           xi   ti   study   group1  group2  group12  offset  intrcpt  mod1
            dat.grp <- cbind(xi=c(rbind(x1i,x2i)))                   ### grp-level outcome data    x1i  t1i  i       1       0       +1/2     t1i     1        x1i
            dat.fam <- poisson                                       ### log(ti) for offset        x2i  t2i  i       0       1       -1/2     t2i     0        0
            dat.off <- log(c(rbind(t1i,t2i)))                        ###
         }

         group1  <- rep(c(1,0), times=k)                             ### group dummy for 1st group (ai,bi for group 1)
         group2  <- rep(c(0,1), times=k)                             ### group dummy for 2nd group (ci,di for group 2) (not really needed)
         group12 <- rep(c(1/2,-1/2), times=k)                        ### group dummy with +- 1/2 coding
         study   <- factor(rep(seq_len(k), each=2))                  ### study factor
         const   <- cbind(rep(1,2*k))                                ### intercept for random study effects model

         X.fit   <- X[rep(seq(k), each=2),,drop=FALSE]               ### duplicate each row in X (drop=FALSE, so column names are preserved)
         X.fit   <- cbind(group1*X.fit[,,drop=FALSE])                ### then multiply by group1 dummy (intercept, if included, becomes the group1 dummy)

         row.names(X.fit) <- seq_len(2*k)

         #return(data.frame(dat.grp, X.fit, study, dat.off=ifelse(!is.null(dat.off), dat.off, NA), const, group1, group2, group12))

         ###################################################################

         ####################################################
         ### unconditional model with fixed study effects ###
         ####################################################

         if (model == "UM.FS") {

            ### fit FE model

            if (verbose)
               message("Fitting FE model ...")

            if (k > 1) {
               res.FE <- try(glm(dat.grp ~ -1 + X.fit + study, offset=dat.off, family=dat.fam, control=glmCtrl), silent=!verbose)
            } else {
               res.FE <- try(glm(dat.grp ~ -1 + X.fit + const, offset=dat.off, family=dat.fam, control=glmCtrl), silent=!verbose)
            }

            if (inherits(res.FE, "try-error"))
               stop("Cannot fit FE model.")

            ### log-likelihood

            #ll.FE <- with(data.frame(dat.grp), sum(dbinom(xi, xi+mi, predict(res.FE, type="response"), log=TRUE))) ### model has a NULL offset
            #ll.FE <- with(data.frame(dat.grp), sum(dpois(xi, predict(res.FE, type="response"), log=TRUE)))         ### offset already incorporated into predict()
            ll.FE <- c(logLik(res.FE)) ### same as above

            ### fit saturated FE model (= QE model)

            if (verbose)
               message("Fitting saturated model ...")

            if (k > 1) {
               X.QE <- model.matrix(~ -1 + X.fit + study + study:group1)
               res.QE <- try(glm(dat.grp ~ -1 + X.QE, offset=dat.off, family=dat.fam, control=glmCtrl), silent=!verbose)
            } else {
               res.QE <- res.FE
            }

            if (inherits(res.QE, "try-error")) {

               warning("Cannot fit saturated model.")
               QEconv <- FALSE
               ll.QE <- NA

            } else {

               QEconv <- TRUE

               ### log-likelihood

               #ll.QE <- with(data.frame(dat.grp), sum(dbinom(xi, xi+mi, xi/(xi+mi), log=TRUE))) ### model has a NULL offset
               #ll.QE <- with(data.frame(dat.grp), sum(dpois(xi, xi, log=TRUE)))                 ### offset not relevant for saturated model
               ll.QE  <- c(logLik(res.QE)) ### same as above

               ### extract coefficients and variance-covariance matrix for Wald-type test for heterogeneity

               b2.QE  <- cbind(na.omit(coef(res.QE)[-seq_len(k+p)]))          ### coef() still includes aliased coefficients as NAs, so have to filter them out
               vb2.QE <- vcov(res.QE)[-seq_len(k+p),-seq_len(k+p),drop=FALSE] ### aliased coefficients are already removed by vcov()

            }

            if (method == "ML") {

               ### fit ML model
               ### notes: 1) not recommended alternative: using group1 instead of group12 for the random effect (since that forces the variance in group 2 to be lower)

               if (verbose)
                  message("Fitting ML model ...")

               res.ML <- try(lme4::glmer(dat.grp ~ -1 + X.fit + study + (group12 - 1 | study), offset=dat.off, family=dat.fam, nAGQ=nAGQ, verbose=verbose, control=do.call(lme4::glmerControl, glmerCtrl)), silent=!verbose)
               #return(res.ML)

               if (inherits(res.ML, "try-error"))
                  stop("Cannot fit ML model.")

               ### log-likelihood

               #ll.ML <- with(data.frame(dat.grp), sum(dbinom(xi, xi+mi, fitted(res.ML), log=TRUE))) ### not correct (since it does not incorporate the random effects; same as ll.FE if tau^2=0)
               #ll.ML <- with(data.frame(dat.grp), sum(dbinom(xi, xi+mi, plogis(qlogis(fitted(res.ML)) + group12*unlist(ranef(res.ML))), log=TRUE))) ### not correct (since one really has to integrate; same as ll.FE if tau^2=0)
               #ll.ML <- c(logLik(res.ML)) ### this is not the same as ll.FE when tau^2 = 0 (not sure why)
               ll.ML <- ll.QE - 1/2 * deviance(res.ML) ### this makes ll.ML comparable to ll.FE (same as ll.FE when tau^2=0)

            }

            #return(list(res.FE, res.QE, ll.FE=ll.FE, ll.QE=ll.QE))
            #return(list(res.FE, res.QE, res.ML, ll.FE=ll.FE, ll.QE=ll.QE, ll.ML=ll.ML))
            #res.FE <- res[[1]]; res.QE <- res[[2]]; res.ML <- res[[3]]

            if (method == "FE") {
               b      <- cbind(coef(res.FE)[seq_len(p)])
               vb     <- vcov(res.FE)[seq_len(p),seq_len(p),drop=FALSE]
               tau2   <- 0
               sigma2 <- NA
               parms  <- p + k
               p.eff  <- p + k
               k.eff  <- 2*k
            }

            if (method == "ML") {
               b      <- cbind(lme4::fixef(res.ML)[seq_len(p)])
               vb     <- as.matrix(vcov(res.ML))[seq_len(p),seq_len(p),drop=FALSE]
               tau2   <- lme4::VarCorr(res.ML)[[1]][1]
               sigma2 <- NA
               parms  <- p + k + 1
               p.eff  <- p + k
               k.eff  <- 2*k
            }

            #return(list(b=b, vb=vb, tau2=tau2, sigma2=sigma2, parms=parms, p.eff=p.eff, k.eff=k.eff, b2.QE=b2.QE, vb2.QE=vb2.QE))

         }

         ###################################################################

         #####################################################
         ### unconditional model with random study effects ###
         #####################################################

         if (model == "UM.RS") {

            ### fit FE model

            if (verbose)
               message("Fitting FE model ...")

            res.FE <- try(lme4::glmer(dat.grp ~ -1 + X.fit + const + (1 | study), offset=dat.off, family=dat.fam, nAGQ=nAGQ, verbose=verbose, control=do.call(lme4::glmerControl, glmerCtrl)), silent=!verbose)

            if (inherits(res.FE, "try-error"))
               stop("Cannot fit FE model.")

            ### log-likelihood

            ll.FE <- c(logLik(res.FE))

            ### fit saturated FE model (= QE model)
            ### notes: 1) must figure out which terms are aliased in saturated model and remove those terms before fitting
            ###        2) sigma^2 for the study random effect must be close to the one from the FE model - so set start value to sigma from that model

            if (verbose)
               message("Fitting saturated model ...")

            if (k > 1) {
               X.QE <- model.matrix(~ -1 + X.fit + const + study:group1)
               res.QE <- try(glm(dat.grp ~ -1 + X.QE, offset=dat.off, family=dat.fam, control=glmCtrl), silent=!verbose)
               X.QE <- X.QE[,!is.na(coef(res.QE)),drop=FALSE]
               res.QE <- try(lme4::glmer(dat.grp ~ -1 + X.QE + (1 | study), offset=dat.off, family=dat.fam, start=c(sqrt(lme4::VarCorr(res.FE)[[1]][1])), nAGQ=nAGQ, verbose=verbose, control=do.call(lme4::glmerControl, glmerCtrl)), silent=!verbose)
            } else {
               res.QE <- res.FE
            }

            if (inherits(res.QE, "try-error")) {

               warning("Cannot fit saturated model.")
               QEconv <- FALSE
               ll.QE <- NA

            } else {

               QEconv <- TRUE

               ### log-likelihood

               ll.QE <- c(logLik(res.QE))

               ### extract coefficients and variance-covariance matrix for Wald-type test for heterogeneity

               b2.QE  <- cbind(lme4::fixef(res.QE)[-seq_len(p+1)])                       ### aliased coefficients are already removed
               vb2.QE <- as.matrix(vcov(res.QE))[-seq_len(p+1),-seq_len(p+1),drop=FALSE] ### aliased coefficients are already removed

            }

            if (method == "ML") {

               ### fit ML model
               ### notes: 1) not recommended alternative: using group1 instead of group12 for the random effect (since that forces the variance in group 2 to be lower)
               ###        2) this approach is okay if we also allow group1 random effect and intercepts to correlate (in fact, this is identical to the bivariate model)
               ###        3) start=c(sqrt(lme4::VarCorr(res.FE)[[1]][1])) has no effect, since the start value for tau^2 is not specified (and using 0 is probably not ideal for that)

               if (verbose)
                  message("Fitting ML model ...")

               res.ML <- try(lme4::glmer(dat.grp ~ -1 + X.fit + const + (1 | study) + (group12 - 1 | study), offset=dat.off, family=dat.fam, nAGQ=nAGQ, verbose=verbose, control=do.call(lme4::glmerControl, glmerCtrl)), silent=!verbose)
               #res.ML <- try(lme4::glmer(dat.grp ~ -1 + X.fit + const + (group1 | study),                   offset=dat.off, family=dat.fam, nAGQ=nAGQ, verbose=verbose, control=do.call(lme4::glmerControl, glmerCtrl)), silent=!verbose)
               #return(res.ML)

               if (inherits(res.ML, "try-error"))
                  stop("Cannot fit ML model.")

               ### log-likelihood

               ll.ML <- c(logLik(res.ML))

            }

            #return(list(res.FE, res.QE, res.ML, ll.FE=ll.FE, ll.QE=ll.QE, ll.ML=ll.ML))
            #res.FE <- res[[1]]; res.QE <- res[[2]]; res.ML <- res[[3]]

            if (method == "FE") {
               b      <- cbind(lme4::fixef(res.FE)[seq_len(p)])
               vb     <- as.matrix(vcov(res.FE))[seq_len(p),seq_len(p),drop=FALSE]
               tau2   <- 0
               sigma2 <- lme4::VarCorr(res.FE)[[1]][1]
               parms  <- p + 1 + 1
               p.eff  <- p + 1
               k.eff  <- 2*k
            }

            if (method == "ML") {
               b      <- cbind(lme4::fixef(res.ML)[seq_len(p)])
               vb     <- as.matrix(vcov(res.ML))[seq_len(p),seq_len(p),drop=FALSE]
               tau2   <- lme4::VarCorr(res.ML)[[2]][1]
               sigma2 <- lme4::VarCorr(res.ML)[[1]][1]
               parms  <- p + 1 + 2
               p.eff  <- p + 1
               k.eff  <- 2*k
            }

            #return(list(b=b, vb=vb, tau2=tau2, sigma2=sigma2, parms=parms, p.eff=p.eff, k.eff=k.eff, b2.QE=b2.QE, vb2.QE=vb2.QE))

         }

         ###################################################################

      }

      ######################################################################

      if ((measure=="IRR" && model == "CM.EL") || (measure=="OR" && model=="CM.AL") || (measure=="OR" && model=="CM.EL")) {

         ### prepare data for the conditional models

         if (measure == "OR") {
            dat.grp <- cbind(xi=ai, mi=ci)   ### conditional outcome data (number of cases in group 1 conditional on total number of cases)
            dat.off <- log((ai+bi)/(ci+di))  ### log(n1i/n2i) for offset
         }

         if (measure == "IRR") {
            dat.grp <- cbind(xi=x1i, mi=x2i) ### conditional outcome data (number of events in group 1 conditional on total number of events)
            dat.off <- log(t1i/t2i)          ### log(t1i/t1i) for offset
         }

         study <- factor(seq_len(k))         ### study factor
         X.fit <- X

         #return(data.frame(dat.grp, X.fit, study, dat.off=ifelse(!is.null(dat.off), dat.off, NA)))

         ###################################################################

         ###############################################################
         ### conditional model (approx. ll for ORs / exact for IRRs) ###
         ###############################################################

         ### fit FE model

         if (verbose)
            message("Fitting FE model ...")

         res.FE <- try(glm(dat.grp ~ -1 + X.fit, offset=dat.off, family=binomial, control=glmCtrl), silent=!verbose)

         if (inherits(res.FE, "try-error"))
            stop("Cannot fit FE model.")

         ### log-likelihood

         #ll.FE <- with(data.frame(dat.grp), sum(dbinom(xi, xi+mi, predict(res.FE, type="response"), log=TRUE))) ### offset already incorporated into predict()
         #ll.FE <- with(data.frame(dat.grp), sum(dpois(xi, predict(res.FE, type="response"), log=TRUE)))         ### offset already incorporated into predict()
         ll.FE <- c(logLik(res.FE)) ### same as above

         ### fit saturated FE model (= QE model)

         if (verbose)
            message("Fitting saturated model ...")

         if (k > 1) {
            X.QE <- model.matrix(~ -1 + X.fit + study)
            res.QE <- try(glm(dat.grp ~ -1 + X.QE, offset=dat.off, family=binomial, control=glmCtrl), silent=!verbose)
         } else {
            res.QE <- res.FE
         }

         if (inherits(res.QE, "try-error")) {

            warning("Cannot fit saturated model.")
            QEconv <- FALSE
            ll.QE <- NA

         } else {

            QEconv <- TRUE

            ### log-likelihood

            #ll.QE  <- with(data.frame(dat.grp), sum(dbinom(xi, xi+mi, xi/(xi+mi), log=TRUE))) ### offset not relevant for saturated model
            #ll.QE  <- with(data.frame(dat.grp), sum(dpois(xi, xi, log=TRUE)))                 ### offset not relevant for saturated model
            ll.QE  <- c(logLik(res.QE)) ### same as above

            ### extract coefficients and variance-covariance matrix for Wald-type test for heterogeneity

            b2.QE  <- cbind(na.omit(coef(res.QE)[-seq_len(p)]))        ### coef() still includes aliased coefficients as NAs, so have to filter them out
            vb2.QE <- vcov(res.QE)[-seq_len(p),-seq_len(p),drop=FALSE] ### aliased coefficients are already removed by vcov()

         }

         #return(list(res.FE, res.QE, ll.FE, ll.QE))
         #res.FE <- res[[1]]; res.QE <- res[[2]]

         if (method == "ML") {

            ### fit ML model
            ### notes: 1) suppressMessages to suppress the 'one random effect per observation' warning

            if (verbose)
               message("Fitting ML model ...")

            if (verbose) {
               res.ML <- try(lme4::glmer(dat.grp ~ -1 + X.fit + (1 | study), offset=dat.off, family=binomial, nAGQ=nAGQ, verbose=verbose, control=do.call(lme4::glmerControl, glmerCtrl)), silent=!verbose)
            } else {
               res.ML <- suppressMessages(try(lme4::glmer(dat.grp ~ -1 + X.fit + (1 | study), offset=dat.off, family=binomial, nAGQ=nAGQ, verbose=verbose, control=do.call(lme4::glmerControl, glmerCtrl)), silent=!verbose))
            }

            if (inherits(res.ML, "try-error"))
               stop("Cannot fit ML model.")

            ### log-likelihood

            ll.ML  <- ll.QE - 1/2 * deviance(res.ML) ### makes ll.ML comparable to ll.FE (same as ll.FE when tau^2=0)

         }

         #return(list(res.FE, res.QE, res.ML, ll.FE=ll.FE, ll.QE=ll.QE, ll.ML=ll.ML))
         #return(list(res.FE, res.QE, res.ML, ll.FE=ll.FE, ll.QE=ll.QE, ll.ML=ll.ML))
         #res.FE <- res[[1]]; res.QE <- res[[2]]; res.ML <- res[[3]]

         if (method == "FE") {
            b      <- cbind(coef(res.FE)[seq_len(p)])
            vb     <- vcov(res.FE)[seq_len(p),seq_len(p),drop=FALSE]
            tau2   <- 0
            sigma2 <- NA
            parms  <- p
            p.eff  <- p
            k.eff  <- k
         }

         if (method == "ML") {
            b      <- cbind(lme4::fixef(res.ML)[seq_len(p)])
            vb     <- as.matrix(vcov(res.ML))[seq_len(p),seq_len(p),drop=FALSE]
            tau2   <- lme4::VarCorr(res.ML)[[1]][1]
            sigma2 <- NA
            parms  <- p + 1
            p.eff  <- p
            k.eff  <- k
         }

         #return(list(b=b, vb=vb, tau2=tau2, sigma2=sigma2, parms=parms, p.eff=p.eff, k.eff=k.eff, b2.QE=b2.QE, vb2.QE=vb2.QE))

         ###################################################################

      }

      if (measure=="OR" && model=="CM.EL") {

         ####################################################
         ### conditional model (exact likelihood for ORs) ###
         ####################################################

         if (verbose)
            message("Fitting FE model ...")

         if (con$optimizer == "optim" || con$optimizer == "nlminb" || con$optimizer == "minqa") {

            ### fit FE model
            ### notes: 1) this routine uses direct optimization over the non-central hypergeometric distribution
            ###        2) start values from CM.AL model and 0 for tau^2 (held at 0 during the optimization since random=FALSE)
            ###        3) no integration for FE model, so intCtrl is not relevant
            ###        4) results can be sensitive to the scaling of moderators

            if (con$optimizer == "optim") {
               res.FE <- try(optim(par=c(coef(res.FE)[seq_len(p)], 0), .dnchg, method=con$optmethod, ai=ai, bi=bi, ci=ci, di=di, X.fit=X.fit, random=FALSE,
                                   verbose=verbose, digits=digits, dnchgcalc=con$dnchgcalc, dnchgprec=con$dnchgprec, control=optCtrl), silent=!verbose)
            }
            if (con$optimizer == "nlminb") {
               res.FE <- try(nlminb(start=c(coef(res.FE)[seq_len(p)], 0), .dnchg, ai=ai, bi=bi, ci=ci, di=di, X.fit=X.fit, random=FALSE,
                                    verbose=verbose, digits=digits, dnchgcalc=con$dnchgcalc, dnchgprec=con$dnchgprec, control=optCtrl), silent=!verbose)
            }
            if (con$optimizer == "minqa") {
               res.FE <- try(minqa(par=c(coef(res.FE)[seq_len(p)], 0), .dnchg, ai=ai, bi=bi, ci=ci, di=di, X.fit=X.fit, random=FALSE,
                                   verbose=verbose, digits=digits, dnchgcalc=con$dnchgcalc, dnchgprec=con$dnchgprec, control=optCtrl), silent=!verbose)
            }

            if (con$optimizer == "optim" || con$optimizer == "nlminb") {
               if (inherits(res.FE, "try-error") || res.FE$convergence != 0)
                  stop("Cannot fit FE model.")
            }
            if (con$optimizer == "minqa") {
               if (inherits(res.FE, "try-error") || res.FE$ierr != 0)
                  stop("Cannot fit FE model.")
            }

            if (verbose > 1)
               message("Computing Hessian ...")

            h.FE <- numDeriv::hessian(.dnchg, x=res.FE$par, method.args=hessianCtrl, ai=ai, bi=bi, ci=ci, di=di, X.fit=X.fit, random=FALSE, verbose=verbose, digits=digits, dnchgcalc=con$dnchgcalc, dnchgprec=con$dnchgprec)
            #return(list(res.FE, h.FE))

            ### log-likelihood

            if (con$optimizer == "optim")
               ll.FE <- -1 * res.FE$value
            if (con$optimizer == "nlminb")
               ll.FE <- -1 * res.FE$objective
            if (con$optimizer == "minqa")
               ll.FE <- -1 * res.FE$fval

            ### fit saturated FE model (= QE model)
            ### notes: 1) must figure out which terms are aliased in saturated model and remove those terms before fitting
            ###        2) start values from CM.AL model and 0 for tau^2 (held at 0 during the optimization since random=FALSE)
            ###        3) therefore only try to fit saturated model if this was possible with CM.AL
            ###        4) no integration for FE model, so intCtrl is not relevant

            if (QEconv) {

               if (verbose)
                  message("Fitting saturated model ...")

               if (k > 1) {

                  is.aliased <- is.na(coef(res.QE))
                  X.QE <- X.QE[,!is.aliased,drop=FALSE] ### res.QE is from CM.AL model

                  if (con$optimizer == "optim") {
                     res.QE <- try(optim(par=c(coef(res.QE)[!is.aliased], 0), .dnchg, method=con$optmethod, ai=ai, bi=bi, ci=ci, di=di, X.fit=X.QE, random=FALSE,
                                         verbose=verbose, digits=digits, dnchgcalc=con$dnchgcalc, dnchgprec=con$dnchgprec, control=optCtrl), silent=!verbose)
                  }
                  if (con$optimizer == "nlminb") {
                     res.QE <- try(nlminb(start=c(coef(res.QE)[!is.aliased], 0), .dnchg, ai=ai, bi=bi, ci=ci, di=di, X.fit=X.QE, random=FALSE,
                                          verbose=verbose, digits=digits, dnchgcalc=con$dnchgcalc, dnchgprec=con$dnchgprec, control=optCtrl), silent=!verbose)
                  }
                  if (con$optimizer == "minqa") {
                     res.QE <- try(minqa(par=c(coef(res.QE)[!is.aliased], 0), .dnchg, ai=ai, bi=bi, ci=ci, di=di, X.fit=X.QE, random=FALSE,
                                         verbose=verbose, digits=digits, dnchgcalc=con$dnchgcalc, dnchgprec=con$dnchgprec, control=optCtrl), silent=!verbose)
                  }

                  if (con$optimizer == "optim" || con$optimizer == "nlminb") {
                     if (inherits(res.QE, "try-error") || res.QE$convergence != 0) {
                        warning("Cannot fit saturated model.")
                        QEconv <- FALSE
                        ll.QE <- NA
                     }
                  }
                  if (con$optimizer == "minqa") {
                     if (inherits(res.QE, "try-error") || res.QE$ierr != 0) {
                        warning("Cannot fit saturated model.")
                        QEconv <- FALSE
                        ll.QE <- NA
                     }
                  }

                  if (QEconv) {
                     if (verbose > 1)
                        message("Computing Hessian ...")
                     h.QE <- numDeriv::hessian(.dnchg, x=res.QE$par, method.args=hessianCtrl, ai=ai, bi=bi, ci=ci, di=di, X.fit=X.QE, random=FALSE, verbose=verbose, digits=digits, dnchgcalc=con$dnchgcalc, dnchgprec=con$dnchgprec)
                  }

               } else {
                  res.QE <- res.FE
                  h.QE   <- h.FE
               }

               #return(list(res.QE, h.QE))

            }

            if (QEconv) {

               ### log-likelihood

               if (con$optimizer == "optim")
                  ll.QE <- -1 * res.QE$value
               if (con$optimizer == "nlminb")
                  ll.QE <- -1 * res.QE$objective
               if (con$optimizer == "minqa")
                  ll.QE <- -1 * res.QE$fval

               ### extract coefficients and variance-covariance matrix for Wald-type test for heterogeneity

               #return(res.QE)
               b2.QE    <- res.QE$par                                  ### recall: aliased coefficients are already removed
               hessian  <- h.QE                                        ### take hessian from hessian() (again, aliased coefs are already removed)
               #hessian <- res.QE$hessian                              ### take hessian from optim() (again, aliased coefs are already removed)
               p.QE     <- length(b2.QE)                               ### how many parameters are left in saturated model?
               b2.QE    <- b2.QE[-p.QE]                                ### remove last element (for tau^2, constrained to 0)
               hessian  <- hessian[-p.QE,-p.QE,drop=FALSE]             ### remove last row/column (for tau^2, constrained to 0)
               p.QE     <- length(b2.QE)                               ### how many parameters are now left?
               is.0     <- colSums(hessian == 0L) == p.QE              ### any columns in hessian entirely composed of 0s?
               b2.QE    <- b2.QE[!is.0]                                ### keep coefficients where this is not the case
               hessian  <- hessian[!is.0,!is.0,drop=FALSE]             ### keep parts of hessian where this is not the case
               b2.QE    <- cbind(b2.QE[-seq_len(p)])                   ### remove first p coefficients
               h.A      <- hessian[seq_len(p),seq_len(p),drop=FALSE]   ### upper left part of hessian
               h.B      <- hessian[seq_len(p),-seq_len(p),drop=FALSE]  ### upper right part of hessian
               h.C      <- hessian[-seq_len(p),seq_len(p),drop=FALSE]  ### lower left part of hessian
               h.D      <- hessian[-seq_len(p),-seq_len(p),drop=FALSE] ### lower right part of hessian (of which we need the inverse)
               chol.h.A <- try(chol(h.A), silent=!verbose)             ### see if h.A can be inverted with chol()
               if (inherits(chol.h.A, "try-error")) {
                  warning("Cannot invert Hessian for saturated model.")
               } else {
                  Ivb2.QE  <- h.D-h.C%*%chol2inv(chol.h.A)%*%h.B       ### inverse of the inverse of the lower right part
                  QE.Wld   <- c(t(b2.QE) %*% Ivb2.QE %*% b2.QE)        ### Wald statistic (note: this approach only requires taking the inverse of h.A)
               }                                                       ### see: http://en.wikipedia.org/wiki/Invertible_matrix#Blockwise_inversion

               #vb2.QE <- chol2inv(chol(hessian))[-seq_len(p),-seq_len(p),drop=FALSE] ### take inverse, then take part relevant for QE test
               #QE.Wld <- c(t(b2.QE) %*% chol2inv(chol(vb2.QE)) %*% b2.QE)

            }

         }

         if (con$optimizer == "clogit" || con$optimizer == "clogistic") {

            ### fit FE model
            ### notes: 1) this routine uses either clogit() from the survival package or clogistic() from the Epi package
            ###        2) the dataset must be in group-level and IPD format (i.e., not in the conditional format)
            ###        3) if the studies are large, the IPD dataset may also be very large, and R may run out of memory
            ###        4) for larger datasets, run time is often excessive (and may essentially freeze R)
            ###        5) suppressMessages for clogit() to suppress the 'beta may be infinite' warning

            ### prepare IPD dataset                                                                                                      ###                   study  event  group1  intrcpt  moderator
                                                                                                                                         ###                   i      1      1       1        x1i       (repeated ai times)
            event   <- unlist(lapply(seq_len(k), function(i) c(rep.int(1,ai[i]), rep.int(0,bi[i]), rep.int(1,ci[i]), rep.int(0,di[i])))) ### event dummy       i      0      1       1        x1i       (repeated bi times)
            group1  <- unlist(lapply(seq_len(k), function(i) c(rep.int(1,ai[i]), rep.int(1,bi[i]), rep.int(0,ci[i]), rep.int(0,di[i])))) ### group1 dummy      i      1      0       0        0         (repeated ci times)
            study.l <- factor(rep(seq_len(k), times=ni))        ### study factor                                                                               i      0      0       0        0         (repeated di times)
            X.fit.l <- X[rep(seq_len(k), times=ni),,drop=FALSE] ### repeat each row in X ni times each
            X.fit.l <- cbind(group1*X.fit.l)                    ### multiply by group1 dummy (including intercept, which becomes the group1 dummy)
            const   <- rep(1,length(event))

            #return(data.frame(event, group1, study.l, X.fit.l, const))

            ### fit FE model

            if (k > 1) {

               if (con$optimizer == "clogit") {
                  args.clogit <- clogitCtrl
                  args.clogit$formula <- event ~ X.fit.l + strata(study.l)
                  res.FE <- try(do.call(survival::clogit, args.clogit), silent=!verbose)
               }
               if (con$optimizer == "clogistic") {
                  args.clogistic <- clogisticCtrl
                  args.clogistic$formula <- event ~ X.fit.l
                  args.clogistic$strata <- study.l
                  res.FE <- try(do.call(Epi::clogistic, args.clogistic), silent=!verbose)
               }

            } else {
               stop(paste("Cannot use '", con$optimizer, "' optimizer when k=1.", sep=""))
            }

            if (inherits(res.FE, "try-error"))
               stop("Cannot fit FE model.")

            #return(res.FE)

            ### fit saturated FE model (= QE model)
            ### notes: 1) must figure out which terms are aliased in saturated model and remove those terms before fitting
            ###        2) fixed effects part does not include 'study' factor, since this is incorporated into the strata
            ###        3) however, for calculating the log likelihood, we need to go back to the conditional data, so we need to reconstruct X.QE (the study.l:group1 coefficients are the study coefficients)

            if (verbose)
               message("Fitting saturated model ...")

            X.QE.l <- model.matrix(~ -1 + X.fit.l + study.l:group1)
            X.QE.l <- X.QE.l[,!is.na(coef(res.QE)),drop=FALSE]
            X.QE   <- X.QE[,!is.na(coef(res.QE)),drop=FALSE]

            if (con$optimizer == "clogit") {
               args.clogit <- clogitCtrl
               args.clogit$formula <- event ~ X.QE.l + strata(study.l)
               if (verbose) {
                  res.QE <- try(do.call(survival::clogit, args.clogit), silent=!verbose)
               } else {
                  res.QE <- try(suppressWarnings(do.call(survival::clogit, args.clogit)), silent=!verbose)
               }
            }

            if (con$optimizer == "clogistic") {
               args.clogistic <- clogisticCtrl
               args.clogistic$formula <- event ~ X.QE.l
               args.clogistic$strata <- study.l
               res.QE <- try(do.call(Epi::clogistic, args.clogistic), silent=!verbose)
            }

            if (inherits(res.QE, "try-error"))
               stop("Cannot fit saturated model.")

            #return(res.QE)

            ### log-likelihood

            ll.FE <- -1 * .dnchg(c(cbind(coef(res.FE)),0), ai=ai, bi=bi, ci=ci, di=di, X.fit=X.fit, random=FALSE, verbose=verbose, digits=digits, dnchgcalc=con$dnchgcalc, dnchgprec=con$dnchgprec)
            ll.QE <- -1 * .dnchg(c(cbind(coef(res.QE)),0), ai=ai, bi=bi, ci=ci, di=di, X.fit=X.QE,  random=FALSE, verbose=verbose, digits=digits, dnchgcalc=con$dnchgcalc, dnchgprec=con$dnchgprec)

            ### extract coefficients and variance-covariance matrix for Wald-type test for heterogeneity

            b2.QE  <- cbind(coef(res.QE)[-seq_len(p)])                 ### aliased coefficients are already removed
            vb2.QE <- vcov(res.QE)[-seq_len(p),-seq_len(p),drop=FALSE] ### aliased coefficients are already removed

         }

         #return(list(res.FE, res.QE, ll.FE=ll.FE, ll.QE=ll.QE))
         #res.FE <- res[[1]]; res.QE <- res[[2]]

         if (method == "ML") {

            ### fit ML model
            ### notes: 1) cannot use clogit() or clogistic() for this (do not allow for the addition of random effects)
            ###        2) mclogit() from mclogit package may be an alternative (but it only provides PQL method)
            ###        3) start values from CM.AL model (add .001 to tau^2 estimate, in case estimate of tau^2 is 0)
            ###        4) optimization involves integration, so intCtrl is relevant
            ###        5) results can be sensitive to the scaling of moderators

            if (verbose)
               message("Fitting ML model ...")

            if (con$optimizer == "optim") {
               res.ML <- try(optim(par=c(b, log(tau2+.001)), .dnchg, method=con$optmethod, ai=ai, bi=bi, ci=ci, di=di, X.fit=X.fit, random=TRUE,
                                   verbose=verbose, digits=digits, dnchgcalc=con$dnchgcalc, dnchgprec=con$dnchgprec, intCtrl=intCtrl, control=optCtrl), silent=!verbose)
            }
            if (con$optimizer == "nlminb") {
               res.ML <- try(nlminb(start=c(b, log(tau2+.001)), .dnchg, ai=ai, bi=bi, ci=ci, di=di, X.fit=X.fit, random=TRUE,
                                    verbose=verbose, digits=digits, dnchgcalc=con$dnchgcalc, dnchgprec=con$dnchgprec, intCtrl=intCtrl, control=optCtrl), silent=!verbose)
            }
            if (con$optimizer == "minqa") {
               res.ML <- try(minqa(par=c(b, log(tau2+.001)), .dnchg, ai=ai, bi=bi, ci=ci, di=di, X.fit=X.fit, random=TRUE,
                                   verbose=verbose, digits=digits, dnchgcalc=con$dnchgcalc, dnchgprec=con$dnchgprec, intCtrl=intCtrl, control=optCtrl), silent=!verbose)
            }

            if (con$optimizer == "optim" || con$optimizer == "nlminb") {
               if (inherits(res.ML, "try-error") || res.ML$convergence != 0)
                  stop("Cannot fit ML model.")
            }
            if (con$optimizer == "minqa") {
               if (inherits(res.ML, "try-error") || res.ML$ierr != 0)
                  stop("Cannot fit ML model.")
            }

            if (verbose > 1)
               message("Computing Hessian ...")

            ### TODO: r=8 seems to give more accurate results, but this needs a whole lot more testing
            h.ML <- numDeriv::hessian(.dnchg, x=res.ML$par, method.args=hessianCtrl, ai=ai, bi=bi, ci=ci, di=di, X.fit=X.fit, random=TRUE, verbose=verbose, digits=digits, dnchgcalc=con$dnchgcalc, dnchgprec=con$dnchgprec, intCtrl=intCtrl)
            #return(list(res.ML, h.ML))

            ### log-likelihood

            if (con$optimizer == "optim")
               ll.ML <- -1 * res.ML$value
            if (con$optimizer == "nlminb")
               ll.ML <- -1 * res.ML$objective
            if (con$optimizer == "minqa")
               ll.ML <- -1 * res.ML$fval

         }

         #return(list(res.FE, res.QE, res.ML, ll.FE=ll.FE, ll.QE=ll.QE, ll.ML=ll.ML))
         #res.FE <- res[[1]]; res.QE <- res[[2]]; res.ML <- res[[3]]

         if (method == "FE") {
            if (con$optimizer == "optim" || con$optimizer == "nlminb" || con$optimizer == "minqa") {
               b <- cbind(res.FE$par[seq_len(p)])
               #chol.h <- try(chol(res.FE$hessian[seq_len(p),seq_len(p)]), silent=!verbose) ### see if Hessian can be inverted with chol()
               chol.h <- try(chol(h.FE[seq_len(p),seq_len(p)]), silent=!verbose) ### see if Hessian can be inverted with chol()
               if (inherits(chol.h, "try-error")) {
                  warning("Choleski factorization of Hessian failed. Trying inversion via QR decomposition.")
                  vb <- try(qr.solve(h.FE[seq_len(p),seq_len(p)]), silent=!verbose) ### see if Hessian can be inverted with qr.solve()
                  if (inherits(vb, "try-error"))
                     stop("Cannot invert Hessian for ML model.")
               } else {
                  vb <- chol2inv(chol.h)
               }
            }
            if (con$optimizer == "clogit" || con$optimizer == "clogistic") {
               b  <- cbind(coef(res.FE)[seq_len(p)])
               vb <- vcov(res.FE)[seq_len(p),seq_len(p),drop=FALSE]
            }
            tau2   <- 0
            sigma2 <- NA
            parms  <- p
            p.eff  <- p
            k.eff  <- k
         }

         if (method == "ML") {
            b <- cbind(res.ML$par[seq_len(p)])
            chol.h <- try(chol(h.ML), silent=!verbose) ### see if Hessian can be inverted with chol()
            if (inherits(chol.h, "try-error")) {
               warning("Choleski factorization of Hessian failed. Trying inversion via QR decomposition.")
               vb.f <- try(qr.solve(h.ML), silent=!verbose) ### see if Hessian can be inverted with qr.solve()
               if (inherits(vb.f, "try-error"))
                  stop("Cannot invert Hessian for ML model.")
            } else {
               vb.f <- chol2inv(chol.h)
            }
            vb     <- vb.f[seq_len(p),seq_len(p),drop=FALSE]
            tau2   <- exp(res.ML$par[p+1])
            sigma2 <- NA
            parms  <- p + 1
            p.eff  <- p
            k.eff  <- k
            if (vb.f[p+1,p+1] >= 0) {
               se.tau2 <- sqrt(vb.f[p+1,p+1]) * tau2 ### delta rule: vb[p+1,p+1] is the variance of log(tau2), so vb[p+1,p+1] * tau2^2 is the variance of exp(log(tau2))
            } else {
               se.tau2 <- NA
            }
         }

         #return(list(b=b, vb=vb, tau2=tau2, sigma2=sigma2, parms=parms, p.eff=p.eff, k.eff=k.eff, b2.QE=b2.QE, vb2.QE=vb2.QE))

      }

   }

   #########################################################################
   #########################################################################
   #########################################################################

   ### one group outcomes (log odds and log transformed rates)

   if (is.element(measure, c("PLO","IRLN"))) {

      ### prepare data

      if (measure == "PLO") {
         dat.grp <- cbind(xi=xi,mi=mi)
         dat.fam <- binomial
         dat.off <- NULL
      }

      if (measure == "IRLN") {
         dat.grp <- xi
         dat.fam <- poisson
         dat.off <- log(ti)
      }

      study <- factor(seq_len(k)) ### study factor
      X.fit <- X

      #return(data.frame(study, dat.grp, dat.off=ifelse(!is.null(dat.off), dat.off, NA), X.fit))

      ### fit FE model

      if (verbose)
         message("Fitting FE model ...")

      res.FE <- try(glm(dat.grp ~ -1 + X.fit, offset=dat.off, family=dat.fam, control=glmCtrl), silent=!verbose)

      if (inherits(res.FE, "try-error"))
         stop("Cannot fit FE model.")

      ### log-likelihood

      #ll.FE <- with(data.frame(dat.grp), sum(dbinom(xi, xi+mi, predict(res.FE, type="response"), log=TRUE))) ### model has a NULL offset
      #ll.FE <- with(data.frame(dat.grp), sum(dpois(xi, predict(res.FE, type="response"), log=TRUE)))         ### offset already incorporated into predict()
      ll.FE <- c(logLik(res.FE)) ### same as above

      ### fit saturated FE model (= QE model)
      ### notes: 1) suppressWarnings() to suppress warning "glm.fit: fitted probabilities numerically 0 or 1 occurred"

      if (verbose)
         message("Fitting saturated model ...")

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

         warning("Cannot fit saturated model.")
         QEconv <- FALSE
         ll.QE <- NA

      } else {

         QEconv <- TRUE

         ### log-likelihood

         #ll.QE <- with(data.frame(dat.grp), sum(dbinom(xi, xi+mi, xi/(xi+mi), log=TRUE))) ### model has a NULL offset
         #ll.QE <- with(data.frame(dat.grp), sum(dpois(xi, xi, log=TRUE)))                 ### offset not relevant for saturated model
         ll.QE <- c(logLik(res.QE)) ### same as above

         ### extract coefficients and variance-covariance matrix for Wald-type test for heterogeneity

         b2.QE  <- cbind(na.omit(coef(res.QE)[-seq_len(p)]))        ### coef() still includes aliased coefficients as NAs, so have to filter them out
         vb2.QE <- vcov(res.QE)[-seq_len(p),-seq_len(p),drop=FALSE] ### aliased coefficients are already removed by vcov()

      }

      if (method == "ML") {

         ### fit ML model
         ### notes: 1) suppressMessages to suppress the 'one random effect per observation' warning

         if (verbose)
            message("Fitting ML model ...")

         if (verbose) {
            res.ML <- try(lme4::glmer(dat.grp ~ -1 + X.fit + (1 | study), offset=dat.off, family=dat.fam, nAGQ=nAGQ, verbose=verbose, control=do.call(lme4::glmerControl, glmerCtrl)), silent=!verbose)
         } else {
            res.ML <- suppressMessages(try(lme4::glmer(dat.grp ~ -1 + X.fit + (1 | study), offset=dat.off, family=dat.fam, nAGQ=nAGQ, verbose=verbose, control=do.call(lme4::glmerControl, glmerCtrl)), silent=!verbose))
         }

         if (inherits(res.ML, "try-error"))
            stop("Cannot fit ML model.")

         ### log-likelihood

         #ll.ML <- with(data.frame(dat.grp), sum(dbinom(xi, xi+mi, fitted(res.ML), log=TRUE))) ### not correct (since it does not incorporate the random effects; same as ll.FE if tau^2=0)
         #ll.ML <- with(data.frame(dat.grp), sum(dbinom(xi, xi+mi, plogis(qlogis(fitted(res.ML)) + group12*unlist(ranef(res.ML))), log=TRUE))) ### not correct (since one really has to integrate; same as ll.FE if tau^2=0)
         #ll.ML <- c(logLik(res.ML)) ### this is not the same as ll.FE when tau^2 = 0 (not sure why)
         #ll.ML <- ll.QE - 1/2 * lme4::deviance(res.ML) ### this makes ll.ML comparable to ll.FE (same as ll.FE when tau^2=0)
         ll.ML <- ll.QE - 1/2 * deviance(res.ML) ### this makes ll.ML comparable to ll.FE (same as ll.FE when tau^2=0)

      }

      #return(list(res.FE, res.QE, res.ML, ll.FE=ll.FE, ll.QE=ll.QE, ll.ML=ll.ML))
      #res.FE <- res[[1]]; res.QE <- res[[2]]; res.ML <- res[[3]]

      if (method == "FE") {
         b      <- cbind(coef(res.FE)[seq_len(p)])
         vb     <- vcov(res.FE)[seq_len(p),seq_len(p),drop=FALSE]
         tau2   <- 0
         sigma2 <- NA
         parms  <- p
         p.eff  <- p
         k.eff  <- k
      }

      if (method == "ML") {
         b      <- cbind(lme4::fixef(res.ML)[seq_len(p)])
         vb     <- as.matrix(vcov(res.ML))[seq_len(p),seq_len(p),drop=FALSE]
         tau2   <- lme4::VarCorr(res.ML)[[1]][1]
         sigma2 <- NA
         parms  <- p + 1
         p.eff  <- p
         k.eff  <- k
      }

      #return(list(b=b, vb=vb, tau2=tau2, sigma2=sigma2, parms=parms, p.eff=p.eff, k.eff=k.eff, b2.QE=b2.QE, vb2.QE=vb2.QE))

   }

   #########################################################################
   #########################################################################
   #########################################################################

   ### heterogeneity tests (Wald-type and likelihood ratio tests of the extra coefficients in the saturated model)

   if (verbose > 1)
      message("Heterogeneity testing ...")

   if (QEconv) {

      ### for OR, CM.EL, & optim/nlminb/minqa, QE.Wld is already calculated, so skip this part then

      if (measure!="OR" || model!="CM.EL" || !is.element(con$optimizer, c("optim", "nlminb", "minqa"))) {

         if (dim(vb2.QE)[1] > 0) {

            chol.h <- try(chol(vb2.QE), silent=!verbose) ### see if Hessian can be inverted with chol()

            if (inherits(chol.h, "try-error")) {
               warning("Cannot invert Hessian for saturated model.")
               QE.Wld <- NA
            } else {
               QE.Wld <- c(t(b2.QE) %*% chol2inv(chol.h) %*% b2.QE)
            }

         } else {

            QE.Wld <- 0 ### if vb2.QE has 0x0 dims, then fitted model is the saturated model and QE.Wld must be 0

         }

      }

      QE.LRT <- -2 * (ll.FE - ll.QE)

      QE.Wld[QE.Wld < 0] <- 0
      QE.LRT[QE.LRT < 0] <- 0

      #QE.df <- length(b2.QE) ### removed coefficients are not counted if dfs are determined like this
      QE.df <- k-p            ### this yields always the same dfs regardless of how many coefficients are removed
      if (QE.df > 0L) {
         QEp.Wld <- pchisq(QE.Wld, df=QE.df, lower.tail=FALSE)
         QEp.LRT <- pchisq(QE.LRT, df=QE.df, lower.tail=FALSE)
      } else {
         QEp.Wld <- 1
         QEp.LRT <- 1
      }

   } else {

      QE.Wld <- NA
      QE.LRT <- NA
      QEp.Wld <- NA
      QEp.LRT <- NA
      QE.df <- NA

   }

   ### calculation of I^2 and H^2

   wi      <- 1/vi
   W       <- diag(wi, nrow=k.yi, ncol=k.yi)
   stXWX   <- .invcalc(X=X.yi, W=W, k=k.yi)
   P       <- W - W %*% X.yi %*% stXWX %*% crossprod(X.yi,W)
   #vi.avg <- (k-1) / (sum(wi) - sum(wi^2)/sum(wi)) ### this only applies to the RE model
   #vi.avg <- 1/mean(wi) ### harmonic mean of vi's (see Takkouche et al., 1999)
   vi.avg  <- (k.yi-p) / .tr(P)
   I2      <- 100 * tau2 / (vi.avg + tau2)
   H2      <- tau2 / vi.avg + 1

   chol.h <- try(chol(vb[btt,btt]), silent=!verbose) ### see if Hessian can be inverted with chol()

   if (inherits(chol.h, "try-error")) {
      warning("Cannot invert Hessian for QM test.")
      QM <- NA
   } else {
      QM <- c(t(b)[btt] %*% chol2inv(chol.h) %*% b[btt])
   }

   ### scale back b and vb

   if (!int.only && int.incl && con$scale) {
      mX <- rbind(c(1, -1*ifelse(is.d[-1], 0, meanX/sdX)), cbind(0, diag(ifelse(is.d[-1], 1, 1/sdX), nrow=length(is.d)-1, ncol=length(is.d)-1)))
      b  <- mX %*% b
      vb <- mX %*% vb %*% t(mX)
      X  <- Xsave
   }

   rownames(b) <- rownames(vb) <- colnames(vb) <- colnames(X)

   ve <- diag(vb)
   se <- ifelse(ve >= 0, sqrt(ve), NA)
   names(se) <- NULL
   zval <- c(b/se)

   if (knha) {
      dfs <- k-p
      QM  <- QM / m
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

   #return(list(b=b, se=se, zval=zval, ci.lb=ci.lb, ci.ub=ci.ub, vb=vb, tau2=tau2, QM=QM, QMp=QMp))

   #########################################################################

   ###### fit statistics

   if (verbose > 1)
      message("Computing fit statistics and log likelihood ...")

   ll.ML     <- ifelse(method == "FE", ll.FE, ll.ML)
   ll.REML   <- NA
   dev.ML    <- -2 * (ll.ML - ll.QE)
   AIC.ML    <- -2 * ll.ML + 2*parms
   BIC.ML    <- -2 * ll.ML +   parms * log(k.eff)
   AICc.ML   <- -2 * ll.ML + 2*parms * max(k.eff, parms+2) / (max(k.eff, parms+2) - parms - 1)
   dev.REML  <- NA
   AIC.REML  <- NA
   BIC.REML  <- NA
   AICc.REML <- NA

   fit.stats <- matrix(c(ll.ML, dev.ML, AIC.ML, BIC.ML, AICc.ML, ll.REML, dev.REML, AIC.REML, BIC.REML, AICc.REML), ncol=2, byrow=FALSE)
   dimnames(fit.stats) <- list(c("ll","dev","AIC","BIC","AICc"), c("ML","REML"))
   fit.stats <- data.frame(fit.stats)

   #########################################################################

   ###### prepare output

   if (verbose > 1)
      message("Preparing output ...")

   weighted  <- TRUE

   res <- list(b=b, se=se, zval=zval, pval=pval, ci.lb=ci.lb, ci.ub=ci.ub, vb=vb,
               tau2=tau2, se.tau2=se.tau2, sigma2=sigma2,
               k=k, k.f=k.f, k.yi=k.yi, k.eff=k.eff, p=p, p.eff=p.eff, parms=parms, m=m,
               QE.Wld=QE.Wld, QEp.Wld=QEp.Wld, QE.LRT=QE.LRT, QEp.LRT=QEp.LRT, QE.df=QE.df, QM=QM, QMp=QMp, I2=I2, H2=H2,
               int.only=int.only, int.incl=int.incl,
               yi=yi, vi=vi, X=X, yi.f=yi.f, vi.f=vi.f, X.f=X.f,
               ai=ai, bi=bi, ci=ci, di=di, ai.f=ai.f, bi.f=bi.f, ci.f=ci.f, di.f=di.f,
               x1i=x1i, x2i=x2i, t1i=t1i, t2i=t2i, x1i.f=x1i.f, x2i.f=x2i.f, t1i.f=t1i.f, t2i.f=t2i.f,
               xi=xi, mi=mi, ti=ti, xi.f=xi.f, mi.f=mi.f, ti.f=ti.f, ni=ni, ni.f=ni.f,
               ids=ids, not.na=not.na, not.na.yivi=not.na.yivi, slab=slab, slab.null=slab.null,
               measure=measure, method=method, model=model, weighted=weighted, knha=knha, dfs=dfs, btt=btt, intercept=intercept, digits=digits, level=level, control=control, verbose=verbose,
               add=add, to=to, drop00=drop00,
               fit.stats=fit.stats, version=packageVersion("metafor"), call=mf)

   class(res) <- c("rma.glmm", "rma")
   return(res)

}
