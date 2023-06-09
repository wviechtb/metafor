# fixed/random/mixed-effects multivariate/multilevel model with:
#    - possibly one or multiple random intercepts (sigma2) with potentially known correlation matrices
#    - possibly correlated random effects for arms/groups/levels within studies (tau2 and rho for 1st term, gamma2 and phi for 2nd term)
# model also allows for correlated sampling errors via non-diagonal V matrix

# V      = variance-covariance matrix of the sampling errors
# sigma2 = (preset) value(s) for the variance of the random intercept(s)
# tau2   = (preset) value(s) for the variance of the random effects
# rho    = (preset) value(s) for the correlation(s) between random effects
# gamma2 = (preset) value(s) for the variance of the random effects
# phi    = (preset) value(s) for the correlation(s) between random effects

# structures when there is an '~ inner | outer' term in the random argument:
# - CS   (compound symmetry)
# - HCS  (heteroscedastic compound symmetry)
# - UN   (general positive-definite matrix with no structure)
# - UNR  (general positive-definite correlation matrix with a single tau2/gamma2 value)
# - AR   (AR1 structure with a single tau2/gamma2 value and autocorrelation rho/phi)
# - HAR  (heteroscedastic AR1 structure with multiple tau2/gamma2 values and autocorrelation rho/phi)
# - CAR  (continuous time AR1 structure)
# - ID   (same as CS but with rho/phi=0)
# - DIAG (same as HCS but with rho/phi=0)
# - SPEXP/SPGAU/SPLIN/SPRAT/SPSPH (spatial structures: exponential, gaussian, linear, rational quadratic, spherical)
# - GEN (general positive-definite matrix for an arbitrary number of predictors)
# - PHYBM/PHYPL/PHYPD (phylogenetic structures: Brownian motion, Pagel's lambda, Pagel's delta)

rma.mv <- function(yi, V, W, mods, random, struct="CS", intercept=TRUE,
data, slab, subset, method="REML",
test="z", dfs="residual", level=95, btt,
R, Rscale="cor", sigma2, tau2, rho, gamma2, phi,
cvvc=FALSE, sparse=FALSE, verbose=FALSE, digits, control, ...) {

# add ni as argument in the future

   #########################################################################

   ###### setup

   ### check argument specifications

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!is.element(method, c("FE","EE","CE","ML","REML")))
      stop(mstyle$stop("Unknown 'method' specified."))

   if (any(!is.element(struct, c("CS","HCS","UN","AR","HAR","CAR","ID","DIAG","SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","GEN","GDIAG")))) # "UNR", "PHYBM","PHYPL","PHYPD"))))
      stop(mstyle$stop("Unknown 'struct' specified."))

   if (length(struct) == 1L)
      struct <- c(struct, struct)

   na.act <- getOption("na.action")
   on.exit(options(na.action=na.act), add=TRUE)

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   if (missing(random))
      random <- NULL

   if (missing(R))
      R <- NULL

   if (missing(sigma2))
      sigma2 <- NULL

   if (missing(tau2))
      tau2 <- NULL

   if (missing(rho))
      rho <- NULL

   if (missing(gamma2))
      gamma2 <- NULL

   if (missing(phi))
      phi <- NULL

   if (missing(control))
      control <- list()

   ### set defaults for digits

   if (missing(digits)) {
      digits <- .set.digits(dmiss=TRUE)
   } else {
      digits <- .set.digits(digits, dmiss=FALSE)
   }

   time.start <- proc.time()

   ### get ... argument and check for extra/superfluous arguments

   ddd <- list(...)

   .chkdots(ddd, c("tdist", "outlist", "time", "dist", "abbrev", "restart", "beta", "vccon"))

   ### handle 'tdist' argument from ... (note: overrides test argument)

   if (.isFALSE(ddd$tdist))
      test <- "z"
   if (.isTRUE(ddd$tdist))
      test <- "t"

   if (!is.element(test, c("z", "t", "knha", "adhoc")))
      stop(mstyle$stop("Invalid option selected for 'test' argument."))

   if (is.character(dfs))
      dfs <- match.arg(dfs, c("residual", "contain"))

   if (is.numeric(dfs) || (dfs == "contain" && test == "z"))
      test <- "t"

   ### handle Rscale argument (either character, logical, or integer)

   if (is.character(Rscale))
      Rscale <- match.arg(Rscale, c("none", "cor", "cor0", "cov0"))

   if (is.logical(Rscale))
      Rscale <- ifelse(Rscale, "cor", "none")

   if (is.numeric(Rscale)) {
      Rscale <- round(Rscale)
      if (Rscale > 3 | Rscale < 0)
         stop(mstyle$stop("Unknown 'Rscale' value specified."))
      Rscale <- switch(as.character(Rscale), "0"="none", "1"="cor", "2"="cor0", "3"="cov0")
   }

   ### handle 'dist' argument from ...

   if (!is.null(ddd$dist)) {

      if (is.data.frame(ddd$dist) || .is.matrix(ddd$dist))
         ddd$dist <- list(ddd$dist)

      if (!inherits(ddd$dist, "list"))
         ddd$dist <- as.list(ddd$dist)

      if (length(ddd$dist) == 1L)
         ddd$dist <- c(ddd$dist, ddd$dist)

      dist.methods <- c("euclidean", "maximum", "manhattan", "gcd")

      for (j in 1:2) {

         if (is.data.frame(ddd$dist[[j]]))
            ddd$dist[[j]] <- as.matrix(ddd$dist[[j]])

         if (!is.function(ddd$dist[[j]]) && !.is.matrix(ddd$dist[[j]])) {
            ddd$dist[[j]] <- charmatch(ddd$dist[[j]], dist.methods, nomatch = 0)
            if (ddd$dist[[j]] == 0) {
               stop(mstyle$stop("Argument 'dist' must be one of 'euclidean', 'maximum', 'manhattan', or 'gcd'."))
            } else {
               ddd$dist[[j]] <- dist.methods[ddd$dist[[j]]]
            }
         }

      }

      if (any(ddd$dist == "gcd")) {
         if (!requireNamespace("sp", quietly=TRUE))
            stop(mstyle$stop("Please install the 'sp' package to compute great-circle distances."))
      }

   } else {
      ddd$dist <- list("euclidean", "euclidean")
   }

   if (is.null(ddd$vccon)) {
      vccon <- NULL
   } else {
      vccon <- ddd$vccon
      sigma2 <- .chkvccon(vccon$sigma2, sigma2)
      tau2   <- .chkvccon(vccon$tau2,   tau2)
      rho    <- .chkvccon(vccon$rho,    rho)
      gamma2 <- .chkvccon(vccon$gamma2, gamma2)
      phi    <- .chkvccon(vccon$phi,    phi)
   }

   ### set defaults for formulas

   formula.yi <- NULL
   formula.mods <- NULL

   ### in case user specifies v (instead of V), verbose is set to v, which is non-sensical
   ### - if v is set to the name of a variable in 'data', it won't be found; can check for
   ###   this with try() and inherits(verbose, "try-error")
   ### - if v is set to vi or var (or anything else that might be interpreted as a function),
   ###   then can catch this by checking if verbose is a function

   verbose <- try(verbose, silent=TRUE)

   if (inherits(verbose, "try-error") || is.function(verbose) || length(verbose) != 1L || !(is.logical(verbose) || is.numeric(verbose)))
      stop(mstyle$stop("Argument 'verbose' must be a scalar (logical or numeric/integer)."))

   ### set options(warn=1) if verbose > 2

   if (verbose > 2) {
      opwarn <- options(warn=1)
      on.exit(options(warn=opwarn$warn), add=TRUE)
   }

   #########################################################################

   if (verbose > 1) .space()

   if (verbose > 1)
      message(mstyle$message("Extracting yi/V values ..."))

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

   ### extract yi, V, W, ni, slab, subset, and mods values, possibly from the data frame specified via data (arguments not specified are NULL)

   yi     <- .getx("yi",     mf=mf, data=data)
   V      <- .getx("V",      mf=mf, data=data)
   W      <- .getx("W",      mf=mf, data=data)
   ni     <- .getx("ni",     mf=mf, data=data) ### not yet possible to specify this
   slab   <- .getx("slab",   mf=mf, data=data)
   subset <- .getx("subset", mf=mf, data=data)
   mods   <- .getx("mods",   mf=mf, data=data)

   ### if yi is a formula, extract yi and X (this overrides anything specified via the mods argument further below)

   if (inherits(yi, "formula")) {
      formula.yi <- yi
      formula.mods <- formula.yi[-2]
      options(na.action = "na.pass")                   ### set na.action to na.pass, so that NAs are not filtered out (we'll do that later)
      mods <- model.matrix(yi, data=data)              ### extract model matrix (now mods is no longer a formula, so [a] further below is skipped)
      attr(mods, "assign") <- NULL                     ### strip assign attribute (not needed at the moment)
      attr(mods, "contrasts") <- NULL                  ### strip contrasts attribute (not needed at the moment)
      yi <- model.response(model.frame(yi, data=data)) ### extract yi values from model frame
      options(na.action = na.act)                      ### set na.action back to na.act
      names(yi) <- NULL                                ### strip names (1:k) from yi (so res$yi is the same whether yi is a formula or not)
      intercept <- FALSE                               ### set to FALSE since formula now controls whether the intercept is included or not
   }                                                   ### note: code further below ([b]) actually checks whether intercept is included or not

   ### in case user passed a data frame to yi, convert it to a vector (if possible)

   if (is.data.frame(yi)) {
      if (ncol(yi) == 1L) {
         yi <- yi[[1]]
      } else {
         stop(mstyle$stop("The object/variable specified for the 'yi' argument is a data frame with multiple columns."))
      }
   }

   ### in case user passed a matrix to yi, convert it to a vector (if possible)

   if (.is.matrix(yi)) {
      if (nrow(yi) == 1L || ncol(yi) == 1L) {
         yi <- as.vector(yi)
      } else {
         stop(mstyle$stop("The object/variable specified for the 'yi' argument is a matrix with multiple rows/columns."))
      }
   }

   ### check if yi is numeric

   if (!is.numeric(yi))
      stop(mstyle$stop("The object/variable specified for the 'yi' argument is not numeric."))

   ### number of outcomes before subsetting

   k <- length(yi)
   k.all <- k

   ### set default measure argument

   measure <- "GEN"

   if (!is.null(attr(yi, "measure"))) ### take 'measure' from yi (if it is there)
      measure <- attr(yi, "measure")

   ### add measure attribute (back) to the yi vector

   attr(yi, "measure") <- measure

   ### some checks on V (and turn V into a diagonal matrix if it is a column/row vector)

   if (is.null(V))
      stop(mstyle$stop("Must specify 'V' argument."))

   ### catch cases where 'V' is the utils::vi() function

   if (identical(V, utils::vi))
      stop(mstyle$stop("Variable specified for 'V' argument cannot be found."))

   if (is.list(V) && !is.data.frame(V)) {

      ### list elements may be data frames (or scalars), so coerce to matrices

      V <- lapply(V, as.matrix)

      ### check that all elements are square

      if (any(!sapply(V, .is.square)))
         stop(mstyle$stop("All list elements in 'V' must be square matrices."))

      ### turn list into block-diagonal (sparse) matrix

      if (sparse) {
         V <- bdiag(V)
      } else {
         V <- bldiag(V)
      }

   }

   ### check if user constrained V to 0 (can skip a lot of the steps below then)

   if ((.is.vector(V) && length(V) == 1L && V == 0) || (.is.vector(V) && length(V) == k && !anyNA(V) && all(V == 0))) {
      V0 <- TRUE
   } else {
      V0 <- FALSE
   }

   ### turn V into a diagonal matrix if it is a column/row vector
   ### note: if V is a scalar (e.g., V=0), then this will turn V into a kxk
   ### matrix with the value of V along the diagonal

   if (V0 || .is.vector(V) || nrow(V) == 1L || ncol(V) == 1L) {
      if (sparse) {
         V <- Diagonal(k, as.vector(V))
      } else {
         V <- diag(as.vector(V), nrow=k, ncol=k)
      }
   }

   ### turn V into a matrix if it is a data frame

   if (is.data.frame(V))
      V <- as.matrix(V)

   ### remove row and column names (important for isSymmetric() function)
   ### (but only do this if V has row/column names to avoid making an unnecessary copy)

   if (!is.null(dimnames(V)))
      V <- unname(V)

   ### check whether V is square and symmetric (can skip when V0)

   if (!V0 && !.is.square(V))
      stop(mstyle$stop("'V' must be a square matrix."))

   if (!V0 && !isSymmetric(V)) ### note: copy of V is made when doing this
      stop(mstyle$stop("'V' must be a symmetric matrix."))

   ### check length of yi and V

   if (nrow(V) != k)
      stop(mstyle$stop(paste0("Length of 'yi' (", k, ") and length/dimensions of 'V' (", nrow(V), ") is not the same.")))

   ### force V to be sparse when sparse=TRUE (and V is not yet sparse)

   if (sparse && inherits(V, "matrix"))
      V <- Matrix(V, sparse=TRUE)

   ### check if V is numeric (but only for 'regular' matrices, since this is always FALSE for sparse matrices)

   if (inherits(V, "matrix") && !is.numeric(V))
      stop(mstyle$stop("The object/variable specified for the 'V' argument is not numeric."))

   ### process W if it was specified

   if (!is.null(W)) {

      ### turn W into a diagonal matrix if it is a column/row vector
      ### in general, turn W into A (arbitrary weight matrix)

      if (.is.vector(W) || nrow(W) == 1L || ncol(W) == 1L) {

         W <- as.vector(W)

         ### allow easy setting of W to a single value

         if (length(W) == 1L)
            W <- rep(W, k)

         A <- diag(W, nrow=length(W), ncol=length(W))

      } else {

         A <- W

      }

      if (is.data.frame(A))
         A <- as.matrix(A)

      ### remove row and column names (important for isSymmetric() function)
      ### (but only do this if A has row/column names to avoid making an unnecessary copy)

      if (!is.null(dimnames(A)))
         A <- unname(A)

      ### check whether A is square and symmetric

      if (!.is.square(A))
         stop(mstyle$stop("'W' must be a square matrix."))

      if (!isSymmetric(A))
         stop(mstyle$stop("'W' must be a symmetric matrix."))

      ### check length of yi and A

      if (nrow(A) != k)
         stop(mstyle$stop(paste0("Length of 'yi' (", k, ") and length/dimensions of 'W' (", nrow(A), ") is not the same.")))

      ### force A to be sparse when sparse=TRUE (and A is not yet sparse)

      if (sparse && inherits(A, "matrix"))
         A <- Matrix(A, sparse=TRUE)

      if (inherits(A, "matrix") && !is.numeric(A))
         stop(mstyle$stop("The object/variable specified for the 'W' argument is not numeric."))

   } else {

      A <- NULL

   }

   ### if ni has not been specified (and hence is NULL), try to get it from the attributes of yi
   ### note: currently ni argument removed, so this is the only way to pass ni to the function

   if (is.null(ni))
      ni <- attr(yi, "ni")

   ### check length of yi and ni
   ### if there is a mismatch, then ni cannot be trusted, so set it to NULL

   if (!is.null(ni) && length(ni) != k)
      ni <- NULL

   ### if ni is now available, add it (back) as an attribute to yi
   ### this is currently pointless, but may be useful if function has an ni argument

   #if (!is.null(ni))
   #   attr(yi, "ni") <- ni

   #########################################################################

   if (verbose > 1)
      message(mstyle$message("Creating model matrix ..."))

   ### convert mods formula to X matrix and set intercept equal to FALSE
   ### skipped if formula has already been specified via yi argument, since mods is then no longer a formula (see [a])

   if (inherits(mods, "formula")) {
      formula.mods <- mods
      if (isTRUE(all.equal(formula.mods, ~ 1))) { # needed so 'mods = ~ 1' without 'data' specified works
         mods <- matrix(1, nrow=k, ncol=1)
         intercept <- FALSE
      } else {
         options(na.action = "na.pass")        ### set na.action to na.pass, so that NAs are not filtered out (we'll do that later)
         mods <- model.matrix(mods, data=data) ### extract model matrix
         attr(mods, "assign") <- NULL          ### strip assign attribute (not needed at the moment)
         attr(mods, "contrasts") <- NULL       ### strip contrasts attribute (not needed at the moment)
         options(na.action = na.act)           ### set na.action back to na.act
         intercept <- FALSE                    ### set to FALSE since formula now controls whether the intercept is included or not
      }                                        ### note: code further below ([b]) actually checks whether intercept is included or not
   }

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

   #########################################################################
   #########################################################################
   #########################################################################

   ### process random argument

   if (!is.element(method, c("FE","EE","CE")) && !is.null(random)) {

      if (verbose > 1)
         message(mstyle$message("Processing 'random' argument ..."))

      ### make sure random argument is always a list (so lapply() below works)

      if (!is.list(random))
         random <- list(random)

      ### check that all elements are formulas

      if (any(sapply(random, function(x) !inherits(x, "formula"))))
         stop(mstyle$stop("All elements of 'random' must be formulas."))

      ### check that all formulas have a vertical bar

      has.vbar <- sapply(random, function(f) grepl("|", paste0(f, collapse=""), fixed=TRUE))

      if (any(!has.vbar))
         stop(mstyle$stop("All formulas in 'random' must contain a grouping variable after the | symbol."))

      ### check if any formula have a $

      has.dollar <- sapply(random, function(f) grepl("$", paste0(f, collapse=""), fixed=TRUE))

      if (any(has.dollar))
         stop(mstyle$stop("Cannot use '$' notation in formulas in the 'random' argument (use the 'data' argument instead)."))

      ### check if any formula have a :

      has.colon <- sapply(random, function(f) grepl(":", paste0(f, collapse=""), fixed=TRUE))

      if (any(has.colon))
         stop(mstyle$stop("Cannot use ':' notation in formulas in the 'random' argument (use 'interaction()' instead)."))

      ### check if any formula have a %in%

      has.in <- sapply(random, function(f) grepl("%in%", paste0(f, collapse=""), fixed=TRUE))

      if (any(has.in))
         stop(mstyle$stop("Cannot use '%in%' notation in formulas in the 'random' argument (use 'interaction()' instead)."))

      ### check which formulas have a ||

      has.dblvbar <- sapply(random, function(f) grepl("||", paste0(f, collapse=""), fixed=TRUE))

      ### replace || with |

      random <- lapply(random, function(f) {
         if (grepl("||", paste0(f, collapse=""), fixed=TRUE)) {
            f <- paste0(f, collapse="")
            f <- gsub("||", "|", f, fixed=TRUE)
            f <- as.formula(f)
         }
         return(f)
      })

      ### check which formulas in random are '~ inner | outer' formulas

      formulas <- list(NULL, NULL)
      split.formulas <- sapply(random, function(f) strsplit(paste0(f, collapse=""), " | ", fixed=TRUE))
      is.inner.outer <- sapply(split.formulas, function(f) f[1] != "~1")

      ### make sure that there are only up to two '~ inner | outer' formulas

      if (sum(is.inner.outer) > 2)
         stop(mstyle$stop("Only up to two '~ inner | outer' formulas allowed in the 'random' argument."))

      ### get '~ inner | outer' formulas

      if (any(is.inner.outer))
         formulas[[1]] <- random[is.inner.outer][1][[1]]
      if (sum(is.inner.outer) == 2)
         formulas[[2]] <- random[is.inner.outer][2][[1]]

      ### figure out if a formulas has a slash (as in '~ 1 | study/id')

      has.slash <- sapply(random, function(f) grepl("/", paste0(f, collapse=""), fixed=TRUE))

      ### check if slash is used in combination with an '~ inner | outer' term

      if (any(is.inner.outer & has.slash))
         stop(mstyle$stop("Cannot use '~ inner | outer1/outer2' type terms in the 'random' argument."))

      ### substitute + for | in all formulas (so that model.frame() below works)

      random.plus <- lapply(random, function(f) formula(sub("\\|", "+", paste0(f, collapse=""))))

      ### get all model frames corresponding to the formulas in the random argument
      ### mf.r <- lapply(random, get_all_vars, data=data)
      ### note: get_all_vars() does not carry out any functions calls within the formula
      ###       so use model.frame(), which allows for things like 'random = ~ factor(arm) | study'
      ###       need to use na.pass so that NAs are passed through (checks for NAs are done later)

      #mf.r <- lapply(random.plus, model.frame, data=data, na.action=na.pass)

      mf.r <- list()

      io <- 0

      for (j in seq_along(is.inner.outer)) {

         if (is.inner.outer[j]) {

            io <- io + 1

            ### for an '~ inner | outer' term with struct="GEN", expand the inner formula to the
            ### model matrix and re-combine this with the outer variable

            if (is.element(struct[io], c("GEN","GDIAG"))) {

               f.inner <- as.formula(strsplit(paste(random[[j]], collapse=""), " | ", fixed=TRUE)[[1]][1])
               f.outer <- as.formula(paste("~", strsplit(paste(random[[j]], collapse=""), " | ", fixed=TRUE)[[1]][2]))
               options(na.action = "na.pass")
               X.inner <- model.matrix(f.inner, data=data)
               options(na.action = na.act)
               is.int <- apply(X.inner, 2, .is.intercept)
               colnames(X.inner)[is.int] <- "intrcpt"
               mf.r[[j]] <- cbind(X.inner, model.frame(f.outer, data=data, na.action=na.pass))

               if (has.dblvbar[j]) # change "GEN" to "GDIAG" if the formula had a ||
                  struct[io] <- "GDIAG"

            } else {

               mf.r[[j]] <- model.frame(random.plus[[j]], data=data, na.action=na.pass)

            }

         } else {

            mf.r[[j]] <- model.frame(random.plus[[j]], data=data, na.action=na.pass)

         }

      }

      ### count number of columns in each model frame

      mf.r.ncols <- sapply(mf.r, ncol)

      ### for formulas with slashes, create interaction terms

      for (j in seq_along(has.slash)) {

         if (!has.slash[j])
            next

         ### need to go backwards; otherwise, with 3 or more terms (e.g., ~ 1 | var1/var2/var3), the third term would be an
         ### interaction between var1, var1:var2, and var3; by going backwards, we get var1, var1:var2, and var1:var2:var3

         for (p in mf.r.ncols[j]:1) {
            mf.r[[j]][,p] <- interaction(mf.r[[j]][1:p], drop=TRUE, lex.order=TRUE, sep = "/")
            colnames(mf.r[[j]])[p] <- paste(colnames(mf.r[[j]])[1:p], collapse="/")
         }

      }

      ### create list where model frames with multiple columns based on slashes are flattened out

      if (any(has.slash)) {

         if (length(mf.r) == 1L) {

            ### if formula only has one element of the form ~ 1 | var1/var2/..., create a list of the data frames (each with one column)

            mf.r <- lapply(seq(ncol(mf.r[[1]])), function(x) mf.r[[1]][x])

         } else {

            ### if there are non-slash elements, then this flattens things out (obviously ...)

            mf.r <- unlist(mapply(function(mf, sl) if (sl) lapply(seq(mf), function(x) mf[x]) else list(mf), mf.r, has.slash, SIMPLIFY=FALSE), recursive=FALSE, use.names=FALSE)

         }

         ### recount number of columns in each model frame

         mf.r.ncols <- sapply(mf.r, ncol)

      }

      #return(mf.r)

      ### separate mf.r into mf.s (~ 1 | id), mf.g (~ inner | outer), and mf.h (~ inner | outer) parts

      mf.s <- mf.r[which(mf.r.ncols == 1)]      ### if there is no '~ 1 | factor' term, this is list() ([] so that we get a list of data frames)
      mf.g <- mf.r[[which(mf.r.ncols >= 2)[1]]] ### if there is no 1st '~ inner | outer' terms, this is NULL ([[]] so that we get a data frame, not a list)
      mf.h <- mf.r[[which(mf.r.ncols >= 2)[2]]] ### if there is no 2nd '~ inner | outer' terms, this is NULL ([[]] so that we get a data frame, not a list)

      ### if there is no (~ 1 | factor) term, then mf.s is list(), so turn that into NULL

      if (length(mf.s) == 0L)
         mf.s <- NULL

      ### does the random argument include at least one (~ 1 | id) term?

      withS <- !is.null(mf.s)

      ### does the random argument include '~ inner | outer' terms?

      withG <- !is.null(mf.g)
      withH <- !is.null(mf.h)

      ### count number of rows in each model frame

      mf.r.nrows <- sapply(mf.r, nrow)

      ### make sure that rows in each model frame match the length of the data

      if (any(mf.r.nrows != k))
         stop(mstyle$stop("Length of variables specified via the 'random' argument does not match length of the data."))

      ### need this for profile(); with things like 'random = ~ factor(arm) | study', 'mf.r' contains variables 'factor(arm)' and 'study'
      ### but the former won't work when using the same formula for the refitting (same when using interaction() in the random formula)
      ### note: with ~ 1 | interaction(var1, var2), mf.r will have 2 columns, but is actually a 'one variable' term
      ###   and with ~ interaction(var1, var2) | var3, mf.r will have 3 columns, but is actually a 'two variable' term
      ### mf.r.ncols above is correct even in these cases (since it is based on the model.frame() results), but need
      ### to be careful that this doesn't screw up anything in other functions (for now, mf.r.ncols is not used in any other function)

      mf.r <- lapply(random.plus, get_all_vars, data=data)

   } else {

      ### set defaults for some elements when method="FE/EE/CE"

      formulas <- list(NULL, NULL)
      mf.r  <- NULL
      mf.s  <- NULL
      mf.g  <- NULL
      mf.h  <- NULL
      withS <- FALSE
      withG <- FALSE
      withH <- FALSE

   }

   ### warn that 'struct' argument is disregarded if it has been specified, but model contains no '~ inner | outer' terms

   if (!withG && "struct" %in% names(mf))
      warning(mstyle$warning("Model does not contain an '~ inner | outer' term, so 'struct' argument is disregaded."), call.=FALSE)

   ### warn that 'random' argument is disregarded if it has been specified, but method="FE/EE/CE"

   if (is.element(method, c("FE","EE","CE")) && "random" %in% names(mf))
      warning(mstyle$warning(paste0("The 'random' argument is disregaded when method=\"", method, "\".")), call.=FALSE)

   #return(list(mf.r=mf.r, mf.s=mf.s, mf.g=mf.g, mf.h=mf.h))

   ### note: checks on NAs in mf.s, mf.g, and mf.h after subsetting (since NAs may be removed by subsetting)

   #########################################################################
   #########################################################################
   #########################################################################

   ### generate study labels if none are specified (or none can be found in yi argument)

   if (verbose > 1)
      message(mstyle$message("Generating/extracting study labels ..."))

   ### study ids (1:k sequence before subsetting)

   ids <- seq_len(k)

   ### if slab has not been specified but is an attribute of yi, get it

   if (is.null(slab)) {

      slab <- attr(yi, "slab") # will be NULL if there is no slab attribute

      ### check length of yi and slab (only if slab is now not NULL)
      ### if there is a mismatch, then slab cannot be trusted, so set it to NULL

      if (!is.null(slab) && length(slab) != k)
         slab <- NULL

   }

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

      if (verbose > 1)
         message(mstyle$message("Subsetting ..."))

      subset <- .chksubset(subset, k)

      yi   <- .getsubset(yi,   subset)
      V    <- .getsubset(V,    subset, col=TRUE)
      A    <- .getsubset(A,    subset, col=TRUE)
      ni   <- .getsubset(ni,   subset)
      mods <- .getsubset(mods, subset)
      slab <- .getsubset(slab, subset)
      mf.r <- lapply(mf.r, .getsubset, subset)
      mf.s <- lapply(mf.s, .getsubset, subset)
      mf.g <- .getsubset(mf.g, subset)
      mf.h <- .getsubset(mf.h, subset)
      ids  <- .getsubset(ids,  subset)

      k <- length(yi)

      attr(yi, "measure") <- measure ### add measure attribute back
      attr(yi, "ni")      <- ni      ### add ni attribute back

   }

   ### check if study labels are unique; if not, make them unique

   if (anyDuplicated(slab))
      slab <- .make.unique(slab)

   ### add slab attribute back

   attr(yi, "slab") <- slab

   ### get the sampling variances from the diagonal of V

   vi <- diag(V)

   ### save full data (including potential NAs in yi/vi/V/W/ni/mods)

   yi.f   <- yi
   vi.f   <- vi
   V.f    <- V
   W.f    <- A
   ni.f   <- ni
   mods.f <- mods
   #mf.g.f <- mf.g ### copied further below
   #mf.h.f <- mf.h ### copied further below
   #mf.s.f <- mf.s ### copied further below

   k.f <- k ### total number of observed outcomes including all NAs

   #########################################################################
   #########################################################################
   #########################################################################

   ### stuff that need to be done after subsetting

   if (withS) {

      if (verbose > 1)
         message(mstyle$message(paste0("Processing '", paste0("~ 1 | ", sapply(mf.s, names), collapse=", "), "' term(s) ...")))

      ### get variables names in mf.s

      s.names <- sapply(mf.s, names) ### one name per term

      ### turn each variable in mf.s into a factor (and turn each column vector into just a vector)
      ### if a variable was a factor to begin with, this drops any unused levels, but order of existing levels is preserved

      mf.s <- lapply(mf.s, function(x) factor(x[[1]]))

      ### check if there are any NAs anywhere in mf.s

      if (any(sapply(mf.s, anyNA)))
         stop(mstyle$stop("No NAs allowed in variables specified in the 'random' argument."))

      ### how many (~ 1 | id) terms does the random argument include? (0 if none, but if withS is TRUE, must be at least 1)

      sigma2s <- length(mf.s)

      ### set default value(s) for sigma2 argument if it is unspecified

      if (is.null(sigma2))
         sigma2 <- rep(NA_real_, sigma2s)

      ### allow quickly setting all sigma2 values to a fixed value

      if (length(sigma2) == 1L)
         sigma2 <- rep(sigma2, sigma2s)

      ### check if sigma2 is of the correct length

      if (length(sigma2) != sigma2s)
         stop(mstyle$stop(paste0("Length of 'sigma2' argument (", length(sigma2), ") does not match actual number of variance components (", sigma2s, ").")))

      ### checks on any fixed values of sigma2 argument

      if (any(sigma2 < 0, na.rm=TRUE))
         stop(mstyle$stop("Specified value(s) of 'sigma2' must be non-negative."))

      ### get number of levels of each variable in mf.s (vector with one value per term)

      s.nlevels <- sapply(mf.s, nlevels)

      ### get levels of each variable in mf.s (list with levels for each variable)

      s.levels <- lapply(mf.s, levels)

      ### checks on R (note: do this after subsetting, so user can filter out ids with no info in R)

      if (is.null(R)) {

         withR <- FALSE
         Rfix  <- rep(FALSE, sigma2s)

      } else {

         if (verbose > 1)
            message(mstyle$message("Processing 'R' argument ..."))

         withR <- TRUE

         ### make sure R is always a list (so lapply() below works)

         if (is.data.frame(R) || !is.list(R))
            R <- list(R)

         ### check if R list has no names at all or some names are missing
         ### (if only some elements of R have names, then names(R) is "" for the unnamed elements, so use nchar()==0 to check for that)

         if (is.null(names(R)) || any(nchar(names(R)) == 0L))
            stop(mstyle$stop("Argument 'R' must be a *named* list."))

         ### remove elements in R that are NULL (not sure why this is needed; why would anybody ever do this?)
         ### maybe this had something to do with functions that repeatedly call rma.mv(); so leave this be for now

         R <- R[!sapply(R, is.null)]

         ### turn all elements in R into matrices (this would fail with a NULL element)

         R <- lapply(R, as.matrix)

         ### match up R matrices based on the s.names (and correct names of R)
         ### so if a particular ~ 1 | id term has a matching id=R element, the corresponding R element is that R matrix
         ###    if a particular ~ 1 | id term does not have a matching id=R element, the corresponding R element is NULL

         R <- R[s.names]

         ### NULL elements in R would have no name, so this makes sure that all R elements have the correct s.names

         names(R) <- s.names

         ### check for which components an R matrix has been specified

         Rfix <- !sapply(R, is.null)

         ### Rfix could be all FALSE (if user has used id names in R that are not actually in 'random')
         ### so only do the rest below if that is *not* the case

         if (any(Rfix)) {

            ### check if given R matrices are square and symmetric

            if (any(!sapply(R[Rfix], .is.square)))
               stop(mstyle$stop("Elements of 'R' must be square matrices."))
            if (any(!sapply(R[Rfix], function(x) isSymmetric(unname(x)))))
               stop(mstyle$stop("Elements of 'R' must be symmetric matrices."))

            for (j in seq_along(R)) {

               if (!Rfix[j])
                  next

               ### even if isSymmetric() is TRUE, there may still be minor numerical differences between the lower and upper triangular
               ### parts that could lead to isSymmetric() being FALSE once we do any potentially rescaling of the R matrices further
               ### below; this ensures strict symmetry to avoid this issue
               #R[[j]][lower.tri(R[[j]])] <- t(R[[j]])[lower.tri(R[[j]])]
               R[[j]] <- symmpart(R[[j]])

               ### if rownames are missing, copy colnames to rownames and vice-versa

               if (is.null(rownames(R[[j]])))
                  rownames(R[[j]]) <- colnames(R[[j]])
               if (is.null(colnames(R[[j]])))
                  colnames(R[[j]]) <- rownames(R[[j]])

               ### if colnames are still missing at this point, R element did not have dimension names to begin with

               if (is.null(colnames(R[[j]])))
                  stop(mstyle$stop("Elements of 'R' must have dimension names."))

            }

            ### if user specifies the entire (k x k) correlation matrix, this removes the duplicate rows/columns

            #R[Rfix] <- lapply(R[Rfix], unique, MARGIN=1)
            #R[Rfix] <- lapply(R[Rfix], unique, MARGIN=2)

            ### no, the user can specify an entire (k x k) matrix; the problem is repeated dimension names
            ### so let's filter out rows/columns with the same dimension names

            R[Rfix] <- lapply(R[Rfix], function(x) x[!duplicated(rownames(x)), !duplicated(colnames(x)), drop=FALSE])

            ### after the two commands above, this should always be FALSE, but leave for now just in case

            if (any(sapply(R[Rfix], function(x) length(colnames(x)) != length(unique(colnames(x))))))
               stop(mstyle$stop("Each element of 'R' must have unique dimension names."))

            ### check for R being positive definite
            ### skipped: even if R is not positive definite, the marginal var-cov matrix can still be; so just check for pd during optimization

            #if (any(sapply(R[Rfix], !.chkpd)))
            #   stop(mstyle$stop("Matrix in R is not positive definite."))

            for (j in seq_along(R)) {

               if (!Rfix[j])
                  next

               ### check if there are NAs in a matrix specified via R

               if (anyNA(R[[j]]))
                  stop(mstyle$stop("No missing values allowed in matrices specified via 'R'."))

               ### check if there are levels in s.levels which are not in R (if yes, issue an error and stop)

               if (any(!is.element(s.levels[[j]], colnames(R[[j]]))))
                  stop(mstyle$stop(paste0("There are levels in '", s.names[j], "' for which there are no matching rows/columns in the corresponding 'R' matrix.")))

               ### check if there are levels in R which are not in s.levels (if yes, issue a warning)

               if (any(!is.element(colnames(R[[j]]), s.levels[[j]])))
                  warning(mstyle$warning(paste0("There are rows/columns in the 'R' matrix for '", s.names[j], "' for which there are no data.")))

            }

         } else {

            warning(mstyle$warning("Argument 'R' specified, but list name(s) not in 'random'."), call.=FALSE)

            withR <- FALSE
            Rfix  <- rep(FALSE, sigma2s)
            R     <- NULL

         }

      }

   } else {

      ### need one fixed sigma2 value for optimization function

      sigma2s <- 1
      sigma2  <- 0

      s.nlevels <- NULL
      s.levels  <- NULL
      s.names   <- NULL

      withR   <- FALSE
      Rfix    <- FALSE
      R       <- NULL

   }

   #mf.s.f <- mf.s ### not needed at the moment

   ### copy s.nlevels and s.levels (needed for ranef())

   s.nlevels.f <- s.nlevels
   s.levels.f  <- s.levels

   #########################################################################

   ### stuff that need to be done after subsetting

   if (withG) {

      tmp <- .process.G.aftersub(mf.g, struct[1], formulas[[1]], tau2, rho, isG=TRUE, k, sparse, verbose)

      mf.g      <- tmp$mf.g
      g.names   <- tmp$g.names
      g.nlevels <- tmp$g.nlevels
      g.levels  <- tmp$g.levels
      g.values  <- tmp$g.values
      tau2s     <- tmp$tau2s
      rhos      <- tmp$rhos
      tau2      <- tmp$tau2
      rho       <- tmp$rho
      Z.G1      <- tmp$Z.G1
      Z.G2      <- tmp$Z.G2

   } else {

      ### need one fixed tau2 and rho value for optimization function

      tau2s <- 1
      rhos  <- 1
      tau2  <- 0
      rho   <- 0

      ### need Z.G1 and Z.G2 to exist further below and for optimization function

      Z.G1 <- NULL
      Z.G2 <- NULL
      g.nlevels <- NULL
      g.levels  <- NULL
      g.values  <- NULL

      g.names <- NULL

   }

   mf.g.f <- mf.g ### needed for predict()

   #########################################################################

   ### stuff that need to be done after subsetting

   if (withH) {

      tmp <- .process.G.aftersub(mf.h, struct[2], formulas[[2]], gamma2, phi, isG=FALSE, k, sparse, verbose)

      mf.h      <- tmp$mf.g
      h.names   <- tmp$g.names
      h.nlevels <- tmp$g.nlevels
      h.levels  <- tmp$g.levels
      h.values  <- tmp$g.values
      gamma2s   <- tmp$tau2s
      phis      <- tmp$rhos
      gamma2    <- tmp$tau2
      phi       <- tmp$rho
      Z.H1      <- tmp$Z.G1
      Z.H2      <- tmp$Z.G2

   } else {

      ### need one fixed gamma2 and phi value for optimization function

      gamma2s <- 1
      phis    <- 1
      gamma2  <- 0
      phi     <- 0

      ### need Z.H1 and Z.H2 to exist further below and for optimization function

      Z.H1 <- NULL
      Z.H2 <- NULL
      h.nlevels <- NULL
      h.levels  <- NULL
      h.values  <- NULL

      h.names <- NULL

   }

   mf.h.f <- mf.h ### needed for predict()

   # return(list(Z.G1=Z.G1, Z.G2=Z.G2, g.nlevels=g.nlevels, g.levels=g.levels, g.values=g.values, tau2=tau2, rho=rho,
   #             Z.H1=Z.H1, Z.H2=Z.H2, h.nlevels=h.nlevels, h.levels=h.levels, h.values=h.values, gamma2=gamma2, phi=phi))

   #########################################################################
   #########################################################################
   #########################################################################

   ### check for NAs and act accordingly

   has.na <- is.na(yi) | (if (is.null(mods)) FALSE else apply(is.na(mods), 1, any)) | (if (V0) FALSE else .anyNAv(V)) | (if (is.null(A)) FALSE else apply(is.na(A), 1, any))
   not.na <- !has.na

   if (any(has.na)) {

      if (verbose > 1)
         message(mstyle$message("Handling NAs ..."))

      if (na.act == "na.omit" || na.act == "na.exclude" || na.act == "na.pass") {

         yi   <- yi[not.na]
         V    <- V[not.na,not.na,drop=FALSE]
         A    <- A[not.na,not.na,drop=FALSE]
         vi   <- vi[not.na]
         ni   <- ni[not.na]
         mods <- mods[not.na,,drop=FALSE]
         mf.r <- lapply(mf.r, function(x) x[not.na,,drop=FALSE])
         mf.s <- lapply(mf.s, function(x) x[not.na]) ### note: mf.s is a list of vectors at this point
         mf.g <- mf.g[not.na,,drop=FALSE]
         mf.h <- mf.h[not.na,,drop=FALSE]
         if (is.element(struct[1], c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYBM","PHYPL","PHYPD"))) {
            Z.G1 <- Z.G1[not.na,not.na,drop=FALSE]
         } else {
            Z.G1 <- Z.G1[not.na,,drop=FALSE]
         }
         Z.G2 <- Z.G2[not.na,,drop=FALSE]
         if (is.element(struct[2], c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYBM","PHYPL","PHYPD"))) {
            Z.H1 <- Z.H1[not.na,not.na,drop=FALSE]
         } else {
            Z.H1 <- Z.H1[not.na,,drop=FALSE]
         }
         Z.H2 <- Z.H2[not.na,,drop=FALSE]
         k    <- length(yi)
         warning(mstyle$warning(paste(sum(has.na), ifelse(sum(has.na) > 1, "rows", "row"), "with NAs omitted from model fitting.")), call.=FALSE)

         attr(yi, "measure") <- measure ### add measure attribute back
         attr(yi, "ni")      <- ni      ### add ni attribute back

         ### note: slab is always of the same length as the full yi vector (after subsetting), so missings are not removed and slab is not added back to yi

      }

      if (na.act == "na.fail")
         stop(mstyle$stop("Missing values in data."))

   }

   ### more than one study left?

   if (k <= 1)
      stop(mstyle$stop("Processing terminated since k <= 1."))

   ### check for non-positive sampling variances (and set negative values to 0)

   if (any(vi <= 0)) {
      allvipos <- FALSE
      if (!V0)
         warning(mstyle$warning("There are outcomes with non-positive sampling variances."), call.=FALSE)
      vi.neg <- vi < 0
      if (any(vi.neg)) {
         V[vi.neg,] <- 0 ### note: entire row set to 0 (so covariances are also 0)
         V[,vi.neg] <- 0 ### note: entire col set to 0 (so covariances are also 0)
         vi[vi.neg] <- 0
         warning(mstyle$warning("Negative sampling variances constrained to zero."), call.=FALSE)
      }
   } else {
      allvipos <- TRUE
   }

   ### check for V being positive definite (this should also cover non-positive variances)
   ### skipped: even if V is not positive definite, the marginal var-cov matrix can still be; so just check for pd during the optimization
   ### but at least issue a warning, since a fixed-effects model can then not be fitted and there is otherwise no indication why this is the case

   if (!V0 && !.chkpd(V))
      warning(mstyle$warning("'V' appears to be not positive definite."), call.=FALSE)

   ### check ratio of largest to smallest sampling variance
   ### note: need to exclude some special cases (0/0 = NaN, max(vi)/0 = Inf)
   ### TODO: use the condition number of V here instead?

   vimaxmin <- max(vi) / min(vi)

   if (is.finite(vimaxmin) && vimaxmin >= 1e7)
      warning(mstyle$warning("Ratio of largest to smallest sampling variance extremely large. May not be able to obtain stable results."), call.=FALSE)

   ### make sure that there is at least one column in X ([b])

   if (is.null(mods) && !intercept) {
      warning(mstyle$warning("Must either include an intercept and/or moderators in model.\nCoerced intercept into the model."), call.=FALSE)
      intercept <- TRUE
   }

   if (!is.null(mods) && ncol(mods) == 0L) {
      warning(mstyle$warning("Cannot fit model with an empty model matrix. Coerced intercept into the model."), call.=FALSE)
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
   if (inherits(tmp, "try-error")) {
      stop(mstyle$stop("Error in check for redundant predictors."))
   } else {
      coef.na <- is.na(coef(tmp))
      if (any(coef.na)) {
         warning(mstyle$warning("Redundant predictors dropped from the model."), call.=FALSE)
         X   <- X[,!coef.na,drop=FALSE]
         X.f <- X.f[,!coef.na,drop=FALSE]
      }
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

   ### number of columns in X (including the intercept if it is included)

   p <- NCOL(X)

   ### make sure variable names in X are unique

   colnames(X) <- colnames(X.f) <- .make.unique(colnames(X))

   ### check whether this is an intercept-only model

   if ((p == 1L) && .is.intercept(X)) {
      int.only <- TRUE
   } else {
      int.only <- FALSE
   }

   ### check if there are too many parameters for given k (currently skipped)

   ### set/check 'btt' argument

   btt <- .set.btt(btt, p, int.incl, colnames(X))
   m <- length(btt) ### number of betas to test (m = p if all betas are tested)

   ### check which beta elements are estimated versus fixed

   if (is.null(ddd$beta)) {
      beta.val <- rep(NA_real_, p)
      beta.est <- rep(TRUE, p)
   } else {
      beta.val <- ddd$beta
      if (length(beta.val) != p)
         stop(mstyle$stop(paste0("Length of 'beta' argument (", length(beta.val), ") does not match actual number of fixed effects (", p, ").")))
      beta.est <- is.na(beta.val)
   }

   #########################################################################
   #########################################################################
   #########################################################################

   ### stuff that need to be done after subsetting and filtering out NAs

   if (withS) {

      ### redo: turn each variable in mf.s into a factor (reevaluates the levels present, but order of existing levels is preserved)

      mf.s <- lapply(mf.s, factor)

      ### redo: get number of levels of each variable in mf.s (vector with one value per term)

      s.nlevels <- sapply(mf.s, nlevels)

      ### redo: get levels of each variable in mf.s

      s.levels <- lapply(mf.s, levels)

      ### for any single-level factor with unfixed sigma2, fix the sigma2 value to 0

      if (any(is.na(sigma2) & s.nlevels == 1)) {
         sigma2[is.na(sigma2) & s.nlevels == 1] <- 0
         warning(mstyle$warning("Single-level factor(s) found in 'random' argument. Corresponding 'sigma2' value(s) fixed to 0."), call.=FALSE)
      }

      ### create model matrix for each element in mf.s

      Z.S <- vector(mode="list", length=sigma2s)

      for (j in seq_len(sigma2s)) {
         if (s.nlevels[j] == 1) {
            Z.S[[j]] <- cbind(rep(1,k))
         } else {
            if (sparse) {
               Z.S[[j]] <- sparse.model.matrix(~ mf.s[[j]] - 1) ### cannot use this for factors with a single level
            } else {
               Z.S[[j]] <- model.matrix(~ mf.s[[j]] - 1) ### cannot use this for factors with a single level
            }
         }
         attr(Z.S[[j]], "assign")    <- NULL
         attr(Z.S[[j]], "contrasts") <- NULL
      }

   } else {

      Z.S <- NULL

   }

   #########################################################################

   ### stuff that need to be done after subsetting and filtering out NAs

   if (withR) {

      ### R may contain levels that are not in ids (that's fine; just filter them out)
      ### also, R may not be in the order that Z.S is in, so this fixes that up

      for (j in seq_along(R)) {
         if (!Rfix[j])
            next
         R[[j]] <- R[[j]][s.levels[[j]], s.levels[[j]]]
      }

      ### TODO: allow Rscale to be a vector so that different Rs can be scaled differently

      ### force each element of R to be a correlation matrix (and do some checks on that)

      if (Rscale=="cor" || Rscale=="cor0") {
         R[Rfix] <- lapply(R[Rfix], function(x) {
            if (any(diag(x) <= 0))
               stop(mstyle$stop("Cannot use Rscale=\"cor\" or Rscale=\"cor0\" with non-positive values on the diagonal of an 'R' matrix."), call.=FALSE)
            tmp <- cov2cor(x)
            if (any(abs(tmp) > 1))
               warning(mstyle$warning("Some values are larger than +-1 in an 'R' matrix after cov2cor() (see 'Rscale' argument)."), call.=FALSE)
            return(tmp)
         })
      }

      ### rescale R so that entries are 0 to (max(R) - min(R)) / (1 - min(R))
      ### this preserves the ultrametric properties of R and makes levels split at the root uncorrelated

      if (Rscale=="cor0")
         R[Rfix] <- lapply(R[Rfix], function(x) (x - min(x)) / (1 - min(x)))

      ### rescale R so that min(R) is zero (this is for the case that R is covariance matrix)

      if (Rscale=="cov0")
         R[Rfix] <- lapply(R[Rfix], function(x) (x - min(x)))

   }

   #########################################################################

   ### create (kxk) indicator/correlation matrices for random intercepts

   if (withS) {

      D.S <- vector(mode="list", length=sigma2s)

      for (j in seq_len(sigma2s)) {
         if (Rfix[j]) {
            if (sparse) {
               D.S[[j]] <- Z.S[[j]] %*% Matrix(R[[j]], sparse=TRUE) %*% t(Z.S[[j]])
            } else {
               D.S[[j]] <- Z.S[[j]] %*% R[[j]] %*% t(Z.S[[j]])
            }
            # D.S[[j]] <- as.matrix(nearPD(D.S[[j]])$mat)
            ### this avoids that the full matrix becomes non-positive definite but adding
            ### a tiny amount to the diagonal of D.S[[j]] is easier and works just as well
            ### TODO: consider doing something like this by default
         } else {
            D.S[[j]] <- tcrossprod(Z.S[[j]])
         }
      }

   } else {

      D.S <- NULL

   }

   #########################################################################

   ### stuff that need to be done after subsetting and filtering out NAs

   if (withG) {

      tmp <- .process.G.afterrmna(mf.g, g.nlevels, g.levels, g.values, struct[1], formulas[[1]], tau2, rho, Z.G1, Z.G2, isG=TRUE, sparse, ddd$dist[[1]], verbose)

      mf.g <- tmp$mf.g

      g.nlevels       <- tmp$g.nlevels
      g.nlevels.f     <- tmp$g.nlevels.f
      g.levels        <- tmp$g.levels
      g.levels.f      <- tmp$g.levels.f
      g.levels.r      <- tmp$g.levels.r
      g.levels.k      <- tmp$g.levels.k
      g.levels.comb.k <- tmp$g.levels.comb.k

      tau2   <- tmp$tau2
      rho    <- tmp$rho
      G      <- tmp$G
      g.Dmat <- tmp$Dmat
      g.rho.init <- tmp$rho.init

   } else {

      g.nlevels.f     <- NULL
      g.levels.f      <- NULL
      g.levels.r      <- NULL
      g.levels.k      <- NULL
      g.levels.comb.k <- NULL

      G <- NULL
      g.Dmat <- NULL
      g.rho.init <- NULL

   }

   #########################################################################

   ### stuff that need to be done after subsetting and filtering out NAs

   if (withH) {

      tmp <- .process.G.afterrmna(mf.h, h.nlevels, h.levels, h.values, struct[2], formulas[[2]], gamma2, phi, Z.H1, Z.H2, isG=FALSE, sparse, ddd$dist[[2]], verbose)

      mf.h <- tmp$mf.g

      h.nlevels       <- tmp$g.nlevels
      h.nlevels.f     <- tmp$g.nlevels.f
      h.levels        <- tmp$g.levels
      h.levels.f      <- tmp$g.levels.f
      h.levels.r      <- tmp$g.levels.r
      h.levels.k      <- tmp$g.levels.k
      h.levels.comb.k <- tmp$g.levels.comb.k

      gamma2 <- tmp$tau2
      phi    <- tmp$rho
      H      <- tmp$G
      h.Dmat <- tmp$Dmat
      h.phi.init <- tmp$rho.init

   } else {

      h.nlevels.f     <- NULL
      h.levels.f      <- NULL
      h.levels.r      <- NULL
      h.levels.k      <- NULL
      h.levels.comb.k <- NULL

      H <- NULL
      h.Dmat <- NULL
      h.phi.init <- NULL

   }

   #########################################################################

   #return(list(Z.S=Z.S, sigma2=sigma2, Z.G1=Z.G1, Z.G2=Z.G2, tau2=tau2, rho=rho, G=G, Z.H1=Z.H1, Z.H2=Z.H2, gamma2=gamma2, phi=phi, H=H, Rfix=Rfix, R=R))

   #########################################################################
   #########################################################################
   #########################################################################

   Y <- as.matrix(yi)

   ### initial values for variance components (need to do something better here in the future; see rma.mv2() and rma.bv() for some general ideas)

   if (verbose > 1)
      message(mstyle$message("Extracting/computing initial values ..."))

   QE <- NA_real_

   if (!V0) { # for V0 case, this always fails, so can skip it

      if (verbose > 1) {
         U <- try(chol(chol2inv(chol(V))), silent=FALSE)
      } else {
         U <- try(suppressWarnings(chol(chol2inv(chol(V)))), silent=TRUE)
      }

   }

   if (V0 || inherits(U, "try-error") || any(is.infinite(U))) {

      ### note: if V is sparse diagonal with 0 along the diagonal, U will not be a 'try-error'
      ### but have Inf along the diagonal, so need to check for this as well

      total <- sigma(lm(Y ~ X - 1))^2

      if (is.na(total)) # if X is a saturated model, then sigma() yields NaN
         total <- var(as.vector(Y)) / 100

   } else {

      sX <- U %*% X
      sY <- U %*% Y
      beta.FE <- try(solve(crossprod(sX), crossprod(sX, sY)), silent=TRUE)

      if (inherits(beta.FE, "try-error")) {

         total <- var(as.vector(Y))

      } else {

         ### TODO: consider a better way to set initial values
         #total <- max(.001*(sigma2s + tau2s + gamma2s), var(c(Y - X %*% res.FE$beta)) - 1/mean(1/diag(V)))
         #total <- max(.001*(sigma2s + tau2s + gamma2s), var(as.vector(sY - sX %*% beta)) - 1/mean(1/diag(V)))
         total  <- max(.001*(sigma2s + tau2s + gamma2s), var(as.vector(Y) - as.vector(X %*% beta.FE)) - 1/mean(1/diag(V)))

         #beta.FE <- ifelse(beta.est, beta.FE, beta.val)
         QE <- sum(as.vector(sY - sX %*% beta.FE)^2)

         ### QEp calculated further below

      }

   }

   sigma2.init <- rep(total / (sigma2s + tau2s + gamma2s), sigma2s)
   tau2.init   <- rep(total / (sigma2s + tau2s + gamma2s), tau2s)
   gamma2.init <- rep(total / (sigma2s + tau2s + gamma2s), gamma2s)

   if (is.null(g.rho.init)) {
      rho.init <- rep(.50, rhos)
   } else {
      rho.init <- g.rho.init
   }

   if (is.null(h.phi.init)) {
      phi.init <- rep(.50, phis)
   } else {
      phi.init <- h.phi.init
   }

   #########################################################################

   ### set default control parameters

   con <- list(verbose = FALSE,
               optimizer = "nlminb",      # optimizer to use ("optim","nlminb","uobyqa","newuoa","bobyqa","nloptr","nlm","hjk","nmk","mads","ucminf","lbfgsb3c","subplex","BBoptim","optimParallel","Rcgmin","Rvmmin")
               optmethod = "BFGS",        # argument 'method' for optim() ("Nelder-Mead" and "BFGS" are sensible options)
               parallel = list(),         # parallel argument for optimParallel() (note: 'cl' argument in parallel is not passed; this is directly specified via 'cl')
               cl = NULL,                 # arguments for optimParallel()
               ncpus = 1L,                # arguments for optimParallel()
               sigma2.init = sigma2.init, # initial value(s) for sigma2
               tau2.init = tau2.init,     # initial value(s) for tau2
               rho.init = rho.init,       # initial value(s) for rho
               gamma2.init = gamma2.init, # initial value(s) for gamma2
               phi.init = phi.init,       # initial value(s) for phi
               REMLf = TRUE,              # full REML likelihood (including all constants)
               evtol = 1e-07,             # lower bound for eigenvalues to determine if model matrix is positive definite
               cholesky = ifelse(is.element(struct, c("UN","UNR","GEN")), TRUE, FALSE), # by default, use Cholesky factorization for G and H matrix for "UN", "UNR", and "GEN" structures
               nearpd = FALSE,            # to force G and H matrix to become positive definite
               hessianCtrl = list(r=8),   # arguments passed on to 'method.args' of hessian()
               hesstol = .Machine$double.eps^0.5, # threshold for detecting fixed elements in Hessian
               hesspack = "numDeriv")     # package for computing the Hessian (numDeriv or pracma)

   ### replace defaults with any user-defined values

   con.pos <- pmatch(names(control), names(con))
   con[c(na.omit(con.pos))] <- control[!is.na(con.pos)]

   if (verbose)
      con$verbose <- verbose

   verbose <- con$verbose

   ### when restart=TRUE, restart at current estimates

   if (isTRUE(ddd$restart)) {

      ### check that the restart is done for a model that has the same type/number of var-cor components as the initial one

      okrestart <- TRUE

      if (withS && (is.null(.getfromenv("rma.mv", "sigma2")) || length(.getfromenv("rma.mv", "sigma2")) != sigma2s))
         okrestart <- FALSE
      if (withG && (is.null(.getfromenv("rma.mv", "tau2")) || length(.getfromenv("rma.mv", "tau2")) != tau2s))
         okrestart <- FALSE
      if (withG && (is.null(.getfromenv("rma.mv", "rho")) || length(.getfromenv("rma.mv", "rho")) != rhos))
         okrestart <- FALSE
      if (withH && (is.null(.getfromenv("rma.mv", "gamma2")) || length(.getfromenv("rma.mv", "gamma2")) != gamma2s))
         okrestart <- FALSE
      if (withH && (is.null(.getfromenv("rma.mv", "phi")) || length(.getfromenv("rma.mv", "phi")) != phis))
         okrestart <- FALSE

      if (!okrestart)
         stop(mstyle$stop(paste0("Restarting for a different model than the initial one.")))

      con$sigma2.init <- .getfromenv("rma.mv", "sigma2", default=con$sigma2.init)
      con$tau2.init   <- .getfromenv("rma.mv", "tau2",   default=con$tau2.init)
      con$rho.init    <- .getfromenv("rma.mv", "rho",    default=con$rho.init)
      con$gamma2.init <- .getfromenv("rma.mv", "gamma2", default=con$gamma2.init)
      con$phi.init    <- .getfromenv("rma.mv", "phi",    default=con$phi.init)

   }

   ### check for missings in initial values

   if (anyNA(con$sigma2.init))
      stop(mstyle$stop(paste0("No missing values allowed in 'sigma2.init'.")))
   if (anyNA(con$tau2.init))
      stop(mstyle$stop(paste0("No missing values allowed in 'tau2.init'.")))
   if (anyNA(con$rho.init))
      stop(mstyle$stop(paste0("No missing values allowed in 'rho.init'.")))
   if (anyNA(con$gamma2.init))
      stop(mstyle$stop(paste0("No missing values allowed in 'gamma2.init'.")))
   if (anyNA(con$phi.init))
      stop(mstyle$stop(paste0("No missing values allowed in 'phi.init'.")))

   ### expand initial values to correct length

   if (length(con$sigma2.init) == 1L)
      con$sigma2.init <- rep(con$sigma2.init, sigma2s)
   if (length(con$tau2.init) == 1L)
      con$tau2.init <- rep(con$tau2.init, tau2s)
   if (length(con$rho.init) == 1L)
      con$rho.init <- rep(con$rho.init, rhos)
   if (length(con$gamma2.init) == 1L)
      con$gamma2.init <- rep(con$gamma2.init, gamma2s)
   if (length(con$phi.init) == 1L)
      con$phi.init <- rep(con$phi.init, phis)

   ### checks on initial values set by the user (the initial values computed by the function are replaced by the user defined ones at this point)

   if (withS && any(con$sigma2.init <= 0))
      stop(mstyle$stop("Value(s) of 'sigma2.init' must be > 0"))

   if (withG && any(con$tau2.init <= 0))
      stop(mstyle$stop("Value(s) of 'tau2.init' must be > 0."))
   if (withG && struct[1]=="CAR" && (con$rho.init <= 0 | con$rho.init >= 1))
      stop(mstyle$stop("Value(s) of 'rho.init' must be in (0,1)."))
   if (withG && is.element(struct[1], c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH")) && any(con$rho.init <= 0))
      stop(mstyle$stop("Value(s) of 'rho.init' must be > 0."))
   if (withG && is.element(struct[1], c("PHYPL","PHYPD")) && con$rho.init < 0)
      stop(mstyle$stop("Value(s) of 'rho.init' must be in >= 0."))
   if (withG && !is.element(struct[1], c("CAR","SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYBM","PHYPL","PHYPD")) && any(con$rho.init <= -1 | con$rho.init >= 1))
      stop(mstyle$stop("Value(s) of 'rho.init' must be in (-1,1)."))

   if (withH && any(con$gamma2.init <= 0))
      stop(mstyle$stop("Value(s) of 'gamma2.init' must be > 0."))
   if (withH && struct[2]=="CAR" && (con$phi.init <= 0 | con$phi.init >= 1))
      stop(mstyle$stop("Value(s) of 'phi.init' must be in (0,1)."))
   if (withH && is.element(struct[2], c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH")) && any(con$phi.init <= 0))
      stop(mstyle$stop("Value(s) of 'phi.init' must be > 0."))
   if (withH && is.element(struct[2], c("PHYPL","PHYPD")) && con$phi.init < 0)
      stop(mstyle$stop("Value(s) of 'phi.init' must be in >= 0."))
   if (withH && !is.element(struct[2], c("CAR","SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYBM","PHYPL","PHYPD")) && any(con$phi.init <= -1 | con$phi.init >= 1))
      stop(mstyle$stop("Value(s) of 'phi.init' must be in (-1,1)."))

   ### in case user manually sets con$cholesky and specifies only a single value

   if (length(con$cholesky) == 1L)
      con$cholesky <- rep(con$cholesky, 2L)

   ### use of Cholesky factorization only applicable for models with "UN", "UNR", and "GEN" structure

   if (!withG) ### in case user sets cholesky=TRUE and struct="UN", struct="UNR", or struct="GEN" even though there is no 1st 'inner | outer' term
      con$cholesky[1] <- FALSE

   if (con$cholesky[1] && !is.element(struct[1], c("UN","UNR","GEN")))
      con$cholesky[1] <- FALSE

   if (!withH) ### in case user sets cholesky=TRUE and struct="UN", struct="UNR", or struct="GEN" even though there is no 2nd 'inner | outer' term
      con$cholesky[2] <- FALSE

   if (con$cholesky[2] && !is.element(struct[2], c("UN","UNR","GEN")))
      con$cholesky[2] <- FALSE

   ### copy initial values back (in case they were replaced by user-defined values); those values are
   ### then shown in the 'Variance Components in Model' table that is given when verbose=TRUE; cannot
   ### replace any fixed values, since that can lead to -Inf/+Inf below when transforming the initial
   ### values and then optim() throws an error and chol(G) and/or chol(H) is then likely to fail

   #sigma2.init <- ifelse(is.na(sigma2), con$sigma2.init, sigma2)
   #tau2.init   <- ifelse(is.na(tau2), con$tau2.init, tau2)
   #rho.init    <- ifelse(is.na(rho), con$rho.init, rho)
   sigma2.init <- con$sigma2.init
   tau2.init   <- con$tau2.init
   rho.init    <- con$rho.init
   gamma2.init <- con$gamma2.init
   phi.init    <- con$phi.init

   ### plug in fixed values for sigma2, tau2, rho, gamma2, and phi and transform initial values

   con$sigma2.init <- log(sigma2.init)

   if (con$cholesky[1]) {
      if (struct[1] == "UNR") {
         G <- .con.vcov.UNR(tau2.init, rho.init)
      } else {
         G <- .con.vcov.UN(tau2.init, rho.init)
      }
      G <- try(chol(G), silent=TRUE)
      if (inherits(G, "try-error") || anyNA(G))
         stop(mstyle$stop("Cannot take Choleski decomposition of initial 'G' matrix."))
      if (struct[1] == "UNR") {
         con$tau2.init <- log(tau2.init)
      } else {
         con$tau2.init <- diag(G)        ### note: con$tau2.init and con$rho.init are the 'choled' values of the initial G matrix, so con$rho.init really
         con$rho.init <- G[lower.tri(G)] ### contains the 'choled' covariances; and these values are also passed on the .ll.rma.mv as the initial values
      }
      if (length(con$rho.init) == 0L)
         con$rho.init <- 0
   } else {
      con$tau2.init <- log(tau2.init)
      if (struct[1] == "CAR")
         con$rho.init  <- qlogis(rho.init)
      if (is.element(struct[1], c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYBM","PHYPL","PHYPD")))
         con$rho.init  <- log(rho.init)
      if (!is.element(struct[1], c("CAR","SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYBM","PHYPL","PHYPD")))
         con$rho.init  <- atanh(rho.init)
   }

   if (con$cholesky[2]) {
      H <- .con.vcov.UN(gamma2.init, phi.init)
      H <- try(chol(H), silent=TRUE)
      if (inherits(H, "try-error") || anyNA(H))
         stop(mstyle$stop("Cannot take Choleski decomposition of initial 'H' matrix."))
      con$gamma2.init <- diag(H)      ### note: con$gamma2.init and con$phi.init are the 'choled' values of the initial H matrix, so con$phi.init really
      con$phi.init <- H[lower.tri(H)] ### contains the 'choled' covariances; and these values are also passed on the .ll.rma.mv as the initial values
      if (length(con$phi.init) == 0L)
         con$phi.init <- 0
   } else {
      con$gamma2.init <- log(gamma2.init)
      if (struct[2] == "CAR")
         con$phi.init  <- qlogis(phi.init)
      if (is.element(struct[2], c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYBM","PHYPL","PHYPD")))
         con$phi.init  <- log(phi.init)
      if (!is.element(struct[2], c("CAR","SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYBM","PHYPL","PHYPD")))
         con$phi.init  <- atanh(phi.init)
   }

   optimizer <- match.arg(con$optimizer, c("optim","nlminb","uobyqa","newuoa","bobyqa","nloptr","nlm","hjk","nmk","mads","ucminf","lbfgsb3c","subplex","BBoptim","optimParallel","Nelder-Mead","BFGS","CG","L-BFGS-B","SANN","Brent","Rcgmin","Rvmmin"))
   optmethod <- match.arg(con$optmethod, c("Nelder-Mead","BFGS","CG","L-BFGS-B","SANN","Brent"))
   if (optimizer %in% c("Nelder-Mead","BFGS","CG","L-BFGS-B","SANN","Brent")) {
      optmethod <- optimizer
      optimizer <- "optim"
   }
   nearpd     <- con$nearpd
   cholesky   <- con$cholesky
   parallel   <- con$parallel
   cl         <- con$cl
   ncpus      <- con$ncpus
   optcontrol <- control[is.na(con.pos)] ### get arguments that are control arguments for optimizer

   if (length(optcontrol) == 0L)
      optcontrol <- list()

   ### if control argument 'ncpus' is larger than 1, automatically switch to optimParallel optimizer

   if (ncpus > 1L)
      optimizer <- "optimParallel"

   reml <- ifelse(method == "REML", TRUE, FALSE)

   con$hesspack <- match.arg(con$hesspack, c("numDeriv","pracma"))

   if ((.isTRUE(cvvc) || cvvc %in% c("varcor","varcov","transf")) && !requireNamespace(con$hesspack, quietly=TRUE))
      stop(mstyle$stop(paste0("Please install the '", con$hesspack, "' package to compute the Hessian.")))

   ### check if length of sigma2.init, tau2.init, rho.init, gamma2.init, and phi.init matches number of variance components
   ### note: if a particular component is not included, reset (transformed) initial values (in case the user still specifies multiple initial values)

   if (withS) {
      if (length(con$sigma2.init) != sigma2s)
         stop(mstyle$stop(paste0("Length of 'sigma2.init' argument (", length(con$sigma2.init), ") does not match actual number of variance components (", sigma2s, ").")))
   } else {
      con$sigma2.init <- 0
   }

   if (withG) {
      if (length(con$tau2.init) != tau2s)
         stop(mstyle$stop(paste0("Length of 'tau2.init' argument (", length(con$tau2.init), ") does not match actual number of variance components (", tau2s, ").")))
   } else {
      con$tau2.init <- 0
   }

   if (withG) {
      if (length(con$rho.init) != rhos)
         stop(mstyle$stop(paste0("Length of 'rho.init' argument (", length(con$rho.init), ") does not match actual number of correlations (", rhos, ").")))
   } else {
      con$rho.init <- 0
   }

   if (withH) {
      if (length(con$gamma2.init) != gamma2s)
         stop(mstyle$stop(paste0("Length of 'gamma2.init' argument (", length(con$gamma2.init), ") does not match actual number of variance components (", gamma2s, ").")))
   } else {
      con$gamma2.init <- 0
   }

   if (withH) {
      if (length(con$phi.init) != phis)
         stop(mstyle$stop(paste0("Length of 'phi.init' argument (", length(con$phi.init), ") does not match actual number of correlations (", phis, ").")))
   } else {
      con$phi.init <- 0
   }

   #########################################################################

   ### check whether model matrix is of full rank

   if (!.chkpd(crossprod(X), tol=con$evtol))
      stop(mstyle$stop("Model matrix not of full rank. Cannot fit model."))

   ### which variance components are fixed? (TRUE/FALSE or NA if not applicable = not included)

   if (withS) {
      sigma2.fix <- !is.na(sigma2)
   } else {
      sigma2.fix <- NA
   }
   if (withG) {
      tau2.fix <- !is.na(tau2)
      rho.fix  <- !is.na(rho)
   } else {
      tau2.fix <- NA
      rho.fix  <- NA
   }
   if (withH) {
      gamma2.fix <- !is.na(gamma2)
      phi.fix    <- !is.na(phi)
   } else {
      gamma2.fix <- NA
      phi.fix    <- NA
   }

   vc.fix <- list(sigma2=sigma2.fix, tau2=tau2.fix, rho=rho.fix, gamma2=gamma2.fix, phi=phi.fix)

   ### show which variance components are included in the model, their initial value, and their specified value (NA if not specified)

   if (verbose) {
      cat("\n")
      cat(mstyle$verbose("Variance Components in Model:"))
      if (!withS && !withG && !withH) {
         cat(mstyle$verbose(" none"))
         cat("\n\n")
      } else {
         cat("\n\n")
         vcs <- rbind(c("sigma2" = if (withS) round(sigma2.init, digits[["var"]]) else NA_real_,
                        "tau2"   = if (withG) round(tau2.init, digits[["var"]]) else NA_real_,
                        "rho"    = if (withG) round(rho.init, digits[["var"]]) else NA_real_,
                        "gamma2" = if (withH) round(gamma2.init, digits[["var"]]) else NA_real_,
                        "phi"    = if (withH) round(phi.init, digits[["var"]]) else NA_real_),
                        round(c(   if (withS) sigma2 else NA_real_,
                                   if (withG) tau2 else NA_real_,
                                   if (withG) rho else NA_real_,
                                   if (withH) gamma2 else NA_real_,
                                   if (withH) phi else NA_real_), digits[["var"]]))
         vcs <- data.frame(vcs, stringsAsFactors=FALSE)
         rownames(vcs) <- c("initial", "specified")
         vcs <- rbind(included=ifelse(c(rep(withS, sigma2s), rep(withG, tau2s), rep(withG, rhos), rep(withH, gamma2s), rep(withH, phis)), "Yes", "No"), fixed=unlist(vc.fix), vcs)
         tmp <- capture.output(print(vcs, na.print="---"))
         .print.output(tmp, mstyle$verbose)
         cat("\n")
      }
   }

   level <- .level(level)

   #return(list(sigma2s, tau2s, rhos, gamma2s, phis))

   #########################################################################
   #########################################################################
   #########################################################################

   ###### model fitting, test statistics, and confidence intervals

   if (verbose > 1)
      message(mstyle$message("Model fitting ...\n"))

   ### estimate sigma2, tau2, rho, gamma2, and phi as needed

   tmp <- .chkopt(optimizer, optcontrol)
   optimizer  <- tmp$optimizer
   optcontrol <- tmp$optcontrol
   par.arg    <- tmp$par.arg
   ctrl.arg   <- tmp$ctrl.arg

   if (optimizer == "optimParallel::optimParallel") {

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

   if (!is.element(method, c("FE","EE","CE")) && !is.null(random)) {

      ### if at least one parameter needs to be estimated

      if (anyNA(c(sigma2, tau2, rho, gamma2, phi))) {

         optcall <- paste(optimizer, "(", par.arg, "=c(con$sigma2.init, con$tau2.init, con$rho.init, con$gamma2.init, con$phi.init),
            .ll.rma.mv, reml=reml, ", ifelse(optimizer=="optim", "method=optmethod, ", ""), "Y=Y, M=V, A=NULL, X=X, k=k, pX=p,
            D.S=D.S, Z.G1=Z.G1, Z.G2=Z.G2, Z.H1=Z.H1, Z.H2=Z.H2, g.Dmat=g.Dmat, h.Dmat=h.Dmat,
            sigma2.val=sigma2, tau2.val=tau2, rho.val=rho, gamma2.val=gamma2, phi.val=phi, beta.val=beta.val,
            sigma2s=sigma2s, tau2s=tau2s, rhos=rhos, gamma2s=gamma2s, phis=phis,
            withS=withS, withG=withG, withH=withH, struct=struct,
            g.levels.r=g.levels.r, h.levels.r=h.levels.r, g.values=g.values, h.values=h.values,
            sparse=sparse, cholesky=cholesky, nearpd=nearpd, vctransf=TRUE, vccov=FALSE, vccon=vccon,
            verbose=verbose, digits=digits, REMLf=con$REMLf, dofit=FALSE, hessian=FALSE", ctrl.arg, ")\n", sep="")

         #return(optcall)

         iteration <- 0
         try(assign("iteration", iteration, envir=.metafor), silent=TRUE)

         if (verbose) {
            opt.res <- try(eval(str2lang(optcall)), silent=!verbose)
         } else {
            opt.res <- try(suppressWarnings(eval(str2lang(optcall))), silent=!verbose)
         }

         #return(opt.res)

         ### convergence checks (if verbose print optimParallel log, if verbose > 2 print opt.res, and unify opt.res$par)

         opt.res$par <- .chkconv(optimizer=optimizer, opt.res=opt.res, optcontrol=optcontrol, fun="rma.mv", verbose=verbose)

         if (p == k) {

            ### when fitting a saturated model (with REML estimation), estimated values of variance components can remain stuck
            ### at their initial values; this ensures that the values are fixed to zero (unless values were fixed by the user)

            sigma2[is.na(sigma2)] <- 0
            tau2[is.na(tau2)]     <- 0
            rho[is.na(rho)]       <- 0
            gamma2[is.na(gamma2)] <- 0
            phi[is.na(phi)]       <- 0

         }

      } else {

         ### if all parameter are fixed to known values, can skip optimization

         opt.res <- list(par=c(sigma2, tau2, rho, gamma2, phi))

      }

      ### save these for Hessian computation

      sigma2.val <- sigma2
      tau2.val   <- tau2
      rho.val    <- rho
      gamma2.val <- gamma2
      phi.val    <- phi

   } else {

      opt.res <- list(par=c(0,0,0,0,0))

   }

   #########################################################################

   ### do the final model fit with estimated variance components

   fitcall <- .ll.rma.mv(opt.res$par, reml=reml, Y=Y, M=V, A=A, X=X, k=k, pX=p,
      D.S=D.S, Z.G1=Z.G1, Z.G2=Z.G2, Z.H1=Z.H1, Z.H2=Z.H2, g.Dmat=g.Dmat, h.Dmat=h.Dmat,
      sigma2.val=sigma2, tau2.val=tau2, rho.val=rho, gamma2.val=gamma2, phi.val=phi, beta.val=beta.val,
      sigma2s=sigma2s, tau2s=tau2s, rhos=rhos, gamma2s=gamma2s, phis=phis,
      withS=withS, withG=withG, withH=withH, struct=struct,
      g.levels.r=g.levels.r, h.levels.r=h.levels.r, g.values=g.values, h.values=h.values,
      sparse=sparse, cholesky=cholesky, nearpd=nearpd, vctransf=TRUE, vccov=FALSE, vccon=vccon,
      verbose=FALSE, digits=digits, REMLf=con$REMLf, dofit=TRUE)

   ### extract elements

   beta <- as.matrix(fitcall$beta)
   vb   <- as.matrix(fitcall$vb)

   vb[!beta.est,] <- NA_real_
   vb[,!beta.est] <- NA_real_

   if (withS)
      sigma2 <- fitcall$sigma2

   if (withG) {
      G <- as.matrix(fitcall$G)
      if (is.element(struct[1], c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYBM","PHYPL","PHYPD")))
         colnames(G) <- rownames(G) <- seq_len(nrow(G))
      if (is.element(struct[1], c("CS","HCS","UN","UNR","AR","HAR","CAR","ID","DIAG")))
         colnames(G) <- rownames(G) <- g.levels.f[[1]]
      if (is.element(struct[1], c("GEN","GDIAG")))
         colnames(G) <- rownames(G) <- g.names[-length(g.names)]
      tau2 <- fitcall$tau2
      rho  <- fitcall$rho
      cov1 <- G[lower.tri(G)]
   } else {
      cov1 <- 0
   }

   if (withH) {
      H <- as.matrix(fitcall$H)
      if (is.element(struct[2], c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYBM","PHYPL","PHYPD")))
         colnames(H) <- rownames(H) <- seq_len(nrow(H))
      if (is.element(struct[2], c("CS","HCS","UN","UNR","AR","HAR","CAR","ID","DIAG")))
         colnames(H) <- rownames(H) <- h.levels.f[[1]]
      if (is.element(struct[2], c("GEN","GDIAG")))
         colnames(H) <- rownames(H) <- h.names[-length(h.names)]
      gamma2 <- fitcall$gamma2
      phi    <- fitcall$phi
      cov2   <- H[lower.tri(H)]
   } else {
      cov2   <- 0
   }

   M <- fitcall$M

   ### remove row and column names of M
   ### (but only do this if M has row/column names)

   if (!is.null(dimnames(M)))
      M <- unname(M)

   #print(M[1:8,1:8])

   if (verbose > 1)
      message(mstyle$message(ifelse(verbose > 2, "", "\n"), "Conducting tests of the fixed effects ..."))

   ### ddf calculation

   if (is.element(test, c("knha","adhoc","t"))) {
      ddf <- .ddf.calc(dfs, X=X, k=k, p=p, mf.s=mf.s, mf.g=mf.g, mf.h=mf.h)
   } else {
      ddf <- rep(NA_integer_, p)
   }

   ### the Knapp & Hartung method (this is experimental)

   s2w <- 1

   if (is.element(test, c("knha","adhoc"))) {
      knha.rma.mv.warn <- .getfromenv("knha.rma.mv.warn", default=TRUE)
      if (knha.rma.mv.warn) {
         warning(mstyle$warning("Use of Knapp and Hartung method for 'rma.mv()' models is experimental.\nNote: This warning is only issued once per session (ignore at your peril)."), call.=FALSE)
         try(assign("knha.rma.mv.warn", FALSE, envir=.metafor), silent=TRUE)
      }
      RSS <- try(as.vector(t(Y - X %*% beta) %*% chol2inv(chol(M)) %*% (Y - X %*% beta)), silent=TRUE)
      if (inherits(RSS, "try-error"))
         stop(mstyle$stop(paste0("Failure when trying to compute adjustment factor for Knapp and Hartung method.")))
      if (RSS <= .Machine$double.eps) {
         s2w <- 0
      } else {
         s2w <- as.vector(RSS / (k - p))
      }
   }

   if (test == "adhoc")
      s2w[s2w < 1] <- 1

   vb <- s2w * vb

   ### QM calculation

   QM <- try(as.vector(t(beta)[btt] %*% chol2inv(chol(vb[btt,btt])) %*% beta[btt]), silent=TRUE)

   if (inherits(QM, "try-error"))
      QM <- NA_real_

   ### abbreviate some types of coefficient names

   if (.isTRUE(ddd$abbrev)) {
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

   se <- sqrt(diag(vb))
   names(se) <- NULL
   zval <- c(beta/se)

   if (is.element(test, c("knha","adhoc","t"))) {
      QM   <- QM / m
      QMdf <- c(m, min(ddf[btt]))
      QMp  <- if (QMdf[2] > 0) pf(QM, df1=QMdf[1], df2=QMdf[2], lower.tail=FALSE) else NA_real_
      pval <- sapply(seq_along(ddf), function(j) if (ddf[j] > 0) 2*pt(abs(zval[j]), df=ddf[j], lower.tail=FALSE) else NA_real_)
      crit <- sapply(seq_along(ddf), function(j) if (ddf[j] > 0) qt(level/2, df=ddf[j], lower.tail=FALSE) else NA_real_)
   } else {
      QMdf <- c(m, NA_integer_)
      QMp  <- pchisq(QM, df=QMdf[1], lower.tail=FALSE)
      pval <- 2*pnorm(abs(zval), lower.tail=FALSE)
      crit <- qnorm(level/2, lower.tail=FALSE)
   }

   ci.lb <- c(beta - crit * se)
   ci.ub <- c(beta + crit * se)

   #########################################################################

   ### heterogeneity test (Wald-type test of the extra coefficients in the saturated model)

   if (verbose > 1)
      message(mstyle$message("Conducting heterogeneity test ..."))

   QEdf <- k-p

   if (QEdf > 0L) {

      ### if V is not positive definite, FE model fit will fail; then QE is NA
      ### otherwise compute the RSS (which is equal to the Q/QE-test statistic)

      QEp <- pchisq(QE, df=QEdf, lower.tail=FALSE)

   } else {

      ### if the user fits a saturated model, then fit must be perfect and QE = 0 and QEp = 1

      QE  <- 0
      QEp <- 1

   }

   ### log-likelihood under a saturated model with ML estimation

   ll.QE <- -1/2 * (k) * log(2*base::pi) - 1/2 * determinant(V, logarithm=TRUE)$modulus

   #########################################################################

   ###### compute Hessian

   hessian <- NA_real_
   vvc <- NA_real_

   if (.isTRUE(cvvc) || cvvc %in% c("varcor","varcov","transf")) {

      if (verbose > 1)
         message(mstyle$message("Computing Hessian ...\n"))

      if (cvvc == "varcov" && (any(sigma2.fix, na.rm=TRUE) || any(tau2.fix, na.rm=TRUE) || any(rho.fix, na.rm=TRUE) || any(gamma2.fix, na.rm=TRUE) || any(phi.fix, na.rm=TRUE))) {
         warning(mstyle$warning("Cannot use cvvc='varcov' when one or more components are fixed. Setting cvvc='varcor'."), call.=FALSE)
         cvvc <- "varcor"
      }

      if (cvvc == "varcov" && any(!is.element(struct, c("UN","GEN")))) {
         warning(mstyle$warning("Cannot use cvvc='varcov' for the specified structure(s). Setting cvvc='varcor'."), call.=FALSE)
         cvvc <- "varcor"
      }

      if (cvvc == "varcov") {

         if (con$hesspack == "numDeriv")
            hessian <- try(numDeriv::hessian(func=.ll.rma.mv, x = c(sigma2, tau2, cov1, gamma2, cov2),
               method.args=con$hessianCtrl, reml=reml, Y=Y, M=V, A=NULL, X=X, k=k, pX=p,
               D.S=D.S, Z.G1=Z.G1, Z.G2=Z.G2, Z.H1=Z.H1, Z.H2=Z.H2, g.Dmat=g.Dmat, h.Dmat=h.Dmat,
               sigma2.val=sigma2.val, tau2.val=tau2.val, rho.val=rho.val, gamma2.val=gamma2.val, phi.val=phi.val, beta.val=beta.val,
               sigma2s=sigma2s, tau2s=tau2s, rhos=rhos, gamma2s=gamma2s, phis=phis,
               withS=withS, withG=withG, withH=withH, struct=struct,
               g.levels.r=g.levels.r, h.levels.r=h.levels.r, g.values=g.values, h.values=h.values,
               sparse=sparse, cholesky=c(FALSE,FALSE), nearpd=nearpd, vctransf=FALSE, vccov=TRUE, vccon=vccon,
               verbose=verbose, digits=digits, REMLf=con$REMLf, hessian=TRUE), silent=TRUE)
         if (con$hesspack == "pracma")
            hessian <- try(pracma::hessian(f=.ll.rma.mv, x0 = c(sigma2, tau2, cov1, gamma2, cov2),
               reml=reml, Y=Y, M=V, A=NULL, X=X, k=k, pX=p,
               D.S=D.S, Z.G1=Z.G1, Z.G2=Z.G2, Z.H1=Z.H1, Z.H2=Z.H2, g.Dmat=g.Dmat, h.Dmat=h.Dmat,
               sigma2.val=sigma2.val, tau2.val=tau2.val, rho.val=rho.val, gamma2.val=gamma2.val, phi.val=phi.val, beta.val=beta.val,
               sigma2s=sigma2s, tau2s=tau2s, rhos=rhos, gamma2s=gamma2s, phis=phis,
               withS=withS, withG=withG, withH=withH, struct=struct,
               g.levels.r=g.levels.r, h.levels.r=h.levels.r, g.values=g.values, h.values=h.values,
               sparse=sparse, cholesky=c(FALSE,FALSE), nearpd=nearpd, vctransf=FALSE, vccov=TRUE, vccon=vccon,
               verbose=verbose, digits=digits, REMLf=con$REMLf, hessian=TRUE), silent=TRUE)

         # note: vctransf=FALSE and cholesky=c(FALSE,FALSE), so we get the Hessian for the untransfored variances and covariances

      } else {

         if (con$hesspack == "numDeriv")
            hessian <- try(numDeriv::hessian(func=.ll.rma.mv, x = if (cvvc=="transf") opt.res$par else c(sigma2, tau2, rho, gamma2, phi),
               method.args=con$hessianCtrl, reml=reml, Y=Y, M=V, A=NULL, X=X, k=k, pX=p,
               D.S=D.S, Z.G1=Z.G1, Z.G2=Z.G2, Z.H1=Z.H1, Z.H2=Z.H2, g.Dmat=g.Dmat, h.Dmat=h.Dmat,
               sigma2.val=sigma2.val, tau2.val=tau2.val, rho.val=rho.val, gamma2.val=gamma2.val, phi.val=phi.val, beta.val=beta.val,
               sigma2s=sigma2s, tau2s=tau2s, rhos=rhos, gamma2s=gamma2s, phis=phis,
               withS=withS, withG=withG, withH=withH, struct=struct,
               g.levels.r=g.levels.r, h.levels.r=h.levels.r, g.values=g.values, h.values=h.values,
               sparse=sparse, cholesky=ifelse(c(cvvc=="transf",cvvc=="transf") & cholesky, TRUE, FALSE),
               nearpd=nearpd, vctransf=cvvc=="transf", vccov=FALSE, vccon=vccon,
               verbose=verbose, digits=digits, REMLf=con$REMLf, hessian=TRUE), silent=TRUE)
         if (con$hesspack == "pracma")
            hessian <- try(pracma::hessian(f=.ll.rma.mv, x0 = if (cvvc=="transf") opt.res$par else c(sigma2, tau2, rho, gamma2, phi),
               reml=reml, Y=Y, M=V, A=NULL, X=X, k=k, pX=p,
               D.S=D.S, Z.G1=Z.G1, Z.G2=Z.G2, Z.H1=Z.H1, Z.H2=Z.H2, g.Dmat=g.Dmat, h.Dmat=h.Dmat,
               sigma2.val=sigma2.val, tau2.val=tau2.val, rho.val=rho.val, gamma2.val=gamma2.val, phi.val=phi.val, beta.val=beta.val,
               sigma2s=sigma2s, tau2s=tau2s, rhos=rhos, gamma2s=gamma2s, phis=phis,
               withS=withS, withG=withG, withH=withH, struct=struct,
               g.levels.r=g.levels.r, h.levels.r=h.levels.r, g.values=g.values, h.values=h.values,
               sparse=sparse, cholesky=ifelse(c(cvvc=="transf",cvvc=="transf") & cholesky, TRUE, FALSE),
               nearpd=nearpd, vctransf=cvvc=="transf", vccov=FALSE, vccon=vccon,
               verbose=verbose, digits=digits, REMLf=con$REMLf, hessian=TRUE), silent=TRUE)

         # note: when cvvc=TRUE/"covcor", get the Hessian for the (untransfored) variances and correlations
         #       when cvvc="transf", get the Hessian for the transformed variances and correlations

      }

      if (inherits(hessian, "try-error")) {

         warning(mstyle$warning("Error when trying to compute the Hessian."), call.=FALSE)
         hessian <- NA_real_

      } else {

         ### row/column names

         colnames(hessian) <- seq_len(ncol(hessian)) ### need to do this, so the subsetting of colnames below works

         if (sigma2s == 1) {
            colnames(hessian)[1] <- "sigma^2"
         } else {
            colnames(hessian)[1:sigma2s] <- paste("sigma^2.", seq_len(sigma2s), sep="")
         }
         if (tau2s == 1) {
            colnames(hessian)[sigma2s+1] <- "tau^2"
         } else {
            colnames(hessian)[(sigma2s+1):(sigma2s+tau2s)] <- paste("tau^2.", seq_len(tau2s), sep="")
         }
         term <- ifelse(cvvc == "varcov", ifelse(withH, "cov1", "cov"), "rho")
         if (rhos == 1) {
            colnames(hessian)[sigma2s+tau2s+1] <- term
         } else {
            colnames(hessian)[(sigma2s+tau2s+1):(sigma2s+tau2s+rhos)] <- paste(term, ".", outer(seq_len(g.nlevels.f[1]), seq_len(g.nlevels.f[1]), paste, sep=".")[lower.tri(matrix(NA,nrow=g.nlevels.f,ncol=g.nlevels.f))], sep="")
            #colnames(hessian)[(sigma2s+tau2s+1):(sigma2s+tau2s+rhos)] <- paste(term, ".", seq_len(rhos), sep="")
         }
         if (gamma2s == 1) {
            colnames(hessian)[sigma2s+tau2s+rhos+1] <- "gamma^2"
         } else {
            colnames(hessian)[(sigma2s+tau2s+rhos+1):(sigma2s+tau2s+rhos+gamma2s)] <- paste("gamma^2.", seq_len(gamma2s), sep="")
         }
         term <- ifelse(cvvc == "varcov", "cov2", "phi")
         if (phis == 1) {
            colnames(hessian)[sigma2s+tau2s+rhos+gamma2s+1] <- term
         } else {
            colnames(hessian)[(sigma2s+tau2s+rhos+gamma2s+1):(sigma2s+tau2s+rhos+gamma2s+phis)] <- paste(term, ".", outer(seq_len(h.nlevels.f[1]), seq_len(h.nlevels.f[1]), paste, sep=".")[lower.tri(matrix(NA,nrow=h.nlevels.f,ncol=h.nlevels.f))], sep="")
            #colnames(hessian)[(sigma2s+tau2s+rhos+gamma2s+1):(sigma2s+tau2s+rhos+gamma2s+phis)] <- paste(term, ".", seq_len(phis), sep="")
         }

         rownames(hessian) <- colnames(hessian)

         ### select correct rows/columns from Hessian depending on components in the model
         ### FIXME: this isn't quite right, since "DIAG" and "ID" have a rho/phi element, but this is fixed at 0, so should also exclude this
         ###        in fact, all fixed elements should be filtered out (done below)

         #if (withS && withG && withH)
            #hessian <- hessian[1:nrow(hessian),1:ncol(hessian), drop=FALSE]
         if (withS && withG && !withH)
            hessian <- hessian[1:(nrow(hessian)-2),1:(ncol(hessian)-2), drop=FALSE]
         if (withS && !withG && !withH)
            hessian <- hessian[1:(nrow(hessian)-4),1:(ncol(hessian)-4), drop=FALSE]
         if (!withS && withG && withH)
            hessian <- hessian[2:nrow(hessian),2:ncol(hessian), drop=FALSE]
         if (!withS && withG && !withH)
            hessian <- hessian[2:(nrow(hessian)-2),2:(ncol(hessian)-2), drop=FALSE]
         if (!withS && !withG && !withH)
            hessian <- NA_real_

         ### reorder hessian when vccov into the order of the lower triangular elements of G/H

         if (cvvc == "varcov" && withG) {
            posG <- matrix(NA_real_, nrow=tau2s, ncol=tau2s)
            diag(posG) <- 1:tau2s
            posG[lower.tri(posG)] <- (tau2s+1):(tau2s*(tau2s+1)/2)
            posG <- posG[lower.tri(posG, diag=TRUE)]
            if (withS) {
               pos <- c(1:sigma2s, sigma2s+posG)
            } else {
               pos <- posG
            }
            if (withH) {
               posH <- matrix(NA_real_, nrow=gamma2s, ncol=gamma2s)
               diag(posH) <- 1:gamma2s
               posH[lower.tri(posH)] <- (gamma2s+1):(gamma2s*(gamma2s+1)/2)
               posH <- posH[lower.tri(posH, diag=TRUE)]
               pos <- c(pos, max(pos)+posH)
            }
            hessian <- hessian[pos,pos]
         }

         ### detect rows/columns that are essentially all equal to 0 (fixed elements) and filter them out

         hest <- !apply(hessian, 1, function(x) all(abs(x) <= con$hesstol))
         hessian <- hessian[hest, hest, drop=FALSE]

         ### try to invert Hessian

         vvc <- try(suppressWarnings(chol2inv(chol(hessian))), silent=TRUE)

         if (inherits(vvc, "try-error") || anyNA(vvc) || any(is.infinite(vvc))) {
            warning(mstyle$warning("Error when trying to invert Hessian."), call.=FALSE)
            vvc <- NA_real_
         } else {
            dimnames(vvc) <- dimnames(hessian)
         }

      }

      if (verbose > 1)
         cat("\n")

   }

   #########################################################################

   ###### fit statistics

   if (verbose > 1)
      message(mstyle$message("Computing fit statistics and log likelihood ..."))

   ### note: this only counts *estimated* variance components and correlations for the total number of parameters

   p <- sum(beta.est)

   if (is.null(vccon)) {

      parms <- p + ifelse(withS, sum(ifelse(sigma2.fix, 0, 1)), 0) +
                   ifelse(withG, sum(ifelse(tau2.fix,   0, 1)), 0) +
                   ifelse(withG, sum(ifelse(rho.fix,    0, 1)), 0) +
                   ifelse(withH, sum(ifelse(gamma2.fix, 0, 1)), 0) +
                   ifelse(withH, sum(ifelse(phi.fix,    0, 1)), 0)

   } else {

      parms <- p + ifelse(withS && !is.null(vccon$sigma2), length(unique(vccon$sigma2)) - sum(sigma2.fix), 0) +
                   ifelse(withG && !is.null(vccon$tau2),   length(unique(vccon$tau2))   - sum(tau2.fix),   0) +
                   ifelse(withG && !is.null(vccon$rho),    length(unique(vccon$rho))    - sum(rho.fix),    0) +
                   ifelse(withH && !is.null(vccon$gamma2), length(unique(vccon$gamma2)) - sum(gamma2.fix), 0) +
                   ifelse(withH && !is.null(vccon$phi),    length(unique(vccon$phi))    - sum(phi.fix),    0)

   }

   ll.ML   <- fitcall$llvals[1]
   ll.REML <- fitcall$llvals[2]

   if (allvipos) {
      dev.ML <- -2 * (ll.ML - ll.QE)
   } else {
      dev.ML <- -2 * ll.ML
   }
   AIC.ML    <- -2 * ll.ML   + 2*parms
   BIC.ML    <- -2 * ll.ML   +   parms * log(k)
   AICc.ML   <- -2 * ll.ML   + 2*parms * max(k, parms+2) / (max(k, parms+2) - parms - 1)
   dev.REML  <- -2 * (ll.REML - 0) # saturated model has ll = 0 when using the full REML likelihood
   AIC.REML  <- -2 * ll.REML + 2*parms
   BIC.REML  <- -2 * ll.REML +   parms * log(k-p)
   AICc.REML <- -2 * ll.REML + 2*parms * max(k-p, parms+2) / (max(k-p, parms+2) - parms - 1)

   fit.stats <- matrix(c(ll.ML, dev.ML, AIC.ML, BIC.ML, AICc.ML, ll.REML, dev.REML, AIC.REML, BIC.REML, AICc.REML), ncol=2, byrow=FALSE)
   dimnames(fit.stats) <- list(c("ll","dev","AIC","BIC","AICc"), c("ML","REML"))
   fit.stats <- data.frame(fit.stats)

   #########################################################################

   ### replace interaction() notation with : notation for nicer output (also for paste() and paste0())

   replfun <- function(x) {
      if (grepl("interaction(", x, fixed=TRUE) || grepl("paste(", x, fixed=TRUE) || grepl("paste0(", x, fixed=TRUE)) {
         #x <- gsub("^interaction\\(", "", x)
         #x <- gsub(", ", ":", x, fixed=TRUE)
         #x <- gsub("\\)$", "", x, fixed=FALSE)
         #x <- gsub("(.*)interaction\\(\\s*(.*)\\s*,\\s*(.*)\\s*\\)(.*)", "\\1(\\2:\\3)\\4", x)
         #x <- gsub("interaction\\((.*)\\s*,\\s*(.*)\\)", "(\\1:\\2)", x)
         x <- gsub("interaction\\((.*)\\)", "(\\1)", x)
         x <- gsub("paste[0]?\\((.*)\\)", "(\\1)", x)
         x <- gsub(",", ":", x, fixed=TRUE)
         x <- gsub(" ", "", x, fixed=TRUE)
         x <- gsub("^\\((.*)\\)$", "\\1", x) # if a name is "(...)", then can remove the ()
      }
      return(x)
   }

   s.names <- sapply(s.names, replfun)
   g.names <- sapply(g.names, replfun)
   h.names <- sapply(h.names, replfun)

   ############################################################################

   ###### prepare output

   if (verbose > 1)
      message(mstyle$message("Preparing output ..."))

   p.eff <- p
   k.eff <- k

   weighted <- TRUE

   if (is.null(ddd$outlist) || ddd$outlist == "nodata") {

      res <- list(b=beta, beta=beta, se=se, zval=zval, pval=pval, ci.lb=ci.lb, ci.ub=ci.ub, vb=vb,
                  sigma2=sigma2, tau2=tau2, rho=rho, gamma2=gamma2, phi=phi,
                  QE=QE, QEdf=QEdf, QEp=QEp, QM=QM, QMdf=QMdf, QMp=QMp,
                  k=k, k.f=k.f, k.eff=k.eff, k.all=k.all, p=p, p.eff=p.eff, parms=parms,
                  int.only=int.only, int.incl=int.incl, intercept=intercept, allvipos=allvipos, coef.na=coef.na,
                  yi=yi, vi=vi, V=V, W=A, X=X, yi.f=yi.f, vi.f=vi.f, V.f=V.f, X.f=X.f, W.f=W.f, ni=ni, ni.f=ni.f,
                  M=M, G=G, H=H, hessian=hessian, vvc=vvc, vccon=vccon,
                  ids=ids, not.na=not.na, subset=subset, slab=slab, slab.null=slab.null,
                  measure=measure, method=method, weighted=weighted,
                  test=test, dfs=dfs, ddf=ddf, s2w=s2w, btt=btt, m=m,
                  digits=digits, level=level, sparse=sparse, dist=ddd$dist, control=control, verbose=verbose,
                  fit.stats=fit.stats,
                  vc.fix=vc.fix,
                  withS=withS, withG=withG, withH=withH, withR=withR, formulas=formulas,
                  sigma2s=sigma2s, tau2s=tau2s, rhos=rhos, gamma2s=gamma2s, phis=phis,
                  s.names=s.names, g.names=g.names, h.names=h.names,
                  s.levels=s.levels, s.levels.f=s.levels.f,
                  s.nlevels=s.nlevels, s.nlevels.f=s.nlevels.f,
                  g.nlevels.f=g.nlevels.f, g.nlevels=g.nlevels,
                  h.nlevels.f=h.nlevels.f, h.nlevels=h.nlevels,
                  g.levels.f=g.levels.f, g.levels.k=g.levels.k, g.levels.comb.k=g.levels.comb.k,
                  h.levels.f=h.levels.f, h.levels.k=h.levels.k, h.levels.comb.k=h.levels.comb.k,
                  struct=struct, Rfix=Rfix, R=R, Rscale=Rscale,
                  mf.r=mf.r, mf.s=mf.s, mf.g=mf.g, mf.g.f=mf.g.f, mf.h=mf.h, mf.h.f=mf.h.f, Z.S=Z.S, Z.G1=Z.G1, Z.G2=Z.G2, Z.H1=Z.H1, Z.H2=Z.H2,
                  formula.yi=formula.yi, formula.mods=formula.mods, random=random,
                  version=packageVersion("metafor"), call=mf)

      if (is.null(ddd$outlist))
         res <- append(res, list(data=data), which(names(res) == "fit.stats"))

   } else {

      if (ddd$outlist == "minimal") {
         res <- list(b=beta, beta=beta, se=se, zval=zval, pval=pval, ci.lb=ci.lb, ci.ub=ci.ub, vb=vb,
                     sigma2=sigma2, tau2=tau2, rho=rho, gamma2=gamma2, phi=phi,
                     QE=QE, QEdf=QEdf, QEp=QEp, QM=QM, QMdf=QMdf, QMp=QMp,
                     k=k, k.eff=k.eff, k.all=k.all, p=p, p.eff=p.eff, parms=parms,
                     int.only=int.only,
                     measure=measure, method=method,
                     test=test, dfs=dfs, ddf=ddf, btt=btt, m=m,
                     digits=digits,
                     fit.stats=fit.stats,
                     vc.fix=vc.fix,
                     withS=withS, withG=withG, withH=withH, withR=withR,
                     s.names=s.names, g.names=g.names, h.names=h.names,
                     s.nlevels=s.nlevels,
                     g.nlevels.f=g.nlevels.f, g.nlevels=g.nlevels,
                     h.nlevels.f=h.nlevels.f, h.nlevels=h.nlevels,
                     g.levels.f=g.levels.f, g.levels.k=g.levels.k, g.levels.comb.k=g.levels.comb.k,
                     h.levels.f=h.levels.f, h.levels.k=h.levels.k, h.levels.comb.k=h.levels.comb.k,
                     struct=struct, Rfix=Rfix)
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

   class(res) <- c("rma.mv", "rma")
   return(res)

}
