selmodel.rma.uni <- function(x, type, alternative="greater", prec, delta, steps, verbose=FALSE, digits, control, ...) {

   # TODO: add a H0 argument? since p-value may not be based on H0: theta_i = 0
   # TODO: argument for which deltas to include in LRT (a delta may also not be constrained under H0, so it should not be included in the LRT then)

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="rma.uni", notav=c("rma.ls", "robust.rma"))

   alternative <- match.arg(alternative, c("two.sided", "greater", "less"))

   if (missing(digits)) {
      digits <- .get.digits(xdigits=x$digits, dmiss=TRUE)
   } else {
      digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
   }

   time.start <- proc.time()

   if (!x$allvipos)
      stop(mstyle$stop("Cannot fit selection model when one or more sampling variances are non-positive."))

   if (!x$weighted || !is.null(x$weights))
      stop(mstyle$stop("Cannot fit selection model for unweighted models or models with custom weights."))

   if (missing(type))
      stop(mstyle$stop("Must choose a specific selection model via the 'type' argument."))

   type.options <- c("beta", "halfnorm", "negexp", "logistic", "power", "negexppow", "halfnorm2", "negexp2", "logistic2", "power2", "stepfun")

   #type <- match.arg(type, type.options)
   type <- type.options[grep(type, type.options)[1]]
   if (is.na(type))
      stop(mstyle$stop("Unknown 'type' specified."))

   if (missing(control))
      control <- list()

   ### refit RE/ME models with ML estimation

   if (!is.element(x$method, c("FE","EE","CE","ML"))) {
      #stop(mstyle$stop("Argument 'x' must either be an equal/fixed-effects model or a model fitted with ML estimation."))
      #x <- try(update(x, method="ML"), silent=TRUE)
      #x <- suppressWarnings(update(x, method="ML"))
      #x <- try(suppressWarnings(rma.uni(x$yi, x$vi, weights=x$weights, mods=x$X, intercept=FALSE, method="ML", weighted=x$weighted, test=x$test, level=x$level, tau2=ifelse(x$tau2.fix, x$tau2, NA), control=x$control, skipr2=TRUE)), silent=TRUE)
      args <- list(yi=x$yi, vi=x$vi, weights=x$weights, mods=x$X, intercept=FALSE, method="ML", weighted=x$weighted, test=x$test, level=x$level, tau2=ifelse(x$tau2.fix, x$tau2, NA), control=x$control, skipr2=TRUE)
      x <- try(suppressWarnings(.do.call(rma.uni, args)), silent=TRUE)
      if (inherits(x, "try-error"))
         stop(mstyle$stop("Could not refit input model using method='ML'."))
   }

   ### get ... argument and check for extra/superfluous arguments

   ddd <- list(...)

   .chkdots(ddd, c("time", "tau2", "beta", "skiphes", "skiphet", "skipintcheck", "scaleprec", "defmap", "mapfun", "mapinvfun"))

   ### handle 'tau2' argument from ...

   if (is.null(ddd$tau2)) {

      if (is.element(x$method, c("FE","EE","CE"))) {
         tau2 <- 0
      } else {
         if (x$tau2.fix) {
            tau2 <- x$tau2
         } else {
            tau2 <- NA
         }
      }

   } else {

      tau2 <- ddd$tau2
      if (!is.na(tau2))
         x$tau2.fix <- TRUE

   }

   ### handle 'beta' argument from ...

   if (is.null(ddd$beta)) {
      beta <- rep(NA, x$p)
      betaspec <- FALSE
   } else {
      beta <- ddd$beta
      betaspec <- TRUE
   }

   yi <- x$yi
   vi <- x$vi
   X  <- x$X
   p  <- x$p
   k  <- x$k

   ### set precision measure

   if (!missing(prec) && !is.null(prec)) {

      precspec <- TRUE

      prec <- match.arg(prec, c("sei", "vi", "ninv", "sqrtninv"))

      ### check if sample size information is available if prec is "ninv" or "sqrtninv"

      if (is.element(prec, c("ninv", "sqrtninv"))) {
         if (is.null(x$ni) || anyNA(x$ni))
            stop(mstyle$stop("No sample size information stored in model object (or sample size information stored in model object contains NAs)."))
      }

      if (prec == "sei")
         preci <- sqrt(vi)
      if (prec == "vi")
         preci <- vi
      if (prec == "ninv")
         preci <- 1/x$ni
      if (prec == "sqrtninv")
         preci <- 1/sqrt(x$ni)

      if (is.null(ddd$scaleprec) || isTRUE(ddd$scaleprec))
         preci <- preci / max(preci)

   } else {

      prec <- NULL
      precspec <- FALSE
      preci <- rep(1, k)

   }

   precis <- c(min = min(preci), max = max(preci), mean = mean(preci), median = median(preci))

   ### compute p-values

   pvals <- .selmodel.pval(yi=c(yi), vi=vi, alternative=alternative)

   ### checks on steps argument

   if (missing(steps) || (length(steps) == 1L && is.na(steps))) {

      stepsspec <- FALSE
      steps <- NA

   } else {

      stepsspec <- TRUE

      if (anyNA(steps))
         stop(mstyle$stop("No missing values allowed in 'steps' argument."))

      if (any(steps < 0 | steps > 1))
         stop(mstyle$stop("Value(s) specified for 'steps' argument must be between 0 and 1."))

      steps <- unique(sort(steps))

      if (steps[length(steps)] != 1)
         steps <- c(steps, 1)

   }

   ############################################################################

   ### set default control parameters

   con <- list(verbose = FALSE,
               delta.init = NULL,     # initial value(s) for selection model parameter(s)
               beta.init = NULL,      # initial value(s) for fixed effect(s)
               tau2.init = NULL,      # initial value for tau^2
               delta.min = NULL,      # min possible value(s) for selection model parameter(s)
               delta.max = NULL,      # max possible value(s) for selection model parameter(s)
               tau2.max = Inf,        # max possible value for tau^2
               pval.min = NULL,       # minimum p-value to intergrate over (for selection models where this matters)
               optimizer = "optim",   # optimizer to use ("optim","nlminb","uobyqa","newuoa","bobyqa","nloptr","nlm","hjk","nmk","mads","ucminf","lbfgsb3c","subplex","BBoptim","optimParallel")
               optmethod = "BFGS",    # argument 'method' for optim() ("Nelder-Mead" and "BFGS" are sensible options)
               parallel = list(),     # parallel argument for optimParallel() (note: 'cl' argument in parallel is not passed; this is directly specified via 'cl')
               cl = NULL,             # arguments for optimParallel()
               ncpus = 1L,            # arguments for optimParallel()
               beta.fix = FALSE,      # fix beta in Hessian computation
               tau2.fix = FALSE,      # fix tau2 in Hessian computation
               delta.fix = FALSE,     # fix delta in Hessian computation
               htransf = FALSE,       # FALSE/TRUE: Hessian is computed for the untransformed/transformed delta and tau^2 estimates
               hessianCtrl=list(r=6), # arguments passed on to 'method.args' of hessian()
               hesspack = "numDeriv", # package for computing the Hessian (numDeriv or pracma)
               scaleX = !betaspec)    # whether non-dummy variables in the X matrix should be rescaled before model fitting

   ### replace defaults with any user-defined values

   con.pos <- pmatch(names(control), names(con))
   con[c(na.omit(con.pos))] <- control[!is.na(con.pos)]

   if (verbose)
      con$verbose <- verbose

   verbose <- con$verbose

   optimizer <- match.arg(con$optimizer, c("optim","nlminb","uobyqa","newuoa","bobyqa","nloptr","nlm","hjk","nmk","mads","ucminf","lbfgsb3c","subplex","BBoptim","optimParallel","Nelder-Mead","BFGS","CG","L-BFGS-B","SANN","Brent"))
   optmethod <- match.arg(con$optmethod, c("Nelder-Mead","BFGS","CG","L-BFGS-B","SANN","Brent"))
   if (optimizer %in% c("Nelder-Mead","BFGS","CG","L-BFGS-B","SANN","Brent")) {
      optmethod <- optimizer
      optimizer <- "optim"
   }
   parallel   <- con$parallel
   cl         <- con$cl
   ncpus      <- con$ncpus
   optcontrol <- control[is.na(con.pos)] # get arguments that are control arguments for optimizer
   optcontrol$intCtrl <- NULL            # but remove intCtrl from this list

   if (length(optcontrol) == 0L)
      optcontrol <- list()

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

   ### if control argument 'ncpus' is larger than 1, automatically switch to optimParallel optimizer

   if (ncpus > 1L)
      optimizer <- "optimParallel"

   ### rescale X matrix (only for models with moderators and models including an intercept term)

   if (!x$int.only && x$int.incl && con$scaleX) {
      Xsave <- X
      meanX <- colMeans(X[, 2:p, drop=FALSE])
      sdX   <- apply(X[, 2:p, drop=FALSE], 2, sd) ### consider using colSds() from matrixStats package
      is.d  <- apply(X, 2, .is.dummy) ### is each column a dummy variable (i.e., only 0s and 1s)?
      mX    <- rbind(c(intrcpt=1, -1*ifelse(is.d[-1], 0, meanX/sdX)), cbind(0, diag(ifelse(is.d[-1], 1, 1/sdX), nrow=length(is.d)-1, ncol=length(is.d)-1)))
      X[,!is.d] <- apply(X[, !is.d, drop=FALSE], 2, scale) ### rescale the non-dummy variables
   }

   ### initial value(s) for beta

   if (is.null(con$beta.init)) {
      beta.init <- c(x$beta)
   } else {
      if (length(con$beta.init) != p)
         stop(mstyle$stop(paste0("Length of 'beta.init' argument (", length(con$beta.init), ") does not match actual number of parameters (", p, ").")))
      beta.init <- con$beta.init
   }

   if (!x$int.only && x$int.incl && con$scaleX) {
      imX <- try(suppressWarnings(solve(mX)), silent=TRUE)
      if (inherits(imX, "try-error"))
         stop(mstyle$stop(paste0("Unable to rescale starting values for fixed effects.")))
      beta.init <- c(imX %*% cbind(beta.init))
   }

   ### check that tau2.max is larger than the tau^2 value

   if (x$tau2 >= con$tau2.max)
      stop(mstyle$stop("Value of 'tau2.max' must be > tau^2 value."))

   tau2.max <- con$tau2.max

   ### initial value for tau^2

   if (is.null(con$tau2.init)) {
      tau2.init <- log(x$tau2 + 0.001)
   } else {
      if (length(con$tau2.init) != 1L)
         stop(mstyle$stop(paste0("Argument 'tau2.init' should specify a single value.")))
      if (con$tau2.init <= 0)
         stop(mstyle$stop("Value of 'tau2.init' must be > 0."))
      if (con$tau2.init >= tau2.max)
         stop(mstyle$stop("Value of 'tau2.init' must be < 'tau2.max'."))
      tau2.init <- log(con$tau2.init)
   }

   ### set NLOPT_LN_BOBYQA as the default algorithm for nloptr optimizer
   ### and by default use a relative convergence criterion of 1e-8 on the function value

   if (optimizer == "nloptr" && !is.element("algorithm", names(optcontrol)))
      optcontrol$algorithm <- "NLOPT_LN_BOBYQA"

   if (optimizer == "nloptr" && !is.element("ftol_rel", names(optcontrol)))
      optcontrol$ftol_rel <- 1e-8

   ### for mads, set trace=FALSE and tol=1e-6 by default

   if (optimizer == "mads" && !is.element("trace", names(optcontrol)))
      optcontrol$trace <- FALSE

   if (optimizer == "mads" && !is.element("tol", names(optcontrol)))
      optcontrol$tol <- 1e-6

   ### for subplex, set reltol=1e-8 by default (the default in subplex() is .Machine$double.eps)

   if (optimizer == "subplex" && !is.element("reltol", names(optcontrol)))
      optcontrol$reltol <- 1e-8

   ### for BBoptim, set trace=FALSE by default

   if (optimizer == "BBoptim" && !is.element("trace", names(optcontrol)))
      optcontrol$trace <- FALSE

   ### check that the required packages are installed

   if (is.element(optimizer, c("uobyqa","newuoa","bobyqa"))) {
      if (!requireNamespace("minqa", quietly=TRUE))
         stop(mstyle$stop("Please install the 'minqa' package to use this optimizer."))
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

   if (!isTRUE(ddd$skiphes) && !requireNamespace(con$hesspack, quietly=TRUE))
      stop(mstyle$stop(paste0("Please install the '", con$hesspack, "' package to compute the Hessian.")))

   ############################################################################

   if (type == "beta") {

      if (stepsspec)
         warning(mstyle$warning("Argument 'steps' ignored (not applicable to this type of selection model)."), call.=FALSE)

      if (precspec)
         warning(mstyle$warning("Argument 'prec' ignored (not applicable to this type of selection model)."), call.=FALSE)

      deltas <- 2L
      delta.transf.fun <- c("exp", "exp")
      delta.transf.fun.inv <- c("log", "log")
      delta.lb <- c(0, 0)
      delta.ub <- c(Inf, Inf)
      delta.lb.excl <- c(TRUE, TRUE)
      delta.ub.excl <- c(FALSE, FALSE)
      delta.init <- c(1, 1)
      delta.min <- c(1e-05, 1e-05)
      delta.max <- c(100, 100)
      H0.delta <- c(1, 1)
      delta.LRT <- c(TRUE, TRUE)
      pval.min <- 1e-5
      wi.fun <- function(x, delta, yi, vi, preci, alternative, steps)
         x^(delta[1]-1) * (1-x)^(delta[2]-1)
      .selmodel.ll <- ".selmodel.ll.cont"

   }

   if (is.element(type, c("halfnorm", "negexp", "logistic", "power"))) {

      if (stepsspec) {
         if (length(steps) != 2L) # steps should be c(alpha,1)
            stop(mstyle$stop("Can only specify a single value for the 'steps' argument for this type of selection model."))
      } else {
         steps <- 0
      }

      deltas <- 1L
      delta.transf.fun <- "exp"
      delta.transf.fun.inv <- "log"
      delta.lb <- 0
      delta.ub <- Inf
      delta.lb.excl <- FALSE
      delta.ub.excl <- FALSE
      delta.init <- 1
      delta.min <- 0
      delta.max <- 100
      H0.delta <- 0
      delta.LRT <- TRUE
      if (type == "power") {
         pval.min <- 1e-5
      } else {
         pval.min <- 0
      }
      if (type == "halfnorm") {
         wi.fun <- function(x, delta, yi, vi, preci, alternative, steps)
            ifelse(x <= steps[1], 1, exp(-delta * preci * x^2) / exp(-delta * preci * steps[1]^2))
      }
      if (type == "negexp") {
         wi.fun <- function(x, delta, yi, vi, preci, alternative, steps)
            ifelse(x <= steps[1], 1, exp(-delta * preci * x) / exp(-delta * preci * steps[1]))
      }
      if (type == "logistic") {
         wi.fun <- function(x, delta, yi, vi, preci, alternative, steps)
            ifelse(x <= steps[1], 1, (2 * exp(-delta * preci * x) / (1 + exp(-delta * preci * x))) / (2 * exp(-delta * preci * steps[1]) / (1 + exp(-delta * preci * steps[1]))))
      }
      if (type == "power") {
         wi.fun <- function(x, delta, yi, vi, preci, alternative, steps)
            ifelse(x <= steps[1], 1, (1-x)^(preci*delta) / (1-steps[1])^(preci*delta))
      }
      .selmodel.ll <- ".selmodel.ll.cont"

   }

   if (type == "negexppow") {

      if (stepsspec) {
         if (length(steps) != 2L) # steps should be c(alpha,1)
            stop(mstyle$stop("Can only specify a single value for the 'steps' argument for this type of selection model."))
      } else {
         steps <- 0
      }

      deltas <- 2L
      delta.transf.fun <- c("exp", "exp")
      delta.transf.fun.inv <- c("log", "log")
      delta.lb <- c(0, 0)
      delta.ub <- c(Inf, Inf)
      delta.lb.excl <- c(FALSE, FALSE)
      delta.ub.excl <- c(FALSE, FALSE)
      delta.init <- c(1, 1)
      delta.min <- c(0, 0)
      delta.max <- c(100, 100)
      H0.delta <- c(0, 0)
      delta.LRT <- c(TRUE, TRUE)
      pval.min <- 0
      wi.fun <- function(x, delta, yi, vi, preci, alternative, steps)
         ifelse(x <= steps[1], 1, exp(-delta[1] * preci * x^(1/delta[2])) / exp(-delta[1] * preci * steps[1]^(1/delta[2])))
      .selmodel.ll <- ".selmodel.ll.cont"

   }

   if (is.element(type, c("halfnorm2", "negexp2", "logistic2", "power2"))) {

      if (stepsspec) {
         if (length(steps) != 2L) # steps should be c(alpha,1)
            stop(mstyle$stop("Can only specify a single value for the 'steps' argument for this type of selection model."))
      } else {
         steps <- 0
      }

      deltas <- 2L
      delta.transf.fun <- c("exp", "exp")
      delta.transf.fun.inv <- c("log", "log")
      delta.lb <- c(0,0)
      delta.ub <- c(Inf, Inf)
      delta.lb.excl <- c(FALSE, FALSE)
      delta.ub.excl <- c(FALSE, FALSE)
      delta.init <- c(1, 0.25)
      delta.min <- c(0, 0)
      delta.max <- c(100, 100)
      H0.delta <- c(0, 0)
      delta.LRT <- c(TRUE, TRUE)
      pval.min <- 0
      if (type == "halfnorm2") {
         wi.fun <- function(x, delta, yi, vi, preci, alternative, steps)
            ifelse(x <= steps[1], 1, (delta[1] + exp(-delta[2] * preci * x^2) / exp(-delta[2] * preci * steps[1]^2)) / (1 + delta[1]))
      }
      if (type == "negexp2") {
         wi.fun <- function(x, delta, yi, vi, preci, alternative, steps)
            ifelse(x <= steps[1], 1, (delta[1] + exp(-delta[2] * preci * x) / exp(-delta[2] * preci * steps[1])) / (1 + delta[1]))
      }
      if (type == "logistic2") {
         wi.fun <- function(x, delta, yi, vi, preci, alternative, steps)
            ifelse(x <= steps[1], 1, (delta[1] + (2 * exp(-delta[2] * preci * x) / (1 + exp(-delta[2] * preci * x))) / (2 * exp(-delta[2] * preci * steps[1]) / (1 + exp(-delta[2] * preci * steps[1])))) / (1 + delta[1]))
      }
      if (type == "power2") {
         wi.fun <- function(x, delta, yi, vi, preci, alternative, steps)
            ifelse(x <= steps[1], 1, (delta[1] + (1-x)^(preci*delta[2]) / (1-steps[1])^(preci*delta[2])) / (1 + delta[1]))
      }
      .selmodel.ll <- ".selmodel.ll.cont"

   }

   if (type == "stepfun") {

      if (!stepsspec)
         stop(mstyle$stop("Must specify the 'steps' argument for this type of selection model."))

      if (precspec)
         warning(mstyle$warning("Adding a precision measure to this selection model is undocumented and experimental."), call.=FALSE)

      deltas <- length(steps)
      delta.transf.fun <- rep("exp", deltas)
      delta.transf.fun.inv <- rep("log", deltas)
      delta.lb <- rep(0, deltas)
      delta.ub <- rep(Inf, deltas)
      delta.lb.excl <- rep(FALSE, deltas)
      delta.ub.excl <- rep(FALSE, deltas)
      delta.init <- rep(1, deltas)
      delta.min <- rep(0, deltas)
      delta.max <- rep(100, deltas)
      H0.delta <- rep(1, deltas)
      delta.LRT <- rep(TRUE, deltas) # note: delta[1] should actually not be included, but this gets constrained to 1 below anyway
      pval.min <- 0
      if (type == "stepfun") {
         wi.fun <- function(x, delta, yi, vi, preci, alternative, steps)
            delta[sapply(x, function(p) which(p <= steps)[1])] / preci
      }
      .selmodel.ll <- ".selmodel.ll.stepfun"

   }

   ############################################################################

   ### checks on delta, delta.init, delta.min, delta.max

   if (missing(delta)) {
      delta <- rep(NA, deltas)
   } else {
      if (length(delta) == 1L)
         delta <- rep(delta, deltas)
      if (length(delta) != deltas)
         stop(mstyle$stop(paste0("Argument 'delta' should be of length ", deltas, " for this type of selection model.")))
      for (j in seq_len(deltas)) {
         if (delta.lb.excl[j] && isTRUE(delta[j] <= delta.lb[j]))
            stop(mstyle$stop(paste0("Value of 'delta[", j, "]' must be > ", delta.lb[j], " for this type of selection model.")))
         if (!delta.lb.excl[j] && isTRUE(delta[j] < delta.lb[j]))
            stop(mstyle$stop(paste0("Value of 'delta[", j, "]' must be >= ", delta.lb[j], " for this type of selection model.")))
      }
      for (j in seq_len(deltas)) {
         if (delta.ub.excl[j] && isTRUE(delta[j] >= delta.ub[j]))
            stop(mstyle$stop(paste0("Value of 'delta[", j, "]' must be < ", delta.ub[j], " for this type of selection model.")))
         if (!delta.ub.excl[j] && isTRUE(delta[j] > delta.ub[j]))
            stop(mstyle$stop(paste0("Value of 'delta[", j, "]' must be <= ", delta.ub[j], " for this type of selection model.")))
      }
   }

   if (type == "stepfun" && is.na(delta[1]))
      delta[1] <- 1

   if (!is.null(con$delta.min))
      delta.min <- con$delta.min

   if (length(delta.min) == 1L)
      delta.min <- rep(delta.min, deltas)
   if (length(delta.min) != deltas)
      stop(mstyle$stop(paste0("Argument 'delta.min' should be of length ", deltas, " for this type of selection model.")))
   if (anyNA(delta.min))
      stop(mstyle$stop("No missing values allowed in 'delta.min'."))
   for (j in seq_len(deltas)) {
      if (delta.lb.excl[j] && delta.min[j] <= delta.lb[j])
         stop(mstyle$stop(paste0("Value of 'delta.min[", j, "]' must be > ", delta.lb[j], " for this type of selection model.")))
      if (!delta.lb.excl[j] && delta.min[j] < delta.lb[j])
         stop(mstyle$stop(paste0("Value of 'delta.min[", j, "]' must be >= ", delta.lb[j], " for this type of selection model.")))
   }
   for (j in seq_len(deltas)) {
      if (delta.ub.excl[j] && delta.min[j] >= delta.ub[j])
         stop(mstyle$stop(paste0("Value of 'delta.min[", j, "]' must be < ", delta.ub[j], " for this type of selection model.")))
      if (!delta.ub.excl[j] && delta.min[j] > delta.ub[j])
         stop(mstyle$stop(paste0("Value of 'delta.min[", j, "]' must be <= ", delta.ub[j], " for this type of selection model.")))
   }

   delta.min <- ifelse(!is.na(delta) & delta.min > delta, delta - .Machine$double.eps^0.2, delta.min)

   if (!is.null(con$delta.max))
      delta.max <- con$delta.max

   if (length(delta.max) == 1L)
      delta.max <- rep(delta.max, deltas)
   if (length(delta.max) != deltas)
      stop(mstyle$stop(paste0("Argument 'delta.max' should be of length ", deltas, " for this type of selection model.")))
   if (anyNA(delta.max))
      stop(mstyle$stop("No missing values allowed in 'delta.max'."))
   for (j in seq_len(deltas)) {
      if (delta.lb.excl[j] && delta.max[j] <= delta.lb[j])
         stop(mstyle$stop(paste0("Value of 'delta.max[", j, "]' must be > ", delta.lb[j], " for this type of selection model.")))
      if (!delta.lb.excl[j] && delta.max[j] < delta.lb[j])
         stop(mstyle$stop(paste0("Value of 'delta.max[", j, "]' must be >= ", delta.lb[j], " for this type of selection model.")))
   }
   for (j in seq_len(deltas)) {
      if (delta.ub.excl[j] && delta.max[j] >= delta.ub[j])
         stop(mstyle$stop(paste0("Value of 'delta.max[", j, "]' must be < ", delta.ub[j], " for this type of selection model.")))
      if (!delta.ub.excl[j] && delta.max[j] > delta.ub[j])
         stop(mstyle$stop(paste0("Value of 'delta.max[", j, "]' must be <= ", delta.ub[j], " for this type of selection model.")))
   }

   if (any(delta.max < delta.min))
      stop(mstyle$stop("Value(s) of 'delta.max' must be >= value(s) of 'delta.min'."))

   delta.max <- ifelse(!is.na(delta) & delta.max < delta, delta + .Machine$double.eps^0.2, delta.max)

   if (!is.null(con$delta.init))
      delta.init <- con$delta.init

   if (length(delta.init) == 1L)
      delta.init <- rep(delta.init, deltas)
   if (length(delta.init) != deltas)
      stop(mstyle$stop(paste0("Argument 'delta.init' should be of length ", deltas, " for this type of selection model.")))
   if (anyNA(delta.init))
      stop(mstyle$stop("No missing values allowed in 'delta.init'."))
   for (j in seq_len(deltas)) {
      if (delta.lb.excl[j] && delta.init[j] <= delta.lb[j])
         stop(mstyle$stop(paste0("Value of 'delta.init[", j, "]' must be > ", delta.lb[j], " for this type of selection model.")))
      if (!delta.lb.excl[j] && delta.init[j] < delta.lb[j])
         stop(mstyle$stop(paste0("Value of 'delta.init[", j, "]' must be >= ", delta.lb[j], " for this type of selection model.")))
   }
   for (j in seq_len(deltas)) {
      if (delta.ub.excl[j] && delta.init[j] >= delta.ub[j])
         stop(mstyle$stop(paste0("Value of 'delta.init[", j, "]' must be < ", delta.ub[j], " for this type of selection model.")))
      if (!delta.ub.excl[j] && delta.init[j] > delta.ub[j])
         stop(mstyle$stop(paste0("Value of 'delta.init[", j, "]' must be <= ", delta.ub[j], " for this type of selection model.")))
   }

   delta.init <- ifelse(!is.na(delta), delta, delta.init)

   if (.isTRUE(ddd$defmap) || any(is.infinite(delta.max))) {
      ddd$mapfun <- delta.transf.fun
      ddd$mapinvfun <- delta.transf.fun.inv
   }

   if (is.null(ddd$mapfun)) {
      mapfun <- rep(NA, deltas)
   } else {
      if (length(ddd$mapfun) == 1L) {
         mapfun <- rep(ddd$mapfun, deltas)
      } else {
         mapfun <- ddd$mapfun
      }
   }
   if (is.null(ddd$mapinvfun)) {
      mapinvfun <- rep(NA, deltas)
   } else {
      if (length(ddd$mapinvfun) == 1L) {
         mapinvfun <- rep(ddd$mapinvfun, deltas)
      } else {
         mapinvfun <- ddd$mapinvfun
      }
   }

   delta.init <- mapply(.mapinvfun, delta.init, delta.min, delta.max, mapinvfun)

   if (!is.null(con$pval.min))
      pval.min <- con$pval.min

   if (k < p + ifelse(is.element(x$method, c("FE","EE","CE")) || x$tau2.fix, 0, 1) + sum(is.na(delta)))
      stop(mstyle$stop(paste0("Number of studies (k=", k, ") is too small to fit the selection model.")))

   ############################################################################

   pvals[pvals < pval.min]     <- pval.min
   pvals[pvals > (1-pval.min)] <- 1-pval.min

   if (stepsspec) {

      pgrp     <- sapply(pvals, function(p) which(p <= steps)[1])
      psteps.l <- as.character(c(0,steps[-length(steps)]))
      psteps.r <- as.character(steps)
      len.l    <- nchar(psteps.l)
      pad.l    <- sapply(max(len.l) - len.l, function(x) paste0(rep(" ", x), collapse=""))
      psteps.l <- paste0(psteps.l, pad.l)
      psteps   <- paste0(psteps.l, " < p <= ", psteps.r)
      ptable   <- table(factor(pgrp, levels=seq_along(steps), labels=psteps))
      ptable   <- data.frame(k=as.vector(ptable), row.names=names(ptable))

      if (any(ptable[["k"]] == 0L)) {
         if (verbose >= 1)
            print(ptable)
         if (!isTRUE(ddd$skipintcheck) && type == "stepfun" && (any(ptable[["k"]] == 0L & is.na(delta))))
            stop(mstyle$stop(paste0("One or more intervals do not contain any observed p-values", if (!verbose) " (use 'verbose=TRUE' to see which)", ".")))
         if (!isTRUE(ddd$skipintcheck) && type != "stepfun")
            stop(mstyle$stop(paste0("One of the intervals does not contain any observed p-values", if (!verbose) " (use 'verbose=TRUE' to see which)", ".")))
      }

   } else {

      pgrp   <- NA
      ptable <- NA

   }

   ############################################################################

   if (optimizer == "optim") {
      par.arg <- "par"
      ctrl.arg <- ", control=optcontrol"
   }

   if (optimizer == "nlminb") {
      par.arg <- "start"
      ctrl.arg <- ", control=optcontrol"
   }

   if (is.element(optimizer, c("uobyqa","newuoa","bobyqa"))) {
      par.arg <- "par"
      optimizer <- paste0("minqa::", optimizer)
      ctrl.arg <- ", control=optcontrol"
   }

   if (optimizer == "nloptr") {
      par.arg <- "x0"
      optimizer <- paste0("nloptr::nloptr")
      ctrl.arg <- ", opts=optcontrol"
   }

   if (optimizer == "nlm") {
      par.arg <- "p" ### because of this, must use argument name pX for p (number of columns in X matrix)
      ctrl.arg <- paste(names(optcontrol), unlist(optcontrol), sep="=", collapse=", ")
      if (nchar(ctrl.arg) != 0L)
         ctrl.arg <- paste0(", ", ctrl.arg)
   }

   if (is.element(optimizer, c("hjk","nmk","mads"))) {
      par.arg <- "par"
      optimizer <- paste0("dfoptim::", optimizer)
      ctrl.arg <- ", control=optcontrol"
   }

   if (is.element(optimizer, c("ucminf","lbfgsb3c","subplex"))) {
      par.arg <- "par"
      optimizer <- paste0(optimizer, "::", optimizer)
      ctrl.arg <- ", control=optcontrol"
   }

   if (optimizer == "BBoptim") {
      par.arg <- "par"
      optimizer <- "BB::BBoptim"
      ctrl.arg <- ", quiet=TRUE, control=optcontrol"
   }

   if (optimizer == "optimParallel") {

      par.arg <- "par"
      optimizer <- paste0("optimParallel::optimParallel")
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

   # note: X.fit due to hessian(); pX due to nlm()

   optcall <- paste(optimizer, "(", par.arg, "=c(beta.init, tau2.init, delta.init), ",
      .selmodel.ll, ", ", ifelse(optimizer=="optim", "method=optmethod, ", ""),
      "yi=yi, vi=vi, X.fit=X, preci=preci, k=k, pX=p, pvals=pvals,
      deltas=deltas, delta.val=delta, delta.transf=TRUE, mapfun=mapfun, delta.min=delta.min, delta.max=delta.max,
      tau2.val=tau2, tau2.transf=TRUE, tau2.max=tau2.max, beta.val=beta,
      wi.fun=wi.fun, steps=steps, pgrp=pgrp,
      alternative=alternative, pval.min=pval.min, intCtrl=intCtrl, verbose=verbose, digits=digits", ctrl.arg, ")\n", sep="")

   if (verbose > 1)
      message(mstyle$message("\nModel fitting ...\n"))

   #return(optcall)

   if (verbose) {
      opt.res <- try(eval(str2lang(optcall)), silent=!verbose)
   } else {
      opt.res <- try(suppressWarnings(eval(str2lang(optcall))), silent=!verbose)
   }

   #return(opt.res)

   if (optimizer == "optimParallel::optimParallel" && verbose) {
      tmp <- capture.output(print(opt.res$loginfo))
      .print.output(tmp, mstyle$verbose)
   }

   if (inherits(opt.res, "try-error"))
      stop(mstyle$stop("Error during the optimization. Use verbose=TRUE and see help(selmodel) for more details on the optimization routines."))

   ### convergence checks

   if (is.element(optimizer, c("optim","nlminb","dfoptim::hjk","dfoptim::nmk","optimParallel::optimParallel","lbfgsb3c::lbfgsb3c","subplex::subplex","BB::BBoptim")) && opt.res$convergence != 0)
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

   ### estimates/values of tau2 and delta on the transformed scale

   tau2.transf  <- opt.res$par[p+1]
   delta.transf <- opt.res$par[(p+2):(p+1+deltas)]

   ### save for Hessian computation

   beta.val  <- beta
   tau2.val  <- tau2
   delta.val <- delta

   ### do the final model fit with estimated values

   fitcall <- paste(.selmodel.ll, "(par=opt.res$par,
      yi=yi, vi=vi, X.fit=X, preci=preci, k=k, pX=p, pvals=pvals,
      deltas=deltas, delta.val=delta, delta.transf=TRUE, mapfun=mapfun, delta.min=delta.min, delta.max=delta.max,
      tau2.val=tau2, tau2.transf=TRUE, tau2.max=tau2.max, beta.val=beta,
      wi.fun=wi.fun, steps=steps, pgrp=pgrp,
      alternative=alternative, pval.min=pval.min, intCtrl=intCtrl, verbose=FALSE, digits=digits, dofit=TRUE)\n", sep="")

   #return(fitcall)
   fitcall <- try(eval(str2lang(fitcall)), silent=!verbose)
   #return(fitcall)

   if (inherits(fitcall, "try-error"))
      stop(mstyle$stop("Error during the optimization. Use verbose=TRUE and see help(selmodel) for more details on the optimization routines."))

   ll    <- fitcall$ll
   beta  <- cbind(fitcall$beta)
   tau2  <- fitcall$tau2
   delta <- fitcall$delta

   if (any(delta <= delta.min + .Machine$double.eps^0.25) || any(delta >= delta.max - 100*.Machine$double.eps^0.25))
      warning(mstyle$warning("One or more 'delta' estimates are (almost) equal to their lower or upper bound.\nTreat results with caution (or consider adjusting 'delta.min' and/or 'delta.max')."), call.=FALSE)

   ############################################################################

   ### computing (inverse) Hessian

   H        <- NA
   vb       <- matrix(NA_real_, nrow=p, ncol=p)
   se.tau2  <- NA
   vb.delta <- matrix(NA_real_, nrow=deltas, ncol=deltas)

   if (con$beta.fix) {
      beta.hes <- c(beta)
   } else {
      beta.hes <- beta.val
   }

   if (con$tau2.fix) {
      tau2.hes <- tau2
   } else {
      tau2.hes <- tau2.val
   }

   if (con$delta.fix) {
      delta.hes <- delta
   } else {
      delta.hes <- delta.val
   }

   hest <- c(is.na(beta.hes), is.na(tau2.hes), is.na(delta.hes))

   if (any(hest) && !isTRUE(ddd$skiphes)) {

      if (verbose > 1)
         message(mstyle$message("\nComputing Hessian ..."))

      if (verbose > 3)
         cat("\n")

      if (con$htransf) {

         if (con$hesspack == "numDeriv")
            hescall <- paste("numDeriv::hessian(", .selmodel.ll, ", x=opt.res$par, method.args=con$hessianCtrl,
               yi=yi, vi=vi, X.fit=X, preci=preci, k=k, pX=p, pvals=pvals,
               deltas=deltas, delta.val=delta.hes, delta.transf=TRUE, mapfun=mapfun, delta.min=delta.min, delta.max=delta.max,
               tau2.val=tau2.hes, tau2.transf=TRUE, tau2.max=tau2.max, beta.val=beta.hes,
               wi.fun=wi.fun, steps=steps, pgrp=pgrp,
               alternative=alternative, pval.min=pval.min, intCtrl=intCtrl, verbose=ifelse(verbose > 3, verbose, 0), digits=digits)\n", sep="")
      if (con$hesspack == "pracma")
            hescall <- paste("pracma::hessian(", .selmodel.ll, ", x0=opt.res$par,
               yi=yi, vi=vi, X.fit=X, preci=preci, k=k, pX=p, pvals=pvals,
               deltas=deltas, delta.val=delta.hes, delta.transf=TRUE, mapfun=mapfun, delta.min=delta.min, delta.max=delta.max,
               tau2.val=tau2.hes, tau2.transf=TRUE, tau2.max=tau2.max, beta.val=beta.hes,
               wi.fun=wi.fun, steps=steps, pgrp=pgrp,
               alternative=alternative, pval.min=pval.min, intCtrl=intCtrl, verbose=ifelse(verbose > 3, verbose, 0), digits=digits)\n", sep="")

      } else {

         if (con$hesspack == "numDeriv")
            hescall <- paste("numDeriv::hessian(", .selmodel.ll, ", x=c(beta, tau2, delta), method.args=con$hessianCtrl,
               yi=yi, vi=vi, X.fit=X, preci=preci, k=k, pX=p, pvals=pvals,
               deltas=deltas, delta.val=delta.hes, delta.transf=FALSE, mapfun=mapfun, delta.min=delta.min, delta.max=delta.max,
               tau2.val=tau2.hes, tau2.transf=FALSE, tau2.max=tau2.max, beta.val=beta.hes,
               wi.fun=wi.fun, steps=steps, pgrp=pgrp,
               alternative=alternative, pval.min=pval.min, intCtrl=intCtrl, verbose=ifelse(verbose > 3, verbose, 0), digits=digits)\n", sep="")
      if (con$hesspack == "pracma")
            hescall <- paste("pracma::hessian(", .selmodel.ll, ", x0=c(beta, tau2, delta),
               yi=yi, vi=vi, X.fit=X, preci=preci, k=k, pX=p, pvals=pvals,
               deltas=deltas, delta.val=delta.hes, delta.transf=FALSE, mapfun=mapfun, delta.min=delta.min, delta.max=delta.max,
               tau2.val=tau2.hes, tau2.transf=FALSE, tau2.max=tau2.max, beta.val=beta.hes,
               wi.fun=wi.fun, steps=steps, pgrp=pgrp,
               alternative=alternative, pval.min=pval.min, intCtrl=intCtrl, verbose=ifelse(verbose > 3, verbose, 0), digits=digits)\n", sep="")

      }

      #return(hescall)

      H <- try(eval(str2lang(hescall)), silent=TRUE)

      #return(H)

      if (verbose > 3)
         cat("\n")

      if (inherits(H, "try-error")) {

         warning(mstyle$warning("Error when trying to compute the Hessian."), call.=FALSE)

      } else {

         if (deltas == 1L) {
            rownames(H) <- colnames(H) <- c(colnames(X), "tau2", "delta")
         } else {
            rownames(H) <- colnames(H) <- c(colnames(X), "tau2", paste0("delta.", seq_len(deltas)))
         }

         H.hest  <- H[hest, hest, drop=FALSE]
         iH.hest <- try(suppressWarnings(chol2inv(chol(H.hest))), silent=TRUE)

         if (inherits(iH.hest, "try-error") || anyNA(iH.hest) || any(is.infinite(iH.hest))) {
            warning(mstyle$warning("Error when trying to invert Hessian."), call.=FALSE)
         } else {
            iH <- matrix(0, nrow=length(hest), ncol=length(hest))
            iH[hest, hest] <- iH.hest
            if (anyNA(beta.hes))
               vb[is.na(beta.hes), is.na(beta.hes)] <- iH[c(is.na(beta.hes),FALSE,rep(FALSE,deltas)), c(is.na(beta.hes),FALSE,rep(FALSE,deltas)), drop=FALSE]
            if (is.na(tau2.hes))
               se.tau2 <- sqrt(iH[c(rep(FALSE,p),TRUE,rep(FALSE,deltas)), c(rep(FALSE,p),TRUE,rep(FALSE,deltas))])
            if (anyNA(delta.hes))
               vb.delta[is.na(delta.hes), is.na(delta.hes)] <- iH[c(rep(FALSE,p),FALSE,is.na(delta.hes)), c(rep(FALSE,p),FALSE,is.na(delta.hes)), drop=FALSE]
         }

      }

   }

   ############################################################################

   ### Wald-type tests of the fixed effects

   if (verbose > 1)
      message(mstyle$message("Conducting tests of the fixed effects ..."))

   ### scale back beta and vb

   if (!x$int.only && x$int.incl && con$scaleX) {
      beta <- mX %*% beta
      vb <- mX %*% vb %*% t(mX)
      X <- Xsave
   }

   ### QM calculation

   QM <- try(as.vector(t(beta)[x$btt] %*% chol2inv(chol(vb[x$btt,x$btt])) %*% beta[x$btt]), silent=TRUE)

   if (inherits(QM, "try-error"))
      QM <- NA

   QMp <- pchisq(QM, df=x$m, lower.tail=FALSE)

   rownames(beta) <- rownames(vb) <- colnames(vb) <- colnames(X)

   se <- sqrt(diag(vb))
   names(se) <- NULL

   ### inference for beta parameters

   zval  <- c(beta/se)
   pval  <- 2*pnorm(abs(zval), lower.tail=FALSE)
   crit  <- qnorm(x$level/2, lower.tail=FALSE)
   ci.lb <- c(beta - crit * se)
   ci.ub <- c(beta + crit * se)

   ### inference for delta parameters

   se.delta <- sqrt(diag(vb.delta))

   if (con$htransf) {

      zval.delta <- rep(NA_real_, deltas)
      pval.delta <- rep(NA_real_, deltas)
      ci.lb.delta <- c(delta.transf - crit * se.delta)
      ci.ub.delta <- c(delta.transf + crit * se.delta)
      ci.lb.delta <- mapply(.mapfun, ci.lb.delta, delta.min, delta.max, mapfun)
      ci.ub.delta <- mapply(.mapfun, ci.ub.delta, delta.min, delta.max, mapfun)
      vb.delta <- matrix(NA_real_, nrow=deltas, ncol=deltas)
      se.delta <- rep(NA_real_, deltas)

   } else {

      zval.delta  <- (delta - H0.delta) / se.delta
      pval.delta  <- 2*pnorm(abs(zval.delta), lower.tail=FALSE)
      ci.lb.delta <- c(delta - crit * se.delta)
      ci.ub.delta <- c(delta + crit * se.delta)

   }

   ci.lb.delta <- ifelse(ci.lb.delta < delta.lb, delta.lb, ci.lb.delta)
   ci.ub.delta <- ifelse(ci.ub.delta > delta.ub, delta.ub, ci.ub.delta)

   ### inference for tau^2 parameter

   if (con$htransf) {
      ci.lb.tau2 <- exp(tau2.transf - crit * se.tau2)
      ci.ub.tau2 <- exp(tau2.transf + crit * se.tau2)
      se.tau2 <- NA
   } else {
      ci.lb.tau2 <- tau2 - crit * se.tau2
      ci.ub.tau2 <- tau2 + crit * se.tau2
   }

   ci.lb.tau2[ci.lb.tau2 < 0] <- 0

   ############################################################################

   ### LRT for H0: tau^2 = 0 (only when NOT fitting a FE model)

   LRT.tau2  <- NA
   LRTp.tau2 <- NA

   if (!x$tau2.fix && !is.element(x$method, c("FE","EE","CE")) && !isTRUE(ddd$skiphet)) {

      if (verbose > 1)
         message(mstyle$message("Conducting heterogeneity test ..."))

      if (verbose > 4)
         cat("\n")

      optcall <- paste(optimizer, "(", par.arg, "=c(beta.init, tau2.init, delta.init), ",
         .selmodel.ll, ", ", ifelse(optimizer=="optim", "method=optmethod, ", ""),
         "yi=yi, vi=vi, X.fit=X, preci=preci, k=k, pX=p, pvals=pvals,
         deltas=deltas, delta.val=delta.val, delta.transf=TRUE, mapfun=mapfun, delta.min=delta.min, delta.max=delta.max,
         tau2.val=0, tau2.transf=FALSE, tau2.max=tau2.max, beta.val=beta.val,
         wi.fun=wi.fun, steps=steps, pgrp=pgrp,
         alternative=alternative, pval.min=pval.min, intCtrl=intCtrl, verbose=ifelse(verbose > 4, verbose, 0), digits=digits", ctrl.arg, ")\n", sep="")

      opt.res <- try(eval(str2lang(optcall)), silent=!verbose)

      if (verbose > 4)
         cat("\n")

      if (!inherits(opt.res, "try-error")) {

         fitcall <- paste(.selmodel.ll, "(par=opt.res$par,
            yi=yi, vi=vi, X.fit=X, preci=preci, k=k, pX=p, pvals=pvals,
            deltas=deltas, delta.val=delta.val, delta.transf=TRUE, mapfun=mapfun, delta.min=delta.min, delta.max=delta.max,
            tau2.val=0, tau2.transf=FALSE, tau2.max=tau2.max, beta.val=beta.val,
            wi.fun=wi.fun, steps=steps, pgrp=pgrp,
            alternative=alternative, pval.min=pval.min, intCtrl=intCtrl, verbose=FALSE, digits=digits, dofit=TRUE)\n", sep="")

         fitcall <- try(eval(str2lang(fitcall)), silent=!verbose)

         if (!inherits(fitcall, "try-error")) {
            ll0 <- fitcall$ll
            LRT.tau2  <- max(0, -2 * (ll0 - ll))
            LRTp.tau2 <- pchisq(LRT.tau2, df=1, lower.tail=FALSE)
         }

      }

   }

   ############################################################################

   ### LRT for selection model parameter(s)

   if (verbose > 1)
      message(mstyle$message("Conducting LRT for selection model parameter(s) ..."))

   ll0   <- c(logLik(x, REML=FALSE))
   LRT   <- max(0, -2 * (ll0 - ll))
   LRTdf <- sum(is.na(delta.val) & delta.LRT)
   LRTp  <- ifelse(LRTdf > 0, pchisq(LRT, df=LRTdf, lower.tail=FALSE), NA_real_)

   ############################################################################

   ### fit statistics

   if (verbose > 1)
      message(mstyle$message("Computing fit statistics and log likelihood ..."))

   ### note: tau2 and delta are not counted as parameters when they were fixed by the user
   parms <- p + ifelse(is.element(x$method, c("FE","EE","CE")) || x$tau2.fix, 0, 1) + sum(is.na(delta.val))

   ll.ML   <- ll
   dev.ML  <- -2 * ll.ML
   AIC.ML  <- -2 * ll.ML + 2*parms
   BIC.ML  <- -2 * ll.ML +   parms * log(k)
   AICc.ML <- -2 * ll.ML + 2*parms * max(k, parms+2) / (max(k, parms+2) - parms - 1)

   fit.stats <- matrix(c(ll.ML, dev.ML, AIC.ML, BIC.ML, AICc.ML, ll.REML=NA, dev.REML=NA, AIC.REML=NA, BIC.REML=NA, AICc.REML=NA), ncol=2, byrow=FALSE)
   dimnames(fit.stats) <- list(c("ll","dev","AIC","BIC","AICc"), c("ML","REML"))
   fit.stats <- data.frame(fit.stats)

   ############################################################################

   ### prepare output

   if (verbose > 1)
      message(mstyle$message("Preparing output ..."))

   res <- x

   res$beta    <- res$b <- beta
   res$se      <- se
   res$zval    <- zval
   res$pval    <- pval
   res$ci.lb   <- ci.lb
   res$ci.ub   <- ci.ub
   res$vb      <- vb

   res$betaspec <- betaspec

   res$tau2       <- res$tau2.f <- tau2
   res$se.tau2    <- se.tau2
   res$ci.lb.tau2 <- ci.lb.tau2
   res$ci.ub.tau2 <- ci.ub.tau2

   res$dfs <- res$ddf <- NA
   res$test <- "z"
   res$s2w  <- 1

   res$QE <- res$QEp <- NA
   res$I2 <- res$H2 <- res$vt <- NA
   res$R2 <- NULL

   res$QM  <- QM
   res$QMp <- QMp

   res$delta       <- delta
   res$vb.delta    <- vb.delta
   res$se.delta    <- se.delta
   res$zval.delta  <- zval.delta
   res$pval.delta  <- pval.delta
   res$ci.lb.delta <- ci.lb.delta
   res$ci.ub.delta <- ci.ub.delta
   res$deltas      <- deltas
   res$delta.fix   <- !is.na(delta.val)

   res$hessian <- H
   res$hest    <- hest
   res$ll      <- ll
   res$ll0     <- ll0
   res$LRT     <- LRT
   res$LRTdf   <- LRTdf
   res$LRTp    <- LRTp

   res$LRT.tau2  <- LRT.tau2
   res$LRTp.tau2 <- LRTp.tau2

   res$M         <- diag(vi + tau2, nrow=k, ncol=k)
   res$model     <- "rma.uni.selmodel"
   res$parms     <- parms
   res$fit.stats <- fit.stats
   res$pvals     <- pvals

   res$digits  <- digits
   res$verbose <- verbose

   res$type        <- type
   res$steps       <- steps
   res$stepsspec   <- stepsspec
   res$pgrp        <- pgrp
   res$ptable      <- ptable
   res$alternative <- alternative
   res$pval.min    <- pval.min
   res$prec        <- prec
   res$precspec    <- precspec
   res$precis      <- precis
   res$scaleprec   <- ddd$scaleprec

   res$wi.fun        <- wi.fun
   res$delta.lb      <- delta.lb
   res$delta.ub      <- delta.ub
   res$delta.lb.excl <- delta.lb.excl
   res$delta.ub.excl <- delta.ub.excl
   res$delta.min     <- delta.min
   res$delta.max     <- delta.max
   res$tau2.max      <- tau2.max

   res$call      <- match.call()
   res$control   <- control
   res$defmap    <- ddd$defmap
   res$mapfun    <- ddd$mapfun
   res$mapinvfun <- ddd$mapinvfun

   time.end <- proc.time()
   res$time <- unname(time.end - time.start)[3]

   if (.isTRUE(ddd$time))
      .print.time(res$time)

   if (verbose || .isTRUE(ddd$time))
      cat("\n")

   class(res) <- c("rma.uni.selmodel", class(res))
   return(res)

}
