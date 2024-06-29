selmodel.rma.uni <- function(x, type, alternative="greater", prec, subset, delta, steps, decreasing=FALSE, verbose=FALSE, digits, control, ...) {

   # TODO: add a H0 argument? since p-value may not be based on H0: theta_i = 0
   # TODO: argument for which deltas to include in LRT (a delta may also not be constrained under H0, so it should not be included in the LRT then)

   mstyle <- .get.mstyle()

   .chkclass(class(x), must="rma.uni", notav=c("rma.ls", "rma.gen", "robust.rma"))

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
      stop(mstyle$stop("Must choose a specific selection model via the 'type' argument (see 'help(selmodel)' for options)."))

   type.options <- c("beta", "halfnorm", "negexp", "logistic", "power", "negexppow", "halfnorm2", "negexp2", "logistic2", "power2", "stepfun", "stepcon", "trunc", "truncest")

   #type <- match.arg(type, type.options)
   type <- type.options[grep(type, type.options)[1]]
   if (is.na(type))
      stop(mstyle$stop("Unknown 'type' specified (see 'help(selmodel)' for options)."))

   if (is.element(type, c("trunc","truncest")) && alternative == "two.sided")
      stop(mstyle$stop("Cannot use alternative='two-sided' with this type of selection model."))

   decreasing <- isTRUE(decreasing)

   if (type != "stepfun" && decreasing) {
      warning(mstyle$warning("Argument 'decreasing' ignored (not applicable to this type of selection model)."), call.=FALSE)
      decreasing <- FALSE
   }

   if (missing(control))
      control <- list()

   ### refit RE/ME models with ML estimation

   if (!is.element(x$method, c("FE","EE","CE","ML"))) {
      #stop(mstyle$stop("Argument 'x' must either be an equal/fixed-effects model or a model fitted with ML estimation."))
      #x <- try(update(x, method="ML"), silent=TRUE)
      #x <- suppressWarnings(update(x, method="ML"))
      #x <- try(suppressWarnings(rma.uni(x$yi, x$vi, weights=x$weights, mods=x$X, intercept=FALSE, method="ML", weighted=x$weighted, test=x$test, level=x$level, tau2=ifelse(x$tau2.fix, x$tau2, NA), control=x$control, skipr2=TRUE)), silent=TRUE)
      args <- list(yi=x$yi, vi=x$vi, weights=x$weights, mods=x$X, intercept=FALSE, method="ML", weighted=x$weighted,
                   test=x$test, level=x$level, tau2=ifelse(x$tau2.fix, x$tau2, NA), control=x$control, skipr2=TRUE)
      x <- try(suppressWarnings(.do.call(rma.uni, args)), silent=TRUE)
      if (inherits(x, "try-error"))
         stop(mstyle$stop("Could not refit input model using method='ML'."))
   }

   ### get ... argument and check for extra/superfluous arguments

   ddd <- list(...)

   .chkdots(ddd, c("time", "tau2", "beta", "skiphes", "skiphet", "skipintcheck", "scaleprec", "defmap", "mapfun", "mapinvfun", "pval", "ptable", "retopt"))

   ### handle 'tau2' argument from ...

   if (is.null(ddd$tau2)) {

      if (is.element(x$method, c("FE","EE","CE"))) {
         tau2 <- 0
      } else {
         if (x$tau2.fix) {
            tau2 <- x$tau2
         } else {
            tau2 <- NA_real_
         }
      }

   } else {

      tau2 <- ddd$tau2

      if (!is.na(tau2))
         x$tau2.fix <- TRUE

   }

   ### handle 'beta' argument from ...

   if (is.null(ddd$beta)) {
      beta <- rep(NA_real_, x$p)
      betaspec <- FALSE # [a] sets con$scaleX=TRUE
   } else {
      beta <- ddd$beta
      betaspec <- TRUE # [a] sets con$scaleX=FALSE
   }

   yi <- c(x$yi)
   vi <- x$vi
   X  <- x$X
   p  <- x$p
   k  <- x$k

   ### set precision measure

   if (!missing(prec) && !is.null(prec)) {

      precspec <- TRUE # used to check if prec is set for certain models where this is not applicable or experimental [b]

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

      precspec <- FALSE
      prec <- NULL
      preci <- rep(1, k)

   }

   precis <- c(min = min(preci), max = max(preci), mean = mean(preci), median = median(preci))

   ### compute p-values

   if (is.null(ddd$pval)) {

      pvals <- .selmodel.pval(yi=yi, vi=vi, alternative=alternative)

   } else {

      # can pass p-values directly to the function via 'pval' argument from ... (this is highly experimental)

      pvals <- ddd$pval

      if (length(pvals) != x$k.all)
         stop(mstyle$stop(paste0("Length of the 'pval' argument (", length(pvals), ") does not correspond to the size of the original dataset (", x$k.all, ").")))

      pvals <- .getsubset(pvals, x$subset)
      pvals <- pvals[x$not.na]

      if (anyNA(pvals))
         stop(mstyle$stop(paste0("No missing values in 'pval' argument allowed.")))
      if (any(pvals <= 0) || any(pvals > 1))
         stop(mstyle$stop(paste0("One or more 'pval' values are <= 0 or > 1.")))

   }

   ### checks on steps argument

   if (missing(steps) || (length(steps) == 1L && is.na(steps))) {

      stepsspec <- FALSE
      steps <- NA_real_

   } else {

      stepsspec <- TRUE

      if (anyNA(steps))
         stop(mstyle$stop("No missing values allowed in 'steps' argument."))

      if (type != "trunc" && any(steps < 0 | steps > 1))
         stop(mstyle$stop("Value(s) specified for 'steps' argument must be between 0 and 1."))

      steps <- unique(sort(steps))

      if (type != "trunc") {
         if (steps[1] == 0)
            stop(mstyle$stop("Lowest 'steps' value must be > 0."))
         if (steps[length(steps)] != 1)
            steps <- c(steps, 1)
      }

   }

   if (type == "trunc" && !stepsspec) {
      stepsspec <- TRUE
      #if (alternative == "greater")
      #   steps <- min(yi)
      #if (alternative == "less")
      #   steps <- max(yi)
      steps <- 0
   }

   if (is.element(type, c("trunc","truncest")) && verbose > 2) {
      warning(mstyle$warning("Cannot use 'verbose > 2' for this type of selection model (setting verbose=2)."), call.=FALSE)
      verbose <- 2
   }


   if (missing(subset)) {
      subset <- rep(TRUE, k)
      subset.spec <- FALSE
   } else {
      mf <- match.call()
      subset <- .getx("subset", mf=mf, data=x$data)
      subset <- .chksubset(subset, x$k.all)
      subset <- .getsubset(subset, x$subset)
      subset <- subset[x$not.na]
      subset.spec <- TRUE
   }

   ############################################################################

   ### set defaults for control parameters

   con <- list(verbose = FALSE,
               delta.init = NULL,     # initial value(s) for selection model parameter(s)
               beta.init = NULL,      # initial value(s) for fixed effect(s)
               tau2.init = NULL,      # initial value for tau^2
               delta.min = NULL,      # min possible value(s) for selection model parameter(s)
               delta.max = NULL,      # max possible value(s) for selection model parameter(s)
               tau2.max = Inf,        # max possible value for tau^2
               tau2tol = min(vi/10, 1e-04), # threshold for treating tau^2 as effectively equal to 0 in the Hessian computation
               deltatol = 1e-04,      # threshold for treating deltas as effectively equal to 0 in the Hessian computation (only for stepfun)
               pval.min = NULL,       # minimum p-value to intergrate over (for selection models where this matters)
               optimizer = "optim",   # optimizer to use ("optim","nlminb","uobyqa","newuoa","bobyqa","nloptr","nlm","hjk","nmk","mads","ucminf","lbfgsb3c","subplex","BBoptim","optimParallel","solnp","alabama"/"constrOptim.nl","Rcgmin","Rvmmin")
               optmethod = "BFGS",    # argument 'method' for optim() ("Nelder-Mead" and "BFGS" are sensible options)
               parallel = list(),     # parallel argument for optimParallel() (note: 'cl' argument in parallel is not passed; this is directly specified via 'cl')
               cl = NULL,             # arguments for optimParallel()
               ncpus = 1L,            # arguments for optimParallel()
               beta.fix = FALSE,      # fix beta in Hessian computation
               tau2.fix = FALSE,      # fix tau2 in Hessian computation
               delta.fix = FALSE,     # fix delta in Hessian computation
               htransf = FALSE,       # when FALSE, Hessian is computed directly for the delta and tau^2 estimates (e.g., we get Var(tau^2)); when TRUE, Hessian is computed for the transformed estimates (e.g., we get Var(log(tau2)))
               hessianCtrl=list(r=6), # arguments passed on to 'method.args' of hessian()
               hesspack = "numDeriv", # package for computing the Hessian (numDeriv or pracma)
               scaleX = !betaspec)    # whether non-dummy variables in the X matrix should be rescaled before model fitting [a]

   ### replace defaults with any user-defined values

   con.pos <- pmatch(names(control), names(con))
   con[c(na.omit(con.pos))] <- control[!is.na(con.pos)]

   if (verbose)
      con$verbose <- verbose

   verbose <- con$verbose

   optimizer <- match.arg(con$optimizer, c("optim","nlminb","uobyqa","newuoa","bobyqa","nloptr","nlm","hjk","nmk","mads","ucminf","lbfgsb3c","subplex","BBoptim","optimParallel","solnp","alabama","constrOptim.nl","Nelder-Mead","BFGS","CG","L-BFGS-B","SANN","Brent","Rcgmin","Rvmmin"))
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

   ### can use optimizer="alabama" as a shortcut for optimizer="constrOptim.nl"

   if (optimizer == "alabama")
      optimizer <- "constrOptim.nl"

   ### when type="stepcon", automatically set solnp as the default optimizer

   if (type == "stepcon") {
      if (optimizer == "optim" && optmethod=="BFGS") { # this is the default
         optimizer <- "solnp"
      } else {
         if (!is.element(optimizer, c("solnp","nloptr","constrOptim.nl"))) {
            optimizer <- "solnp"
            warning(mstyle$warning(paste0("Can only use optimizers 'solnp', 'nloptr', or 'constrOptim.nl' when type='stepcon' (resetting to '", optimizer, "').")), call.=FALSE)
         }
      }
   }

   if (type != "stepcon" && optimizer == "constrOptim.nl") { # but can use solnp and nloptr
      optimizer <- "optim"
      warning(mstyle$warning(paste0("Cannot use 'constrOptim.nl' optimizer to fit this model (resetting to '", optimizer, "').")), call.=FALSE)
   }

   ### rescale X matrix (only for models with moderators and models including an intercept term)

   if (!x$int.only && x$int.incl && con$scaleX) {
      Xsave <- X
      meanX <- colMeans(X[, 2:p, drop=FALSE])
      sdX   <- apply(X[, 2:p, drop=FALSE], 2, sd) # consider using colSds() from matrixStats package
      is.d  <- apply(X, 2, .is.dummy) # is each column a dummy variable (i.e., only 0s and 1s)?
      mX    <- rbind(c(intrcpt=1, -1*ifelse(is.d[-1], 0, meanX/sdX)), cbind(0, diag(ifelse(is.d[-1], 1, 1/sdX), nrow=length(is.d)-1, ncol=length(is.d)-1)))
      X[,!is.d] <- apply(X[, !is.d, drop=FALSE], 2, scale) # rescale the non-dummy variables
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
         stop(mstyle$stop("Unable to rescale starting values for the fixed effects."))
      beta.init <- c(imX %*% cbind(beta.init))
   }

   ### check that tau2.max (Inf by default) is larger than the tau^2 value

   tau2.max <- con$tau2.max

   if (x$tau2 >= con$tau2.max)
      stop(mstyle$stop("Value of 'tau2.max' must be > tau^2 value."))

   ### initial value for tau^2

   if (is.null(con$tau2.init)) {
      tau2.init <- log(x$tau2 + 1e-3)
   } else {
      if (length(con$tau2.init) != 1L)
         stop(mstyle$stop("Argument 'tau2.init' should specify a single value."))
      if (con$tau2.init <= 0)
         stop(mstyle$stop("Value of 'tau2.init' must be > 0."))
      if (con$tau2.init >= tau2.max)
         stop(mstyle$stop("Value of 'tau2.init' must be < 'tau2.max'."))
      tau2.init <- log(con$tau2.init)
   }

   con$hesspack <- match.arg(con$hesspack, c("numDeriv","pracma"))

   if (!isTRUE(ddd$skiphes) && !requireNamespace(con$hesspack, quietly=TRUE))
      stop(mstyle$stop(paste0("Please install the '", con$hesspack, "' package to compute the Hessian.")))

   ############################################################################

   ### definition of the various selection model types

   # delta.lb / delta.ub: parameter space of the delta value(s)
   # delta.lb.excl / delta.ub.excl: whether delta must be >/< or can be >=/<=
   # delta.min / delta.max: limits imposed on delta for numerical reasons

   delta.min.check <- TRUE
   delta.max.check <- TRUE

   if (type == "beta") {

      if (stepsspec)
         warning(mstyle$warning("Argument 'steps' ignored (not applicable to this type of selection model)."), call.=FALSE)

      stepsspec <- FALSE
      steps <- NA_real_

      if (precspec) # [b]
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
         stop(mstyle$stop("Must specify 'steps' argument for this type of selection model."))

      if (precspec) { # [b]
         if (decreasing) {
            warning(mstyle$warning("Argument 'prec' ignored (not applicable to this type of selection model)."), call.=FALSE)
            preci <- rep(1, k)
         } else {
            warning(mstyle$warning("Adding a precision measure to this selection model is undocumented and experimental."), call.=FALSE)
         }
      }

      deltas <- length(steps)
      if (decreasing) {
         delta.transf.fun <- rep("I", deltas)
         delta.transf.fun.inv <- rep("I", deltas)
         ddd$defmap <- TRUE # actual mapping is defined directly in .selmodel.ll.stepfun() for this special case
         if (con$htransf)
            stop(mstyle$stop("Cannot use 'htransf=TRUE' for this type of selection model."))
         #delta.lb <- rep(0, deltas)
         #delta.ub <- rep(1, deltas)
         delta.lb <- c(0, rep(-Inf, deltas-1))
         delta.ub <- c(1, rep( Inf, deltas-1))
         delta.lb.excl <- rep(FALSE, deltas)
         delta.ub.excl <- rep(FALSE, deltas)
         #delta.init <- rep(1, deltas)
         delta.init <- c(1, rep(-2, deltas-1))
         delta.min <- rep(0, deltas)
         delta.max <- rep(1, deltas)
         delta.max.check <- FALSE
      } else {
         delta.transf.fun <- rep("exp", deltas)
         delta.transf.fun.inv <- rep("log", deltas)
         delta.lb <- rep(0, deltas)
         delta.ub <- rep(Inf, deltas)
         delta.lb.excl <- rep(FALSE, deltas)
         delta.ub.excl <- rep(FALSE, deltas)
         delta.init <- seq(1, 0.8, length.out=deltas)
         delta.min <- rep(0, deltas)
         delta.max <- rep(100, deltas)
      }
      H0.delta <- rep(1, deltas)
      delta.LRT <- rep(TRUE, deltas) # note: delta[1] should actually not be included in the LRT, but gets constrained to 1 below anyway
      pval.min <- 0
      wi.fun <- function(x, delta, yi, vi, preci, alternative, steps)
         delta[sapply(x, function(p) which(p <= steps)[1])] / preci
      .selmodel.ll <- ".selmodel.ll.stepfun"

   }

   if (type == "stepcon") {

      if (!stepsspec)
         stop(mstyle$stop("Must specify 'steps' argument for this type of selection model."))

      if (precspec) { # [b]
         warning(mstyle$warning("Argument 'prec' ignored (not applicable to this type of selection model)."), call.=FALSE)
         preci <- rep(1, k)
      }

      deltas <- length(steps)
      delta.transf.fun <- rep("plogis", deltas)
      delta.transf.fun.inv <- rep("qlogis", deltas)
      delta.lb <- rep(0, deltas)
      delta.ub <- rep(1, deltas)
      delta.lb.excl <- rep(FALSE, deltas)
      delta.ub.excl <- rep(FALSE, deltas)
      delta.init <- seq(1, 0.5, length.out=deltas)
      delta.min <- rep(0, deltas)
      delta.max <- rep(1, deltas)
      delta.max.check <- FALSE
      H0.delta <- rep(1, deltas)
      delta.LRT <- rep(TRUE, deltas) # note: delta[1] should actually not be included in the LRT, but gets constrained to 1 below anyway
      pval.min <- 0
      wi.fun <- function(x, delta, yi, vi, preci, alternative, steps)
         delta[sapply(x, function(p) which(p <= steps)[1])] / preci
      .selmodel.ll <- ".selmodel.ll.stepfun"

   }

   if (type == "trunc") {

      if (length(steps) != 1L) # steps should be a single value
         stop(mstyle$stop("Can only specify a single value for the 'steps' argument for this type of selection model."))

      if (precspec) # [b]
         warning(mstyle$warning("Argument 'prec' ignored (not applicable to this type of selection model)."), call.=FALSE)

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
      H0.delta <- 1
      delta.LRT <- TRUE
      pval.min <- 0
      wi.fun <- function(x, delta, yi, vi, preci, alternative, steps) {
         if (alternative == "less") {
            yival <- qnorm(x, sd=sqrt(vi), lower.tail=TRUE)
            ifelse(yival < steps[1], 1, delta)
         } else {
            yival <- qnorm(x, sd=sqrt(vi), lower.tail=FALSE)
            ifelse(yival > steps[1], 1, delta)
         }
      }
      #.selmodel.ll <- ".selmodel.ll.cont"
      .selmodel.ll <- ".selmodel.ll.trunc"

   }

   if (type == "truncest") {

      if (stepsspec)
         warning(mstyle$warning("Argument 'steps' ignored (not applicable to this type of selection model)."), call.=FALSE)

      stepsspec <- FALSE
      steps <- NA_real_

      if (precspec) # [b]
         warning(mstyle$warning("Argument 'prec' ignored (not applicable to this type of selection model)."), call.=FALSE)

      deltas <- 2L
      delta.transf.fun <- c("exp", "I")
      delta.transf.fun.inv <- c("log", "I")
      delta.lb <- c(0, -Inf)
      delta.ub <- c(Inf, Inf)
      delta.lb.excl <- c(FALSE, FALSE)
      delta.ub.excl <- c(FALSE, FALSE)
      delta.init <- c(1, mean(yi))
      delta.min <- c(0,   ifelse(alternative=="greater", min(yi)-sd(yi), min(yi)))
      delta.max <- c(100, ifelse(alternative=="greater", max(yi), max(yi)+sd(yi)))
      H0.delta <- c(1, 0)
      delta.LRT <- c(TRUE, FALSE)
      pval.min <- 0
      wi.fun <- function(x, delta, yi, vi, preci, alternative, steps) {
         if (alternative == "less") {
            yival <- qnorm(x, sd=sqrt(vi), lower.tail=TRUE)
            ifelse(yival < delta[2], 1, delta[1])
         } else {
            yival <- qnorm(x, sd=sqrt(vi), lower.tail=FALSE)
            ifelse(yival > delta[2], 1, delta[1])
         }
      }
      #.selmodel.ll <- ".selmodel.ll.cont"
      .selmodel.ll <- ".selmodel.ll.trunc"

   }

   ############################################################################

   ### checks on delta, delta.min, delta.max, and delta.init

   if (missing(delta)) {
      delta <- rep(NA_real_, deltas)
   } else {
      delta <- .expand1(delta, deltas)
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

   if (type == "stepfun") {
      if (decreasing) {
         delta[1] <- 1
      } else if (is.na(delta[1])) {
         delta[1] <- 1
      }
   }

   if (type == "stepcon")
      delta[1] <- 1

   if (!is.null(con$delta.min))
      delta.min <- con$delta.min

   delta.min <- .expand1(delta.min, deltas)
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

   delta.max <- .expand1(delta.max, deltas)
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

   delta.init <- .expand1(delta.init, deltas)
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

   # when ddd$defmap=TRUE or any delta.max value is infinity, use the default mapping functions defined
   # above for the various models (note that this will not be the case with the default settings);
   # otherwise use .mapfun() / .mapinvfun() or the functions passed via ddd$mapfun / ddd$mapinvfun

   if (.isTRUE(ddd$defmap) || any(is.infinite(delta.max))) {
      ddd$mapfun <- delta.transf.fun
      ddd$mapinvfun <- delta.transf.fun.inv
   }

   if (is.null(ddd$mapfun)) {
      mapfun <- rep(NA, deltas)
   } else {
      if (length(ddd$mapfun) == 1L) { # note: mapfun must be given as character string
         mapfun <- rep(ddd$mapfun, deltas)
      } else {
         mapfun <- ddd$mapfun
      }
   }

   if (is.null(ddd$mapinvfun)) {
      mapinvfun <- rep(NA, deltas)
   } else {
      if (length(ddd$mapinvfun) == 1L) { # note: mapinvfun must be given as character string
         mapinvfun <- rep(ddd$mapinvfun, deltas)
      } else {
         mapinvfun <- ddd$mapinvfun
      }
   }

   ### force use of certain transformation functions for mapfunv / mapinvfun for some special cases

   if (type == "truncest") {
      mapfun[2] <- "I"
      mapinvfun[2] <- "I"
   }

   ### remap initial delta values (except for the fixed ones)

   delta.init <- mapply(.mapinvfun, delta.init, delta.min, delta.max, mapinvfun)
   delta.init <- ifelse(is.na(delta), delta.init, delta)

   if (!is.null(con$pval.min))
      pval.min <- con$pval.min

   if (subset.spec) {
      if (sum(subset) < p + ifelse(is.element(x$method, c("FE","EE","CE")) || x$tau2.fix, 0, 1) + sum(is.na(delta)))
         stop(mstyle$stop(paste0("Number of studies (k_subset=", sum(subset), ") is too small to fit the selection model.")))
   } else {
      if (k < p + ifelse(is.element(x$method, c("FE","EE","CE")) || x$tau2.fix, 0, 1) + sum(is.na(delta)))
         stop(mstyle$stop(paste0("Number of studies (k=", k, ") is too small to fit the selection model.")))
   }

   ############################################################################

   pvals[pvals <   pval.min] <-   pval.min
   pvals[pvals > 1-pval.min] <- 1-pval.min

   if (type != "trunc" && stepsspec) {

      tmp <- .ptable(pvals, steps, subset)
      pgrp <- tmp$pgrp
      ptable <- tmp$ptable

      if (isTRUE(ddd$ptable))
         return(ptable)

      if (any(ptable[["k"]] == 0L)) {
         if (!isTRUE(ddd$skipintcheck) && type == "stepfun" && any(is.na(delta[-1])))
            warning(mstyle$warning(paste0("One or more intervals do not contain any observed p-values.")), call.=FALSE)
         if (!isTRUE(ddd$skipintcheck) && type != "stepfun")
            warning(mstyle$warning(paste0("One of the intervals does not contain any observed p-values.")), call.=FALSE)
      }

   } else {

      pgrp   <- NA
      ptable <- NA

   }

   ############################################################################

   ### model fitting

   if (verbose > 1)
      message(mstyle$message("\nModel fitting ...\n"))

   tmp <- .chkopt(optimizer, optcontrol, ineq=type=="stepcon")
   optimizer  <- tmp$optimizer
   optcontrol <- tmp$optcontrol
   par.arg    <- tmp$par.arg
   ctrl.arg   <- tmp$ctrl.arg

   ### set up default cluster when using optimParallel

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

   if (type == "stepcon") {

      if (optimizer == "Rsolnp::solnp")
         optcall <- paste0("Rsolnp::solnp(pars=c(beta.init, tau2.init, delta.init), fun=.selmodel.ll.stepfun,
            ineqfun=.rma.selmodel.ineqfun.pos, ineqLB=rep(0,deltas-1), ineqUB=rep(1,deltas-1),
            yi=yi, vi=vi, X=X, preci=preci, subset=subset, k=k, pX=p, pvals=pvals,
            deltas=deltas, delta.arg=delta, delta.transf=TRUE, mapfun=mapfun, delta.min=delta.min, delta.max=delta.max, decreasing=decreasing,
            tau2.arg=tau2, tau2.transf=TRUE, tau2.max=tau2.max, beta.arg=beta,
            wi.fun=wi.fun, steps=steps, pgrp=pgrp,
            alternative=alternative, pval.min=pval.min, intCtrl=intCtrl, verbose=verbose, digits=digits, dofit=FALSE", ctrl.arg, ")\n")

      if (optimizer == "nloptr::nloptr")
         optcall <- paste0("nloptr::nloptr(x0=c(beta.init, tau2.init, delta.init), eval_f=.selmodel.ll.stepfun,
            eval_g_ineq=.rma.selmodel.ineqfun.neg,
            yi=yi, vi=vi, X=X, preci=preci, subset=subset, k=k, pX=p, pvals=pvals,
            deltas=deltas, delta.arg=delta, delta.transf=TRUE, mapfun=mapfun, delta.min=delta.min, delta.max=delta.max, decreasing=decreasing,
            tau2.arg=tau2, tau2.transf=TRUE, tau2.max=tau2.max, beta.arg=beta,
            wi.fun=wi.fun, steps=steps, pgrp=pgrp,
            alternative=alternative, pval.min=pval.min, intCtrl=intCtrl, verbose=verbose, digits=digits, dofit=FALSE", ctrl.arg, ")\n")

      if (optimizer == "alabama::constrOptim.nl")
         optcall <- paste0("alabama::constrOptim.nl(par=c(beta.init, tau2.init, delta.init), fn=.selmodel.ll.stepfun,
            hin=.rma.selmodel.ineqfun.pos,
            yi=yi, vi=vi, X=X, preci=preci, subset=subset, k=k, pX=p, pvals=pvals,
            deltas=deltas, delta.arg=delta, delta.transf=TRUE, mapfun=mapfun, delta.min=delta.min, delta.max=delta.max, decreasing=decreasing,
            tau2.arg=tau2, tau2.transf=TRUE, tau2.max=tau2.max, beta.arg=beta,
            wi.fun=wi.fun, steps=steps, pgrp=pgrp,
            alternative=alternative, pval.min=pval.min, intCtrl=intCtrl, verbose=verbose, digits=digits, dofit=FALSE", ctrl.arg, ")\n")

   } else {

      optcall <- paste0(optimizer, "(", par.arg, "=c(beta.init, tau2.init, delta.init), ",
         .selmodel.ll, ", ", ifelse(optimizer=="optim", "method=optmethod, ", ""),
         "yi=yi, vi=vi, X=X, preci=preci, subset=subset, k=k, pX=p, pvals=pvals,
         deltas=deltas, delta.arg=delta, delta.transf=TRUE, mapfun=mapfun, delta.min=delta.min, delta.max=delta.max, decreasing=decreasing,
         tau2.arg=tau2, tau2.transf=TRUE, tau2.max=tau2.max, beta.arg=beta,
         wi.fun=wi.fun, steps=steps, pgrp=pgrp,
         alternative=alternative, pval.min=pval.min, intCtrl=intCtrl, verbose=verbose, digits=digits, dofit=FALSE", ctrl.arg, ")\n")

   }

   #return(optcall)

   .start.plot(verbose > 2)

   if (verbose) {
      opt.res <- try(eval(str2lang(optcall)), silent=!verbose)
   } else {
      opt.res <- try(suppressWarnings(eval(str2lang(optcall))), silent=!verbose)
   }

   if (isTRUE(ddd$retopt))
      return(opt.res)

   ### convergence checks (if verbose print optimParallel log, if verbose > 2 print opt.res, and unify opt.res$par)

   opt.res$par <- .chkconv(optimizer=optimizer, opt.res=opt.res, optcontrol=optcontrol, fun="selmodel", verbose=verbose)

   ### estimates/values of tau2 and delta on the transformed scale

   tau2.transf  <- opt.res$par[p+1]
   delta.transf <- opt.res$par[(p+2):(p+1+deltas)]

   ### save for Hessian computation

   beta.arg  <- beta
   tau2.arg  <- tau2
   delta.arg <- delta

   ### do the final model fit with estimated values

   fitcall <- paste0(.selmodel.ll, "(par=opt.res$par,
      yi=yi, vi=vi, X=X, preci=preci, subset=subset, k=k, pX=p, pvals=pvals,
      deltas=deltas, delta.arg=delta, delta.transf=TRUE, mapfun=mapfun, delta.min=delta.min, delta.max=delta.max, decreasing=decreasing,
      tau2.arg=tau2, tau2.transf=TRUE, tau2.max=tau2.max, beta.arg=beta,
      wi.fun=wi.fun, steps=steps, pgrp=pgrp,
      alternative=alternative, pval.min=pval.min, intCtrl=intCtrl, verbose=FALSE, digits=digits, dofit=TRUE)\n")

   #return(fitcall)
   fitcall <- try(eval(str2lang(fitcall)), silent=!verbose)
   #return(fitcall)

   if (inherits(fitcall, "try-error"))
      stop(mstyle$stop("Error during the optimization. Use verbose=TRUE and see help(selmodel) for more details on the optimization routines."))

   ll    <- fitcall$ll
   beta  <- cbind(fitcall$beta)
   tau2  <- fitcall$tau2
   delta <- fitcall$delta

   if ((delta.min.check && any(is.na(delta.arg) & delta <= delta.min + .Machine$double.eps^0.25)) || (delta.max.check && any(is.na(delta.arg) & delta >= delta.max - 100*.Machine$double.eps^0.25)))
      warning(mstyle$warning("One or more 'delta' estimates are (almost) equal to their lower or upper bound.\nTreat results with caution (or consider adjusting 'delta.min' and/or 'delta.max')."), call.=FALSE)

   ############################################################################

   ### computing (inverse) Hessian

   H        <- NA_real_
   vb       <- matrix(NA_real_, nrow=p, ncol=p)
   se.tau2  <- NA_real_
   vd       <- matrix(NA_real_, nrow=deltas, ncol=deltas)

   if (con$beta.fix) {
      beta.hes <- c(beta)
   } else {
      beta.hes <- beta.arg
   }

   if (con$tau2.fix || tau2 < con$tau2tol) {
      tau2.hes <- tau2
   } else {
      tau2.hes <- tau2.arg
   }

   if (con$delta.fix) {
      delta.hes <- delta
   } else {
      if (type == "stepfun") {
         delta.hes <- ifelse(delta < con$deltatol, delta, delta.arg)
      } else {
         delta.hes <- delta.arg
      }
   }

   hest <- c(is.na(beta.hes), is.na(tau2.hes), is.na(delta.hes))

   if (any(hest) && !isTRUE(ddd$skiphes)) {

      if (verbose > 1)
         message(mstyle$message("\nComputing Hessian ..."))

      if (verbose > 3)
         cat("\n")

      if (con$htransf) { # TODO: document these two possibilities?

         if (con$hesspack == "numDeriv")
            hescall <- paste0("numDeriv::hessian(", .selmodel.ll, ", x=opt.res$par, method.args=con$hessianCtrl,
               yi=yi, vi=vi, X=X, preci=preci, subset=subset, k=k, pX=p, pvals=pvals,
               deltas=deltas, delta.arg=delta.hes, delta.transf=TRUE, mapfun=mapfun, delta.min=delta.min, delta.max=delta.max, decreasing=decreasing,
               tau2.arg=tau2.hes, tau2.transf=TRUE, tau2.max=tau2.max, beta.arg=beta.hes,
               wi.fun=wi.fun, steps=steps, pgrp=pgrp,
               alternative=alternative, pval.min=pval.min, intCtrl=intCtrl, verbose=ifelse(verbose > 3, verbose, 0), digits=digits)\n")

         if (con$hesspack == "pracma")
            hescall <- paste0("pracma::hessian(", .selmodel.ll, ", x0=opt.res$par,
               yi=yi, vi=vi, X=X, preci=preci, subset=subset, k=k, pX=p, pvals=pvals,
               deltas=deltas, delta.arg=delta.hes, delta.transf=TRUE, mapfun=mapfun, delta.min=delta.min, delta.max=delta.max, decreasing=decreasing,
               tau2.arg=tau2.hes, tau2.transf=TRUE, tau2.max=tau2.max, beta.arg=beta.hes,
               wi.fun=wi.fun, steps=steps, pgrp=pgrp,
               alternative=alternative, pval.min=pval.min, intCtrl=intCtrl, verbose=ifelse(verbose > 3, verbose, 0), digits=digits)\n")

      } else {

         ### this is the default

         if (con$hesspack == "numDeriv")
            hescall <- paste0("numDeriv::hessian(", .selmodel.ll, ", x=c(beta, tau2, delta), method.args=con$hessianCtrl,
               yi=yi, vi=vi, X=X, preci=preci, subset=subset, k=k, pX=p, pvals=pvals,
               deltas=deltas, delta.arg=delta.hes, delta.transf=FALSE, mapfun=mapfun, delta.min=delta.min, delta.max=delta.max, decreasing=decreasing,
               tau2.arg=tau2.hes, tau2.transf=FALSE, tau2.max=tau2.max, beta.arg=beta.hes,
               wi.fun=wi.fun, steps=steps, pgrp=pgrp,
               alternative=alternative, pval.min=pval.min, intCtrl=intCtrl, verbose=ifelse(verbose > 3, verbose, 0), digits=digits)\n")

         if (con$hesspack == "pracma")
            hescall <- paste0("pracma::hessian(", .selmodel.ll, ", x0=c(beta, tau2, delta),
               yi=yi, vi=vi, X=X, preci=preci, subset=subset, k=k, pX=p, pvals=pvals,
               deltas=deltas, delta.arg=delta.hes, delta.transf=FALSE, mapfun=mapfun, delta.min=delta.min, delta.max=delta.max, decreasing=decreasing,
               tau2.arg=tau2.hes, tau2.transf=FALSE, tau2.max=tau2.max, beta.arg=beta.hes,
               wi.fun=wi.fun, steps=steps, pgrp=pgrp,
               alternative=alternative, pval.min=pval.min, intCtrl=intCtrl, verbose=ifelse(verbose > 3, verbose, 0), digits=digits)\n")

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
               vd[is.na(delta.hes), is.na(delta.hes)] <- iH[c(rep(FALSE,p),FALSE,is.na(delta.hes)), c(rep(FALSE,p),FALSE,is.na(delta.hes)), drop=FALSE]
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
      QM <- NA_real_

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

   se.delta <- sqrt(diag(vd))

   if (con$htransf) {

      zval.delta <- rep(NA_real_, deltas)
      pval.delta <- rep(NA_real_, deltas)
      ci.lb.delta <- c(delta.transf - crit * se.delta)
      ci.ub.delta <- c(delta.transf + crit * se.delta)
      ci.lb.delta <- mapply(.mapfun, ci.lb.delta, delta.min, delta.max, mapfun)
      ci.ub.delta <- mapply(.mapfun, ci.ub.delta, delta.min, delta.max, mapfun)
      vd <- matrix(NA_real_, nrow=deltas, ncol=deltas)
      se.delta <- rep(NA_real_, deltas)

   } else {

      zval.delta  <- (delta - H0.delta) / se.delta
      pval.delta  <- 2*pnorm(abs(zval.delta), lower.tail=FALSE)
      ci.lb.delta <- c(delta - crit * se.delta)
      ci.ub.delta <- c(delta + crit * se.delta)

   }

   ### impose constraints on the CI bounds for the delta value(s)

   ci.lb.delta <- ifelse(ci.lb.delta < delta.lb, delta.lb, ci.lb.delta)
   ci.ub.delta <- ifelse(ci.ub.delta > delta.ub, delta.ub, ci.ub.delta)

   ci.lb.delta <- ifelse(ci.lb.delta < delta.min, delta.min, ci.lb.delta)
   ci.ub.delta <- ifelse(ci.ub.delta > delta.max, delta.max, ci.ub.delta)

   ### inference for tau^2 parameter

   if (con$htransf) {
      ci.lb.tau2 <- exp(tau2.transf - crit * se.tau2) # tau2.transf = log(tau^2) and se.tau2 = SE[log(tau^2)]
      ci.ub.tau2 <- exp(tau2.transf + crit * se.tau2)
      se.tau2 <- se.tau2 * exp(tau2.transf) # delta method
   } else {
      ci.lb.tau2 <- tau2 - crit * se.tau2
      ci.ub.tau2 <- tau2 + crit * se.tau2
   }

   ci.lb.tau2[ci.lb.tau2 < 0] <- 0

   ############################################################################

   ### LRT for H0: tau^2 = 0 (only when NOT fitting a FE model)

   LRT.tau2  <- NA_real_
   LRTp.tau2 <- NA_real_

   if (!x$tau2.fix && !is.element(x$method, c("FE","EE","CE")) && !isTRUE(ddd$skiphet)) {

      if (verbose > 1)
         message(mstyle$message("Conducting heterogeneity test ..."))

      if (verbose > 4)
         cat("\n")

      optcall <- paste0(optimizer, "(", par.arg, "=c(beta.init, tau2.init, delta.init), ",
         .selmodel.ll, ", ", ifelse(optimizer=="optim", "method=optmethod, ", ""),
         "yi=yi, vi=vi, X=X, preci=preci, subset=subset, k=k, pX=p, pvals=pvals,
         deltas=deltas, delta.arg=delta.arg, delta.transf=TRUE, mapfun=mapfun, delta.min=delta.min, delta.max=delta.max, decreasing=decreasing,
         tau2.arg=0, tau2.transf=FALSE, tau2.max=tau2.max, beta.arg=beta.arg,
         wi.fun=wi.fun, steps=steps, pgrp=pgrp,
         alternative=alternative, pval.min=pval.min, intCtrl=intCtrl, verbose=ifelse(verbose > 4, verbose, 0), digits=digits", ctrl.arg, ")\n")

      opt.res <- try(eval(str2lang(optcall)), silent=!verbose)

      if (verbose > 4)
         cat("\n")

      if (!inherits(opt.res, "try-error")) {

         fitcall <- paste0(.selmodel.ll, "(par=opt.res$par,
            yi=yi, vi=vi, X=X, preci=preci, subset=subset, k=k, pX=p, pvals=pvals,
            deltas=deltas, delta.arg=delta.arg, delta.transf=TRUE, mapfun=mapfun, delta.min=delta.min, delta.max=delta.max, decreasing=decreasing,
            tau2.arg=0, tau2.transf=FALSE, tau2.max=tau2.max, beta.arg=beta.arg,
            wi.fun=wi.fun, steps=steps, pgrp=pgrp,
            alternative=alternative, pval.min=pval.min, intCtrl=intCtrl, verbose=FALSE, digits=digits, dofit=TRUE)\n")

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
   LRTdf <- sum(is.na(delta.arg) & delta.LRT)
   LRTp  <- ifelse(LRTdf > 0, pchisq(LRT, df=LRTdf, lower.tail=FALSE), NA_real_)

   ############################################################################

   ### fit statistics

   if (verbose > 1)
      message(mstyle$message("Computing fit statistics and log-likelihood ..."))

   ### note: tau2 and delta are not counted as parameters when they were fixed by the user
   parms <- p + ifelse(is.element(x$method, c("FE","EE","CE")) || x$tau2.fix, 0, 1) + sum(is.na(delta.arg))

   ll.ML   <- ll
   dev.ML  <- -2 * ll.ML
   AIC.ML  <- -2 * ll.ML + 2*parms
   BIC.ML  <- -2 * ll.ML +   parms * log(k)
   AICc.ML <- -2 * ll.ML + 2*parms * max(k, parms+2) / (max(k, parms+2) - parms - 1)

   fit.stats <- matrix(c(ll.ML, dev.ML, AIC.ML, BIC.ML, AICc.ML, ll.REML=NA_real_, dev.REML=NA_real_, AIC.REML=NA_real_, BIC.REML=NA_real_, AICc.REML=NA_real_), ncol=2, byrow=FALSE)
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

   res$dfs <- res$ddf <- NA_integer_
   res$test <- "z"
   res$s2w  <- 1

   res$QE <- res$QEp <- NA_real_
   res$I2 <- res$H2 <- res$vt <- NA_real_
   res$R2 <- NULL

   res$QM  <- QM
   res$QMp <- QMp

   res$delta       <- delta
   res$vd          <- vd
   res$se.delta    <- se.delta
   res$zval.delta  <- zval.delta
   res$pval.delta  <- pval.delta
   res$ci.lb.delta <- ci.lb.delta
   res$ci.ub.delta <- ci.ub.delta
   res$deltas      <- deltas
   res$delta.fix   <- !is.na(delta.arg)

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
   res$decreasing  <- decreasing
   res$stepsspec   <- stepsspec
   res$pgrp        <- pgrp
   res$ptable      <- ptable
   res$k0          <- sum(!subset)
   res$k1          <- sum(subset)
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
