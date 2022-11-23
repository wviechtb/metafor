conv.fivenum <- function(min, q1, median, q3, max, n, data, dist="norm", transf=FALSE,
                         include, test=TRUE, var.names=c("mean","sd"), append=TRUE, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (missing(min) && missing(q1) && missing(median) && missing(q3) && missing(max))
      stop(mstyle$stop("Must specify at least some of these arguments: 'min', 'q1', 'median', 'q3', 'max'."))

   ### get ... argument and check for extra/superfluous arguments

   ddd <- list(...)

   .chkdots(ddd, c("method"))

   if (missing(data))
      data <- NULL

   has.data <- !is.null(data)

   if (is.null(data)) {
      data <- sys.frame(sys.parent())
   } else {
      if (!is.data.frame(data))
         data <- data.frame(data)
   }

   ### checks on var.names argument

   if (length(var.names) != 2L)
      stop(mstyle$stop("Argument 'var.names' must be of length 2."))

   if (any(var.names != make.names(var.names, unique=TRUE))) {
      var.names <- make.names(var.names, unique=TRUE)
      warning(mstyle$warning(paste0("Argument 'var.names' does not contain syntactically valid variable names.\nVariable names adjusted to: var.names = c('", var.names[1], "', '", var.names[2], "').")), call.=FALSE)
   }

   mf <- match.call()

   min     <- .getx("min",     mf=mf, data=data, checknumeric=TRUE)
   q1      <- .getx("q1",      mf=mf, data=data, checknumeric=TRUE)
   median  <- .getx("median",  mf=mf, data=data, checknumeric=TRUE)
   q3      <- .getx("q3",      mf=mf, data=data, checknumeric=TRUE)
   max     <- .getx("max",     mf=mf, data=data, checknumeric=TRUE)
   n       <- .getx("n",       mf=mf, data=data, checknumeric=TRUE)
   include <- .getx("include", mf=mf, data=data)
   dist    <- .getx("dist",    mf=mf, data=data)

   if (is.null(dist)) # if 'dist' is not specified, it is NULL
      dist <- "norm"

   if (!.equal.length(min, q1, median, q3, max, n))
      stop(mstyle$stop("Supplied data vectors are not all of the same length."))

   k <- max(length(min), length(q1), length(median), length(q3), length(max), length(n))

   if (is.null(min))
      min <- rep(NA_real_, k)
   if (is.null(q1))
      q1 <- rep(NA_real_, k)
   if (is.null(median))
      median <- rep(NA_real_, k)
   if (is.null(q3))
      q3 <- rep(NA_real_, k)
   if (is.null(max))
      max <- rep(NA_real_, k)
   if (is.null(n))
      n <- rep(NA_real_, k)

   ### handle dist argument

   if (length(dist) == 1L)
      dist <- rep(dist, k)

   if (length(dist) != k)
      stop(mstyle$stop(paste0("Length of 'dist' argument (", length(dist), ") does not match length of data (", k, ").")))

   dist <- c("norm","lnorm")[pmatch(dist, c("norm","lnorm"), duplicates.ok=TRUE)]

   if (anyNA(dist))
      stop(mstyle$stop("Unknown 'dist' value specified (should either be 'norm' or 'lnorm')."))

   ### if include is NULL, set to TRUE vector

   if (is.null(include))
      include <- rep(TRUE, k)

   ### turn numeric include vector into logical vector

   include <- .chksubset(include, k, stoponk0=FALSE)

   ### exclude rows where n < 5

   include[which(n < 5)] <- FALSE

   ### set inputs to NA for rows not to be included

   min[!include]    <- NA_real_
   q1[!include]     <- NA_real_
   median[!include] <- NA_real_
   q3[!include]     <- NA_real_
   max[!include]    <- NA_real_
   n[!include]      <- NA_real_

   #########################################################################

   # determine cases and check methods

   case1 <- !is.na(min) &  is.na(q1) &  is.na(q3) & !is.na(max)
   case2 <-  is.na(min) & !is.na(q1) & !is.na(q3) &  is.na(max)
   case3 <- !is.na(min) & !is.na(q1) & !is.na(q3) & !is.na(max)

   if (is.null(ddd$method)) {
      method <- "rec"
   } else {
      method <- ddd$method
   }

   if (any(dist == "lnorm")) # if any dist value is 'lnorm', force use of 'recommended' methods
      method <- "rec"

   if (length(method) == 1L)
      method <- c(method, method)

   method <- tolower(method)

   method1.options <- c("rec", "hozo2005", "wan2014", "bland2015", "luo2016")
   method[1] <- method1.options[grep(method[1], method1.options)[1]]
   if (is.na(method[1]))
      stop(mstyle$stop("Unknown 'method[1]' specified."))

   method2.options <- c("rec", "hozo2005", "wan2014", "bland2015", "shi2020")
   method[2] <- method2.options[grep(method[2], method2.options)[1]]
   if (is.na(method[2]))
      stop(mstyle$stop("Unknown 'method[2]' specified."))

   #########################################################################

   means <- rep(NA_real_, k)
   sds   <- rep(NA_real_, k)
   tval  <- rep(NA_real_, k)
   crit  <- rep(NA_real_, k)
   sig   <- rep(NA,       k)

   for (i in 1:k) {

      ### check min <= q1 <= median <= q3 <= max

      if (is.unsorted(c(min[i], q1[i], median[i], q3[i], max[i]), na.rm=TRUE))
         stop(mstyle$stop(paste0("Found 'min <= q1 <= median <= q3 <= max' not true for row ", i, ".")))

      if (dist[i] == "lnorm") {

         ### check that min, q1, median, q3, and max are all > 0 when assuming a log-normal distribution

         if (any(min[i] <= 0 || q1[i] <= 0 || median[i] <= 0 || q3[i] <= 0 || max[i] <= 0, na.rm=TRUE))
            stop(mstyle$stop(paste0("Found a non-positive value of 'min', 'q1', 'median', 'q3', or 'max' for row ", i, ".")))

         ### log-transform inputs

         min[i]    <- log(min[i])
         q1[i]     <- log(q1[i])
         median[i] <- log(median[i])
         q3[i]     <- log(q3[i])
         max[i]    <- log(max[i])

      }

      if (case1[i]) {

         ### case 1: min, median, and max are given

         # test for skewness

         tval[i] <- abs((min[i] + max[i] - 2*median[i]) / (max[i] - min[i]))
         #crit[i] <- 1.01 / log(n[i] + 9) + 2.43 / (n[i] + 1) # Shi et al. (2020b)
         crit[i] <- 1 / log(n[i] + 9) + 2.5 / (n[i] + 1)      # Shi et al. (under review)
         sig[i]  <- isTRUE(tval[i] >= crit[i])

         if (test && sig[i])
            next

         # mean estimation

         if (is.element(method[1], c("rec", "luo2016"))) {
            # Luo et al. (2016), equation (7)
            weight <- 4 / (4 + n[i]^0.75)
            means[i] <- weight * (min[i] + max[i]) / 2 + (1 - weight) * median[i]
         }
         if (method[1] == "hozo2005") {
            if (is.na(n[i])) {
               means[i] <- NA
            } else if (n[i] <= 25) {
               means[i] <- (min[i] + 2*median[i] + max[i]) / 4
            } else {
               means[i] <- median[i]
            }
         }
         if (method[1] == "wan2014")
            means[i] <- (min[i] + 2*median[i] + max[i]) / 4

         # sd estimation

         if (is.element(method[2], c("rec", "wan2014"))) {
            # Wan et al. (2014), equation (9)
            xi <- 2 * qnorm((n[i] - 0.375) / (n[i] + 0.25))
            z1 <- ifelse(dist[i] == "norm", 1, 1.01 + 0.25 / log(n[i])^2)
            #z1 <- 1
            sds[i] <- (max[i] - min[i]) / xi * (1/sqrt(z1))
         }
         if (method[2] == "hozo2005") {
            if (is.na(n[i])) {
               sds[i] <- NA
            } else if (n[i] <= 15) {
               sds[i] <- 1/sqrt(12) * sqrt((min[i] - 2*median[i] + max[i])^2 / 4 + (max[i]-min[i])^2)
            } else if (n[i] <= 70) {
               sds[i] <- (max[i] - min[i]) / 4
            } else {
               sds[i] <- (max[i] - min[i]) / 6
            }
         }

         if (dist[i] == "lnorm" && transf) {
            s41 <- ((max[i] - min[i]) / xi)^4 / (1 + 2.23 / log(n[i])^2)
            phi1 <- 1 + 0.565 * sds[i]^2 / n[i] + 0.37 * s41 / n[i]
            btmean <- exp(means[i] + sds[i]^2 / 2) * (1 / phi1)
            phi11 <- 1 + 2.26 * sds[i]^2 / n[i] + 5.92 * s41 / n[i]
            phi12 <- 1 + 2.26 * sds[i]^2 / n[i] + 1.48 * s41 / n[i]
            btsd   <- sqrt(exp(2*means[i] + 2*sds[i]^2) * (1 / phi11) - exp(2*means[i] + sds[i]^2) * (1 / phi12))
            means[i] <- btmean
            sds[i]   <- btsd
         }

      }

      if (case2[i]) {

         ### case 2: q1, median, and q3 are given

         # test for skewness

         tval[i] <- abs((q1[i] + q3[i] - 2*median[i]) / (q3[i] - q1[i]))
         #crit[i] <- 2.66 / sqrt(n[i]) - 5.92 / n[i]^2 # Shi et al. (2020b)
         crit[i] <- 2.65 / sqrt(n[i]) - 6 / n[i]^2     # Shi et al. (under review)
         sig[i] <- isTRUE(tval[i] >= crit[i])

         if (test && sig[i])
            next

         # mean estimation

         if (is.element(method[1], c("rec", "luo2016"))) {
            # Luo et al. (2016), equation (11)
            weight <- 0.7 + 0.39 / n[i]
            #weight <- 0.699 + 0.4 / n[i]
            means[i] <- weight * (q1[i] + q3[i]) / 2 + (1 - weight) * median[i]
         }
         if (method[1] == "wan2014")
            means[i] <- (q1[i] + median[i] + q3[i]) / 3

         # sd estimation

         if (is.element(method[2], c("rec", "wan2014"))) {
            # Wan et al. (2014), equation (16)
            eta <- 2 * qnorm((0.75*n[i] - 0.125) / (n[i] + 0.25))
            z2 <- ifelse(dist[i] == "norm", 1, 1 + 1.58 / n[i])
            #z2 <- 1
            sds[i] <- (q3[i] - q1[i]) / eta * (1/sqrt(z2))
         }

         if (dist[i] == "lnorm" && transf) {
            s42 <- ((q3[i] - q1[i]) / eta)^4 / (1 + 19.2 / n[i]^1.2)
            phi2 <- 1 + 0.57 * sds[i]^2 / n[i] + 0.75 * s42 / n[i]
            btmean <- exp(means[i] + sds[i]^2 / 2) * (1 / phi2)
            phi21 <- 1 + 2.28 * sds[i]^2 / n[i] + 12 * s42 / n[i]
            phi22 <- 1 + 2.28 * sds[i]^2 / n[i] +  3 * s42 / n[i]
            btsd   <- sqrt(exp(2*means[i] + 2*sds[i]^2) * (1 / phi21) - exp(2*means[i] + sds[i]^2) * (1 / phi22))
            means[i] <- btmean
            sds[i]   <- btsd
         }

      }

      if (case3[i]) {

         ### case 3: min, q1, median, q3, and max are given

         # test for skewness

         tval[i] <- max(2.65 * log(0.6 * n[i]) / sqrt(n[i]) * abs((min[i] + max[i] - 2*median[i]) / (max[i] - min[i])), abs((q1[i] + q3[i] - 2*median[i]) / (q3[i] - q1[i])))
         #crit[i] <- 2.97 / sqrt(n[i]) - 39.1 / n[i]^3 # Shi et al. (2020b)
         crit[i] <- 3 / sqrt(n[i]) - 40 / n[i]^3       # Shi et al. (under review)
         sig[i] <- isTRUE(tval[i] >= crit[i])

         if (test && sig[i])
            next

         # mean estimation

         if (is.element(method[1], c("rec", "luo2016"))) {
            # Luo et al. (2016), equation (15)
            weight1 <- 2.2 / (2.2 + n[i]^0.75)
            weight2 <- 0.7 - 0.72 / n[i]^0.55
            means[i] <- weight1 * (min[i] + max[i]) / 2 + weight2 * (q1[i] + q3[i]) / 2 + (1 - weight1 - weight2) * median[i]
         }
         if (is.element(method[1], c("wan2014", "bland2015")))
            means[i] <- (min[i] + 2*q1[i] + 2*median[i] + 2*q3[i] + max[i]) / 8

         # sd estimation

         if (is.element(method[2], c("rec", "shi2020", "wan2014"))) {
            xi <- 2 * qnorm((n[i] - 0.375) / (n[i] + 0.25))
            eta <- 2 * qnorm((0.75*n[i] - 0.125) / (n[i] + 0.25))
         }
         if (is.element(method[2], c("rec", "shi2020"))) {
            # Shi et al. (2020), equation (10)
            weight <- 1 / (1 + 0.07 * n[i]^0.6)
            z3 <- ifelse(dist[i] == "norm", 1, 1 + 0.28 / log(n[i])^2)
            #z3 <- 1
            sds[i] <- (weight * (max[i] - min[i]) / xi + (1 - weight) * (q3[i] - q1[i]) / eta) * (1/sqrt(z3))
         }
         if (method[2] == "wan2014")
            sds[i] <- 1/2 * ((max[i] - min[i]) / xi + (q3[i] - q1[i]) / eta)
         if (method[2] == "bland2015")
            sds[i] <- sqrt((min[i]^2 + 2*q1[i]^2 + 2*median[i]^2 + 2*q3[i]^2 + max[i]^2) / 16 + (min[i]*q1[i] + q1[i]*median[i] + median[i]*q3[i] + q3[i]*max[i]) / 8 - (min[i] + 2*q1[i] + 2*median[i] + 2*q3[i] + max[i])^2 / 64)

         if (dist[i] == "lnorm" && transf) {
            s43 <- (weight * (max[i] - min[i]) / xi + (1 - weight) * (q3[i] - q1[i]) / eta)^4 / (1 + 3.93 / n[i])
            phi3 <- 1 + 0.405 * sds[i]^2 / n[i] + 0.315 * s43 / n[i]
            btmean <- exp(means[i] + sds[i]^2 / 2) * (1 / phi3)
            phi31 <- 1 + 1.62 * sds[i]^2 / n[i] + 5.04 * s43 / n[i]
            phi32 <- 1 + 1.62 * sds[i]^2 / n[i] + 1.26 * s43 / n[i]
            btsd   <- sqrt(exp(2*means[i] + 2*sds[i]^2) * (1 / phi31) - exp(2*means[i] + sds[i]^2) * (1 / phi32))
            means[i] <- btmean
            sds[i]   <- btsd
         }

      }

   }

   #########################################################################

   if (has.data && append) {

      if (is.element(var.names[1], names(data))) {
         attr(data[[var.names[1]]], "est") <- is.na(data[[var.names[1]]]) & !is.na(means)
         data[[var.names[1]]] <- replmiss(data[[var.names[1]]], means)
      } else {
         data <- cbind(data, means)
         names(data)[length(names(data))] <- var.names[1]
      }

      if (is.element(var.names[2], names(data))) {
         attr(data[[var.names[2]]], "est") <- is.na(data[[var.names[2]]]) & !is.na(sds)
         data[[var.names[2]]] <- replmiss(data[[var.names[2]]], sds)
      } else {
         data <- cbind(data, sds)
         names(data)[length(names(data))] <- var.names[2]
      }

   } else {

      data <- data.frame(means, sds)
      names(data) <- var.names

   }

   attr(data[[var.names[1]]], "tval") <- tval
   attr(data[[var.names[1]]], "crit") <- crit
   attr(data[[var.names[1]]], "sig")  <- sig
   attr(data[[var.names[2]]], "tval") <- tval
   attr(data[[var.names[2]]], "crit") <- crit
   attr(data[[var.names[2]]], "sig")  <- sig

   return(data)

}
