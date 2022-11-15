# check for skewness according to Shi

conv.fivenum <- function(min, q1, median, q3, max, n, data, include, var.names=c("mean","sd"), append=TRUE, method) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (missing(min) && missing(q1) && missing(median) && missing(q3) && missing(max))
      stop(mstyle$stop("Must specify at least some of these arguments: 'min', 'q1', 'median', 'q3', 'max'."))

   if (missing(data))
      data <- NULL

   if (is.null(data)) {
      data <- sys.frame(sys.parent())
   } else {
      if (!is.data.frame(data))
         data <- data.frame(data)
   }

   has.data <- !is.null(data)

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
   include <- .getx("include", mf=mf, data=data, checknumeric=TRUE)

   if (!.equal.length(min, q1, median, q3, max, n))
      stop(mstyle$stop("Supplied data vectors are not all of the same length."))

   k <- max(length(min), length(q1), length(median), length(q3), length(max), length(n))

   means <- rep(NA_real_, k)
   sds   <- rep(NA_real_, k)

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

   ### if include is NULL, set to TRUE vector

   if (is.null(include))
      include <- rep(TRUE, k)

   ### turn numeric include vector into logical vector

   include <- .chksubset(include, k)

   ### exclude rows where n < 5

   include[which(n < 5)] <- FALSE

   ### set inputs to NA for rows not to be included

   min[!include]    <- NA_real_
   q1[!include]     <- NA_real_
   median[!include] <- NA_real_
   q3[!include]     <- NA_real_
   max[!include]    <- NA_real_
   n[!include]      <- NA_real_

   ### check min <= q1 <= median <= q3 <= max

   incorder <- apply(cbind(min, q1, median, q3, max), 1, is.unsorted, na.rm=TRUE)

   if (any(incorder)) {
      incrows <- which(incorder)
      if (length(incrows) > 5L) {
         stop(mstyle$stop(paste0("Found 'min <= q1 <= median <= q3 <= max' not true for row(s): ", paste0(incrows[1:5], collapse=","), ",...")))
      } else {
         stop(mstyle$stop(paste0("Found 'min <= q1 <= median <= q3 <= max' not true for row(s): ", paste0(incrows, collapse=","))))
      }
   }

   #########################################################################

   # determine cases and check methods

   case1 <- !is.na(min) &  is.na(q1) & !is.na(median) &  is.na(q3) & !is.na(max)
   case2 <-  is.na(min) & !is.na(q1) & !is.na(median) & !is.na(q3) &  is.na(max)
   case3 <- !is.na(min) & !is.na(q1) & !is.na(median) & !is.na(q3) & !is.na(max)

   if (missing(method))
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

   for (i in 1:k) {

      if (case1[i]) {

         ### case 1: min, median, and max are given

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
            sds[i] <- (max[i] - min[i]) / xi
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

      }

      if (case2[i]) {

         ### case 2: q1, median, and q3 are given

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
            sds[i] <- (q3[i] - q1[i]) / eta
         }

      }

      if (case3[i]) {

         ### case 3: min, q1, median, q3, and max are given

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
            sds[i] <- weight * (max[i] - min[i]) / xi + (1 - weight) * (q3[i] - q1[i]) / eta
         }
         if (method[2] == "wan2014")
            sds[i] <- 1/2 * ((max[i] - min[i]) / xi + (q3[i] - q1[i]) / eta)
         if (method[2] == "bland2015")
            sds[i] <- sqrt((min[i]^2 + 2*q1[i]^2 + 2*median[i]^2 + 2*q3[i]^2 + max[i]^2) / 16 + (min[i]*q1[i] + q1[i]*median[i] + median[i]*q3[i] + q3[i]*max[i]) / 8 - (min[i] + 2*q1[i] + 2*median[i] + 2*q3[i] + max[i])^2 / 64)

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
         attr(data[[var.names[2]]], "est") <- is.na(data[[var.names[2]]]) & !is.na(means)
         data[[var.names[2]]] <- replmiss(data[[var.names[2]]], sds)
      } else {
         data <- cbind(data, sds)
         names(data)[length(names(data))] <- var.names[2]
      }

   } else {

      data <- data.frame(means, sds)
      names(data) <- var.names

   }

   return(data)

}
