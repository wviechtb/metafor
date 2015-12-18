escalc.data.frame <- function(measure, formula, slab, subset, var.names=c("yi","vi"), add.measure=FALSE, digits=4, ...) {

   ### check for some incorrect argument specifications

   if (!is.element(measure, c("RR","OR","PETO","RD","AS","PHI","YUQ","YUY","RTET", ### 2x2 table measures
                              "PBIT","OR2D","OR2DN","OR2DL",                       ### - transformations to SMD
                              "IRR","IRD","IRSD",                                  ### two-group person-time data measures
                              "MD","SMD","SMDH","ROM",                             ### two-group mean/SD measures
                              "RPB","RBIS","D2OR","D2ORN","D2ORL",                 ### - transformations to r_PB, r_BIS, and log(OR)
                              "COR","UCOR","ZCOR",                                 ### correlations (raw and r-to-z transformed)
                              "PR","PLN","PLO","PAS","PFT",                        ### single proportions (and transformations thereof)
                              "IR","IRLN","IRS","IRFT",                            ### single-group person-time data (and transformations thereof)
                              "MN","MC","SMCC","SMCR","SMCRH","ROMC",              ### raw/standardized mean change and log(ROM) for dependent samples
                              "ARAW","AHW","ABT",                                  ### alpha (and transformations thereof)
                              "GEN")))
      stop("Unknown 'measure' specified.")

   if (add.measure) {

      if (length(var.names) == 2)
         var.names <- c(var.names, "measure")

      if (length(var.names) != 3)
         stop("Argument 'var.names' must be of length 2 or 3.")

   } else {

      if (length(var.names) == 3)
         var.names <- var.names[1:2]

      if (length(var.names) != 2)
         stop("Argument 'var.names' must be of length 2.")

   }

   data <- formula

   ### extract yi/vi variables

   yi <- data[[var.names[1]]]
   vi <- data[[var.names[2]]]

   ### check if resulting vectors are not NULL

   if (is.null(yi))
      stop(paste0("Could not find variable '", var.names[1], "' in data frame."))

   if (is.null(vi))
      stop(paste0("Could not find variable '", var.names[2], "' in data frame."))

   mf <- match.call()

   ### get slab and subset arguments (will be NULL when unspecified)

   mf.slab   <- mf[[match("slab",   names(mf))]]
   mf.subset <- mf[[match("subset", names(mf))]]
   slab      <- eval(mf.slab,   data, enclos=sys.frame(sys.parent()))
   subset    <- eval(mf.subset, data, enclos=sys.frame(sys.parent()))

   ### subsetting of yi/vi values

   if (!is.null(subset)) {
      yi <- yi[subset]
      vi <- vi[subset]
   }

   #########################################################################
   #########################################################################
   #########################################################################

   ### check for infinite values and set them to NA

   is.inf <- is.infinite(yi) | is.infinite(vi)

   if (any(is.inf)) {
      warning("Some 'yi' and/or 'vi' values equal to +-Inf. Recoded to NAs.")
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

      if (!is.null(subset))
         slab <- slab[subset]

      if (anyNA(slab))
         stop("NAs in study labels.")

      ### check if study labels are unique; if not, make them unique

      if (anyDuplicated(slab))
         slab <- make.unique(as.character(slab)) ### make.unique() only works with character vectors

      if (length(slab) != length(yi))
         stop("Study labels not of same length as data.")

      ### add slab attribute to the yi vector
      attr(yi, "slab") <- slab

   }

   ### if a subset of studies is specified (note: subsetting of other parts already done above, so yi/vi/slab are already subsetted)

   if (!is.null(subset))
      data <- data[subset,,drop=FALSE]

   ### add measure attribute to the yi vector

   attr(yi, "measure") <- measure

   ### put together dataset

   dat <- data.frame(data)

   dat[[var.names[1]]] <- yi
   dat[[var.names[2]]] <- vi

   if (add.measure) {
      dat[[var.names[3]]] <- ""
      dat[[var.names[3]]][!is.na(yi)] <- measure
   }

   attr(dat, "digits") <- digits

   ### add 'yi.names' and 'vi.names' to the first position of the corresponding attributes
   attr(dat, "yi.names") <- unique(c(var.names[1], attr(data, "yi.names"))) ### if 'yi.names' is not an attribute, attr() returns NULL, so this works fine
   attr(dat, "vi.names") <- unique(c(var.names[2], attr(data, "vi.names"))) ### if 'vi.names' is not an attribute, attr() returns NULL, so this works fine

   class(dat) <- c("escalc", "data.frame")
   return(dat)

}
