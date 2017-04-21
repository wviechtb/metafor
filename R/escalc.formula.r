### 2x2 table data:  escalc(outcome ~ group | study, weights = freq) # outcome is a factor with two levels, group is a factor with two levels
### IRD/IRR/IRSD:    escalc(events/times ~ group | study)            # group is a factor with two levels
### MD and others:   escalc(mean/sd ~ group | study, weights = ni)   # group is a factor with two levels
### COR and others:  escalc(cor ~ 1 | study, weights = ni)           #
### PR and others:   escalc(outcome ~ 1 | study, weights = freq)     # outcome is a factor with two levels
### IR and others:   escalc(events/times ~ 1 | study)                #
### MN:              escalc(mean/sd ~ 1 | study, weights = ni)       #
### mean change:                                                     # not currently implemented
### ARAW and others: escalc(alpha/items ~ 1 | study, weights=ni)     #

escalc.formula <- function(measure, formula, weights, data, # slab, subset,
add=1/2, to="only0", drop00=FALSE, vtype="LS", var.names=c("yi","vi"), digits=4, ...) {

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
                              "MN","MNLN","CVLN","SDLN",                           ### mean, log(mean), log(CV), log(SD)
                              "MC","SMCC","SMCR","SMCRH","ROMC",                   ### raw/standardized mean change and log(ROM) for dependent samples
                              "ARAW","AHW","ABT")))                                ### alpha (and transformations thereof)
      stop("Unknown 'measure' specified.")

   if (is.element(measure, c("MPRD","MPRR","MPOR","MPORC","MPPETO","CVR","VR","PCOR","ZPCOR","SPCOR","CVLN","SDLN","MC","SMCC","SMCR","SMCRH","ROMC")))
      stop("Formula interface (currently) not implemented for this outcome measure.")

   if (!requireNamespace("Formula", quietly=TRUE))
      stop("Please install the 'Formula' package to use the formula interface.")

   ### make sure that NAs are not dropped from model frame

   na.act <- getOption("na.action")
   options(na.action="na.pass")
   on.exit(options(na.action = na.act))

   ### get position of arguments and put formula first

   mf <- match.call(expand.dots = FALSE)
   #m  <- match(c("formula", "data", "weights", "subset"), names(mf), 0L)
   m  <- match(c("formula", "data", "weights"), names(mf), 0L)
   mf <- mf[c(1L, m)]

   mf$drop.unused.levels <- TRUE
   formula <- Formula::as.Formula(formula) ### make sure that formula is really a multi-part formula

   if (length(formula)[2] < 2L)
      stop("Right-hand side of formula must specify both a grouping and a study factor (i.e., ~ group | study).")

   mf$formula <- formula
   mf[[1L]]   <- as.name("model.frame")
   mf         <- eval(mf, parent.frame()) ### create model frame (subsetting is done automatically; but not currently implemented)

   ### extract sample sizes / frequencies

   weights <- model.weights(mf)

   ### extract the lhs (the "outcomes") of the data frame (may be one or two variables)

   lhs <- Formula::model.part(formula, data = mf, lhs = 1)

   ### extract the 1st part (i.e., before the |) of the rhs (group and possibly outcome factor)

   rhs1 <- Formula::model.part(formula, data = mf, rhs = 1)

   ### extract the 2nd part (i.e., after the |) of the rhs (study factor)

   study <- Formula::model.part(formula, data = mf, rhs = 2)

   if (length(study) != 1)
      stop("A single study factor must be specified.")

   if (!is.factor(study[[1]]))
      stop("Study variable must be a factor.")

   ### get study variable from the data.frame

   study <- study[[1]]

   if (anyNA(study))
      stop("Study factor must not contain NAs.")

   #########################################################################

   if (is.element(measure, c("RR","OR","RD","AS","PETO","PHI","YUQ","YUY","RTET","PBIT","OR2D","OR2DN","OR2DL"))) {

      if (is.null(weights))
         stop("Must specify the 'weights' argument.")

      if (length(lhs) != 1)
         stop("Left-hand side of formula must be a single outcome factor.")

      outcome <- lhs[[1]]

      if (!is.factor(outcome))
         stop("Left-hand side of formula must be a factor.")

      if(nlevels(outcome) != 2)
         stop("Outcome factor on left-hand side of formula should have two levels.")

      if (length(rhs1) != 1)
         stop("A single grouping factor must be specified.")

      if (!is.factor(rhs1[[1]]))
         stop("Grouping variable must be a factor.")

      group <- rhs1[[1]]

      if(nlevels(group) != 2)
         stop("Grouping factor should have two levels.")

      if (anyNA(group) || anyNA(outcome))
         stop("Grouping and outcome factors must not contain NAs.")

      ai <- weights[group == levels(group)[1] & outcome == levels(outcome)[1]]
      bi <- weights[group == levels(group)[1] & outcome == levels(outcome)[2]]
      ci <- weights[group == levels(group)[2] & outcome == levels(outcome)[1]]
      di <- weights[group == levels(group)[2] & outcome == levels(outcome)[2]]

      names(ai) <- mf$study[group == levels(group)[1] & outcome == levels(outcome)[1]] ### add study names to ai
      names(bi) <- mf$study[group == levels(group)[1] & outcome == levels(outcome)[2]] ### add study names to bi
      names(ci) <- mf$study[group == levels(group)[2] & outcome == levels(outcome)[1]] ### add study names to ci
      names(di) <- mf$study[group == levels(group)[2] & outcome == levels(outcome)[2]] ### add study names to di

      #return(cbind(ai,bi,ci,di))
      return(escalc(measure=measure, ai=ai, bi=bi, ci=ci, di=di, add=add, to=to, drop00=drop00, vtype=vtype, var.names=var.names, append="FALSE", digits=digits))

   }

   #########################################################################

   if (is.element(measure, c("IRR","IRD","IRSD"))) {

      if (length(lhs) != 2)
         stop("Left-hand side of formula must specify the number of events and the total person-time at risk (i.e., events/times ~).")

      events <- lhs[,1]
      times  <- lhs[,2]

      if (!is.vector(events) || !is.vector(times))
         stop("The events and person-time at risk variables should be vectors.")

      if (length(rhs1) != 1)
         stop("A single grouping factor must be specified.")

      if (!is.factor(rhs1[[1]]))
         stop("Grouping variable must be a factor.")

      group <- rhs1[[1]]

      if(nlevels(group) != 2)
         stop("Grouping factor should have two levels.")

      if (anyNA(group))
         stop("Grouping factor must not contain NAs.")

      x1i <- events[group == levels(group)[1]]
      x2i <- events[group == levels(group)[2]]
      t1i <- times[group == levels(group)[1]]
      t2i <- times[group == levels(group)[2]]

      names(x1i) <- mf$study[group == levels(group)[1]] ### add study names to x1i
      names(x2i) <- mf$study[group == levels(group)[2]] ### add study names to x2i

      #return(cbind(x1i,x2i,t1i,t2i))
      return(escalc(measure=measure, x1i=x1i, x2i=x2i, t1i=t1i, t2i=t2i, add=add, to=to, drop00=drop00, vtype=vtype, var.names=var.names, append="FALSE", digits=digits))

   }

   #########################################################################

   if (is.element(measure, c("MD","SMD","SMDH","ROM","RPB","RBIS","D2OR","D2ORN","D2ORL"))) {

      if (is.null(weights))
         stop("Must specify the 'weights' argument.")

      if (length(lhs) != 2)
         stop("Left-hand side of formula must specify the means and standard devations (i.e., means/sds ~).")

      means <- lhs[,1]
      sds   <- lhs[,2]

      if (!is.vector(means) || !is.vector(sds))
         stop("The mean and standard devation variables should be vectors.")

      if (length(rhs1) != 1)
         stop("A single grouping factor must be specified.")

      if (!is.factor(rhs1[[1]]))
         stop("Grouping variable must be a factor.")

      group <- rhs1[[1]]

      if(nlevels(group) != 2)
         stop("Grouping factor should have two levels.")

      if (anyNA(group))
         stop("Grouping factor must not contain NAs.")

      m1i  <- means[group == levels(group)[1]]
      m2i  <- means[group == levels(group)[2]]
      sd1i <- sds[group == levels(group)[1]]
      sd2i <- sds[group == levels(group)[2]]
      n1i  <- weights[group == levels(group)[1]]
      n2i  <- weights[group == levels(group)[2]]

      names(m1i) <- mf$study[group == levels(group)[1]] ### add study names to m1i
      names(m2i) <- mf$study[group == levels(group)[2]] ### add study names to m2i

      #return(cbind(m1i, m2i, sd1i, sd2i, n1i, n2i))
      return(escalc(measure=measure, m1i=m1i, m2i=m2i, sd1i=sd1i, sd2i=sd2i, n1i=n1i, n2i=n2i, vtype=vtype, var.names=var.names, append="FALSE", digits=digits))

   }

   #########################################################################

   if (is.element(measure, c("COR","UCOR","ZCOR"))) {

      if (is.null(weights))
         stop("Must specify the 'weights' argument.")

      if (length(lhs) != 1)
         stop("Left-hand side of formula must specify the correlations (i.e., cors ~).")

      ri <- lhs[[1]]

      if (!is.vector(ri))
         stop("The variable specifying the correlation should be a vector.")

      #if (length(rhs1) != 1)
      #   stop("A single grouping factor must be specified.")
      #if (!is.factor(rhs1[[1]]))
      #   stop("Grouping variable must be a factor.")

      #group <- rhs1[[1]]

      #if(nlevels(group) != 1)
      #   stop("Grouping factor should have only one level.")

      #if (anyNA(group))
      #   stop("Grouping factor must not contain NAs.")

      ni <- weights

      names(ri) <- mf$study ### add study names to ri

      #return(cbind(ri, ni))
      return(escalc(measure=measure, ri=ri, ni=ni, vtype=vtype, var.names=var.names, append="FALSE", digits=digits))

   }

   #########################################################################

   if (is.element(measure, c("PR","PLN","PLO","PAS","PFT"))) {

      if (length(lhs) != 1)
         stop("Left-hand side of formula must be a single outcome factor.")

      outcome <- lhs[[1]]

      if (!is.factor(outcome))
         stop("Left-hand side of formula must be a factor.")

      if(nlevels(outcome) != 2)
         stop("Outcome factor on left-hand side of formula should have two levels.")

      #if (length(rhs1) != 1)
      #   stop("A single grouping factor must be specified.")
      #if (!is.factor(rhs1[[1]]))
      #   stop("Grouping variable must be a factor.")

      #group <- rhs1[[1]]

      #if(nlevels(group) != 1)
      #   stop("Grouping factor should have only one level.")

      #if (anyNA(group))
      #   stop("Grouping factor must not contain NAs.")

      if (anyNA(outcome))
         stop("Outcome factor must not contain NAs.")

      xi <- weights[outcome == levels(outcome)[1]]
      mi <- weights[outcome == levels(outcome)[2]]

      names(xi) <- mf$study[outcome == levels(outcome)[1]] ### add study names to xi
      names(mi) <- mf$study[outcome == levels(outcome)[2]] ### add study names to mi

      #return(cbind(xi,mi))
      return(escalc(measure=measure, xi=xi, mi=mi, add=add, to=to, vtype=vtype, var.names=var.names, append="FALSE", digits=digits))

   }

   #########################################################################

   if (is.element(measure, c("IR","IRLN","IRS","IRFT"))) {

      if (length(lhs) != 2)
         stop("Left-hand side of formula must specify the number of cases and the total person-time at risk (i.e., cases/times ~).")

      events <- lhs[,1]
      times  <- lhs[,2]

      if (!is.vector(events) || !is.vector(times))
         stop("The events and person-time at risk variables should be vectors.")

      #if (length(rhs1) != 1)
      #   stop("A single grouping factor must be specified.")
      #if (!is.factor(rhs1[[1]]))
      #   stop("Grouping variable must be a factor.")

      #group <- rhs1[[1]]

      #if(nlevels(group) != 1)
      #   stop("Grouping factor should have only one level.")

      #if (anyNA(group))
      #   stop("Grouping factor must not contain NAs.")

      xi <- events
      ti <- times

      names(xi) <- mf$study ### add study names to xi

      #return(cbind(xi,ti))
      return(escalc(measure=measure, xi=xi, ti=ti, add=add, to=to, vtype=vtype, var.names=var.names, append="FALSE", digits=digits))

   }

   #########################################################################

   if (is.element(measure, c("MN","MNLN"))) {

      if (is.null(weights))
         stop("Must specify the 'weights' argument.")

      if (length(lhs) != 2)
         stop("Left-hand side of formula must specify the means and standard devations (i.e., means/sds ~).")

      means <- lhs[,1]
      sds   <- lhs[,2]

      if (!is.vector(means) || !is.vector(sds))
         stop("The mean and standard devation variables should be vectors.")

      mi   <- means
      sdi  <- sds
      ni   <- weights

      names(mi) <- mf$study ### add study names to mi

      #return(cbind(mi, sdi, ni))
      return(escalc(measure=measure, mi=mi, sdi=sdi, ni=ni, vtype=vtype, var.names=var.names, append="FALSE", digits=digits))

   }

   #########################################################################

   #if (is.element(measure, c("MC","SMCC","SMCR","SMCRH","ROMC"))) {
   #}

   #########################################################################

   #if (is.element(measure, c("MPRD","MPRR","MPOR","MPORC","MPPETO"))) {
   #}

   #########################################################################

   if (is.element(measure, c("ARAW","AHW","ABT"))) {

      if (is.null(weights))
         stop("Must specify the 'weights' argument.")

      if (length(lhs) != 2)
         stop("Left-hand side of formula must specify the alpha values and number of items (i.e., alphas/items ~).")

      alphas <- lhs[,1]
      items  <- lhs[,2]

      if (!is.vector(alphas) || !is.vector(items))
         stop("The alpha and item variables should be vectors.")

      ai <- alphas
      mi <- items
      ni <- weights

      names(ai) <- mf$study ### add study names to ai

      #return(cbind(ai, mi, ni))
      return(escalc(measure=measure, ai=ai, mi=mi, ni=ni, vtype=vtype, var.names=var.names, append="FALSE", digits=digits))

   }

   #########################################################################

}
