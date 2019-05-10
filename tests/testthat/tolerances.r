.tol <- c(est   = .01, # effect size estimates
          coef  = .01, # model coefficients
          pred  = .01, # predicted values, BLUPs, also residuals
          se    = .01, # standard errors
          test  = .01, # test statistics, standardized residuals
          pval  = .01, # p-values
          ci    = .01, # confidence/credibility/prediction interval bounds, CI for effects
          var   = .01, # variance components (and CIs thereof), also if sqrt(), var-cov matrices, sampling variances
          cor   = .01, # correlations, ICCs
          cov   = .01, # covariances
          sevar = .01, # SEs of variance components
          fit   = .01, # fit statistics
          r2    = .01, # R^2 type values
          het   = .01, # heterogeneity statistics (and CIs thereof)
          inf   = .01, # influence statistics, hat values
          den   = .01, # density
          misc  = .01) # miscellanoue, mix of values

# to quickly set all tolerances to a common value
#.tol[1:length(.tol)] <- .001

# note to self: search for "]]/10 or "]]*10 to find adjusted tolerances in tests
# some hardcoded tolerances; search for: tolerance=.0
