fit_ipp <- function(sim.list, scenario.str = scenario) {
  start.time <- Sys.time()
  # get the required elements from the simulation information
  dat <- sim.list$dframe
  lambda.truth <- exp(dat$loglamb[dat$pres == 0])
  # set up TMB data list
  dat.list <- list(
    X_pres = as.matrix(cbind(rep(1, sum(dat$pres == 1)), dat$X[dat$pres == 1])),
    X_quad = as.matrix(cbind(rep(1, sum(dat$pres == 0)), dat$X[dat$pres == 0])),
    quad_size = dat$quad.size[dat$pres == 0]
  )
  # clear space
  rm(dat)
  gc()
  # set the starting parameters
  start.pars <- list(b0 = 0, b1 = 0)
  # load the c++ likelihood
  dyn.load(dynlib("ipp"))
  # create the TMB likelihood function
  obj <- MakeADFun(data = dat.list, parameters = start.pars, DLL = "ipp", silent = T)
  # increase the maximum number of iterations
  obj$control = list(maxit = 10000)
  # fit the model with error catching
  try(assign("fit.time", system.time(assign("res", do.call("optim", obj)))), silent = T)
  # check if the model fit and fill results with NA if not
  if (!exists("res")) {
    # return(c(ftime = NA, sqerr = NA, coverage = NA, dkl = NA, ll = NA, slope_est = NA, ciw_lower = NA, ciw_upper = NA, k = 0, flagfit = TRUE))
    # collate results
    sim.res <- data.frame(ftime = NA,
                          sqerr = NA,
                          coverage = NA,
                          sqerr_median = NA,
                          coverage_cred = NA,
                          dkl = NA,
                          ll = NA,
                          slope_est = NA,
                          ciw = NA,
                          ciw_lower = NA,
                          ciw_upper = NA,
                          k = 0,
                          flagfit = TRUE,
                          expected_n = as.numeric(sim.list$info["Expected_N"]),
                          n = as.numeric(sim.list$info["N"]),
                          covar_range = as.numeric(sim.list$info["covariate_range"]),
                          latent_range = as.numeric(sim.list$info["latent.range"]),
                          latent_variance = as.numeric(sim.list$info["latent.variance"]),
                          intercept = as.numeric(sim.list$info["Intercept"]),
                          slope = as.numeric(sim.list$info["Slope"]),
                          sim = as.numeric(sim.list$info["rseed"]),
                          scenario = scenario.str,
                          atime = NA,
                          model = "ipp" 
                          )
  } else {
    # extract TMB results info with SE
    try(assign("res.estimates", summary(sdreport(obj))))
    # need additional check for failure here when the gradient evaluates to smaller dim in LPRK
    if (!exists("res.estimates")) {
      # the following are mostly legacy code no longer used (to extract more information)
      tmp.loglik <- -res$value
      # tmp.ll1 <- obj$report()$ll1
      # tmp.ll2 <- obj$report()$ll2
      # tmp.ll3 <- NA
      tmp.b0 <- c(Estimate = as.numeric(res$par["b0"]), `Std. Error` = NA)
      tmp.b1 <- c(Estimate = as.numeric(res$par["b1"]), `Std. Error` = NA)
      tmp.betas <- rbind(tmp.b0, tmp.b1)
      # tmp.mus <- NA
      # tmp.post_vars <- NA
      # names(tmp.post_vars) <- rep("sig2", length(tmp.post_vars))
      # tmp.prior_var <- NA
      # tmp.range_par <- NA
      # tmp.aic <- (-2 * tmp.loglik) + (2 * nrow(res.estimates))
      # tmp.nodes <- NA
      tmp.flag_fit <- TRUE # flag this as a potential poor fit
    } else {
      # the following are mostly legacy code no longer used (to extract more information)
      tmp.loglik <- -res$value
      # tmp.ll1 <- obj$report()$ll1
      # tmp.ll2 <- obj$report()$ll2
      # tmp.ll3 <- NA
      tmp.b0 <- res.estimates[rownames(res.estimates) == "b0", ]
      tmp.b1 <- res.estimates[rownames(res.estimates) == "b1", ]
      tmp.betas <- rbind(tmp.b0, tmp.b1)
      # tmp.mus <- NA
      # tmp.post_vars <- NA
      # names(tmp.post_vars) <- rep("sig2", length(tmp.post_vars))
      # tmp.prior_var <- NA
      # tmp.range_par <- NA
      # tmp.aic <- (-2 * tmp.loglik) + (2 * nrow(res.estimates))
      # tmp.nodes <- NA
      if (is.na(tmp.b1["Std. Error"])) {
        tmp.flag_fit <- TRUE # flag this as a potential poor fit
      } else if (res$convergence != 0) {
        tmp.flag_fit <- TRUE # flag this as a potential poor fit
      } else {
        tmp.flag_fit <- FALSE # do not flag this as a potential poor fit
      }
    }
    # get 95% CI
    tmp.pars <- as.data.frame(cbind(
      tmp.betas[ , 1],
      tmp.betas[ , 1] - qnorm(1 - (0.05 / 2)) * tmp.betas[ , 2],
      tmp.betas[ , 1] + qnorm(1 - (0.05 / 2)) * tmp.betas[ , 2]
    ))
    dimnames(tmp.pars) <- list(c("Int", "Slope"), c("est", "lower", "upper"))
    tmp.coverage <- as.logical(as.numeric((tmp.pars["Slope", "lower"] <= sim.list$info['Slope']) & (tmp.pars["Slope", "upper"] >= sim.list$info['Slope'])))
    # get the squared error of the estimate
    tmp.sqerr <- as.numeric((sim.list$info['Slope'] - tmp.pars["Slope", "est"])^2)
    # get fitted loglambdas
    tmp.pred.quad <- dat.list$X_quad %*% tmp.betas[ , 1]
    tmp.pred.pres <- dat.list$X_pres %*% tmp.betas[ , 1]
    # get the Kullback-Leibler Divergence
    tmp.lambda <- exp(as.vector(tmp.pred.quad))
    tmp.dkl <- as.numeric(
      dat.list$quad_size %*% (lambda.truth * log(lambda.truth / tmp.lambda))) -
      as.numeric(dat.list$quad_size %*% (lambda.truth - tmp.lambda))

    # collate results
    end.time <- Sys.time()
    sim.res <- data.frame(ftime = as.numeric(fit.time[3]),
                          sqerr = tmp.sqerr,
                          coverage = tmp.coverage,
                          sqerr_median = NA,
                          coverage_cred = NA,
                          dkl = tmp.dkl,
                          ll = tmp.loglik,
                          slope_est = tmp.pars["Slope", "est"],
                          ciw = tmp.pars["Slope", "est"] - tmp.pars["Slope", "lower"],
                          ciw_lower = tmp.pars["Slope", "est"] - tmp.pars["Slope", "lower"],
                          ciw_upper = tmp.pars["Slope", "upper"] - tmp.pars["Slope", "est"],
                          k = 0,
                          flagfit = tmp.flag_fit,
                          expected_n = as.numeric(sim.list$info["Expected_N"]),
                          n = as.numeric(sim.list$info["N"]),
                          covar_range = as.numeric(sim.list$info["covariate_range"]),
                          latent_range = as.numeric(sim.list$info["latent.range"]),
                          latent_variance = as.numeric(sim.list$info["latent.variance"]),
                          intercept = as.numeric(sim.list$info["Intercept"]),
                          slope = as.numeric(sim.list$info["Slope"]),
                          sim = as.numeric(sim.list$info["rseed"]),
                          scenario = scenario.str,
                          atime = as.numeric(difftime(end.time, start.time, units = "secs")),
                          model = "ipp" 
    )
  }
  # save the result
  save(list = "sim.res", file = paste0("Results/", scenario.str, "/res_ipp_sim_", job, ".RDATA"))
}