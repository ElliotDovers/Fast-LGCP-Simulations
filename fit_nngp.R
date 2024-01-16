fit_nngp <- function(sim.list, n_samples = 1000, n_neighbours = 15, burnin = 500, thin_by = 1,
                     starting = list("phi"=1/10, "sigma.sq"=1),
                     tuning = list("phi"=0.5, "sigma.sq"=0.5),
                     priors = list("phi.Unif"=c(1/50, 1/5), "sigma.sq.IG"=c(1, 1)),
                     cov.fun = "exponential", scenario.str = scenario,
                     IWLR.wt = 4600
                     ) {
  start.time <- Sys.time()
  # get the required elements from the simulation information
  dat <- sim.list$dframe
  lambda.truth <- exp(dat$loglamb[dat$pres == 0])
  
  # set the point weights as per Fithian and Hastie (2013)
  p.wt = (IWLR.wt)^(1 - dat$pres)
  
  # fit the model with error catching
  try(assign("res",spNNGP(pres ~ X, data = dat, coords=c("x", "y"), starting=starting, method="latent", # requires latent for binomial
               family = "binomial", weights = p.wt, n.neighbors=n_neighbours,
               tuning=tuning, priors=priors, cov.model = cov.fun,
               n.samples=n_samples, n.omp.threads=1, fit.rep = T, verbose = F, sub.sample = list(start = burnin, end = n_samples, thin = thin_by))))
  
  # check if the model fit and fill results with NA if not
  if (!exists("res")) {
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
                          k = NA,
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
                          model = "nngp" 
    )
  } else {
    # using credible intervals
    tmp.bayes <- t(apply(res$p.beta.samples[seq(as.integer(res$sub.sample$start), as.integer(res$sub.sample$end), by=as.integer(res$sub.sample$thin)), ], 2, function(x){c(quantile(x, probs = c(0.5, 0.025, 0.975)))}))
    # using Wald intervals
    tmp.betas <- t(apply(res$p.beta.samples[seq(as.integer(res$sub.sample$start), as.integer(res$sub.sample$end), by=as.integer(res$sub.sample$thin)), ], 2, function(x){c(mean(x), sd(x))}))
    
    # the following are mostly legacy code no longer used (to extract more information)
    tmp.loglik <- NA
    
    # get 95% CI
    tmp.pars <- as.data.frame(cbind(
      tmp.betas[ , 1],
      tmp.betas[ , 1] - qnorm(1 - (0.05 / 2)) * tmp.betas[ , 2],
      tmp.betas[ , 1] + qnorm(1 - (0.05 / 2)) * tmp.betas[ , 2],
      tmp.bayes
    ))
    dimnames(tmp.pars) <- list(c("Int", "Slope"), c("est", "lower", "upper", "median", "cred lower", "cred upper"))
    
    tmp.coverage <- as.logical((tmp.pars["Slope", "lower"] <= sim.list$info['Slope']) & (tmp.pars["Slope", "upper"] >= sim.list$info['Slope']))
    tmp.coverage_cred <- as.logical((tmp.pars["Slope", "cred lower"] <= sim.list$info['Slope']) & (tmp.pars["Slope", "cred upper"] >= sim.list$info['Slope']))
    # get the squared error of the estimate
    tmp.sqerr <- as.numeric((sim.list$info['Slope'] - tmp.pars["Slope", "est"])^2)
    tmp.sqerr_median <- as.numeric((sim.list$info['Slope'] - tmp.pars["Slope", "median"])^2)
    # get fitted loglambdas
    tmp.pred.quad <- as.vector(res$y.hat.quants[dat$pres == 0, "50%"])
    # we need to perform the intercept correction
    delta_hat <- log(sum(dat$pres)) - log(as.numeric(dat$quad.size[dat$pres == 0] %*% exp(as.vector(tmp.pred.quad))))
    tmp.lambda <- exp(tmp.pred.quad + delta_hat)
    # get the Kullback-Leibler Divergence
    tmp.lambda <- exp(tmp.pred.quad + delta_hat)
    # tmp.lambda <- exp(as.vector(tmp.pred.quad))
    tmp.dkl <- as.numeric(
      dat$quad.size[dat$pres == 0] %*% (lambda.truth * log(lambda.truth / tmp.lambda))) -
      as.numeric(dat$quad.size[dat$pres == 0] %*% (lambda.truth - tmp.lambda))
    
    # collate results
    end.time <- Sys.time()
    sim.res <- data.frame(ftime = as.numeric(res$run.time[3]),
                          sqerr = tmp.sqerr,
                          coverage = tmp.coverage,
                          sqerr_median = tmp.sqerr_median,
                          coverage_cred = tmp.coverage_cred,
                          dkl = tmp.dkl,
                          ll = tmp.loglik,
                          slope_est = tmp.pars["Slope", "est"],
                          ciw = tmp.pars["Slope", "est"] - tmp.pars["Slope", "lower"],
                          ciw_lower = tmp.pars["Slope", "est"] - tmp.pars["Slope", "cred lower"],
                          ciw_upper = tmp.pars["Slope", "cred upper"] - tmp.pars["Slope", "est"],
                          k = NA,
                          flagfit = FALSE,
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
                          model = "nngp" 
    )
  }
  # save the result
  save(list = "sim.res", file = paste0("Results/", scenario.str, "/res_nngp_sim_", job, ".RDATA"))
}