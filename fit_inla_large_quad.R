fit_inla <- function(sim.list, scenario.str = scenario, mesh_from = c("pp", "quad")) {
  mesh_from = match.arg(mesh_from) # build mesh from point pattern or quadrature
  start.time <- Sys.time()
  # get the required elements from the simulation information
  dat <- sim.list$dframe
  lambda.truth <- exp(dat$loglamb[dat$pres == 0])
  
  # create a mesh
  quad <- expand.grid(list(x= seq(0, 100, by = 1), y= seq(0, 100, by = 1)))
  mesh <- inla.mesh.create(loc = as.matrix(quad[, c("x", "y")]))

  # create the SPDE
  spde.matern <- inla.spde2.pcmatern(mesh,
                                     prior.sigma = c(10, 0.01),
                                     prior.range = c(5, 0.01))
  # set up the INLA-specific objects (as in Illian paper)
  y.pp <- rep(0:1, c(sum(dat$pres == 0), sum(dat$pres == 1)))
  e.pp <- c(dat$quad.size[dat$pres == 0], dat$quad.size[dat$pres == 1])
  mesh2pp <- inla.spde.make.A(mesh, as.matrix(dat[dat$pres == 1, c("x", "y")]))
  mesh2quad <- inla.spde.make.A(mesh, as.matrix(dat[dat$pres == 0, c("x", "y")]))
  Amat <- rbind(mesh2quad, mesh2pp)
  stk.pp <- inla.stack(
    data = list(y = y.pp, e = e.pp),
    A = list(1, Amat),
    effects = list(list(Intercept = 1, Slope = c(dat$X[dat$pres == 0], dat$X[dat$pres == 1])),
                   list(i = 1:mesh$n)),
    tag = 'pp')
  # fit the model with error catching
  #   note: this is without computing the fitted values which yields are far longer comp. time)
  try(assign("inla.fit", inla(y ~ 0 + Intercept + Slope + f(i, model = spde.matern),
                   family = 'poisson', data = inla.stack.data(stk.pp),
                   control.predictor = list(A = inla.stack.A(stk.pp)),
                   E = inla.stack.data(stk.pp)$e)),
             silent = T)
  if (!exists("inla.fit")) {
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
                          k = mesh$n,
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
                          model = "inla" 
    )
  } else {
    # fit the model
    model.fit <- inla(y ~ 0 + Intercept + Slope + f(i, model = spde.matern),
                      family = 'poisson', data = inla.stack.data(stk.pp),
                      control.predictor = list(A = inla.stack.A(stk.pp), compute = T),
                      E = inla.stack.data(stk.pp)$e)
    # extract the fitted values
    dnames <- unlist(lapply(strsplit(dimnames(model.fit$summary.linear.predictor)[[1]], ".", fixed = T), function(x){x[1]}))
    fitted.mat <- model.fit$summary.linear.predictor[dnames == "APredictor", ]
    # the following are mostly legacy code no longer used (to extract more information)
    tmp.loglik <- as.numeric(inla.fit$mlik[1, 1])
    # tmp.ll1 <- NA
    # tmp.ll2 <- NA
    # tmp.ll3 <- NA
    tmp.b0 <- inla.fit$summary.fixed[1, c("mean", "sd")]
    tmp.b1 <- inla.fit$summary.fixed[2, c("mean", "sd")]
    tmp.betas <- rbind(tmp.b0, tmp.b1)
    colnames(tmp.betas) <- c("Estimate", "Std.Error")
    # tmp.mus <- inla.fit$summary.random$i$mean
    # tmp.post_vars <- inla.fit$summary.random$i$sd^2
    # names(tmp.post_vars) <- rep("sig2", length(tmp.post_vars))
    # tmp.prior_var <- inla.fit$summary.hyperpar[2, "mean"]
    # tmp.range_par <- inla.fit$summary.hyperpar[1, "mean"]
    # tmp.aic <- NA
    # tmp.nodes <- inla.fit$size.random[[1]]$N
    
    # get 95% CI
    tmp.pars <- as.data.frame(cbind(
      tmp.betas[ , 1],
      tmp.betas[ , 1] - qnorm(1 - (0.05 / 2)) * tmp.betas[ , 2],
      tmp.betas[ , 1] + qnorm(1 - (0.05 / 2)) * tmp.betas[ , 2],
      inla.fit$summary.fixed[ , c("0.5quant", "0.025quant", "0.975quant")]
    ))
    dimnames(tmp.pars) <- list(c("Int", "Slope"), c("est", "lower", "upper", "median", "cred lower", "cred upper"))
    tmp.coverage <- as.logical((tmp.pars["Slope", "lower"] <= sim.list$info['Slope']) & (tmp.pars["Slope", "upper"] >= sim.list$info['Slope']))
    tmp.coverage_cred <- as.logical((tmp.pars["Slope", "cred lower"] <= sim.list$info['Slope']) & (tmp.pars["Slope", "cred upper"] >= sim.list$info['Slope']))
    # get the squared error of the estimate
    tmp.sqerr <- as.numeric((sim.list$info['Slope'] - tmp.pars["Slope", "est"])^2)
    tmp.sqerr_median <- as.numeric((sim.list$info['Slope'] - tmp.pars["Slope", "median"])^2)
    # get fitted loglambdas
    tmp.pred.quad <- as.vector(fitted.mat[1:sum(dat$pres == 0), "mean"])
    # tmp.pred.pres <- as.vector(fitted.mat[(sum(dat$pres == 0)+1):nrow(fitted.mat), "mean"])
    # get the Kullback-Leibler Divergence
    tmp.lambda <- exp(as.vector(tmp.pred.quad))
    tmp.dkl <- as.numeric(
      dat$quad.size[dat$pres == 0] %*% (lambda.truth * log(lambda.truth / tmp.lambda))) -
      as.numeric(dat$quad.size[dat$pres == 0] %*% (lambda.truth - tmp.lambda))

    # collate results
    end.time <- Sys.time()
    sim.res <- data.frame(ftime = as.numeric(inla.fit$cpu.used[4]),
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
                          k = mesh$n,
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
                          model = "inla" 
    )
  }
  # save the result
  save(list = "sim.res", file = paste0("Results/", scenario.str, "/res_inla_", mesh_from, "_sim_", job, ".RDATA"))
}