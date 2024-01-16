fit_vark <- function(sim.list, bfs, scenario.str = scenario, radius_type = "limiting") {
  start.time <- Sys.time()
  # get the required elements from the simulation information
  dat <- sim.list$dframe
  lambda.truth <- exp(dat$loglamb[dat$pres == 0])
  # get function to set up the basis node matrix
  source("make_basis.R")
  if (bfs == "search") {
    # get the search function for the maximum likelihood basis configuration
    source("basis_search.R")
    search.res <- basis_search(sim.list, approx_type = "VA")
    sqrt_bf_nodes <- search.res[which.max(search.res[,2]), 1]
  } else {
    sqrt_bf_nodes <- as.numeric(strsplit(bfs, "x", fixed = T)[[1]][1])
  }
  bfns <- make_basis(sqrt_bf_nodes, data = dat, radius.type = radius_type)
  # bfns <- switch(bfs,
  #                 `7x7` = make_basis(7, data = dat, radius.type = "limiting"),
  #                 `10x10` = make_basis(10, data = dat, radius.type = "limiting"),
  #                 `14x14` = make_basis(14, data = dat, radius.type = "limiting"),
  #                 `25x25` = make_basis(25, data = dat, radius.type = "limiting")
  # )
  nodes <- bfns[ , c("x", "y")]
  Radius <- unique(bfns$scale)

  if (!library(fields, logical.return = T)) {
    stop("'fields' package needed for to set up basis function matrix")
  }
  # calculate the euclidean distances from presence points to basis nodes
  Dist.Mat.pp <- fields::rdist(dat[dat$pres == 1 , c("x", "y")], nodes)
  # calculate the euclidean distances from quadrature points to basis nodes
  Dist.Mat.quad <- fields::rdist(dat[dat$pres == 0 , c("x", "y")], nodes)
  # set distances beyond the radii to zero
  Dist.Mat.pp[Dist.Mat.pp > Radius] <- 0
  Dist.Mat.quad[Dist.Mat.quad > Radius] <- 0
  # convert to sparse matrices
  Dist.Mat.pp <- as(Dist.Mat.pp, "sparseMatrix")
  Dist.Mat.quad <- as(Dist.Mat.quad, "sparseMatrix")
  # set up the sparse bi-square basis function matrix (logical acts as sparse matrix with 1's)
  Z.pp <- ((Dist.Mat.pp != 0) - (Dist.Mat.pp / Radius)^2)^2
  Z.quad <- ((Dist.Mat.quad != 0) - (Dist.Mat.quad / Radius)^2)^2
  # clear space
  rm(Dist.Mat.pp, Dist.Mat.quad)
  gc()

  # set up TMB data list
  dat.list <- list(
    X_pres = as.matrix(cbind(rep(1, sum(dat$pres == 1)), dat$X[dat$pres == 1])),
    Z_pres = Z.pp,
    X_quad = as.matrix(cbind(rep(1, sum(dat$pres == 0)), dat$X[dat$pres == 0])),
    Z_quad = Z.quad,
    quad_size = dat$quad.size[dat$pres == 0]
  )
  # clear space
  rm(dat, Z.pp, Z.quad)
  gc()
  
  ## get ipp warm starts ##
  # load the c++ likelihood
  dyn.load(dynlib("ipp"))
  # create the TMB likelihood function
  obj0 <- MakeADFun(data = dat.list, parameters = list(b0 = 0, b1 = 0), DLL = "ipp", silent = T)
  ipp <- do.call("optim", obj0)
  #####
  
  # set the starting parameters
  # start.pars <- list(b0 = 0, b1 = 0, mu = rep(0, ncol(dat.list$Z_pres)), logsig2 = rep(0, ncol(dat.list$Z_pres)))
  start.pars <- list(b0 = ipp$par[1], b1 = ipp$par[2], mu = rep(0, ncol(dat.list$Z_pres)), logsig2 = rep(0, ncol(dat.list$Z_pres)))
  # load the c++ likelihood
  dyn.load(dynlib("vark"))
  # create the TMB likelihood function
  obj <- MakeADFun(data = dat.list, parameters = start.pars, DLL = "vark", silent = T)
  # increase the maximum number of iterations
  obj$control = list(maxit = 10000)
  # fit the model with error catching
  try(assign("fit.time", system.time(assign("res", do.call("optim", obj)))), silent = T)
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
                 k = nrow(nodes),
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
                 model = paste0("vark_", bfs)
                 )
    # return(c(ftime = NA, sqerr = NA, coverage = NA, dkl = NA, ll = NA, slope_est = NA, ciw_lower = NA, ciw_upper = NA, k = nrow(nodes), flagfit = TRUE))
  } else {
    # extract TMB results info with SE
    try(assign("res.estimates", summary(sdreport(obj))))
    # need additional check for failure here when the gradient evaluates to smaller dim in LPRK
    if (!exists("res.estimates")) {
      # the following are mostly legacy code no longer used (to extract more information)
      tmp.loglik <- -res$value
      # tmp.ll1 <- obj$report()$ll1
      # tmp.ll2 <- obj$report()$ll2
      # tmp.ll3 <- obj$report()$ll3
      tmp.b0 <- c(Estimate = as.numeric(res$par["b0"]), `Std. Error` = NA)
      tmp.b1 <- c(Estimate = as.numeric(res$par["b1"]), `Std. Error` = NA)
      # When res.estimates fails, it means that LP spiraled off to near Inf logLik
      # This gives a horrible point estimate so it may as well be a fit failure, i.e. NA
      # tmp.b0 <- c(Estimate = NA, `Std. Error` = NA)
      # tmp.b1 <- c(Estimate = NA, `Std. Error` = NA)
      tmp.betas <- rbind(tmp.b0, tmp.b1)
      tmp.mus <- res$par[names(res$par) == "mu"]
      # tmp.post_vars <- res.estimates[rownames(res.estimates) == "sig2", 1]
      # names(tmp.post_vars) <- rep("sig2", length(tmp.post_vars))
      # tmp.prior_var <- res.estimates[rownames(res.estimates) == "del2", 1]
      # tmp.range_par <- Radius
      # tmp.aic <- (-2 * tmp.loglik) + 2 * ((2 * length(tmp.mus)) + nrow(tmp.betas))
      # tmp.nodes <- length(tmp.mus)
      tmp.flag_fit <- TRUE # flag this as a potential poor fit
    } else {
      # the following are mostly legacy code no longer used (to extract more information)
      tmp.loglik <- -res$value
      # tmp.ll1 <- obj$report()$ll1
      # tmp.ll2 <- obj$report()$ll2
      # tmp.ll3 <- obj$report()$ll3
      tmp.b0 <- res.estimates[rownames(res.estimates) == "b0", ]
      tmp.b1 <- res.estimates[rownames(res.estimates) == "b1", ]
      tmp.betas <- rbind(tmp.b0, tmp.b1)
      tmp.mus <- res.estimates[rownames(res.estimates) == "mu", 1]
      # tmp.post_vars <- res.estimates[rownames(res.estimates) == "sig2", 1]
      # names(tmp.post_vars) <- rep("sig2", length(tmp.post_vars))
      # tmp.prior_var <- res.estimates[rownames(res.estimates) == "del2", 1]
      # tmp.range_par <- Radius
      # tmp.aic <- (-2 * tmp.loglik) + 2 * ((2 * length(tmp.mus)) + nrow(tmp.betas))
      # tmp.nodes <- length(tmp.mus)
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
    tmp.coverage <- as.logical((tmp.pars["Slope", "lower"] <= sim.list$info['Slope']) & (tmp.pars["Slope", "upper"] >= sim.list$info['Slope']))
    # get the squared error of the estimate
    tmp.sqerr <- as.numeric((sim.list$info['Slope'] - tmp.pars["Slope", "est"])^2)
    # get fitted loglambdas
    tmp.pred.quad <- as.vector(dat.list$X_quad %*% tmp.betas[ , 1] + dat.list$Z_quad %*% tmp.mus)
    # tmp.pred.pres <- as.vector(dat.list$X_pres %*% tmp.betas[ , 1] + dat.list$Z_pres %*% tmp.mus)
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
                       k = nrow(nodes),
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
                       model = paste0("vark_", bfs) 
                       )
  }
  # save the result
  save(list = "sim.res", file = paste0("Results/", scenario.str, "/res_vark_", bfs, "_sim_", job, ".RDATA"))
}