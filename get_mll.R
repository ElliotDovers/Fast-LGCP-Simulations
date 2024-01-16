get_mll <- function(sim.list, sqrt_nodes = 7, approx_type = c("VA", "Lp")) {
  # function for making the FRK basis configuration
  source("make_basis.R")
  # requires the library fields
  if (!library(fields, logical.return = T)) {
    stop("'fields' package needed for to set up basis function matrix")
  }
  approx_type <- match.arg(approx_type)
  
  dat <- sim.list$dframe
  bfns <- make_basis(sqrt_nodes, data = dat, radius.type = "limiting")
  nodes <- bfns[ , c("x", "y")]
  Radius <- unique(bfns$scale)
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
  rm(Z.pp, Z.quad)
  gc()
  
  ## get ipp warm starts ##
  # load the c++ likelihood
  dyn.load(dynlib("ipp"))
  # create the TMB likelihood function
  obj0 <- MakeADFun(data = dat.list, parameters = list(b0 = 0, b1 = 0), DLL = "ipp", silent = T)
  ipp <- do.call("optim", obj0)
  #####
  
  # change approach depending on approx_type
  if (approx_type == "VA") {
    # set the starting parameters
    start.pars <- list(b0 = ipp$par[1], b1 = ipp$par[2], mu = rep(0, ncol(dat.list$Z_pres)), logsig2 = rep(0, ncol(dat.list$Z_pres)))
    # load the c++ likelihood
    dyn.load(dynlib("vark"))
    # create the TMB likelihood function
    obj <- MakeADFun(data = dat.list, parameters = start.pars, DLL = "vark", silent = T)
  } else {
    # set the starting parameters
    start.pars <- list(b0 = ipp$par[1], b1 = ipp$par[2], u = rep(0, ncol(dat.list$Z_pres)), logdel = 0)
    # load the c++ likelihood
    dyn.load(dynlib("lprk"))
    # create the TMB likelihood function
    obj <- MakeADFun(data = dat.list, parameters = start.pars, random = "u", DLL = "lprk", silent = T)
  }
  # increase the maximum number of iterations
  obj$control = list(maxit = 10000)
  try(assign("fit.time", system.time(assign("res", do.call("optim", obj)))), silent = T)
  # res <- do.call("optim", obj)
  if (!exists("res")) {
    return(c(sqrt_nodes = sqrt_nodes, mll = NA))
  } else {
    return(c(sqrt_nodes = sqrt_nodes, mll = -res$value))
  }
}