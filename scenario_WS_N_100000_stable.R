# Control simulation #

############
# Preamble #
############

home.wd <- getwd()

# library(lgcp)
library(spNNGP)
library(TMB)
library(INLA)

#####################
# Simulate the Data #
#####################

source("sim_lgcp_pp.R")
source("fit_ipp.R")
source("fit_vark.R")
source("fit_lprk.R")
source("fit_inla.R")
# source("fit_mala.R")
source("fit_nngp.R")

job = as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))

########### Elements To Change ##############
Int = 1.675054 # Expected N = 100000
cov.par = 5 # Wiggly
range.par = 30 # Smooth
Slope = 1.25
cor.fn = "stable"
############################################

### This defines scenario:
scenario <- "WS_N_100000_stable"

# Data Simulation
sim.list <- sim_lgcp_pp(Intercept = Int, Slope = Slope, latent.range = range.par, rseed = job,
                        cov.function = cor.fn, covariate_range = cov.par)

############################
# Fit the competing models #
############################

# ## IPP:
# if (!paste0("res_ipp_sim_", job, ".RDATA") %in% list.files(paste0("Results/", scenario))) {
#   fit_ipp(sim.list)
# }
# 
# ## VARK/LPRK with 7x7 basis functions:
# # fit VARK (variational approx. with fixed rank kriging)
# if (!paste0("res_vark_7x7_sim_", job, ".RDATA") %in% list.files(paste0("Results/", scenario))) {
#   fit_vark(sim.list, bfs = "7x7")
# }
# # fit LPRK (Laplace approx. with fixed rank kriging)
# if (!paste0("res_lprk_7x7_sim_", job, ".RDATA") %in% list.files(paste0("Results/", scenario))) {
#   fit_lprk(sim.list, bfs = "7x7")
# }
# 
# ## VARK/LPRK with 14x14 basis functions:
# # fit VARK (variational approx. with fixed rank kriging)
# if (!paste0("res_vark_14x14_sim_", job, ".RDATA") %in% list.files(paste0("Results/", scenario))) {
#   fit_vark(sim.list, bfs = "14x14")
# }
# # fit LPRK (Laplace approx. with fixed rank kriging)
# if (!paste0("res_lprk_14x14_sim_", job, ".RDATA") %in% list.files(paste0("Results/", scenario))) {
#   fit_lprk(sim.list, bfs = "14x14")
# }
# 
# ## VARK/LPRK with 21x21 basis functions:
# # fit VARK (variational approx. with fixed rank kriging)
# if (!paste0("res_vark_21x21_sim_", job, ".RDATA") %in% list.files(paste0("Results/", scenario))) {
#   fit_vark(sim.list, bfs = "21x21")
# }
# # fit LPRK (Laplace approx. with fixed rank kriging)
# if (!paste0("res_lprk_21x21_sim_", job, ".RDATA") %in% list.files(paste0("Results/", scenario))) {
#   fit_lprk(sim.list, bfs = "21x21")
# }
# 
# ## VARK/LPRK with 28x28 basis functions:
# # fit VARK (variational approx. with fixed rank kriging)
# if (!paste0("res_vark_28x28_sim_", job, ".RDATA") %in% list.files(paste0("Results/", scenario))) {
#   fit_vark(sim.list, bfs = "28x28")
# }
# # fit LPRK (Laplace approx. with fixed rank kriging)
# if (!paste0("res_lprk_28x28_sim_", job, ".RDATA") %in% list.files(paste0("Results/", scenario))) {
#   fit_lprk(sim.list, bfs = "28x28")
# }
# 
# ## VARK/LPRK with 35x35 basis functions:
# # fit VARK (variational approx. with fixed rank kriging)
# if (!paste0("res_vark_35x35_sim_", job, ".RDATA") %in% list.files(paste0("Results/", scenario))) {
#   fit_vark(sim.list, bfs = "35x35")
# }
# # fit LPRK (Laplace approx. with fixed rank kriging)
# if (!paste0("res_lprk_35x35_sim_", job, ".RDATA") %in% list.files(paste0("Results/", scenario))) {
#   fit_lprk(sim.list, bfs = "35x35")
# }
# 
# ## VARK/LPRK with basis search:
# # fit VARK (variational approx. with fixed rank kriging)
# if (!paste0("res_vark_search_sim_", job, ".RDATA") %in% list.files(paste0("Results/", scenario))) {
#   fit_vark(sim.list, bfs = "search")
# }
# # fit LPRK (Laplace approx. with fixed rank kriging)
# if (!paste0("res_lprk_search_sim_", job, ".RDATA") %in% list.files(paste0("Results/", scenario))) {
#   fit_lprk(sim.list, bfs = "search")
# }
# 
# ## INLA:
# # with a mesh from the point pattern
# if (!paste0("res_inla_pp_sim_", job, ".RDATA") %in% list.files(paste0("Results/", scenario))) {
#   fit_inla(sim.list)
# }
# with a mesh from the quadrature/background points (since n > 10,000) tends crash
if (!paste0("res_inla_quad_sim_", job, ".RDATA") %in% list.files(paste0("Results/", scenario))) {
  fit_inla(sim.list, mesh_from = "quad")
}

## spNNGP:
if (!paste0("res_nngp_sim_", job, ".RDATA") %in% list.files(paste0("Results/", scenario))) {
  fit_nngp(sim.list,
           n_samples = 30000, n_neighbours = 15, burnin = 20000, thin_by = 10,
           starting = list("phi"=1/10, "sigma.sq"=1),
           tuning = list("phi"=0.09, "sigma.sq"=0.5),
           priors = list("phi.Unif"=c(1/50, 1/5), "sigma.sq.IG"=c(1, 2)),
           cov.fun = "exponential", IWLR.wt = 1
  )
}

# # MALA via lgcp package:
# if (!paste0("res_mala_sim_", job, ".RDATA") %in% list.files(paste0("Results/", scenario))) {
#   fit_mala(sim.list, scenario.str = scenario, mcmc.iters = 110000, mcmc.burnin = 10000, mcmc.retain = 100,
#            latent.prior.means = c(1, 10),
#            latent.prior.vars = diag(c(1, 0.5), 2),
#            init.prior = c(exp(1), 10)
#   )
# }

########################
########################

