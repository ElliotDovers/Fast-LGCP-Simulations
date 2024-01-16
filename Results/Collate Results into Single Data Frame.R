## Analyse Scenario ##

home.wd <- getwd()
scenarios <- c(
  "SS_N_1000_exp",
  "SS_N_1000_stable",
  "SW_N_1000_exp",
  "SW_N_1000_stable",
  "SW_N_10000_stable",
  "SW_N_100000_stable",
  "WS_N_1000_exp",
  "WS_N_1000_stable",
  "WS_N_10000_stable",
  "WS_N_100000_stable",
  "WW_N_1000_exp",
  "WW_N_1000_stable"
  )
# scenarios <- c('SS_N_100_stable', 'SW_N_100_stable', 'WS_N_100_stable', 'WW_N_100_stable')

# initialise the result storage
dat <- NULL

# outer loop through scenarios
for (scen in scenarios) {
  # change into appropriate result folder
  setwd(paste(home.wd, scen, sep = "/"))
  # get a list of all the individual result files
  res.list <- list.files()[grepl("res_", list.files(), fixed = T)]
  # inner loop through individual sim-by-model files
  for (sim in res.list) {
    sim.no <- as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(sim, ".", fixed = T), function(x){x[1]})), "_", fixed = T), function(x){x[length(x)]})))
    # load in the individual simulation results
    load(sim)

    # have to distinguish between INLA fits (pp or quad means mesh built from point pattern or quadrature)
    if (grepl("inla_pp", sim, fixed = T)) {
      sim.res$model <- "inla_pp"
    }
    if (grepl("inla_quad", sim, fixed = T)) {
      sim.res$model <- "inla_quad"
    }
    
    # # have to add metrics that are missing from some results
    # if (!"sqerr_median" %in% names(sim.res)) {
    #   sim.res$sqerr_median <- NA
    #   sim.res$coverage_cred <- NA
    #   sim.res$ciw <- NA
    # }
    # have to add metrics that are missing from some results
    if (!"atime" %in% names(sim.res)) {
      sim.res$atime <- NA
    }
    
    # append the result data
    dat <- rbind(dat, sim.res)
    rm(sim.res)
  }
  rm(res.list)
}
setwd(home.wd)
# save(list = c("dat"), file = "raw results in single frame.RDATA")
load("raw results in single frame.RDATA")
# REMOVE INDIVIDUAL SCENARIOS so they can be re-loaded
# dat <- dat[dat$scenario != "SW_N_100000_stable", ]

# Initial checking of results
library(dplyr)
# look at the number of flagged (potential failed) fits for each model within each scenario
View(dat %>% group_by(scenario, model) %>% summarise(flagged = sum(flagfit),
                                                     fitted = length(unique(sim)))
)
# look at the number of erroneous metrics for each model within each scenario
View(dat %>% group_by(scenario, model) %>% summarise(sqerr = sum(is.na(sqerr)),
                                                     coverage = sum(is.na(coverage)),
                                                     dkl = sum(is.na(dkl)),
                                                     ll = sum(is.na(ll)))
)
# look at the number of infinite metrics for each model within each scenario
View(dat %>% group_by(scenario, model) %>% summarise(sqerr = sum(is.infinite(sqerr)),
                                      coverage = sum(is.nan(coverage)),
                                      dkl = sum(is.infinite(dkl)),
                                      ll = sum(is.infinite(ll)))
)

# LAPLACE APPROACH: There are some issues with what are ostensibly fit failures (infinite DKL) for the Laplace version:
data.frame(dat %>% filter(is.infinite(dkl)) %>% group_by(scenario, model) %>% summarise(nsim = length(unique(sim))))
# Most of these are flagged as having a convergence failure, e.g.
data.frame(dat %>% filter(is.infinite(dkl)) %>% group_by(scenario, model) %>% summarise(nsim =  length(unique(sim)), flagged = sum(flagfit)))
# And occur exclusively in the scenarios in which E[N]=1000.
# These are most likely a case of over parametrising since they occur more frequently with greater basis functions.

# models for which these are flagged as poor optimisation probably warrant removal
dat[dat$flagfit, c("sqerr", "dkl", "coverage", "ll")] <- NA

# Even after removing the bad convergence fits there is still:
dat[is.infinite(dat$dkl) & !is.na(dat$dkl) & !dat$flagfit, ]
# These coincide with massive log-likelihoods which also indicate a fit failure - convergence to an incorrect area of the parameter space
dat[dat$ll > 1e8 & !is.na(dat$ll), c("sqerr", "dkl", "coverage", "ll", "scenario", "model")]
dat[dat$ll > 1e8 & !is.na(dat$ll), c("sqerr", "dkl", "coverage", "ll")] <- NA
# NOTE THESE ARE ALL MOSTLY RESULTS IN E[N]=1000 FOR TOO MANY BASIS FUNCTIONS

# The expected number of points has rounding issues:
unique(dat$expected_n)
# So instead we can set the expected number of points via the scenario:
dat$expected_n <- as.numeric(unlist(lapply(strsplit(as.character(dat$scenario), "_", fixed = T), function(x){x[3]})))


## THERE ALSO APPEARS TO BE SOME TERRIBLE CONVERGENCES (BOTH FOR INLA AND VA/Lp) THAT SKEW THE RESULTS:
## Many of the above are INLA n (so not of interest) and VA/Lp 7x7 in the "wrong" scenario.
View(dat %>% filter(expected_n != 1000) %>% group_by(scenario, expected_n, model) %>% summarise(est_more_than_10 = sum(abs(slope_est) > 10 & !is.na(slope_est) & !is.na(sqerr)),
                                                                                                est_more_than_7 = sum(abs(slope_est) > 7 & !is.na(slope_est) & !is.na(sqerr)),
                                                                                                est_more_than_5 = sum(abs(slope_est) > 5 & !is.na(slope_est) & !is.na(sqerr))))
# Most crucially, INLA doesn't have as easy a mechanism for identifying convergence to a poor result
# These results are contributing to INLA's poor performance, particularly in the comp. heavy settings:
table(dat[dat$expected_n != 1000 & abs(dat$slope_est) > 7 & !is.na(dat$slope_est), "model"])
dat[dat$expected_n != 1000 & abs(dat$slope_est) > 7 & !is.na(dat$slope_est) & dat$model == "inla_quad", ]
# I will generously assume these are failed fits for INLA and remove them:
dat[dat$expected_n != 1000 & abs(dat$slope_est) > 7 & !is.na(dat$slope_est) & dat$model == "inla_quad", c("sqerr", "dkl", "coverage", "ll")] <- NA

## Want to generously remove sims for which INLA (inla_quad as reported) fits have appeared to converge terribly
# bad_inla_sims_SW_10000 <- as.vector(dat %>%
#   filter(scenario == "SW_N_10000_stable" & model == "inla_quad") %>%
#   filter(abs(slope_est) > 5) %>% select(sim))
# bad_inla_sims_WS_10000 <- as.vector(dat %>%
#   filter(scenario == "WS_N_10000_stable" & model == "inla_quad") %>%
#   filter(abs(slope_est) > 5) %>% select(sim))
# bad_inla_sims_SW_100000 <- (dat %>%
#   filter(scenario == "SW_N_100000_stable" & model == "inla_quad") %>%
#   filter(abs(slope_est) > 7))$sim
# bad_inla_sims_WS_100000 <- (dat %>%
#   filter(scenario == "WS_N_100000_stable" & model == "inla_quad") %>%
#   filter(abs(slope_est) > 7))$sim
# # Remove all these sims
# dat[(dat$scenario == "WS_N_100000_stable") & (dat$sim %in% bad_inla_sims_WS_100000), c("sqerr", "dkl", "coverage", "ll")] <- NA
# dat[(dat$scenario == "SW_N_100000_stable") & (dat$sim %in% bad_inla_sims_SW_100000), c("sqerr", "dkl", "coverage", "ll")] <- NA

# Re-level models factor
dat$model[dat$model == "vark_5x5_diag"] <- "vark_5x5"
dat$model[dat$model == "vark_7x7_diag"] <- "vark_7x7"
dat$model[dat$model == "vark_14x14_diag"] <- "vark_14x14"
raw.model.names <- c('ipp', 'vark_7x7', 'lprk_7x7', 'vark_14x14', 'lprk_14x14', 'vark_21x21', 'lprk_21x21', 'vark_28x28', 'lprk_28x28', 'vark_35x35', 'lprk_35x35', 'vark_search', 'lprk_search', 'inla_pp', 'inla_quad', 'inla', 'nngp')
dat <- dat[dat$model %in% raw.model.names, ]
dat$model <- factor(dat$model, levels = raw.model.names,
                    labels = c("IPP", "VA 7x7", "Lp 7x7", "VA 14x14", "Lp 14x14", "VA 21x21", "Lp 21x21", "VA 28x28", "Lp 28x28", "VA 35x35", "Lp 35x35", "VA adapt.", "Lp adapt.", "INLA n", "INLA", "INLA old", "spNNGP"))
# raw.model.names <- c('ipp', 'vark_7x7', 'lprk_7x7', 'vark_14x14', 'lprk_14x14', 'vark_21x21', 'lprk_21x21', 'vark_28x28', 'lprk_28x28', 'vark_35x35', 'lprk_35x35', 'vark_search', 'lprk_search', 'inla_pp', 'inla_quad', 'inla', 'nngp')
# dat <- dat[dat$model %in% raw.model.names, ]
# dat$model <- factor(dat$model, levels = raw.model.names,
#                     labels = c("IPP", "VA 7x7", "Lp 7x7", "VA 14x14", "Lp 14x14", "VA 21x21", "Lp 21x21", "VA 28x28", "Lp 28x28", "VA 35x35", "Lp 35x35", "VA adapt.", "Lp adapt.", "INLA n", "INLA", "INLA old", "spNNGP"))

# get the interval width
dat$int_wid <- dat$ciw * 2
# Set the scenario as a factor
dat$full_scenario <- dat$scenario
dat$scenario <- as.factor(unlist(lapply(strsplit(as.character(dat$scenario), "_", fixed = T), function(x){x[1]})))
# Set the latent covariance function type
dat$latent_covar_fn <- as.factor(unlist(lapply(strsplit(as.character(dat$full_scenario), "_", fixed = T), function(x){x[4]})))

save(list = "dat", file = "cleaned simulation results.RDATA")
################################################################################
