# Load the cleaned single data frame of results
load("cleaned simulation results.RDATA")

# Set the plot resolution
plot.res <- 200

library(ggplot2)
library(dplyr)
library(xtable)

###############################
## Plot 1: Computation Times ##
###############################

tdat <- dat %>%
  filter(model %in% c("IPP", "VA 7x7", "Lp 7x7", "VA 14x14", "Lp 14x14", "INLA", "spNNGP") &
           latent_covar_fn == "stable" &
           expected_n != 100 &
           scenario == "WS") %>%
  group_by(model, expected_n) %>%
  summarise(time = mean(ftime, na.rm = T), time2 = median(ftime, na.rm = T),
            n = sum(!is.na(ftime)))

# add in missing
# tdat <- rbind(tdat, data.frame(model = c("spNNGP"), expected_n = c(100000), time = c(NA), n = c(0)))

tdat$expected_n <- factor(tdat$expected_n, levels = c(1000, 10000, 1e+05), labels = c("1,000", "10,000", "100,000"))
# re-level the model factor
tdat$model <- factor(tdat$model,
                     levels = c("IPP", "VA 7x7", "Lp 7x7", "VA 14x14", "Lp 14x14", "INLA", "spNNGP"),
                     labels = c("IPP", "VA 7x7", "Lp 7x7", "VA 14x14", "Lp 14x14", "INLA", "spNNGP"))

# add in some text labels for the timings
tdat$x.dodge <- as.numeric(tdat$model)
tdat$x.dodge[tdat$expected_n == "100,000"] <- tdat$x.dodge[tdat$expected_n == "100,000"] + 1/4
tdat$x.dodge[tdat$expected_n == "1,000"] <- tdat$x.dodge[tdat$expected_n == "1,000"] - 1/4
tdat$time.labs <- paste(round(tdat$time, 0), "secs")
tdat$time.labs[tdat$time < 1 & !is.na(tdat$time)] <- "<1 sec"
tdat$time.labs[tdat$time < 1000 & tdat$time > 100 & !is.na(tdat$time)] <- paste(round(tdat$time[tdat$time < 1000 & tdat$time > 100 & !is.na(tdat$time)] / 60, 0), "mins")
tdat$time.labs[tdat$time > 1000 & !is.na(tdat$time)] <- paste(round(tdat$time[tdat$time > 1000 & !is.na(tdat$time)] / 60^2, 0), "hrs")
tdat$time.labs[is.na(tdat$time)] <- ""
tdat$y.dodge <- rep(0.8, nrow(tdat))
tdat$y.dodge[tdat$time < 1 & !is.na(tdat$time)] <- 5.5

# Main Timing Plot
png(filename = paste0(getwd(), "/Figures/Figure 3.png"), width = 5.8 * plot.res, height = 4.5 * plot.res, res = plot.res)
ggplot(tdat, aes(fill = expected_n, y = time, x = model, width=.75)) +
  geom_bar(position="dodge", stat="identity") + coord_flip() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), axis.text.x = element_text(size=10),
        legend.text=element_text(size=10)) +
  scale_fill_manual(values = c("royalblue", "blue", "navy"), name = "Expected\n number\n of points") +
  scale_y_continuous(name = "Average Computation Time",
                     breaks = c(1, 60,60*10,1*(60^2), 4*(60^2)),
                     labels = c("1sec", "1min", "10mins", "1hrs", "4hrs"),
                     trans = "log") +
  scale_x_discrete(name = "") +
  annotate("text", x = tdat$x.dodge, y = tdat$y.dodge, label = tdat$time.labs, cex = 3.2, hjust = 1)
dev.off()

###################################
## Plot : Accuracy and Inference ##
###################################

# Get the result summaries for simulations E[N]=1000 and stable latent covariance function
res_1000_stable <- dat %>%
  filter(latent_covar_fn == "stable" & expected_n == 1000) %>%
  group_by(scenario, model) %>%
  summarise(rmse = sqrt(mean(sqerr, na.rm = T)),
            rmedse = sqrt(median(sqerr, na.rm = T)),
            dkl_mean = mean(dkl, na.rm = T),
            dkl_median = median(dkl, na.rm = T),
            coverage = mean(coverage, na.rm = T),
            ci_width_mean = mean(int_wid, na.rm = T),
            ci_width_median = median(int_wid, na.rm = T),
            ll_mean = mean(ll, na.rm = T),
            ll_median = median(ll, na.rm = T),
            fit_failures = 1000 - length(unique(sim)),
            convergence_failure = sum(flagfit)
  )

# Select out the models we want for the comparisons here:
reduced.models <- c("IPP", "VA 7x7", "Lp 7x7", "VA 14x14", "Lp 14x14", "INLA", "spNNGP")
dat4plot <- res_1000_stable %>% filter(model %in% reduced.models)
dat4plot$model <- factor(dat4plot$model, levels = reduced.models,
                         labels = c("IPP", "VA 7x7", "Lp 7x7", "VA 14x14", "Lp 14x14", "INLA", "spNNGP"))
dat4plot$scenario <- factor(dat4plot$scenario, levels = c("SS", "SW","WS", "WW"),
                            labels = c("S,S", "S,W", "W,S", "W,W"))
plot.res <- 200
################################################################################
## Plotting single result points for each model in each scenario
png(filename = paste0(getwd(), "/Figures/Figure 4.png"), width = 4.8 * plot.res, height = 5.7 * plot.res, res = plot.res)
layout(mat = matrix(c(1,2,3), nrow = 1, ncol = 3, byrow = TRUE), heights = 1, widths = c(0.35,0.25,0.4))#rep(0.25,4))
jitteredy <- as.numeric(dat4plot$scenario)
jit_unit <- 1/(length(levels(dat4plot$model)) + 1)
if (length(levels(dat4plot$model)) %% 2 == 0) {
  tmp.seq <- seq(1, length(levels(dat4plot$model)), by = 2)
  jit_multiplier <- c(-sort(tmp.seq, decreasing = T), tmp.seq) / 2
} else {
  tmp.seq <- seq(0, length(levels(dat4plot$model)), by = 2)
  jit_multiplier <- c(-sort(tmp.seq, decreasing = T), tmp.seq[-1]) / 2
}
for (i in 1:length(levels(dat4plot$model))) {
  jitteredy[dat4plot$model == levels(dat4plot$model)[i]] <- jitteredy[dat4plot$model == levels(dat4plot$model)[i]] + (jit_multiplier[i] * jit_unit)
}

# set the ylims for the plot
ylims <- range(jitteredy)
# ylims[1] <- ylims[1] - (diff(sort(jitteredy))[1] * 0.60) # the limits seem to add a little
# ylims[2] <- ylims[2] + (diff(sort(jitteredy))[1] * 0.60)

# fix outer margins
par(mgp = c(1.75,0.75,0), xpd = F)

# Plot 1
par(mar = c(4.1, 4.1, 1.2, 0))
plot.default(dat4plot$rmse, jitteredy, type = "p", yaxt = "n", ylim = ylims, log = "x", 
             col = "grey25",
             pch = 16,
             ylab = "",
             xlab = "", xaxt = "n",
             main = "A"
)
axis(2, at = 1:length(levels(dat4plot$scenario)), labels = levels(dat4plot$scenario), cex.axis = 1.2)
axis(1, at = c(0.05, 0.2, 0.8), labels = c(0.05, 0.2, 0.8), cex.axis = 1.2)
abline(h = (1:length(levels(dat4plot$scenario)))[-length(levels(dat4plot$scenario))] + diff(1:length(levels(dat4plot$scenario))) / 2)
abline(h = jitteredy, col = "grey25", lty = "dashed", lwd = 0.5)
title(ylab = "Covariate and Latent Field Types", cex.lab = 1.5, line = 3)
title(xlab = expression(paste("RMSE: ", hat(beta[1]))), cex.lab = 1.3, line = 2.5)

# Plot 2
par(mar = c(4.1, 0.25, 1.2, 0.25))
plot.default(dat4plot$dkl_mean, jitteredy, type = "p", yaxt = "n", ylim = ylims, log = "x", 
             col = "grey25",
             pch = 16,
             ylab = "",
             xlab = "",
             main = "B", cex.axis = 1.2
)
abline(h = (1:length(levels(dat4plot$scenario)))[-length(levels(dat4plot$scenario))] + diff(1:length(levels(dat4plot$scenario))) / 2)
abline(h = jitteredy, col = "grey25", lty = "dashed", lwd = 0.5)
title(xlab = expression(paste("KL Divergence ", hat(lambda))), cex.lab = 1.3, line = 2.5)

# Plot 3
par(mar = c(4.1, 0, 1.2, 6.1))
plot.default(dat4plot$coverage, jitteredy, type = "p", yaxt = "n", ylim = ylims, #log = "x", 
             col = "grey25",
             pch = 16,
             ylab = "",
             xlab = "",
             main = "C", cex.axis = 1.2
)
abline(v = 0.95, lty = "dotted", col = "navyblue")
abline(h = (1:length(levels(dat4plot$scenario)))[-length(levels(dat4plot$scenario))] + diff(1:length(levels(dat4plot$scenario))) / 2)
abline(h = jitteredy, col = "grey25", lty = "dashed", lwd = 0.5)
title(xlab = expression(paste(beta[1], " CI Coverage")), cex.lab = 1.3, line = 2.5)
axis(4, at = jitteredy, labels = substr(dat4plot$model, 1, 9), las = 2, cex.axis = 1.4)

dev.off()
################################################################################

##############################################
## Tables for basis configuration selection ##
##############################################

library(tidyr)

# for stable covariance fn sims
tmp <- dat %>%
  filter(expected_n == 1000 &
           model %in% c("VA 7x7", "VA 14x14", "Lp 7x7", "Lp 14x14") &
           latent_covar_fn == "stable") %>%
  select(scenario, sim, model, ll) %>%
  spread(model, ll)
# tmp[is.na(tmp$`VA 7x7`) | is.na(tmp$`VA 14x14`), c("VA 7x7", "VA 14x14")] <- NA
# tmp[is.na(tmp$`Lp 7x7`) | is.na(tmp$`Lp 14x14`), c("Lp 7x7", "Lp 14x14")] <- NA
# tmp$Lp <- apply(tmp[ , c("Lp 7x7", "Lp 14x14")], 1, which.max)
# tmp$VA <- apply(tmp[ , c("VA 7x7", "VA 14x14")], 1, which.max)
tmp$Lp <- NULL
tmp$VA <- NULL
for (i in 1:nrow(tmp)) {
  if (any(is.na(tmp[i, c("Lp 7x7", "Lp 14x14")]))) {
    tmp$Lp[i] <- NA
  } else {
    tmp$Lp[i] <- which.max(tmp[i , c("Lp 7x7", "Lp 14x14")])
  }
  if (any(is.na(tmp[i, c("VA 7x7", "VA 14x14")]))) {
    tmp$VA[i] <- NA
  } else {
    tmp$VA[i] <- which.max(tmp[i , c("VA 7x7", "VA 14x14")])
  }
}

stable_ll_comp <- tmp %>%
  group_by(scenario) %>%
  summarise(`Lp 7x7 (fits failed)` = paste0(sum(Lp == 1, na.rm = T), " (", sum(is.na(`Lp 7x7`)), ")"),
            `Lp 14x14 (fits failed)` = paste0(sum(Lp == 2, na.rm = T), " (", sum(is.na(`Lp 14x14`)), ")"),
            `VA 7x7 (fits failed)` = paste0(sum(VA == 1, na.rm = T), " (", sum(is.na(`VA 7x7`)), ")"),
            `VA 14x14 (fits failed)` = paste0(sum(VA == 2, na.rm = T), " (", sum(is.na(`VA 14x14`)), ")"))

# for exponential covariance fn sims
tmp <- dat %>%
  filter(expected_n == 1000 &
           model %in% c("VA 7x7", "VA 14x14", "Lp 7x7", "Lp 14x14") &
           latent_covar_fn == "exp") %>%
  select(scenario, sim, model, ll) %>%
  spread(model, ll)
# tmp[is.na(tmp$`VA 7x7`) | is.na(tmp$`VA 14x14`), c("VA 7x7", "VA 14x14")] <- NA
# tmp[is.na(tmp$`Lp 7x7`) | is.na(tmp$`Lp 14x14`), c("Lp 7x7", "Lp 14x14")] <- NA
# tmp$Lp <- apply(tmp[ , c("Lp 7x7", "Lp 14x14")], 1, which.max)
# tmp$VA <- apply(tmp[ , c("VA 7x7", "VA 14x14")], 1, which.max)
tmp$Lp <- NULL
tmp$VA <- NULL
for (i in 1:nrow(tmp)) {
  if (any(is.na(tmp[i, c("Lp 7x7", "Lp 14x14")]))) {
    tmp$Lp[i] <- NA
  } else {
    tmp$Lp[i] <- which.max(tmp[i , c("Lp 7x7", "Lp 14x14")])
  }
  if (any(is.na(tmp[i, c("VA 7x7", "VA 14x14")]))) {
    tmp$VA[i] <- NA
  } else {
    tmp$VA[i] <- which.max(tmp[i , c("VA 7x7", "VA 14x14")])
  }
}

exp_ll_comp <- tmp %>%
  group_by(scenario) %>%
  summarise(`Lp 7x7 (fits failed)` = paste0(sum(Lp == 1, na.rm = T), " (", sum(is.na(`Lp 7x7`)), ")"),
            `Lp 14x14 (fits failed)` = paste0(sum(Lp == 2, na.rm = T), " (", sum(is.na(`Lp 14x14`)), ")"),
            `VA 7x7 (fits failed)` = paste0(sum(VA == 1, na.rm = T), " (", sum(is.na(`VA 7x7`)), ")"),
            `VA 14x14 (fits failed)` = paste0(sum(VA == 2, na.rm = T), " (", sum(is.na(`VA 14x14`)), ")"))

# put into latex table form:
library(xtable)
print(xtable(stable_ll_comp, label = "tab:append_sim_data_bf_choice_stable", digits = 0,
             caption = c("The number of simulations for which each basis function configuration achieves the higher approximate, marginal log-likelihood within each scenario.", "Simulation results: basis function configurations")),
      include.rownames=F, caption.placement = "top")

print(xtable(exp_ll_comp, label = "tab:append_sim_data_bf_choice_exp", digits = 0,
             caption = c("The number of simulations for which each basis function configuration achieves the higher approximate, marginal log-likelihood within each scenario.", "Simulation results: basis function configurations")),
      include.rownames=F, caption.placement = "top")

##############################################################################

############################
## Plots for the appendix ##
############################

## SECTION S2.1

## Comparing results with those for the exponential covariance function

# Get the result summaries for simulations E[N]=1000 and both latent covariance functions
res_1000 <- dat %>%
  filter(expected_n == 1000) %>%
  group_by(scenario, model, latent_covar_fn) %>%
  summarise(rmse = sqrt(mean(sqerr, na.rm = T)),
            rmedse = sqrt(median(sqerr, na.rm = T)),
            dkl_mean = mean(dkl, na.rm = T),
            dkl_median = median(dkl, na.rm = T),
            coverage = mean(coverage, na.rm = T),
            ci_width_mean = mean(int_wid, na.rm = T),
            ci_width_median = median(int_wid, na.rm = T),
            ll_mean = mean(ll, na.rm = T),
            ll_median = median(ll, na.rm = T),
            fit_failures = 1000 - length(unique(sim)),
            convergence_failure = sum(flagfit)
  )

# Select out the models we want for the comparisons here:
reduced.models <- c("IPP", "VA 7x7", "Lp 7x7", "VA 14x14", "Lp 14x14", "INLA", "spNNGP")
dat4plot <- res_1000 %>% filter(model %in% reduced.models)
dat4plot$model <- factor(dat4plot$model, levels = reduced.models,
                         labels = c("IPP", "VA 7x7", "Lp 7x7", "VA 14x14", "Lp 14x14", "INLA", "spNNGP"))
dat4plot$scenario <- factor(dat4plot$scenario, levels = c("SS","SW","WS","WW"),
                            labels = c("Long,Long", "Long,Short", "Short,Long", "Short,Short"))
plot.res <- 200

png(filename = paste0(getwd(), "/Figures/Appendix Figure 2.png"), width = 5 * plot.res, height = 5.5 * plot.res, res = plot.res)
layout(mat = matrix(c(1,2,3), nrow = 1, ncol = 3, byrow = TRUE), heights = 1, widths = c(0.375,0.25,0.375))#rep(0.25,4))
jitteredy <- as.numeric(dat4plot$scenario)
jit_unit <- 1/(length(levels(dat4plot$model)) + 1)
if (length(levels(dat4plot$model)) %% 2 == 0) {
  tmp.seq <- seq(1, length(levels(dat4plot$model)), by = 2)
  jit_multiplier <- c(-sort(tmp.seq, decreasing = T), tmp.seq) / 2
} else {
  tmp.seq <- seq(0, length(levels(dat4plot$model)), by = 2)
  jit_multiplier <- c(-sort(tmp.seq, decreasing = T), tmp.seq[-1]) / 2
}
for (i in 1:length(levels(dat4plot$model))) {
  jitteredy[dat4plot$model == levels(dat4plot$model)[i]] <- jitteredy[dat4plot$model == levels(dat4plot$model)[i]] + (jit_multiplier[i] * jit_unit)
}

# set the ylims for the plot
ylims <- range(jitteredy)
ylims[1] <- ylims[1] - (diff(sort(jitteredy))[1] * 0.60) # the limits seem to add a little
ylims[2] <- ylims[2] + (diff(sort(jitteredy))[1] * 0.60)

# fix outer margins
par(mgp = c(1.75,0.75,0), xpd = F)

# Plot 1
par(mar = c(4.1, 5.1, 1.2, 0))
plot.default(dat4plot$rmse[dat4plot$latent_covar_fn == "stable"], jitteredy[dat4plot$latent_covar_fn == "stable"], type = "p", yaxt = "n", ylim = ylims, log = "x", 
             col = "grey25",
             pch = 1,
             ylab = "",
             xlab = "", xaxt = "n",
             main = "A"
)
points(dat4plot$rmse[dat4plot$latent_covar_fn == "exp"], jitteredy[dat4plot$latent_covar_fn == "exp"], pch = 4)
axis(2, at = 1:length(levels(dat4plot$scenario)), labels = levels(dat4plot$scenario), cex.axis = 1.2)
axis(1, at = c(0.05, 0.2, 0.8), labels = c(0.05, 0.2, 0.8), cex.axis = 1.2)
abline(h = (1:length(levels(dat4plot$scenario)))[-length(levels(dat4plot$scenario))] + diff(1:length(levels(dat4plot$scenario))) / 2)
abline(h = jitteredy, col = "grey25", lty = "dashed", lwd = 0.5)
title(ylab = "Covariate and Latent Field Correlation Ranges", cex.lab = 1.5, line = 3.5)
title(xlab = expression(paste("RMSE: ", hat(beta[1]))), cex.lab = 1.2, line = 2.5)

# Plot 2
par(mar = c(4.1, 0.25, 1.2, 0.25))
plot.default(dat4plot$dkl_mean[dat4plot$latent_covar_fn == "stable"], jitteredy[dat4plot$latent_covar_fn == "stable"], type = "p", yaxt = "n", ylim = ylims, log = "x", 
             col = "grey25",
             pch = 1,
             ylab = "",
             xlab = "",
             main = "B", cex.axis = 1.2
)
points(dat4plot$dkl_mean[dat4plot$latent_covar_fn == "exp"], jitteredy[dat4plot$latent_covar_fn == "exp"], pch = 4)
abline(h = (1:length(levels(dat4plot$scenario)))[-length(levels(dat4plot$scenario))] + diff(1:length(levels(dat4plot$scenario))) / 2)
abline(h = jitteredy, col = "grey25", lty = "dashed", lwd = 0.5)
title(xlab = expression(paste("KL Divergence ", hat(lambda))), cex.lab = 1.2, line = 2.5)

# Plot 3
par(mar = c(4.1, 0, 1.2, 5.1))
plot.default(dat4plot$coverage[dat4plot$latent_covar_fn == "stable"], jitteredy[dat4plot$latent_covar_fn == "stable"], type = "p", yaxt = "n", ylim = ylims, #log = "x", 
             col = "grey25",
             pch = 1,
             ylab = "",
             xlab = "",
             main = "C", cex.axis = 1.2
)
points(dat4plot$coverage[dat4plot$latent_covar_fn == "exp"], jitteredy[dat4plot$latent_covar_fn == "exp"], pch = 4)
abline(v = 0.95, lty = "dotted", col = "navyblue")
abline(h = (1:length(levels(dat4plot$scenario)))[-length(levels(dat4plot$scenario))] + diff(1:length(levels(dat4plot$scenario))) / 2)
abline(h = jitteredy, col = "grey25", lty = "dashed", lwd = 0.5)
title(xlab = expression(paste(beta[1], " CI Coverage")), cex.lab = 1.2, line = 2.5)
axis(4, at = jitteredy, labels = substr(dat4plot$model, 1, 9), las = 2, cex.axis = 1.2)

dev.off()
################################################################################

## SECTION S2.2 ##

####################################
## Table for number of models fit ##
####################################

# vector of models to focus on:
reduced.models <- c("IPP", "VA 7x7", "Lp 7x7", "VA 14x14", "Lp 14x14", "VA 21x21", "Lp 21x21", "VA 28x28", "Lp 28x28", "VA 35x35", "Lp 35x35", 
                    "INLA", "spNNGP")

fit_fails <- dat %>%
  filter(latent_covar_fn == "stable" &
           scenario %in% c("SW", "WS") &
           model %in% reduced.models &
           expected_n != 100)
fit_fails$scenario <- factor(fit_fails$scenario, levels = c("SW","WS"),
                             labels = c("S,W", "W,S"))
fit_fails$expected_n <- factor(fit_fails$expected_n, levels = c(1e+03, 1e+04, 1e+05),
                               labels = c("1,000", "10,000", "100,000"))
fit_fails$model <- factor(fit_fails$model, levels = reduced.models,
                          labels = reduced.models)
fit_fails$`Fit Status` <- "Successful"
fit_fails$`Fit Status`[is.na(fit_fails$sqerr)] <- "Failed Convergence"
fit_fails <- fit_fails %>% select(scenario, expected_n, model, sim, `Fit Status`)
fails4plot <- tidyr::complete(fit_fails, scenario, expected_n, model, sim = 1:1000, fill = list(`Fit Status` = "Exceeded 48hrs"))
fails4plot$`Fit Status` <- factor(fails4plot$`Fit Status`, levels = c("Exceeded 48hrs", "Failed Convergence", "Successful"),
                                  labels = c("Exceeded 48hrs", "Failed Convergence", "Successful"))

png(filename = paste0(home.wd, "/Figures/Appendix Figure 3.png"), width = 5.5 * plot.res, height = 5.5 * plot.res, res = plot.res)
ggplot(fails4plot, aes(x=model, y=sim, fill=`Fit Status`)) +
  geom_bar(position="stack", stat="identity") +
  scale_y_continuous(name = "1000 Simulated Datasets",
                     breaks = NULL,
                     labels = NULL) +
  scale_x_discrete(name = "") +
  scale_fill_manual(values=c("#CC79A7", "#E69F00", "#56B4E9")) +
  facet_wrap(~ scenario+expected_n) + coord_flip()
dev.off()

# Table for exact numbers of fitting failures
fail_tab <- fails4plot %>% group_by(scenario, expected_n, model) %>%
  summarise(Success = sum(`Fit Status` == "Successful"),
            Failed = sum(`Fit Status` == "Failed Convergence"),
            Elapsed = sum(`Fit Status` == "Exceeded 48hrs"))

fail_tab[fail_tab$model == "spNNGP", ]
fail_tab[fail_tab$model == "INLA", ]
fail_tab[fail_tab$model == "VA 7x7", ]
fail_tab[fail_tab$model == "VA 14x14", ]
fail_tab[fail_tab$model == "VA 21x21", ]

## Separate the data into comparable simulations ##

# get the lists of successful INLA fits in each scenario
inla_sims <- aggregate(fails4plot$sim[fails4plot$model == "INLA" & fails4plot$`Fit Status` == "Successful"],
          by = list(fails4plot$scenario[fails4plot$model == "INLA" & fails4plot$`Fit Status` == "Successful"], fails4plot$expected_n[fails4plot$model == "INLA" & fails4plot$`Fit Status` == "Successful"]),
          function(x){unique(x)})
# New results data frame comparing only those sims for which INLA did not fail
newdat <- rbind(
  dat[dat$expected_n == 1000 & dat$scenario == "SW" & dat$sim %in% inla_sims$x[[which(inla_sims$Group.1 == "S,W" & inla_sims$Group.2 == "1,000")]], ],
  dat[dat$expected_n == 1000 & dat$scenario == "WS" & dat$sim %in% inla_sims$x[[which(inla_sims$Group.1 == "W,S" & inla_sims$Group.2 == "1,000")]], ],
  dat[dat$expected_n == 10000 & dat$scenario == "SW" & dat$sim %in% inla_sims$x[[which(inla_sims$Group.1 == "S,W" & inla_sims$Group.2 == "10,000")]], ],
  dat[dat$expected_n == 10000 & dat$scenario == "WS" & dat$sim %in% inla_sims$x[[which(inla_sims$Group.1 == "W,S" & inla_sims$Group.2 == "10,000")]], ],
  dat[dat$expected_n == 100000 & dat$scenario == "SW" & dat$sim %in% inla_sims$x[[which(inla_sims$Group.1 == "S,W" & inla_sims$Group.2 == "100,000")]], ],
  dat[dat$expected_n == 100000 & dat$scenario == "WS" & dat$sim %in% inla_sims$x[[which(inla_sims$Group.1 == "W,S" & inla_sims$Group.2 == "100,000")]], ]
)

# Get the results for heavy computations (these sims use stable covariance function)
res_heavy_computations <- newdat %>%
  filter(latent_covar_fn == "stable" & scenario %in% c("SW", "WS")) %>%
  group_by(scenario, expected_n, model) %>%
  summarise(rmse = sqrt(mean(sqerr, na.rm = T)),
            rmedse = sqrt(median(sqerr, na.rm = T)),
            dkl_mean = mean(dkl, na.rm = T),
            dkl_median = median(dkl, na.rm = T),
            coverage = mean(coverage, na.rm = T),
            ci_width_mean = mean(int_wid, na.rm = T),
            ci_width_median = median(int_wid, na.rm = T),
            ll_mean = mean(ll, na.rm = T),
            ll_median = median(ll, na.rm = T),
            fit_failures = 1000 - length(unique(sim)),
            convergence_failure = sum(flagfit)
  )

# Select out the models we want for the comparisons here:
reduced.models <- c("IPP", "VA 7x7", "Lp 7x7", "VA 14x14", "Lp 14x14", "VA 21x21", "Lp 21x21", "VA 28x28", "Lp 28x28", "VA 35x35", "Lp 35x35", #"VA 35x35 diag", "VA 35x35 limit", "VA 42x42 diag", "VA 42x42 limit", "Lp 42x42",
                    "INLA", "spNNGP")
# reduced.models <- c("IPP", "VA 7x7", "Lp 7x7", "VA 14x14", "Lp 14x14", "INLA q", "spNNGP")
dat4plot <- res_heavy_computations %>% filter(model %in% reduced.models)
dat4plot$model <- factor(dat4plot$model, levels = reduced.models, labels = reduced.models)
dat4plot$scenario <- factor(dat4plot$scenario, levels = c("SW","WS"),
                            labels = c("S,W", "W,S"))
plot.res <- 200
robust <- TRUE
if (robust) {
  dat4plot$met1 <- dat4plot$rmedse
  dat4plot$met2 <- dat4plot$dkl_median
} else {
  dat4plot$met1 <- dat4plot$rmse
  dat4plot$met2 <- dat4plot$dkl_mean
}

# Plot of the heavy computation settings
png(filename = paste0(getwd(), "/Figures/Appendix Figure 4 (up to 35x35) w VA Testing.png"), width = 5 * plot.res, height = 5 * plot.res, res = plot.res)
layout(mat = matrix(c(1,2,3), nrow = 1, ncol = 3, byrow = TRUE), heights = 1, widths = c(0.375,0.25,0.375))#rep(0.25,4))
jitteredy <- as.numeric(dat4plot$scenario)
jit_unit <- 1/(length(levels(dat4plot$model)) + 1)
if (length(levels(dat4plot$model)) %% 2 == 0) {
  tmp.seq <- seq(1, length(levels(dat4plot$model)), by = 2)
  jit_multiplier <- c(-sort(tmp.seq, decreasing = T), tmp.seq) / 2
} else {
  tmp.seq <- seq(0, length(levels(dat4plot$model)), by = 2)
  jit_multiplier <- c(-sort(tmp.seq, decreasing = T), tmp.seq[-1]) / 2
}
for (i in 1:length(levels(dat4plot$model))) {
  jitteredy[dat4plot$model == levels(dat4plot$model)[i]] <- jitteredy[dat4plot$model == levels(dat4plot$model)[i]] + (jit_multiplier[i] * jit_unit)
}

# set the ylims for the plot
ylims <- range(jitteredy)
ylims[1] <- ylims[1] - (diff(sort(jitteredy))[1] * 0.60) # the limits seem to add a little
ylims[2] <- ylims[2] + (diff(sort(jitteredy))[1] * 0.60)
# set the xlims for the plot
xlims_rmse <- range(dat4plot$met1, na.rm = T)
xlims_dkl <- range(dat4plot$met2, na.rm = T)

# fix outer margins
par(mgp = c(1.75,0.75,0), xpd = F)

# Plot 1
par(mar = c(4.1, 5.1, 1.2, 0))
plot.default(dat4plot$met1[dat4plot$expected_n == 100000], jitteredy[dat4plot$expected_n == 100000], type = "p", yaxt = "n", ylim = ylims, log = "x",
             col = adjustcolor("grey25",alpha.f=0.4),
             pch = 16,
             ylab = "",
             xlab = "",
             xaxt = "n", main = "A", xlim = xlims_rmse, cex = 3
)
points(dat4plot$met1[dat4plot$expected_n == 10000], jitteredy[dat4plot$expected_n == 10000], pch = 16, cex = 2, col = adjustcolor("grey25",alpha.f=0.6))
points(dat4plot$met1[dat4plot$expected_n == 1000], jitteredy[dat4plot$expected_n == 1000], pch = 16, cex = 1, col = adjustcolor("grey25",alpha.f=1))
axis(2, at = 1:length(levels(dat4plot$scenario)), labels = levels(dat4plot$scenario), cex.axis = 1.2)
axis(1, at = c(0.05, 0.2, 0.8), labels = c(0.05, 0.2, 0.8), cex.axis = 1.2)
abline(h = (1:length(levels(dat4plot$scenario)))[-length(levels(dat4plot$scenario))] + diff(1:length(levels(dat4plot$scenario))) / 2)
title(ylab = "Covariate and Latent Field Types", cex.lab = 1.5, line = 3.5)
title(xlab = if (robust) {
  expression(paste("RMedSE: ", hat(beta[1])))
} else {
  expression(paste("RMSE: ", hat(beta[1])))
}, cex.lab = 1.2, line = 2.5)
abline(h = jitteredy, col = "grey25", lty = "dashed", lwd = 0.5)

# Plot 2
par(mar = c(4.1, 0.25, 1.2, 0.25))
plot.default(dat4plot$met2[dat4plot$expected_n == 100000], jitteredy[dat4plot$expected_n == 100000], type = "p", yaxt = "n", ylim = ylims, log = "x",
             col = adjustcolor("grey25",alpha.f=0.4),
             pch = 16,
             ylab = "",
             xlab = "", xaxt = "n",
             main = "B", xlim = xlims_dkl, cex = 3, cex.axis = 1.2
)
points(dat4plot$met2[dat4plot$expected_n == 10000], jitteredy[dat4plot$expected_n == 10000], pch = 16, cex = 2, col = adjustcolor("grey25",alpha.f=0.6))
points(dat4plot$met2[dat4plot$expected_n == 1000], jitteredy[dat4plot$expected_n == 1000], pch = 16, cex = 1, col = adjustcolor("grey25",alpha.f=1))
axis(1, at = c(50, 500, 5000), labels = c(50, 500, 5000), cex.axis = 1.2)
abline(h = (1:length(levels(dat4plot$scenario)))[-length(levels(dat4plot$scenario))] + diff(1:length(levels(dat4plot$scenario))) / 2)
abline(h = jitteredy, col = "grey25", lty = "dashed", lwd = 0.5)
title(xlab = expression(paste("KL Divergence ", hat(lambda))), cex.lab = 1.2, line = 2.5)

# Plot 3
par(mar = c(4.1, 0, 1.2, 5.1))
plot.default(dat4plot$coverage[dat4plot$expected_n == 100000], jitteredy[dat4plot$expected_n == 100000], type = "p", yaxt = "n", ylim = ylims, #log = "x",
             col = adjustcolor("grey25",alpha.f=0.4),
             pch = 16,
             ylab = "",
             xlab = "",
             main = "C", cex = 3, cex.axis = 1.2
)
points(dat4plot$coverage[dat4plot$expected_n == 10000], jitteredy[dat4plot$expected_n == 10000], pch = 16, cex = 2, col = adjustcolor("grey25",alpha.f=0.6))
points(dat4plot$coverage[dat4plot$expected_n == 1000], jitteredy[dat4plot$expected_n == 1000], pch = 16, cex = 1, col = adjustcolor("grey25",alpha.f=1))
abline(v = 0.95, lty = "dotted", col = "navyblue")
abline(h = (1:length(levels(dat4plot$scenario)))[-length(levels(dat4plot$scenario))] + diff(1:length(levels(dat4plot$scenario))) / 2)
abline(h = jitteredy, col = "grey25", lty = "dashed", lwd = 0.5)
title(xlab = expression(paste(beta[1], " CI Coverage")), cex.lab = 1.2, line = 2.5)
axis(4, at = jitteredy, labels = substr(dat4plot$model, 1, 9), las = 2, cex.axis = 1.2)

dev.off()
################################################################################

## DATA POOR SETTINGS
# Get the result summaries for simulations E[N]=1000 and stable latent covariance function
res_100_stable <- dat %>%
  filter(latent_covar_fn == "stable" & expected_n == 100) %>%
  group_by(scenario, model) %>%
  summarise(rmse = sqrt(mean(sqerr, na.rm = T)),
            rmedse = sqrt(median(sqerr, na.rm = T)),
            dkl_mean = mean(dkl, na.rm = T),
            dkl_median = median(dkl, na.rm = T),
            coverage = mean(coverage, na.rm = T),
            ci_width_mean = mean(int_wid, na.rm = T),
            ci_width_median = median(int_wid, na.rm = T),
            ll_mean = mean(ll, na.rm = T),
            ll_median = median(ll, na.rm = T),
            fit_failures = 1000 - length(unique(sim)),
            convergence_failure = sum(flagfit)
  )

# Select out the models we want for the comparisons here:
reduced.models <- c("VA 7x7", "Lp 7x7", "VA 14x14", "Lp 14x14")
dat4plot <- res_100_stable %>% filter(model %in% reduced.models)
dat4plot$model <- factor(dat4plot$model, levels = reduced.models,
                         labels = reduced.models)
dat4plot$scenario <- factor(dat4plot$scenario, levels = c("SS", "SW","WS", "WW"),
                            labels = c("S,S", "S,W", "W,S", "W,W"))
plot.res <- 200
################################################################################
## Plotting single result points for each model in each scenario
png(filename = paste0(getwd(), "/Figures/Appendix Figure 5.png"), width = 4.8 * plot.res, height = 5.7 * plot.res, res = plot.res)
layout(mat = matrix(c(1,2,3), nrow = 1, ncol = 3, byrow = TRUE), heights = 1, widths = c(0.35,0.25,0.4))#rep(0.25,4))
jitteredy <- as.numeric(dat4plot$scenario)
jit_unit <- 1/(length(levels(dat4plot$model)) + 1)
if (length(levels(dat4plot$model)) %% 2 == 0) {
  tmp.seq <- seq(1, length(levels(dat4plot$model)), by = 2)
  jit_multiplier <- c(-sort(tmp.seq, decreasing = T), tmp.seq) / 2
} else {
  tmp.seq <- seq(0, length(levels(dat4plot$model)), by = 2)
  jit_multiplier <- c(-sort(tmp.seq, decreasing = T), tmp.seq[-1]) / 2
}
for (i in 1:length(levels(dat4plot$model))) {
  jitteredy[dat4plot$model == levels(dat4plot$model)[i]] <- jitteredy[dat4plot$model == levels(dat4plot$model)[i]] + (jit_multiplier[i] * jit_unit)
}

# set the ylims for the plot
ylims <- range(jitteredy)
# ylims[1] <- ylims[1] - (diff(sort(jitteredy))[1] * 0.60) # the limits seem to add a little
# ylims[2] <- ylims[2] + (diff(sort(jitteredy))[1] * 0.60)

# fix outer margins
par(mgp = c(1.75,0.75,0), xpd = F)

# Plot 1
par(mar = c(4.1, 4.1, 1.2, 0))
plot.default(dat4plot$rmse, jitteredy, type = "p", yaxt = "n", ylim = ylims, log = "x", 
             col = "grey25",
             pch = 16,
             ylab = "",
             xlab = "", xaxt = "n",
             main = "A"
)
axis(2, at = 1:length(levels(dat4plot$scenario)), labels = levels(dat4plot$scenario), cex.axis = 1.2)
axis(1, at = c(0.05, 0.2, 0.8), labels = c(0.05, 0.2, 0.8), cex.axis = 1.2)
abline(h = (1:length(levels(dat4plot$scenario)))[-length(levels(dat4plot$scenario))] + diff(1:length(levels(dat4plot$scenario))) / 2)
abline(h = jitteredy, col = "grey25", lty = "dashed", lwd = 0.5)
title(ylab = "Covariate and Latent Field Types", cex.lab = 1.5, line = 3)
title(xlab = expression(paste("RMSE: ", hat(beta[1]))), cex.lab = 1.3, line = 2.5)

# Plot 2
par(mar = c(4.1, 0.25, 1.2, 0.25))
plot.default(dat4plot$dkl_mean, jitteredy, type = "p", yaxt = "n", ylim = ylims, log = "x", 
             col = "grey25",
             pch = 16,
             ylab = "",
             xlab = "",
             main = "B", cex.axis = 1.2
)
abline(h = (1:length(levels(dat4plot$scenario)))[-length(levels(dat4plot$scenario))] + diff(1:length(levels(dat4plot$scenario))) / 2)
abline(h = jitteredy, col = "grey25", lty = "dashed", lwd = 0.5)
title(xlab = expression(paste("KL Divergence ", hat(lambda))), cex.lab = 1.3, line = 2.5)

# Plot 3
par(mar = c(4.1, 0, 1.2, 6.1))
plot.default(dat4plot$coverage, jitteredy, type = "p", yaxt = "n", ylim = ylims, #log = "x", 
             col = "grey25",
             pch = 16,
             ylab = "",
             xlab = "",
             main = "C", cex.axis = 1.2
)
abline(v = 0.95, lty = "dotted", col = "navyblue")
abline(h = (1:length(levels(dat4plot$scenario)))[-length(levels(dat4plot$scenario))] + diff(1:length(levels(dat4plot$scenario))) / 2)
abline(h = jitteredy, col = "grey25", lty = "dashed", lwd = 0.5)
title(xlab = expression(paste(beta[1], " CI Coverage")), cex.lab = 1.3, line = 2.5)
axis(4, at = jitteredy, labels = substr(dat4plot$model, 1, 9), las = 2, cex.axis = 1.4)

dev.off()
################################################################################

#######################################
## Complete Simulation Result Tables ##
#######################################

tab1 <- dat %>%
  filter(latent_covar_fn == "stable" & expected_n == 1000) %>%
  group_by(scenario, model) %>%
  summarise(comptime = mean(ftime, na.rm = T),
            rmse = sqrt(mean(sqerr, na.rm = T)),
            dkl = mean(dkl, na.rm = T),
            coverage = mean(coverage, na.rm = T),
            ci_width = mean(int_wid, na.rm = T),
            ll = mean(ll, na.rm = T),
            fit_failures = 1000 - length(unique(sim)),
            convergence_failure = sum(is.na(sqerr))
  ) %>%
  filter(model %in% c("IPP", "VA 7x7", "Lp 7x7", "VA 14x14", "Lp 14x14", "INLA", "spNNGP"))
tab1$scenario <- factor(tab1$scenario, levels = c("SS", "SW", "WS", "WW"),
                        labels = c("S,S", "S,W", "W,S", "W,W"))
tab1$fit_failures <- as.integer(tab1$fit_failures)
print(xtable(tab1, label = "tab:append_all_res_1000_stable", digits = 2,
             caption = c("Complete simulation result summaries for point patterns with E[N]=1000 and latent fields with a squared-exponential (Gaussian) covariance function.")),
      include.rownames=F, caption.placement = "top")


tab2 <- dat %>%
  filter(latent_covar_fn == "exp" & expected_n == 1000) %>%
  group_by(scenario, model) %>%
  summarise(comptime = mean(ftime, na.rm = T),
            rmse = sqrt(mean(sqerr, na.rm = T)),
            dkl = mean(dkl, na.rm = T),
            coverage = mean(coverage, na.rm = T),
            ci_width = mean(int_wid, na.rm = T),
            ll = mean(ll, na.rm = T),
            fit_failures = 1000 - length(unique(sim)),
            convergence_failure = sum(is.na(sqerr))
  ) %>%
  filter(model %in% c("IPP", "VA 7x7", "Lp 7x7", "VA 14x14", "Lp 14x14", "INLA", "spNNGP"))
tab2$scenario <- factor(tab2$scenario, levels = c("SS", "SW", "WS", "WW"),
                        labels = c("S,S", "S,W", "W,S", "W,W"))
tab2$fit_failures <- as.integer(tab2$fit_failures)
print(xtable(tab2, label = "tab:append_all_res_1000_exp", digits = 2,
             caption = c("Complete simulation result summaries for point patterns with E[N]=1000 and latent fields with an exponential covariance function.")),
      include.rownames=F, caption.placement = "top")

tab3 <- dat %>%
  filter(expected_n %in% c(10000, 100000)) %>%
  group_by(scenario, model, expected_n) %>%
  summarise(comptime = mean(ftime, na.rm = T),
            rmse = sqrt(median(sqerr, na.rm = T)),
            dkl = median(dkl, na.rm = T),
            coverage = mean(coverage, na.rm = T),
            ci_width = median(int_wid, na.rm = T),
            ll = median(ll, na.rm = T),
            fit_failures = 1000 - length(unique(sim)),
            convergence_failure = sum(is.na(sqerr))
  ) %>%
  filter(model %in% c("IPP", "VA 7x7", "Lp 7x7", "VA 14x14", "Lp 14x14", "VA 21x21", "Lp 21x21", "INLA", "spNNGP"))
tab3$scenario <- factor(tab3$scenario, levels = c("SW", "WS"),
                        labels = c("S,W", "W,S"))
tab3$expected_n <- factor(tab3$expected_n, levels = c(1e+04, 1e+05),
                        labels = c("10,000", "100,000"))
tab3$fit_failures <- as.integer(tab3$fit_failures)
print(xtable(tab3, label = "tab:append_all_res_heavy_comp_settings", digits = 2,
             caption = c("Complete simulation result summaries for point patterns with E[N]=10000 and 100000.")),
      include.rownames=F, caption.placement = "top")

tab4 <- dat %>%
  filter(latent_covar_fn == "stable" & expected_n == 100) %>%
  group_by(scenario, model) %>%
  summarise(rmse = sqrt(mean(sqerr, na.rm = T)),
            dkl_mean = mean(dkl, na.rm = T),
            coverage = mean(coverage, na.rm = T),
            ci_width_mean = mean(int_wid, na.rm = T),
            ll_mean = mean(ll, na.rm = T),
            fit_failures = 1000 - length(unique(sim)),
            convergence_failure = sum(flagfit)
  ) %>%
  filter(model %in% c("VA 7x7", "Lp 7x7", "VA 14x14", "Lp 14x14"))
tab4$scenario <- factor(tab4$scenario, levels = c("SS", "SW","WS", "WW"),
                            labels = c("S,S", "S,W", "W,S", "W,W"))

print(xtable(tab4, label = "tab:append_all_res_data_poor_settings", digits = 2,
             caption = c("Complete simulation result summaries for point patterns with E[N]=100.")),
      include.rownames=F, caption.placement = "top")


################################################################################
## Data Poor (E[N]=200) Simulation Settings ## OLD RESULTS
load("data poor results n200.RDATA")

dat_poor <- rawdat %>%
  filter(expected_n == 200 &
           !scenario %in% c("S,M", "W,M") &
           model %in% c("ipp.res", "vark_S.res", "lprk_S.res", "vark_L.res", "lprk_L.res", "inla.res"))
dat_poor$model <- factor(dat_poor$model,
                         levels = c("ipp.res", "vark_S.res", "lprk_S.res", "vark_L.res", "lprk_L.res", "inla.res"),
                         labels = c("IPP", "VA 7x7", "Lp 7x7", "VA 14x14", "Lp 14x14", "INLA"))
dat_poor$scenario <- as.factor(as.character(dat_poor$scenario))
dat_poor[dat_poor$mll > 0 & !is.na(dat_poor$mll), c("slope_err", "dkl", "cp", "int_wid", "mll")] <- NA
dat_poor[is.infinite(dat_poor$dkl) & !is.na(dat_poor$dkl), c("slope_err", "dkl", "cp", "int_wid", "mll")] <- NA
tab4 <- dat_poor %>%
  group_by(scenario, model, expected_n) %>%
  summarise(comptime = mean(time, na.rm = T),
            rmse = sqrt(mean(slope_err, na.rm = T)),
            dkl = mean(dkl, na.rm = T),
            coverage = mean(cp, na.rm = T),
            ci_width = mean(int_wid, na.rm = T),
            ll = mean(mll, na.rm = T),
            fit_failures = 1000 - length(unique(sim)),
            convergence_failure = sum(is.na(slope_err))
  )

print(xtable(tab4, label = "tab:append_res_data_poor_settings", digits = 2,
             caption = c("Complete simulation result summaries for point patterns with E[N]=200.")),
      include.rownames=F, caption.placement = "top")

# Need to count the failed fits #

# Check all the metrics that are MA
apply(rawdat[ , c("slope_est", "slope_err", "dkl", "cp", "time", "mll", "int_wid")], 2,
      function(x){
        sum(is.na(x))
      }
)
# Everywhere that slope_est (slope estimates) are NA (e.g. failed) all other metrics are all NA too
apply(rawdat[is.na(rawdat$slope_est) , c("slope_err", "dkl", "cp", "time", "mll", "int_wid")], 2,
      function(x){
        sum(is.na(x))
      }
)
# These can be denoted Fit Fail
rawdat$`Fit Fail` <- rep(F, nrow(rawdat))
rawdat$`Fit Fail`[is.na(rawdat$slope_est)] <- T

# Check how many of the remaining metrics are NA
apply(rawdat[!rawdat$`Fit Fail`, c("slope_err", "dkl", "cp", "time", "mll", "int_wid")], 2,
      function(x) {
        sum(is.na(x))
      }
)

tmp <- rawdat[!rawdat$`Fit Fail`, c("slope_err", "dkl", "cp", "time", "mll", "int_wid")]
apply(tmp[is.na(tmp$int_wid), ], 2,
      function(x) {
        sum(is.na(x))
      }
)
# All the NAs in DKL and mLL are found within these
# Are all the instances where the rest are not NA still wild answers? e.g. >1e10
tmp <- tmp[is.na(tmp$int_wid), ]
tmp[!is.na(tmp$dkl), ]
# Some of these seem OK as in the SE fails but point estimate is good
tmp[!is.nan(tmp$int_wid), ] # Also seems to be a difference between NaN and NA int_wid
# What are the LLs when SE fails?
summary(tmp$mll[!is.nan(tmp$int_wid)]) # All are > 1e8
summary(tmp$dkl[!is.nan(tmp$int_wid)]) # All are NA
# What are these when SE is just an NA?
summary(tmp$mll[is.nan(tmp$int_wid)]) # Most are reasonable (max is 4e13)
summary(tmp$dkl[is.nan(tmp$int_wid)]) # There are also some infinite values here
# There are also these problems for Inf DKL
tmp[!is.finite(tmp$dkl) & !is.na(tmp$dkl), ]

# Need to label these failures
rawdat$`SE Failed` <- F
rawdat$`SE Failed`[!rawdat$`Fit Fail` & is.nan(rawdat$int_wid)] <- T
rawdat$`Bad Convergence` <- F
rawdat$`Bad Convergence`[!rawdat$`Fit Fail` & !rawdat$`SE Failed` & !is.finite(rawdat$dkl)] <- T

# Now we can aggregate the metrics based off these and report errors #
rdat <- rawdat[!rawdat$`Fit Fail` & !rawdat$`SE Failed` & !rawdat$`Bad Convergence`, ]

# Get the RMSE
tmp1 <- aggregate(rdat$slope_err,
                  by = list(rdat$model, rdat$scenario, rdat$expected_n, rdat$cov_range, rdat$lat_range),
                  FUN = function(x){sqrt(mean(x, na.rm = T))}
)
# Get the averages of the remaining metrics
tmp2 <- aggregate(rdat[,c("dkl", "cp", "int_wid", "mll", "time")],
                  by = list(rdat$model, rdat$scenario, rdat$expected_n, rdat$cov_range, rdat$lat_range),
                  FUN = function(x){mean(x, na.rm = T)}
)
# Count the error rates
tmp3 <- aggregate(rawdat[,c("Fit Fail", "SE Failed", "Bad Convergence")],
                  by = list(rawdat$model, rawdat$scenario, rawdat$expected_n, rawdat$cov_range, rawdat$lat_range),
                  FUN = function(x){sum(x)}
)

# Improve vis. for models
model <- NULL
model[grepl("lprk", tmp1$Group.1, fixed = T)] <- "Lp"
model[grepl("scampr_lp", tmp1$Group.1, fixed = T)] <- "Lp"
model[grepl("vark", tmp1$Group.1, fixed = T)] <- "VA"
model[grepl("scampr_va", tmp1$Group.1, fixed = T)] <- "VA"
model[grepl("ipp.res", tmp1$Group.1, fixed = T)] <- "IPP"
model[grepl("inla.res", tmp1$Group.1, fixed = T)] <- "INLA"

# Separately indicate basis functions
bfs <- rep(0, nrow(tmp1))
bfs[grepl("_L.res", tmp1$Group.1, fixed = T)] <- "14x14"
bfs[grepl("_M.res", tmp1$Group.1, fixed = T)] <- "10x10"
bfs[grepl("_S.res", tmp1$Group.1, fixed = T)] <- "7x7"
bfs[grepl("inla.res", tmp1$Group.1, fixed = T)] <- "k(n)"
bfs[grepl("scampr_", tmp1$Group.1, fixed = T)] <- "Dual Res."

# check that the aggregations match up
identical(1:nrow(tmp1),
          match(paste0(tmp1$Group.1, tmp1$Group.2, tmp1$Group.3, tmp1$Group.4, tmp1$Group.5),
                paste0(tmp2$Group.1, tmp2$Group.2, tmp2$Group.3, tmp2$Group.4, tmp2$Group.5)
          )
)
identical(1:nrow(tmp1),
          match(paste0(tmp1$Group.1, tmp1$Group.2, tmp1$Group.3, tmp1$Group.4, tmp1$Group.5),
                paste0(tmp3$Group.1, tmp3$Group.2, tmp3$Group.3, tmp3$Group.4, tmp3$Group.5)
          )
)

tab_append <- cbind.data.frame(tmp1$Group.2, tmp1$Group.3, tmp1$Group.4, tmp1$Group.5, tmp1$x, tmp2$dkl, tmp2$cp, tmp2$int_wid, tmp2$mll, tmp2$time, model, bfs, tmp3$`Fit Fail`, tmp3$`SE Failed`, tmp3$`Bad Convergence`)
# colnames(tab_append) <- c("Scenario", "E[N]", "Covariate Range", "Latent Range", "RMSE b", "KL Div.", "Coverage Prob b", "Interval Width", "lnL", "Comp. Time", "Model", "k", "Fit Failure", "SE Failure", "Poor Convergence")
colnames(tab_append) <- c("Scenario", "$\\mathbb{E}\\left[N\\right]$", "Covariate Range", "Latent Range", "RMSE $\\beta_{1}$", "KL Div. $\\lambda||\\hat{\\lambda}$", "$95\\%$ Coverage Prob. $\\hat{\\beta_{1}}$", "CI $\\beta_{1}$ Width", "$\\ell\\left(\\boldsymbol{\\beta}\\right)$", "Comp. Time", "Model", "$k$", "Fit Failure", "$\\hat{\\textrm{SE}}$ Failure", "Poor Convergence")
# library(xtable)
# print(xtable(tab_append, caption = c("Complete simulation result summaries for Chapter 2.", "Simulation Results Summaries"), label = "tab:append_chpt2_all_results_summary"), caption.placement = "top", include.rownames=FALSE)

# Does scenario breakup of the table fit a single page? YES

## Supplementary Material - Table 1:
tab_append_SS <- tab_append[tab_append$Scenario == "S,S", c(-1, -3, -4)]
print(xtable(tab_append_SS, caption = c("Complete simulation result summaries for scenario S,S --- i.e. both the covariate and latent field range of effect are approximately 30 units.",
                                        "Simulation Results Summaries: Scenario S,S"),
             label = "tab:append_chpt2_sim_res_SS", digits = c(0,0,2,2,2,2,2,2,0,0,0,0,0)), caption.placement = "top", include.rownames=FALSE, sanitize.colnames.function = identity)

## Supplementary Material - Table 2:
tab_append_SW <- tab_append[tab_append$Scenario == "S,W", c(-1, -3, -4)]
print(xtable(tab_append_SW, caption = c("Complete simulation result summaries for scenario S,W --- i.e. the range of effects are approximately 30 units for the covariate and approximately 5 units for the latent field.",
                                        "Simulation Results Summaries: Scenario S,W"),
             label = "tab:append_chpt2_sim_res_SW", digits = c(0,0,2,2,2,2,2,2,0,0,0,0,0)), caption.placement = "top", include.rownames=FALSE, sanitize.colnames.function = identity)

## Supplementary Material - Table 3:
tab_append_WS <- tab_append[tab_append$Scenario == "W,S", c(-1, -3, -4)]
print(xtable(tab_append_SW, caption = c("Complete simulation result summaries for scenario W,S --- i.e. the range of effects are approximately 5 units for the covariate and approximately 30 units for the latent field.",
                                        "Simulation Results Summaries: Scenario W,S"),
             label = "tab:append_chpt2_sim_res_WS", digits = c(0,0,2,2,2,2,2,2,0,0,0,0,0)), caption.placement = "top", include.rownames=FALSE, sanitize.colnames.function = identity)

## Supplementary Material - Table 4:
tab_append_WW <- tab_append[tab_append$Scenario == "W,W", c(-1, -3, -4)]
print(xtable(tab_append_SW, caption = c("Complete simulation result summaries for scenario W,W --- i.e. both the covariate and latent field range of effect are approximately 5 units.",
                                        "Simulation Results Summaries: Scenario W,W"),
             label = "tab:append_chpt2_sim_res_WS", digits = c(0,0,2,2,2,2,2,2,0,0,0,0,0)), caption.placement = "top", include.rownames=FALSE, sanitize.colnames.function = identity)
