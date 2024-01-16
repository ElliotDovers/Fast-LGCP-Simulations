home.wd <- getwd()

# library(scampr)
# library(TMB)
library(spatstat)
library(viridis)

source("sim_lgcp_pp.R")
source("vec2mat.r")

job = 27#sample(1:100, 1)
plot.res <- 500

# Data Simulation
sim.list <- sim_lgcp_pp(Intercept = -4.867439, Slope = 1.25, latent.range = 10, rseed = job,
            cov.function = "stable", covariate_range = 5)
wiggly_sim <- sim.list$dframe
sim.list <-sim_lgcp_pp(Intercept = -4.867439, Slope = 1.25, latent.range = 30, rseed = job,
                       cov.function = "stable", covariate_range = 20)
smooth_sim <- sim.list$dframe

job = 47
sim.list <- sim_lgcp_pp(Intercept = -4.867439, Slope = 1.25, latent.range = 10, rseed = job,
                        cov.function = "exp", covariate_range = 5)
wiggly_sim_exp <- sim.list$dframe
sim.list <-sim_lgcp_pp(Intercept = -4.867439, Slope = 1.25, latent.range = 30, rseed = job,
                       cov.function = "exp", covariate_range = 20)
smooth_sim_exp <- sim.list$dframe

quad <- smooth_sim[smooth_sim$pres == 0, c("x", "y")]
cov_zlims <- range(c(smooth_sim$X[smooth_sim$pres == 0], wiggly_sim$X[wiggly_sim$pres == 0]))
lat_zlims <- range(c(smooth_sim$xi[smooth_sim$pres == 0], wiggly_sim$xi[wiggly_sim$pres == 0]))
lat_zlims_exp <- range(c(smooth_sim_exp$xi[smooth_sim_exp$pres == 0], wiggly_sim_exp$xi[wiggly_sim_exp$pres == 0]))

#######################################################
# Simulation Set Up Plot - stable covariance function #
#######################################################

# 2 x 2 Simulation Design
png(filename = paste0(getwd(), "/Results/Figures/Figure 2.png"), res = plot.res, width = 4.5 * plot.res, height = 4.5 * plot.res)
layout(mat = matrix(c(1,2,3,4,5,6,7,8,9), nrow = 3, ncol = 3, byrow = TRUE), heights = rep(1/3,3), widths = rep(1/3,3))
par(mar = c(0, 2.1, 2.1, 0))
# Dummy plot
plot(1, type = "n", axes=FALSE, xlab="", ylab="")

# Smooth Latent Field
image(x = sort(unique(quad$x)),
      y = sort(unique(quad$y)),
      z = vec2mat(smooth_sim$xi[smooth_sim$pres == 0], quad$x, quad$y),
      col = plasma(100), xlab = "", ylab = "", main = "", bty = 'n',
      axes = F,
)
# title(main = expression(paste(xi, "(s) : Smooth")))
par(xpd=T)
text(x = mean(smooth_sim$x[smooth_sim$pres == 0]), y = max(smooth_sim$y[smooth_sim$pres == 0]) + 10, labels = expression(paste(xi, "(s) : Smooth")), cex = 1.5)
par(xpd=F)

# Wiggly Latent Field
image(x = sort(unique(quad$x)),
      y = sort(unique(quad$y)),
      z = vec2mat(wiggly_sim$xi[wiggly_sim$pres == 0], quad$x, quad$y),
      col = plasma(100), xlab = "", ylab = "", main = "", bty = 'n',
      axes = F,
)
# title(main = expression(paste(xi, "(s) : Wiggly")))
par(xpd=T)
text(x = mean(wiggly_sim$x[wiggly_sim$pres == 0]), y = max(wiggly_sim$y[wiggly_sim$pres == 0]) + 10, labels = expression(paste(xi, "(s) : Wiggly")), cex = 1.5)
par(xpd=F)

# Smooth Covariate Field
image(x = sort(unique(quad$x)),
      y = sort(unique(quad$y)),
      z = vec2mat(smooth_sim$X[smooth_sim$pres == 0], quad$x, quad$y),
      col = topo.colors(100), xlab = "", ylab = "", main = "", bty = 'n',
      axes = F,
)
par(xpd=T)
text(x = min(smooth_sim$x[smooth_sim$pres == 0]) - 10, y = mean(smooth_sim$y[smooth_sim$pres == 0]), labels = "X(s) : Smooth", srt = 90, cex = 1.5)
par(xpd=F)

# Smooth, Smooth Scenario Text
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
text(x = 1, labels = "S,S", cex = 2.5)

# Smooth, Wiggly Scenario Text
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
text(x = 1, labels = "S,W", cex = 2.5)

# Wiggly Covariate Field
image(x = sort(unique(quad$x)),
      y = sort(unique(quad$y)),
      z = vec2mat(wiggly_sim$X[wiggly_sim$pres == 0], quad$x, quad$y),
      col = topo.colors(100), xlab = "", ylab = "", main = "", bty = 'n',
      axes = F, #zlim = cov_zlims
)
par(xpd=T)
text(x = min(wiggly_sim$x[wiggly_sim$pres == 0]) - 10, y = mean(wiggly_sim$y[wiggly_sim$pres == 0]), labels = "X(s) : Wiggly", srt = 90, cex = 1.5)
par(xpd=F)

# Wiggly, Smooth Scenario Text
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
text(x = 1, labels = "W,S", cex = 2.5)

# Wiggly, Wiggly Scenario Text
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
text(x = 1, labels = "W,W", cex = 2.5)
dev.off()

############################################################
# Simulation Set Up Plot - exponential covariance function #
############################################################

# 2 x 2 Simulation Design
png(filename = paste0(getwd(), "/Results/Figures/Appendix Figure Ax.png"), res = plot.res, width = 4.5 * plot.res, height = 4.5 * plot.res)
layout(mat = matrix(c(1,2,3,4,5,6,7,8,9), nrow = 3, ncol = 3, byrow = TRUE), heights = rep(1/3,3), widths = rep(1/3,3))
par(mar = c(0, 2.1, 2.1, 0))
# Dummy plot
plot(1, type = "n", axes=FALSE, xlab="", ylab="")

# Smooth Latent Field
image(x = sort(unique(quad$x)),
      y = sort(unique(quad$y)),
      z = vec2mat(smooth_sim_exp$xi[smooth_sim_exp$pres == 0], quad$x, quad$y),
      col = plasma(100), xlab = "", ylab = "", main = "", bty = 'n',
      axes = F,
)
# title(main = expression(paste(xi, "(s) : Smooth")))
par(xpd=T)
text(x = mean(smooth_sim_exp$x[smooth_sim_exp$pres == 0]), y = max(smooth_sim_exp$y[smooth_sim_exp$pres == 0]) + 10, labels = expression(paste(xi[exp], "(s) : Long")), cex = 1.5)
par(xpd=F)

# Wiggly Latent Field
image(x = sort(unique(quad$x)),
      y = sort(unique(quad$y)),
      z = vec2mat(wiggly_sim_exp$xi[wiggly_sim_exp$pres == 0], quad$x, quad$y),
      col = plasma(100), xlab = "", ylab = "", main = "", bty = 'n',
      axes = F,
)
# title(main = expression(paste(xi, "(s) : Wiggly")))
par(xpd=T)
text(x = mean(wiggly_sim_exp$x[wiggly_sim_exp$pres == 0]), y = max(wiggly_sim_exp$y[wiggly_sim_exp$pres == 0]) + 10, labels = expression(paste(xi[exp], "(s) : Short")), cex = 1.5)
par(xpd=F)

# Smooth Covariate Field
image(x = sort(unique(quad$x)),
      y = sort(unique(quad$y)),
      z = vec2mat(smooth_sim_exp$X[smooth_sim_exp$pres == 0], quad$x, quad$y),
      col = topo.colors(100), xlab = "", ylab = "", main = "", bty = 'n',
      axes = F,
)
par(xpd=T)
text(x = min(smooth_sim_exp$x[smooth_sim_exp$pres == 0]) - 10, y = mean(smooth_sim_exp$y[smooth_sim_exp$pres == 0]), labels = "X(s) : Long", srt = 90, cex = 1.5)
par(xpd=F)

# Smooth, Smooth Scenario Text
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
text(x = 1, labels = "Lg,Lg", cex = 2.5)

# Smooth, Wiggly Scenario Text
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
text(x = 1, labels = "Lg,St", cex = 2.5)

# Wiggly Covariate Field
image(x = sort(unique(quad$x)),
      y = sort(unique(quad$y)),
      z = vec2mat(wiggly_sim_exp$X[wiggly_sim_exp$pres == 0], quad$x, quad$y),
      col = topo.colors(100), xlab = "", ylab = "", main = "", bty = 'n',
      axes = F, #zlim = cov_zlims
)
par(xpd=T)
text(x = min(wiggly_sim_exp$x[wiggly_sim_exp$pres == 0]) - 10, y = mean(wiggly_sim_exp$y[wiggly_sim_exp$pres == 0]), labels = "X(s) : Short", srt = 90, cex = 1.5)
par(xpd=F)

# Wiggly, Smooth Scenario Text
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
text(x = 1, labels = "St,Lg", cex = 2.5)

# Wiggly, Wiggly Scenario Text
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
text(x = 1, labels = "St,St", cex = 2.5)
dev.off()


######################################################
# Simulation Set Up Plot - Both covariance functions #
######################################################

# 2 x 2 Simulation Design
png(filename = paste0(getwd(), "/Results/Figures/Appendix Figure 1.png"), res = plot.res, width = 5.5 * plot.res, height = 5.5 * plot.res)
layout(mat = matrix(c(1,2,3,4,5,6,7,8,9,10,11,12), nrow = 4, ncol = 3, byrow = TRUE), heights = rep(1/4,4), widths = rep(1/3,3))
par(mar = c(0, 2.1, 2.1, 0))
# Dummy plot
plot(1, type = "n", axes=FALSE, xlab="", ylab="")

# Smooth Latent Field
image(x = sort(unique(quad$x)),
      y = sort(unique(quad$y)),
      z = vec2mat(smooth_sim$xi[smooth_sim$pres == 0], quad$x, quad$y),
      col = plasma(100), xlab = "", ylab = "", main = "", bty = 'n',
      axes = F,
)
# title(main = expression(paste(xi, "(s) : Smooth")))
par(xpd=T)
text(x = mean(smooth_sim$x[smooth_sim$pres == 0]), y = max(smooth_sim$y[smooth_sim$pres == 0]) + 10, labels = expression(paste(xi, "(s) : Long")), cex = 1.5)
par(xpd=F)

# Wiggly Latent Field
image(x = sort(unique(quad$x)),
      y = sort(unique(quad$y)),
      z = vec2mat(wiggly_sim$xi[wiggly_sim$pres == 0], quad$x, quad$y),
      col = plasma(100), xlab = "", ylab = "", main = "", bty = 'n',
      axes = F,
)
# title(main = expression(paste(xi, "(s) : Wiggly")))
par(xpd=T)
text(x = mean(wiggly_sim$x[wiggly_sim$pres == 0]), y = max(wiggly_sim$y[wiggly_sim$pres == 0]) + 10, labels = expression(paste(xi, "(s) : Short")), cex = 1.5)
par(xpd=F)

# Dummy plot
plot(1, type = "n", axes=FALSE, xlab="", ylab="")

# Smooth Latent Field - Exponential Covariance
image(x = sort(unique(quad$x)),
      y = sort(unique(quad$y)),
      z = vec2mat(smooth_sim_exp$xi[smooth_sim_exp$pres == 0], quad$x, quad$y),
      col = plasma(100), xlab = "", ylab = "", main = "", bty = 'n',
      axes = F,
)
# title(main = expression(paste(xi, "(s) : Smooth")))
par(xpd=T)
text(x = mean(smooth_sim_exp$x[smooth_sim_exp$pres == 0]), y = max(smooth_sim_exp$y[smooth_sim_exp$pres == 0]) + 10, labels = expression(paste(xi[exp], "(s) : Long")), cex = 1.5)
par(xpd=F)

# Wiggly Latent Field - Exponential Covariance
image(x = sort(unique(quad$x)),
      y = sort(unique(quad$y)),
      z = vec2mat(wiggly_sim_exp$xi[wiggly_sim_exp$pres == 0], quad$x, quad$y),
      col = plasma(100), xlab = "", ylab = "", main = "", bty = 'n',
      axes = F,
)
# title(main = expression(paste(xi, "(s) : Wiggly")))
par(xpd=T)
text(x = mean(wiggly_sim_exp$x[wiggly_sim_exp$pres == 0]), y = max(wiggly_sim_exp$y[wiggly_sim_exp$pres == 0]) + 10, labels = expression(paste(xi[exp], "(s) : Short")), cex = 1.5)
par(xpd=F)

# Smooth Covariate Field
image(x = sort(unique(quad$x)),
      y = sort(unique(quad$y)),
      z = vec2mat(smooth_sim$X[smooth_sim$pres == 0], quad$x, quad$y),
      col = topo.colors(100), xlab = "", ylab = "", main = "", bty = 'n',
      axes = F,
)
par(xpd=T)
text(x = min(smooth_sim$x[smooth_sim$pres == 0]) - 10, y = mean(smooth_sim$y[smooth_sim$pres == 0]), labels = "X(s) : Long", srt = 90, cex = 1.5)
par(xpd=F)

# Smooth, Smooth Scenario Text
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
text(x = 1, labels = "Long,Long", cex = 2.5)

# Smooth, Wiggly Scenario Text
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
text(x = 1, labels = "Long,Short", cex = 2.5)

# # Smooth, Smooth Scenario Text
# plot(1, type = "n", axes=FALSE, xlab="", ylab="")
# text(x = 1, labels = expression('S,S'[exp]), cex = 2.5)
# 
# # Smooth, Wiggly Scenario Text
# plot(1, type = "n", axes=FALSE, xlab="", ylab="")
# text(x = 1, labels = expression('S,W'[exp]), cex = 2.5)

# Wiggly Covariate Field
image(x = sort(unique(quad$x)),
      y = sort(unique(quad$y)),
      z = vec2mat(wiggly_sim$X[wiggly_sim$pres == 0], quad$x, quad$y),
      col = topo.colors(100), xlab = "", ylab = "", main = "", bty = 'n',
      axes = F, #zlim = cov_zlims
)
par(xpd=T)
text(x = min(wiggly_sim$x[wiggly_sim$pres == 0]) - 10, y = mean(wiggly_sim$y[wiggly_sim$pres == 0]), labels = "X(s) : Short", srt = 90, cex = 1.5)
par(xpd=F)

# Wiggly, Smooth Scenario Text
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
text(x = 1, labels = "Short,Long", cex = 2.5)

# Wiggly, Wiggly Scenario Text
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
text(x = 1, labels = "Short,Short", cex = 2.5)

# # Wiggly, Smooth Scenario Text
# plot(1, type = "n", axes=FALSE, xlab="", ylab="")
# text(x = 1, labels = expression('W,S'[exp]), cex = 2.5)
# 
# # Wiggly, Wiggly Scenario Text
# plot(1, type = "n", axes=FALSE, xlab="", ylab="")
# text(x = 1, labels = expression('W,W'[exp]), cex = 2.5)
dev.off()

### Example 2D Basis Function Plots
library(fields)
s <- seq(0, 100, by = 0.1)
bf.nodes <- seq(0, 100, by = 25)
radius <- diff(bf.nodes)[1] * sqrt(2)
dist2bfs <- rdist(s, bf.nodes)
bfs <- (1 - (dist2bfs / radius)^2)^2
bfs[dist2bfs > radius] <- 0
us <- c(0.5, 1.2, 0.1, -3, 0.2)

png(filename = paste0(getwd(), "/Results/Figures/Chpt2_basis_fun_eg.png"), res = 500, width = 5 * 500, height = 5 * 500)
# par(mfrow = c(2, 2))
par(mfrow=c(2,2),mar=c(3,3,2.1,1),mgp=c(1.75,0.75,0))
# a)
plot(s, bfs[,3], col = 1, ylab = "Z(s)", xlab = "s", type = "l", main = "a)", bty = "n")
# lines(c(50, 50 + radius), c(1, 1), lty = "dashed")
arrows(x0 = 50, y0 = 1, x1 = 50 + radius, y1 = 1, code = 3, angle = 45, length = 0.05, col = 1, lwd = 1.5)
text(50 + (radius/2), y = 0.9, labels = expression(varphi), col = 1)
# b)
# par(mar = c(5.1, 0, 4.1, 2.1))
par(mar = c(3, 0, 2, 1))
plot(s, bfs[,3], col = 1, ylab = "", xlab = "s'", type = "l", main = "b)", bty = "n", xaxt = "n", yaxt = "n")
lines(s, bfs[,2], col = 2)
lines(c(50, 50), c(1, 0), lty = "dashed")
lines(c(25, 25), c(1, 0), lty = "dashed", col = 2)
# axis(1, at = c(25, 50), labels = c(expression("s'"[2]*" = 25"), expression("s'"[3]*" = 50")), tick = F)
axis(1, at = c(25, 50), labels = c("25", "50"), tick = F, col = "red")
# c)
# par(mar = c(5.1, 0, 4.1, 0.1))
plot(s, bfs[,3], col = 1, ylab = "", xlab = "s", type = "l", main = "c)", bty = "n", xlim = c(-10, 110), ylim = c(0, 1.5), yaxt = "n")
for (i in 1:ncol(bfs)) {
   lines(s, bfs[,i], col = i)
}
# text(bf.nodes, y = 1.2, labels = c(expression("u"[1]*" = 0.5"), expression("u"[2]*" = 1.2"), expression("u"[3]*" = 0.1"), expression("u"[4]*" = -3"), expression("u"[5]*" = 0.2")))
text(bf.nodes[1], y = 1.2, labels = expression("u"[1]*" = 0.5"), col = 1, cex = 0.8)
text(bf.nodes[2], y = 1.2, labels = expression("u"[2]*" = 1.2"), col = 2, cex = 0.8)
text(bf.nodes[3], y = 1.2, labels = expression("u"[3]*" = 0.1"), col = 3, cex = 0.8)
text(bf.nodes[4], y = 1.2, labels = expression("u"[4]*" = -3"), col = 4, cex = 0.8)
text(bf.nodes[5], y = 1.2, labels = expression("u"[5]*" = 0.2"), col = 5, cex = 0.8)
# d)
# par(mar = c(5.1, 4.1, 4.1, 2.1))
par(mar = c(3,3,2.1,1))
plot(s, bfs %*% us, bty = "n", type = "l", ylab = expression(xi%~~%bold(Z)*"(s)"*bold(u)), main = "d)")
par(mfrow = c(1, 1))
dev.off()

#################################################################################################

# WIDE FORMAT FOR PAINT EDITING
png(filename = paste0(getwd(), "/Results/Figures/Chpt2_sim_setup_fig1.png"), res = plot.res, width = 6 * plot.res, height = 4 * plot.res)
par(oma = c(0, 0, 0, 0))
layout(mat = matrix(c(1,2,5,3,4,5), nrow = 2, ncol = 3, byrow = TRUE), heights = c(0.5,0.5), widths = c(0.3,0.3,0.4))
par(mar = c(0.1, 2.1, 7.1, 0.1))
image(x = sort(unique(quad$x)),
      y = sort(unique(quad$y)),
      z = vec2mat(smooth_sim$X[smooth_sim$pres == 0], quad$x, quad$y),
      col = topo.colors(100), xlab = "", ylab = "", main = "", bty = 'n',
      axes = F, #zlim = cov_zlims
)
par(xpd=T)
text(x = min(wiggly_sim$x[wiggly_sim$pres == 0]) - 10, y = mean(wiggly_sim$y[wiggly_sim$pres == 0]), labels = "Smooth", srt = 90, cex = 2)
text(x = mean(wiggly_sim$x[wiggly_sim$pres == 0]), y = max(wiggly_sim$y[wiggly_sim$pres == 0]) + 10, labels = "Covariate", cex = 2)
par(xpd=F)
par(mar = c(0.1, 0.1, 7.1, 2.1))
image(x = sort(unique(quad$x)),
      y = sort(unique(quad$y)),
      z = vec2mat(smooth_sim$xi[smooth_sim$pres == 0], quad$x, quad$y),
      col = topo.colors(100), xlab = "", ylab = "", main = "", bty = 'n',
      axes = F, #zlim = lat_zlims
)
par(xpd=T)
text(x = mean(wiggly_sim$x[wiggly_sim$pres == 0]), y = max(wiggly_sim$y[wiggly_sim$pres == 0]) + 10, labels = "Latent Field", cex = 2)
par(xpd=F)
par(mar = c(7.1, 2.1, 0.1, 0))
image(x = sort(unique(quad$x)),
      y = sort(unique(quad$y)),
      z = vec2mat(wiggly_sim$X[wiggly_sim$pres == 0], quad$x, quad$y),
      col = topo.colors(100), xlab = "", ylab = "", main = "", bty = 'n',
      axes = F, #zlim = cov_zlims
)
par(xpd=T)
text(x = min(wiggly_sim$x[wiggly_sim$pres == 0]) - 10, y = mean(wiggly_sim$y[wiggly_sim$pres == 0]), labels = "Wiggly", srt = 90, cex = 2)
par(xpd=F)
par(mar = c(7.1, 0.1, 0.1, 2.1))
image(x = sort(unique(quad$x)),
      y = sort(unique(quad$y)),
      z = vec2mat(wiggly_sim$xi[wiggly_sim$pres == 0], quad$x, quad$y),
      col = topo.colors(100), xlab = "", ylab = "", main = "", bty = 'n',
      axes = F, #zlim = lat_zlims
)
dev.off()

# NON-WIDE FORMAT
png(filename = paste0(getwd(), "/Results/Figures/Chpt2_sim_setup_fig1.png"), res = plot.res, width = 6 * plot.res, height = 4 * plot.res)
par(oma = c(0, 0, 0, 0))
layout(mat = matrix(c(1,2,3,4), nrow = 2, ncol = 2, byrow = TRUE), heights = c(0.5,0.5), widths = c(0.575,0.425))
par(mar = c(0, 7.2, 2.1, 0))
image(x = sort(unique(quad$x)),
      y = sort(unique(quad$y)),
      z = vec2mat(smooth_sim$X[smooth_sim$pres == 0], quad$x, quad$y),
      col = topo.colors(100), xlab = "", ylab = "", main = "", bty = 'n',
      axes = F, #zlim = cov_zlims
)
par(xpd=T)
text(x = min(wiggly_sim$x[wiggly_sim$pres == 0]) - 10, y = mean(wiggly_sim$y[wiggly_sim$pres == 0]), labels = "Smooth", srt = 90, cex = 2)
text(x = mean(wiggly_sim$x[wiggly_sim$pres == 0]), y = max(wiggly_sim$y[wiggly_sim$pres == 0]) + 10, labels = "Covariate", cex = 2)
par(xpd=F)
par(mar = c(0, 2.1, 2.1, 0))
image(x = sort(unique(quad$x)),
      y = sort(unique(quad$y)),
      z = vec2mat(smooth_sim$xi[smooth_sim$pres == 0], quad$x, quad$y),
      col = topo.colors(100), xlab = "", ylab = "", main = "", bty = 'n',
      axes = F, #zlim = lat_zlims
      )
par(xpd=T)
text(x = mean(wiggly_sim$x[wiggly_sim$pres == 0]), y = max(wiggly_sim$y[wiggly_sim$pres == 0]) + 10, labels = "Latent Field", cex = 2)
par(xpd=F)
par(mar = c(0, 7.2, 2.1, 0))
image(x = sort(unique(quad$x)),
      y = sort(unique(quad$y)),
      z = vec2mat(wiggly_sim$X[wiggly_sim$pres == 0], quad$x, quad$y),
      col = topo.colors(100), xlab = "", ylab = "", main = "", bty = 'n',
      axes = F, #zlim = cov_zlims
)
par(xpd=T)
text(x = min(wiggly_sim$x[wiggly_sim$pres == 0]) - 10, y = mean(wiggly_sim$y[wiggly_sim$pres == 0]), labels = "Wiggly", srt = 90, cex = 2)
par(xpd=F)
par(mar = c(0, 2.1, 2.1, 0))
image(x = sort(unique(quad$x)),
      y = sort(unique(quad$y)),
      z = vec2mat(wiggly_sim$xi[wiggly_sim$pres == 0], quad$x, quad$y),
      col = topo.colors(100), xlab = "", ylab = "", main = "", bty = 'n',
      axes = F, #zlim = lat_zlims
)
dev.off()
