# 1D explainer of why we introduce bias
int <- c(1,2,3)
s <- seq(1,100, by = 0.1)
X <- exp(-((s-50)/25)^2)
xi <- cos(s/5)
lnlambda1 <- int[1] + X + xi
lnlambda2 <- int[2] + X + xi
lnlambda3 <- int[3] + X + xi
# quad weights
w <- rep(0.1, length(s))
w[c(1, length(w))] <- 0.05
lims <- range(lnlambda1,lnlambda2,lnlambda3)

par(mfrow=c(2,2))
plot(s,X,type = "l")
plot(s,xi,type = "l")
plot(s,lnlambda1,type = "l",ylim = lims)
lines(s,lnlambda2,col = "grey")
lines(s,lnlambda3,col = "lightgrey")
plot(s,exp(lnlambda1),type = "l",ylim = exp(lims))
lines(s,exp(lnlambda2),col = "grey")
lines(s,exp(lnlambda3),col = "lightgrey")
par(mfrow=c(1,1))

sim1D <- function(lambda.int) {
  s <- seq(1,100, by = 0.1)
  lambda.at.quad <- exp(lambda.int + exp(-((s-50)/25)^2) + cos(s/5))
  pts.max <- runif(sum(w) * max(lambda), 1, 100)
  lambda.at.pts <- exp(lambda.int + exp(-((pts.max-50)/25)^2) + cos(pts.max/5))
  thinning.prob <- lambda.at.pts / max(lambda)
  thinning.prob[thinning.prob > 1] <- 1
  keep <- rbinom(length(thinning.prob), 1, thinning.prob)
  pts <- pts.max[keep == 1]
  return(pts)
}

p1 <- list(sim1D(1), sim1D(1), sim1D(1))

newlims <- exp(lims)
newlims[2] <- newlims[2] + 3*newlims[2]
par(mfrow=c(2,2))
plot(s,exp(lnlambda1),type = "l",ylim = newlims)
points(p1[[1]], rep(exp(lims)[2], length(p1[[1]])), pch = "1", col = "black")
points(p1[[2]], rep(exp(lims)[2], length(p1[[2]])), pch = "2", col = "black")
lines(s,exp(lnlambda3),col = "lightgrey")
par(mfrow=c(1,1))

w * exp(lnlambda1)
spatstat.random::rpoispp(exp(lnlambda1))
