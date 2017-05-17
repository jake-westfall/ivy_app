library("lavaan") # for SEMs
library("MASS") # for mvrnorm()
library("corpcor") # for pcor2cor()
library("mgcv") # for gam()
library("mvtnorm") # for rmvnorm() and pmvt()
library("zoo") # for rollapply()

# inverse logit function
il <- function(x) 1/(1 + exp(-x))

path <- "/Users/Jake/Desktop/Dropbox/Dropbox/SIC/"

# using lavaan’s simulateData() -------------------------------------------

n <- 100000
alpha <- .4
delta <- .5
bx <- .5
bt <- .5

cor(simulateData(paste0("
                        y ~ ",bx,"*x + ",bt,"*t
                        x ~~ ",delta,"*t
                        t =~ 1*z
                        z ~~ (1/",alpha," - 1)*z
                        "),
                 std.lv=TRUE, sample.nobs=n))

# seems to work, but inputs must be in terms of slopes. want correlations

# simulation function -----------------------------------------------------

semPower <- function(n, alpha1, alpha2, delta, rho1_p, rho2_p){
  # compute the implied set of 3 simple correlations
  delta_p <- delta*sqrt(prod(c(rho1_p, rho2_p)^2 - 1)) - prod(c(rho1_p, rho2_p))
  rp <- rbind(c(1, rho1_p, rho2_p),
              c(rho1_p, 1, delta_p),
              c(rho2_p, delta_p, 1))
  r <- suppressWarnings(pcor2cor(rp))
  # throw error if illegal partial correlation matrix supplied
  if(any(is.nan(r)))
    stop(paste0("Illegal partial correlation matrix:\n
                rho1_p = ",rho1_p,", rho2_p = ",rho2_p,", delta_p = ",delta_p))
  
  # make data
  dat <- data.frame(mvrnorm(n, mu=numeric(3), Sigma=r))
  names(dat) <- c("y","t1","t2")
  dat <- within(dat, {
    x1 <- t1 + rnorm(n, sd=sqrt((1-alpha1)/alpha1))
    x2 <- t2 + rnorm(n, sd=sqrt((1-alpha2)/alpha2))
  })
  
  # fit model
  z <- "NULL"
  mod <- suppressWarnings(sem(paste0("
                                     y ~ t1 + t2
                                     t1 =~ 1*x1
                                     x1 ~~ ",(1-alpha1)*var(dat$x1),"*x1
                                     t2 =~ 1*x2
                                     x2 ~~ ",(1-alpha2)*var(dat$x2),"*x2
                                     "), std.lv=TRUE, data=dat))
  
  # return stats
  if(lavInspect(mod, "converged")){
    z <- try(lavInspect(mod, "est")$beta["y","t1"]/sqrt(diag(vcov(mod)))["y~t1"],
             silent=TRUE)
    if(class(z) == "try-error") z <- NA
    return(as.numeric(abs(z) > qnorm(.975)))
  } else return(NA)
}

# test sim ----------------------------------------------------------------

x <- seq(5,9,.1)

# 1000 sims takes about 54 seconds

png(paste0(path,"power_test.png"),
    units="in", height=5, width=6, res=200, pointsize=15)

# rho1_p = .5
n <- runif(2000, min=log2(32), max=log2(512))
results <- sapply(n, function(x){
  semPower(n=round(2^x), alpha1=1, alpha2=.4, delta=.5, rho1_p=.5, rho2_p=.5)
})
fit <- gam(results ~ s(n), family=binomial)
pred <- predict(fit, se.fit=TRUE, newdata=data.frame(n=x))
plot(y=0:1, x=c(5,9), cex=0, xlab="Sample size", xaxt="n",
     ylab="Probability of rejecting null")
axis(side=1, at=5:9, labels=2^(5:9))
grid(lty=1, col="gray")
polygon(x=c(x, rev(x), x[1]), col=rgb(0,0,0,.25), border=FALSE,
        y=il(c(pred$fit+pred$se.fit, rev(pred$fit-pred$se.fit), pred$fit[1]+pred$se.fit[1])))
lines(x=x, y=il(pred$fit), lwd=2)

# rho1_p = .3
n <- runif(2000, min=log2(32), max=log2(512))
results <- sapply(n, function(x){
  semPower(n=round(2^x), alpha1=1, alpha2=.4, delta=.5, rho1_p=.3, rho2_p=.5)
})
fit <- gam(results ~ s(n), family=binomial)
pred <- predict(fit, se.fit=TRUE, newdata=data.frame(n=x))
polygon(x=c(x, rev(x), x[1]), col=rgb(0,0,0,.25), border=FALSE,
        y=il(c(pred$fit+pred$se.fit, rev(pred$fit-pred$se.fit), pred$fit[1]+pred$se.fit[1])))
lines(x=x, y=il(pred$fit), lwd=2)

# rho1_p = .1
n <- runif(2000, min=log2(32), max=log2(512))
results <- sapply(n, function(x){
  semPower(n=round(2^x), alpha1=1, alpha2=.4, delta=.5, rho1_p=.1, rho2_p=.5)
})
fit <- gam(results ~ s(n), family=binomial)
pred <- predict(fit, se.fit=TRUE, newdata=data.frame(n=x))
polygon(x=c(x, rev(x), x[1]), col=rgb(0,0,0,.25), border=FALSE,
        y=il(c(pred$fit+pred$se.fit, rev(pred$fit-pred$se.fit), pred$fit[1]+pred$se.fit[1])))
lines(x=x, y=il(pred$fit), lwd=2)

# rho1_p = 0
n <- runif(2000, min=log2(32), max=log2(512))
results <- sapply(n, function(x){
  semPower(n=round(2^x), alpha1=1, alpha2=.4, delta=.5, rho1_p=0, rho2_p=.5)
})
fit <- gam(results ~ s(n), family=binomial)
pred <- predict(fit, se.fit=TRUE, newdata=data.frame(n=x))
polygon(x=c(x, rev(x), x[1]), col=rgb(0,0,0,.25), border=FALSE,
        y=il(c(pred$fit+pred$se.fit, rev(pred$fit-pred$se.fit), pred$fit[1]+pred$se.fit[1])))
lines(x=x, y=il(pred$fit), lwd=2)

dev.off()

# full simulation ---------------------------------------------------------

# length of n = number of iterations per parameter combination
n <- seq(xaxis[1], xaxis[7], length=2000)
xaxis <- log(c(30,50,75,120,200,300,500))
xpred <- seq(xaxis[1], xaxis[7], length=40)

# sim parameters to iterate over
params <- within(expand.grid(
  rho1_p=c(0,.1,.3,.5),
  alpha2=c(1, .8, .4),
  delta=c(.3, .5, .7)), {
    rho2_p <- delta
  })

# function to collect and name results
iterate <- function(alpha2, delta, rho1_p, rho2_p){
  results <- list(sapply(n, function(x){
    semPower(n=round(exp(x)), alpha1=1, alpha2=alpha2, delta=delta,
             rho1_p=rho1_p, rho2_p=rho2_p)
  }))
  names(results) <- paste0(colnames(params), "=",
    c(rho1_p, alpha2, delta, rho2_p), collapse=" ")
  return(results)
}

# run sim and save results
system.time({
  powSim1 <- mapply(iterate, alpha2=params$alpha2, delta=params$delta,
                   rho1_p=params$rho1_p, rho2_p=params$rho2_p)
}) # took 70 minutes
saveRDS(powSim1, file=paste0(path,"powSimResults1.rds"))

# do 13000 more, for 15000 in total
n <- seq(xaxis[1], xaxis[7], length=13000)
system.time({
  powSim2 <- mapply(iterate, alpha2=params$alpha2, delta=params$delta,
                   rho1_p=params$rho1_p, rho2_p=params$rho2_p)
}) # took 5.5 hours
saveRDS(powSim2, file=paste0(path,"powSimResults2.rds"))

# combine results
powSim1 <- readRDS(paste0(path,"powSimResults1.rds"))
powSim2 <- readRDS(paste0(path,"powSimResults2.rds"))
powSim <- lapply(seq(length(powSim1)), function(i){
  c(powSim1[[i]], powSim2[[i]])
})
names(powSim) <- names(powSim1)
n <- c(seq(xaxis[1], xaxis[7], length=2000), n)

# plot results

plots <- within(expand.grid(
  alpha2=c(1, .8, .4),
  delta=c(.3, .5, .7)), {
    rho2_p <- delta
  })
gamma <- log(length(n))/2 # BIC-like penalization

png(paste0(path,"semPowerSim.png"),
    units="in", height=7, width=8.8, res=200, pointsize=15)
layout(matrix(1:9, nrow=3))
par(mar=c(3,3,1,1)+.1)
par(oma=c(0,2,2,0))
mapply(function(alpha2, delta, rho2_p){
  select <- paste0(colnames(plots), "=", c(alpha2, delta, rho2_p), collapse=" ")
  
  large <- paste("rho1_p=0.5",select)
  fit <- gam(powSim[[large]] ~ s(n), family=binomial, gamma=gamma)
  pred <- predict(fit, se.fit=TRUE, newdata=data.frame(n=xpred))
  plot(y=0:1, x=range(xaxis), cex=0, xlab="Sample size (log scale)",
       yaxt="n", xaxt="n", ylab="Probability of rejecting null", mgp=2:0)
  axis(side=1, at=xaxis, labels=exp(xaxis), cex.axis=.82)
  axis(side=2, at=seq(0,1,.2), cex.axis=.9)
  abline(h=seq(0,1,.2), v=xaxis, col="gray")
  polygon(x=c(xpred, rev(xpred), xpred[1]), col=rgb(0,0,0,.25), border=FALSE,
          y=il(c(pred$fit+2*pred$se.fit, rev(pred$fit-2*pred$se.fit), pred$fit[1]+2*pred$se.fit[1])))
  lines(x=xpred, y=il(pred$fit), lwd=2)
  
  medium <- paste("rho1_p=0.3",select)
  fit <- gam(powSim[[medium]] ~ s(n), family=binomial, gamma=gamma)
  pred <- predict(fit, se.fit=TRUE, newdata=data.frame(n=xpred))
  polygon(x=c(xpred, rev(xpred), xpred[1]), col=rgb(0,0,0,.25), border=FALSE,
          y=il(c(pred$fit+2*pred$se.fit, rev(pred$fit-2*pred$se.fit), pred$fit[1]+2*pred$se.fit[1])))
  lines(x=xpred, y=il(pred$fit), lwd=2)
  
  small <- paste("rho1_p=0.1",select)
  fit <- gam(powSim[[small]] ~ s(n), family=binomial, gamma=gamma)
  pred <- predict(fit, se.fit=TRUE, newdata=data.frame(n=xpred))
  polygon(x=c(xpred, rev(xpred), xpred[1]), col=rgb(0,0,0,.25), border=FALSE,
          y=il(c(pred$fit+2*pred$se.fit, rev(pred$fit-2*pred$se.fit), pred$fit[1]+2*pred$se.fit[1])))
  lines(x=xpred, y=il(pred$fit), lwd=2)
  
  none <- paste("rho1_p=0",select)
  fit <- gam(powSim[[none]] ~ s(n), family=binomial, gamma=gamma)
  pred <- predict(fit, se.fit=TRUE, newdata=data.frame(n=xpred))
  polygon(x=c(xpred, rev(xpred), xpred[1]), col=rgb(0,0,0,.25), border=FALSE,
          y=il(c(pred$fit+2*pred$se.fit, rev(pred$fit-2*pred$se.fit), pred$fit[1]+2*pred$se.fit[1])))
  lines(x=xpred, y=il(pred$fit), lwd=2)
}, alpha2=plots$alpha2, delta=plots$delta, rho2_p=plots$rho2_p)
mtext(c("Large indirect effect","Medium indirect effect","Small indirect effect"),
      side=2, line=0, outer=T, at=c(.19,.525,.86))
mtext(c("Perfect reliability","High reliability","Low reliability"),
      side=3, line=0, outer=T, at=c(.19,.525,.86))
dev.off()

# simulation for large indirect effect only -------------------------------

exp(seq(log(50), log(5000), length=7))
# 50.0000  107.7217  232.0794  500.0000 1077.2173 2320.7944 5000.0000
# 50, 100, 1000, 5000

# length of n = number of iterations per parameter combination
xaxis <- log(c(50, 100, 1000, 5000))
xpred <- seq(xaxis[1], xaxis[4], length=40)
n <- seq(xaxis[1], xaxis[4], length=30000)

# sim parameters to iterate over
params <- expand.grid(
  rho1_p=c(0,.1,.2,.3),
  alpha2=c(1, .8, .4),
  delta=.7,
  rho2_p=.7)

# function to collect and name results
iterate <- function(alpha2, delta, rho1_p, rho2_p){
  results <- list(sapply(n, function(x){
    semPower(n=round(exp(x)), alpha1=1, alpha2=alpha2, delta=delta,
             rho1_p=rho1_p, rho2_p=rho2_p)
  }))
  names(results) <- paste0(colnames(params), "=",
                           c(rho1_p, alpha2, delta, rho2_p), collapse=" ")
  return(results)
}

# run sim and save results
system.time({
  powSim2 <- mapply(iterate, alpha2=params$alpha2, delta=params$delta,
                    rho1_p=params$rho1_p, rho2_p=params$rho2_p)
}) # elapsed time = 4.5 hours
saveRDS(powSim2, file=paste0(path,"powSimResults2.rds"))
powSim2 <- readRDS(file=paste0(path,"powSimResults2.rds"))

plots <- within(expand.grid(
  alpha2=c(1, .8, .4), delta=.7, rho2_p=.7), {
    x3 <- log(c(80, 200, 1200))
    x2 <- log(c(500, 1000, 3500))
    x1 <- log(c(900, 2000, 3500))
    x0 <- log(c(3000,3000,3000))
    y3 <- c(.975, .975, .95)
    y2 <- c(.9, .9, .75)
    y1 <- c(.65, .65, .35)
    y0 <- c(.1,.1,.125)
  })

tiff(paste0(path,"Fig12.tif"),
    units="in", height=5, width=15, res=150, pointsize=23)
layout(rbind(1:3))
par(mar=c(3,3,.5,1)+.1)
par(oma=c(0,0,2,1))
mapply(function(alpha2, delta, rho2_p, x3,x2,x1,x0, y3,y2,y1,y0){
  select <- paste0(colnames(plots)[1:3], "=", c(alpha2, delta, rho2_p), collapse=" ")
  
  large <- paste("rho1_p=0.3",select)
  fit <- gam(powSim2[[large]] ~ s(n), family=binomial)
  pred <- predict(fit, newdata=data.frame(n=xpred))
  plot(y=0:1, x=range(xaxis), cex=0, xlab="Sample size (log scale)",
       yaxt="n", xaxt="n", ylab="Probability of rejecting null", mgp=2:0)
  axis(side=1, at=xaxis, labels=exp(xaxis), cex.axis=.82)
  axis(side=2, at=seq(0,1,.2), cex.axis=.9)
  abline(h=seq(0,1,.2), col="gray",
         v=log(c(seq(50,100,10), seq(200,1000,100), seq(2000,5000,1000))))
  lines(x=xpred, y=il(pred), lwd=2)
  
  medium <- paste("rho1_p=0.2",select)
  fit <- gam(powSim2[[medium]] ~ s(n), family=binomial)
  pred <- predict(fit, newdata=data.frame(n=xpred))
  lines(x=xpred, y=il(pred), lwd=2)
  
  small <- paste("rho1_p=0.1",select)
  fit <- gam(powSim2[[small]] ~ s(n), family=binomial)
  pred <- predict(fit, newdata=data.frame(n=xpred))
  lines(x=xpred, y=il(pred), lwd=2)
  
  none <- paste("rho1_p=0",select)
  fit <- gam(powSim2[[none]] ~ s(n), family=binomial)
  pred <- predict(fit, newdata=data.frame(n=xpred))
  lines(x=xpred, y=il(pred), lwd=2)
  
  text(x=c(x3,x2,x1,x0), y=c(y3,y2,y1,y0), labels=c(
    expression(rho[1.2]==.3), expression(rho[1.2]==.2),
    expression(rho[1.2]==.1), expression(rho[1.2]==0)
  ), cex=.8)
}, alpha2=plots$alpha2, delta=plots$delta, rho2_p=plots$rho2_p,x3=plots$x3,
x2=plots$x2,x1=plots$x1,x0=plots$x0,y3=plots$y3,y2=plots$y2,y1=plots$y1,y0=plots$y0)
mtext(expression(paste("A: Perfect reliability ",(alpha==1.0))),
                       side=3, line=0, outer=T, at=.19)
mtext(expression(paste("B: High reliability ",(alpha==0.8))),
                       side=3, line=0, outer=T, at=.515)
mtext(expression(paste("C: Low reliability ",(alpha==0.4))),
                       side=3, line=0, outer=T, at=.85)
dev.off()

# comparison of SEM vs. MR ------------------------------------------------

# function for analytic power
pow <- function(rho1=NA, rho2=NA, delta=NA, rho1_p=NA,
                rho2_p=NA, delta_p=NA, alpha1, alpha2, n){
  # check what form of input is given
  corrs <- rbind(c(rho1, rho2, delta),
                 c(rho1_p, rho2_p, delta_p))
  if(any(colSums(apply(corrs, 2, is.na)) != 1)){
    stop("Supply exactly one each of {rho1, rho1_p}, {rho2, rho2_p}, {delta, delta_p}")
  }
  numSimples <- sum(!is.na(corrs[1,]))
  
  # find implied simple correlations. each case must be handled
  # separately, except numSimples==3, which requires no action
  if(numSimples==2){
    miss <- which(is.na(corrs[1,]))
    nonmiss <- which(!is.na(corrs[1,]))
    corrs[1,miss] <- corrs[2,miss]*prod(sqrt(1-corrs[1,nonmiss]^2)) + prod(corrs[1,nonmiss])
  }
  if(numSimples==1){
    miss <- which(is.na(corrs[2,]))
    nonmiss <- which(!is.na(corrs[2,]))
    corrs[2,miss] <- corrs[1,miss]*sqrt(prod(corrs[2,nonmiss]^2 - 1)) - prod(corrs[2,nonmiss])
    rp <- rbind(c(1, corrs[2,1], corrs[2,2]),
                c(corrs[2,1], 1, corrs[2,3]),
                c(corrs[2,2], corrs[2,3], 1))
    r <- suppressWarnings(pcor2cor(rp))
    # throw error if illegal partial correlation matrix supplied
    if(any(is.nan(r)))
      stop(paste0("Illegal partial correlation matrix:\n
                  rho1_p = ",corrs[2,1],", rho2_p = ",corrs[2,2],", delta_p = ",corrs[2,3]))
    corrs[1,] <- r[c(2,3,6)]
  }
  if(numSimples==0){
    rp <- rbind(c(1, rho1_p, rho2_p),
                c(rho1_p, 1, delta_p),
                c(rho2_p, delta_p, 1))
    r <- suppressWarnings(pcor2cor(rp))
    # throw error if illegal partial correlation matrix supplied
    if(any(is.nan(r)))
      stop(paste0("Illegal partial correlation matrix:\n
                  rho1_p = ",rho1_p,", rho2_p = ",rho2_p,", delta_p = ",delta_p))
    corrs[1,] <- r[c(2,3,6)]
  }
  
  # build full simple correlation matrix
  r <- rbind(c(1, corrs[1,1], corrs[1,2], corrs[1,1]*sqrt(alpha1), corrs[1,2]*sqrt(alpha2)),
             c(corrs[1,1], 1, corrs[1,3], sqrt(alpha1), corrs[1,3]*sqrt(alpha2)),
             c(corrs[1,2], corrs[1,3], 1, corrs[1,3]*sqrt(alpha1), sqrt(alpha2)),
             c(corrs[1,1]*sqrt(alpha1), sqrt(alpha1), corrs[1,3]*sqrt(alpha1), 1, corrs[1,3]*sqrt(alpha1*alpha2)),
             c(corrs[1,2]*sqrt(alpha2), corrs[1,3]*sqrt(alpha2), sqrt(alpha2), corrs[1,3]*sqrt(alpha1*alpha2), 1))
  # set names
  dimnames(r) <- list(c("y", "t1", "t2", "x1", "x2"),
                      c("y", "t1", "t2", "x1", "x2"))
  # throw error if illegal correlation matrix supplied
  if(det(r[1:3, 1:3]) < -.00001)
    stop(paste0("Illegal correlation matrix (determinant < 0):\n
                rho1 = ",r["y","t1"],", rho2 = ",r["y","t2"],", delta = ",r["t1","t2"]))
  
  ### get elements of mu
  # degrees of freedom
  nu <- n-3
  # partial corr for x1
  pcorr1 <- (r["y","x1"] - r["y","x2"]*r["x1","x2"])/
    sqrt(1 - r["y","x2"]^2)/sqrt(1 - r["x1","x2"]^2)
  # partial corr for x2
  pcorr2 <- (r["y","x2"] - r["y","x1"]*r["x1","x2"])/
    sqrt(1 - r["y","x1"]^2)/sqrt(1 - r["x1","x2"]^2)
  # ncp for testing partial correlation between y and x1
  t1 <- sqrt(nu*pcorr1^2/(1 - pcorr1^2))
  # ncp for testing partial correlation between y and x1
  t2 <- sqrt(nu*pcorr2^2/(1 - pcorr2^2))
  
  ### get elements of sigma
  # residual variance
  err <- 1 - (r["y","x1"]^2 + r["y","x2"]^2 - 2*r["y","x1"]*r["y","x2"]*r["x1","x2"])/
    (1 - r["t1","t2"]^2*alpha1*alpha2)
  # variances of b1 and b2
  vb <- err/nu/(1 - r["x1","x2"]^2)
  # cov(b1,b2)
  covb <- err*r["x1","x2"]/n/(r["x1","x2"]^2 - 1)
  
  # return probability of observing significant x1 slope
  # p1.1: probability that x1 significant & x2 n.s.
  # p1.2: probability that x1 n.s. & x2 significant
  # p2: probability that x1 significant & x2 significant
  # p1: probability that only 1 predictor is significant (p1.1 + p1.2)
  # p: probability that x1 is significant (p1.1 + p2)
  bCovMat <- cbind(c(vb, covb), c(covb, vb))
  p1.1 <- pmvt(lower=c(qt(.975, nu), qt(.025, nu)),
               upper=c(Inf, qt(.975, nu)),
               df=nu, delta=c(t1, t2), corr=cov2cor(bCovMat)) + 
    pmvt(lower=c(-Inf, qt(.025, nu)),
         upper=c(qt(.025, nu), qt(.975, nu)),
         df=nu, delta=c(t1, t2), corr=cov2cor(bCovMat))
  p1.2 <- pmvt(lower=c(qt(.025, nu), qt(.975, nu)),
               upper=c(qt(.975, nu), Inf),
               df=nu, delta=c(t1, t2), corr=cov2cor(bCovMat)) + 
    pmvt(lower=c(qt(.025, nu), -Inf),
         upper=c(qt(.975, nu), qt(.025, nu)),
         df=nu, delta=c(t1, t2), corr=cov2cor(bCovMat))
  p2 <- pmvt(lower=c(qt(.975, nu), qt(.975, nu)),
             upper=c(Inf, Inf),
             df=nu, delta=c(t1, t2), corr=cov2cor(bCovMat)) + 
    pmvt(lower=c(-Inf, qt(.975, nu)),
         upper=c(qt(.025, nu), Inf),
         df=nu, delta=c(t1, t2), corr=cov2cor(bCovMat)) +
    pmvt(lower=c(qt(.975, nu), -Inf),
         upper=c(Inf, qt(.025, nu)),
         df=nu, delta=c(t1, t2), corr=cov2cor(bCovMat)) + 
    pmvt(lower=c(-Inf, -Inf),
         upper=c(qt(.025, nu), qt(.025, nu)),
         df=nu, delta=c(t1, t2), corr=cov2cor(bCovMat))
  rp <- cor2pcor(r[1:3,1:3])
  colnames(rp) <- colnames(r)[1:3]
  rownames(rp) <- rownames(r)[1:3]
  r2true <- (r["y","t1"]^2 + r["y","t2"]^2 - 2*r["y","t1"]*r["y","t2"]*r["t1","t2"])/
    (1 - r["t1","t2"]^2)
  return(c(p1.1=p1.1, p1.2=p1.2, p1=p1.1+p1.2, p2=p2, p=p1.1+p2,
           rho1=r["y","t1"], rho2=r["y","t2"], delta=r["t1","t2"],
           rho1_p=rp["y","t1"], rho2_p=rp["y","t2"], delta_p=rp["t1","t2"],
           pcorr1=pcorr1, scorr1=r["y","x1"], pcorr2=pcorr2, scorr2=r["y","x2"],
           r_x1x2=r["x1","x2"], r2true=r2true, err=err))
}

# get analytic and simulation results
# estimate time for 300K sims = 5.5 hours
alpha <- seq(.01, 1, length=300000)
rng <- seq(0,1,.01)
anResults <- sapply(rng, function(x){
  pow(n=100, alpha1=1, alpha2=x, delta=.5, rho1_p=.2, rho2_p=.5)["p"]
})
system.time({
  SEMvsMR_results <- sapply(alpha, function(x){
    semPower(n=100, alpha1=1, alpha2=x, delta=.5, rho1_p=.2, rho2_p=.5)
  })
}) # time = 5 hours
saveRDS(SEMvsMR_results, file=paste0(path,"SEMvsMR_results.rds"))
# SEMvsMR_results <- readRDS(paste0(path,"SEMvsMR_results.rds"))
fit <- gam(SEMvsMR_results ~ s(alpha), family=binomial)
pred <- predict(fit, newdata=data.frame(alpha=rng))

# plot
tiff(paste0(path,"Fig11.tif"),
    units="in", height=6, width=11, res=125, pointsize=15)
layout(rbind(1:2))
par(mar=c(4,4,1,1)+.1)
par(oma=c(0,0,2,0))
plot(y=0:1, x=range(rng), cex=0, xlab="Reliability of covariate",
     ylab="Probability of rejecting null", mgp=2:0)
grid(lty=1, col="gray")
# polygon(x=c(rng, rev(rng), rng[1]), col=rgb(0,0,0,.25), border=FALSE,
#         y=il(c(pred$fit+2*pred$se.fit, rev(pred$fit-2*pred$se.fit), pred$fit[1]+2*pred$se.fit[1])))
lines(x=rng, y=il(pred), lwd=2)
lines(x=rng, y=anResults, lwd=2)
text(x=c(.4,.6), y=c(.9,.3), labels=c("Regression","SEM"))

# beta +/- SE for regression vs. SEM --------------------------------------

# analytic b + SE for regression
bSE <- function(n, alpha1, alpha2, delta, rho1_p, rho2_p){
  # compute the implied set of 3 simple correlations
  delta_p <- delta*sqrt(prod(c(rho1_p, rho2_p)^2 - 1)) - prod(c(rho1_p, rho2_p))
  rp <- rbind(c(1, rho1_p, rho2_p),
              c(rho1_p, 1, delta_p),
              c(rho2_p, delta_p, 1))
  r <- suppressWarnings(pcor2cor(rp))
  # throw error if illegal partial correlation matrix supplied
  if(any(is.nan(r)))
    stop(paste0("Illegal partial correlation matrix:\n
                rho1_p = ",rho1_p,", rho2_p = ",rho2_p,", delta_p = ",delta_p))
  
  # get b
  dimnames(r) <- list(c("y","t1","t2"), c("y","t1","t2"))
  r_y1 <- r["y","t1"]*sqrt(alpha1)
  r_y2 <- r["y","t2"]*sqrt(alpha2)
  r_12 <- r["t1","t2"]*sqrt(alpha1*alpha2)
  b <- (r_y1 - r_y2*r_12)/(1-r_12^2)
  
  # get SE
  err <- 1 - (r_y1^2 + r_y2^2 - 2*r_y1*r_y2*r_12)/(1 - r_12^2)
  SE <- sqrt(err/(n-3)/(1 - r_12^2))
  
  return(c(b=b, SE=SE))
}

# modify semPower() to return b + SE
semPower_bSE <- function(n, alpha1, alpha2, delta, rho1_p, rho2_p){
  # compute the implied set of 3 simple correlations
  delta_p <- delta*sqrt(prod(c(rho1_p, rho2_p)^2 - 1)) - prod(c(rho1_p, rho2_p))
  rp <- rbind(c(1, rho1_p, rho2_p),
              c(rho1_p, 1, delta_p),
              c(rho2_p, delta_p, 1))
  r <- suppressWarnings(pcor2cor(rp))
  # throw error if illegal partial correlation matrix supplied
  if(any(is.nan(r)))
    stop(paste0("Illegal partial correlation matrix:\n
                rho1_p = ",rho1_p,", rho2_p = ",rho2_p,", delta_p = ",delta_p))
  
  # make data
  dat <- data.frame(mvrnorm(n, mu=numeric(3), Sigma=r))
  names(dat) <- c("y","t1","t2")
  dat <- within(dat, {
    x1 <- t1 + rnorm(n, sd=sqrt((1-alpha1)/alpha1))
    x2 <- t2 + rnorm(n, sd=sqrt((1-alpha2)/alpha2))
  })

  # fit model
  z <- "NULL"
  mod <- suppressWarnings(sem(paste0("
                                     y ~ t1 + t2
                                     t1 =~ 1*x1
                                     x1 ~~ ",(1-alpha1)*var(dat$x1),"*x1
                                     t2 =~ 1*x2
                                     x2 ~~ ",(1-alpha2)*var(dat$x2),"*x2
                                     "), std.lv=TRUE, data=dat))
  
  # return stats
  if(lavInspect(mod, "converged")){
#     z <- try(lavInspect(mod, "est")$beta["y","t1"]/sqrt(diag(vcov(mod)))["y~t1"],
#              silent=TRUE)
    z <- try(c(b=lavInspect(mod, "est")$beta["y","t1"], SE=sqrt(diag(vcov(mod)))["y~t1"]),
             silent=TRUE)
    if(class(z) == "try-error") z <- NA
    return(z)
  } else return(NA)
}

# analytic results
xpred <- seq(0,1,.01)
est <- sapply(xpred, function(x){
  bSE(n=100, alpha1=1, alpha2=x, delta=.5, rho1_p=.2, rho2_p=.5)
})

# simulation results
# 100K iterations should take ~90 minutes
alpha <- seq(.01, 1, length=100000)
system.time({
  results <- sapply(alpha, function(x){
    semPower_bSE(n=100, alpha1=1, alpha2=x, delta=.5, rho1_p=.2, rho2_p=.5)
  })
})
# saveRDS(results, file=paste0(path,"bSE_MR_vs_SEM.rds"))
results <- readRDS(paste0(path,"bSE_MR_vs_SEM.rds"))
alpha_trim <- tail(head(alpha, -50), -50)
fit1 <- gam(rollapply(sapply(results, "[", "b"), width=101, FUN=median, na.rm=TRUE) ~ s(alpha_trim))
b <- predict(fit1, newdata=data.frame(alpha_trim=xpred))
fit2 <- gam(log(rollapply(sapply(results, "[", "SE.y~t1"), width=101, FUN=median, na.rm=TRUE)) ~ s(alpha_trim))
SE <- exp(predict(fit2, newdata=data.frame(alpha_trim=xpred)))

# plot!
plot(y=c(-.25, .6), x=0:1, cex=0, xlab="Reliability of covariate", mgp=2:0,
     ylab="Expected regression coefficient\n+/- expected standard error")
polygon(x=c(xpred, rev(xpred), xpred[1]), col=rgb(1,0,0,.25), border=FALSE,
        y=c(est["b",]+est["SE",], rev(est["b",]-est["SE",]), est["b",1]+est["SE",1]))
polygon(x=c(xpred, rev(xpred), xpred[1]), col=rgb(0,0,1,.25), border=FALSE,
        y=c(b+SE, rev(b-SE), b[1]+SE[1]))
lines(y=est["b",], x=xpred, lwd=3, col="red")
lines(y=b, x=xpred, lwd=3, col="blue")
abline(h=c(0, tail(est["b",],1)), lty=1:2)
text(x=c(.6,.6), y=c(.375,.125), labels=c("Regression","SEM"),
     col=c("red","blue"))

mtext(expression(paste("N=100, medium partial correlation (",rho[1.2]," = ",.2,")",
                       ", medium indirect effect (",rho[2.1]," = ",delta," = ",.5,")")),
      outer=TRUE)
dev.off()

# simulation for small indirect effect only -------------------------------

exp(seq(log(50), log(5000), length=7))
# 50.0000  107.7217  232.0794  500.0000 1077.2173 2320.7944 5000.0000
# 50, 100, 1000, 5000

# length of n = number of iterations per parameter combination
xaxis <- log(c(50, 100, 1000, 5000))
xpred <- seq(xaxis[1], xaxis[4], length=40)
n <- seq(xaxis[1], xaxis[4], length=30000)

# sim parameters to iterate over
params <- expand.grid(
  rho1_p=c(0,.1,.2,.3),
  alpha2=c(1, .8, .4),
  delta=.3,
  rho2_p=.3)

# function to collect and name results
iterate <- function(alpha2, delta, rho1_p, rho2_p){
  results <- list(sapply(n, function(x){
    semPower(n=round(exp(x)), alpha1=1, alpha2=alpha2, delta=delta,
             rho1_p=rho1_p, rho2_p=rho2_p)
  }))
  names(results) <- paste0(colnames(params), "=",
                           c(rho1_p, alpha2, delta, rho2_p), collapse=" ")
  return(results)
}

# run sim and save results
system.time({
  powSim3 <- mapply(iterate, alpha2=params$alpha2, delta=params$delta,
                    rho1_p=params$rho1_p, rho2_p=params$rho2_p)
}) # elapsed time = 4.5 hours
saveRDS(powSim3, file=paste0(path,"powSimResults3.rds"))
powSim3 <- readRDS(file=paste0(path,"powSimResults3.rds"))

plots <- within(expand.grid(
  alpha2=c(1, .8, .4), delta=.3, rho2_p=.3), {
    x3 <- log(c(80, 200, 1200))
    x2 <- log(c(500, 1000, 3500))
    x1 <- log(c(900, 2000, 3500))
    x0 <- log(c(3000,3000,3000))
    y3 <- c(.975, .975, .95)
    y2 <- c(.9, .9, .75)
    y1 <- c(.65, .65, .35)
    y0 <- c(.1,.1,.125)
  })

png(paste0(path,"semPowerSim_smallAB.png"),
    units="in", height=5, width=15, res=200, pointsize=23)
layout(rbind(1:3))
par(mar=c(3,3,.5,1)+.1)
par(oma=c(0,0,2,1))
mapply(function(alpha2, delta, rho2_p, x3,x2,x1,x0, y3,y2,y1,y0){
  select <- paste0(colnames(plots)[1:3], "=", c(alpha2, delta, rho2_p), collapse=" ")
  
  large <- paste("rho1_p=0.3",select)
  fit <- gam(powSim3[[large]] ~ s(n), family=binomial)
  pred <- predict(fit, newdata=data.frame(n=xpred))
  plot(y=0:1, x=range(xaxis), cex=0, xlab="Sample size (log scale)",
       yaxt="n", xaxt="n", ylab="Probability of rejecting null", mgp=2:0)
  axis(side=1, at=xaxis, labels=exp(xaxis), cex.axis=.82)
  axis(side=2, at=seq(0,1,.2), cex.axis=.9)
  abline(h=seq(0,1,.2), col="gray",
         v=log(c(seq(50,100,10), seq(200,1000,100), seq(2000,5000,1000))))
  lines(x=xpred, y=il(pred), lwd=2)
  
  medium <- paste("rho1_p=0.2",select)
  fit <- gam(powSim3[[medium]] ~ s(n), family=binomial)
  pred <- predict(fit, newdata=data.frame(n=xpred))
  lines(x=xpred, y=il(pred), lwd=2)
  
  small <- paste("rho1_p=0.1",select)
  fit <- gam(powSim3[[small]] ~ s(n), family=binomial)
  pred <- predict(fit, newdata=data.frame(n=xpred))
  lines(x=xpred, y=il(pred), lwd=2)
  
  none <- paste("rho1_p=0",select)
  fit <- gam(powSim3[[none]] ~ s(n), family=binomial)
  pred <- predict(fit, newdata=data.frame(n=xpred))
  lines(x=xpred, y=il(pred), lwd=2)
  
  text(x=c(x3,x2,x1,x0), y=c(y3,y2,y1,y0), labels=c(
    expression(rho[1.2]==.3), expression(rho[1.2]==.2),
    expression(rho[1.2]==.1), expression(rho[1.2]==0)
  ), cex=.8)
}, alpha2=plots$alpha2, delta=plots$delta, rho2_p=plots$rho2_p,x3=plots$x3,
x2=plots$x2,x1=plots$x1,x0=plots$x0,y3=plots$y3,y2=plots$y2,y1=plots$y1,y0=plots$y0)
mtext(expression(paste("A: Perfect reliability ",(alpha==1.0))),
      side=3, line=0, outer=T, at=.19)
mtext(expression(paste("B: High reliability ",(alpha==0.8))),
      side=3, line=0, outer=T, at=.515)
mtext(expression(paste("C: Low reliability ",(alpha==0.4))),
      side=3, line=0, outer=T, at=.85)
dev.off()
