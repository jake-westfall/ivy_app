library(mvtnorm) # for rmvnorm() and pmvt()
library(MBESS) # for cor2cov()
library(ppcor) # for pcor
library(corpcor) # for pcor2cor

# path <- "/Users/Jake/Desktop/Dropbox/Dropbox/SIC/"
path <- "/Users/Jake/Desktop/figs/"

# simulation --------------------------------------------------------------


# rho1: true correlation between y and latent variable 1
# rho2: true correlation between y and latent variable 2
# delta: correlation between latent variables 1 and 2
# alpha1: reliability of x1
# alpha2: reliability of x2
# n: sample size

sim <- function(rho1, rho2, delta, alpha1, alpha2, n){
  # build correlation matrix
  mat <- rbind(c(1, rho1, rho2),
               c(rho1, 1, delta),
               c(rho2, delta, 1))
  
  # throw error if illegal correlation matrix supplied
  if(det(mat) < -.00001)
    stop(paste0("Illegal correlation matrix (determinant < 0):\n
                rho1 = ",rho1,", rho2 = ",rho2,", delta = ",delta))
  
  # generate data
  dat <- rmvnorm(n, mean=numeric(3), sigma=cor2cov(mat, sd=c(1,1,1)))
  dat <- data.frame(y=dat[,1],
                    x1=sqrt(alpha1)*dat[,2] + sqrt(1-alpha1)*rnorm(n),
                    x2=sqrt(alpha2)*dat[,3] + sqrt(1-alpha2)*rnorm(n))
  
  # return stats
  mod <- lm(y ~ x1 + x2, data=dat)
  summ <- summary(mod)
  return(c(r_yx1=cor(dat$y,dat$x1), rp_yx1=pcor(dat)$estimate[2,1], r_yx2=cor(dat$y,dat$x2),
           r_x1x2=cor(dat$x1,dat$x2), r2=summ$r.squared, coef(mod), 
           t_x1=coef(summ)["x1","t value"], p_x1=coef(summ)["x1","Pr(>|t|)"]<.05,
           t_x2=coef(summ)["x2","t value"], p_x2=coef(summ)["x2","Pr(>|t|)"]<.05))
}
dat <- sim(rho1=.5, rho2=.5, alpha1=.8, alpha2=.8, delta=.5, n=100)
dat

# test sim for 1 set of values
result <- replicate(10000, {
  sim(rho1=.5, rho2=.5, alpha1=.8, alpha2=.8, delta=.5, n=100)
})
round(rowMeans(result), digits=3)
round(c(r_yx1=rho1*sqrt(alpha1), r_yx2=rho2*sqrt(alpha2),
  r_x1x2=delta*sqrt(alpha1*alpha2)), digits=3)

### run full sim -- takes around 70 minutes ###

Sys.time()

# vary delta and effect sizes
params1 <- expand.grid(rho1=c(.1,.3,.5), n=2^(4:8), delta=c(.2,.5,.8))
params1 <- within(params1, {
  rho2 <- rho2
  alpha1 <- .7
  alpha2 <- .7
})
results1 <- with(params1, replicate(10000, mapply(sim, rho1=rho1,
  rho2=rho2, alpha1=alpha1, alpha2=alpha2, delta=delta, n=n)))
Sys.time()

# fix delta=1, vary reliabilities
params2 <- expand.grid(rho1=c(.1,.3,.5), n=2^(4:8), alpha1=c(.4,.6,.8))
params2 <- within(params2, {
  rho2 <- rho1
  delta <- 1
  alpha2 <- alpha2
})
results2 <- with(params2, replicate(10000, mapply(sim, rho1=rho1,
  rho2=rho2, alpha1=alpha1, alpha2=alpha2, delta=delta, n=n)))
Sys.time()

# # test that all correlations are legal
# test <- expand.grid(rho1=c(.1,.3,.5),rho2=c(.1,.3,.5))
# apply(test, 1, function(x) {
#   c(cos(acos(x[1]) + acos(x[2])), cos(acos(x[1]) - acos(x[2])))
# })

# vary rho1 and rho2
params3 <- expand.grid(rho1=c(.1,.3,.5), rho2=c(.1, .3, .5), n=2^(4:8))
params3 <- within(params3, {
  delta <- .5
  alpha1 <- .7
  alpha2 <- .7
})
results3 <- with(params3, replicate(10000, mapply(sim, rho1=rho1,
  rho2=rho2, alpha1=alpha1, alpha2=alpha2, delta=delta, n=n)))
Sys.time()

### save results
results <- list(results1, results2, results3)
saveRDS(results, file="/Users/Jake/Desktop/Dropbox/SIC/simResults.rds")


# analytic results --------------------------------------------------------


### note: specify either 3 simples OR 3 partials--don't mix!
# rho1: true correlation between y and latent variable 1
# rho2: true correlation between y and latent variable 2
# delta: correlation between latent variables 1 and 2
# rho1_p: true *partial* correlation between y and latent variable 1
# rho2_p: true *partial* correlation between y and latent variable 2
# delta_p: *partial* correlation between latent variables 1 and 2
# alpha1: reliability of x1
# alpha2: reliability of x2
# n: sample size

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

# compare 3 simples to 3 partials
pow(rho1=.3, rho2=.3, delta=.5, alpha1=.7, alpha2=.7, n=50)
pow(rho1_p=0.1815683, rho2_p=0.1815683, delta_p=0.4505495, alpha1=.7, alpha2=.7, n=50)

# compare 3 partials to 3 simples
pow(rho1_p=.3, rho2_p=.3, delta_p=.5, alpha1=.7, alpha2=.7, n=50)
pow(rho1=0.5447048, rho2=0.5447048, delta=0.6483516, alpha1=.7, alpha2=.7, n=50)

# compare 1 partial and 2 simples to 3 simples
pow(rho1_p=.3, rho2=.3, delta=.5, alpha1=.7, alpha2=.7, n=50)
pow(rho1=0.3978407, rho2=.3, delta=.5, alpha1=.7, alpha2=.7, n=50)

# compare 2 partials and 1 simple to 3 simples
pow(rho1=.3, rho2=.3, delta=.5, alpha1=.7, alpha2=.7, n=50)
pow(rho1_p=0.18156826, rho2_p=0.18156826, delta=.5, alpha1=.7, alpha2=.7, n=50)


# test predictions against sims -------------------------------------------


### p1: probability that x1 significant & x2 n.s.

# results1
obs1 <- array(rowMeans(results1["p_x1",,] & !results1["p_x2",,]), dim=c(3,5,3),
              dimnames=list(rho1=c(.1,.3,.5), n=2^(4:8), delta=c(.2,.5,.8)))
pred1 <- with(params1, array(mapply(pow, rho1=rho1, rho2=rho2, alpha1=alpha1,
                                    alpha2=alpha2, delta=delta, n=n)["p1",], dim=c(3,5,3),
                             dimnames=list(rho1=c(.1,.3,.5), n=2^(4:8), delta=c(.2,.5,.8))))
obs1
round(pred1, 4)
obs1 - round(pred1, 4)

# results2
obs2 <- array(rowMeans(results2["p_x1",,] & !results2["p_x2",,]), dim=c(3,5,3),
              dimnames=list(rho=c(.1,.3,.5), n=2^(4:8), alpha1=c(.4,.6,.8)))
pred2 <- with(params2, array(mapply(pow, rho1=rho1, rho2=rho2, alpha1=alpha1,
                                    alpha2=alpha2, delta=delta, n=n)["p1",], dim=c(3,5,3),
                             dimnames=list(rho=c(.1,.3,.5), n=2^(4:8), alpha1=c(.4,.6,.8))))
obs2
round(pred2, 4)
obs2 - round(pred2, 4)

# results3
obs3 <- array(rowMeans(results3["p_x1",,] & !results3["p_x2",,]), dim=c(3,3,5),
              dimnames=list(rho1=c(.1,.3,.5), rho2=c(.1, .3, .5), n=2^(4:8)))
pred3 <- with(params3, array(mapply(pow, rho1=rho1, rho2=rho2, alpha1=alpha1,
                                    alpha2=alpha2, delta=delta, n=n)["p1",], dim=c(3,3,5),
                             dimnames=list(rho1=c(.1,.3,.5), rho2=c(.1, .3, .5), n=2^(4:8))))
obs3
round(pred3, 4)
obs3 - round(pred3, 4)

# plot
plot(y=c(obs1, obs2, obs3), x=c(pred1, pred2, pred3), pch=20,
     ylab="Simulation results", xlab="Analytic predictions",
     main="Probability of rejecting x1 but not x2")
abline(0, 1)


### p2: probability that x1 significant & x2 significant

# results1
obs1 <- array(rowMeans(results1["p_x1",,] & results1["p_x2",,]), dim=c(3,5,3),
              dimnames=list(rho1=c(.1,.3,.5), n=2^(4:8), delta=c(.2,.5,.8)))
pred1 <- with(params1, array(mapply(pow, rho1=rho1, rho2=rho2, alpha1=alpha1,
                                    alpha2=alpha2, delta=delta, n=n)["p2",], dim=c(3,5,3),
                             dimnames=list(rho1=c(.1,.3,.5), n=2^(4:8), delta=c(.2,.5,.8))))
obs1
round(pred1, 4)
obs1 - round(pred1, 4)

# results2
obs2 <- array(rowMeans(results2["p_x1",,] & results2["p_x2",,]), dim=c(3,5,3),
              dimnames=list(rho=c(.1,.3,.5), n=2^(4:8), alpha1=c(.4,.6,.8)))
pred2 <- with(params2, array(mapply(pow, rho1=rho1, rho2=rho2, alpha1=alpha1,
                                    alpha2=alpha2, delta=delta, n=n)["p2",], dim=c(3,5,3),
                             dimnames=list(rho=c(.1,.3,.5), n=2^(4:8), alpha1=c(.4,.6,.8))))
obs2
round(pred2, 4)
obs2 - round(pred2, 4)

# results3
obs3 <- array(rowMeans(results3["p_x1",,] & results3["p_x2",,]), dim=c(3,3,5),
              dimnames=list(rho1=c(.1,.3,.5), rho2=c(.1, .3, .5), n=2^(4:8)))
pred3 <- with(params3, array(mapply(pow, rho1=rho1, rho2=rho2, alpha1=alpha1,
                                    alpha2=alpha2, delta=delta, n=n)["p2",], dim=c(3,3,5),
                             dimnames=list(rho1=c(.1,.3,.5), rho2=c(.1, .3, .5), n=2^(4:8))))
obs3
round(pred3, 4)
obs3 - round(pred3, 4) 

# plot
plot(y=c(obs1, obs2, obs3), x=c(pred1, pred2, pred3), pch=20,
     ylab="Simulation results", xlab="Analytic predictions",
     main="Probability of rejecting both x1 and x2")
abline(0, 1)


# LRT stuff ---------------------------------------------------------------


# critical LR (asymptotic)
exp(qchisq(.95, df=1)/2)

# critical LR as function of sample size
LRcrit1 <- function(n){
  df <- n-2
  (1 + qf(.95, 1, df)/df)^(n/2)
}
curve(LRcrit1, from=10, to=200)
LRcrit1(999999) # matches asymptotic, as expected

# critical LR as function of sample size (log2 scale)
LRcrit2 <- function(n){
  n <- 2^n
  df <- n-2
  (1 + qf(1-alpha, 1, df)/df)^(n/2)
}
alpha <- .05
curve(LRcrit2, xlim=c(4,10), ylim=c(1,20), xaxt="n", lwd=2,
      main="Critical value of likelihood ratio for\na t-test or simple regression model",
#      main=expression(paste("     Critical value of likelihood ratio\nas a function of sample size and ",alpha,"-level")),
      xlab="Sample size", ylab="Likelihood ratio critical value")
axis(side=1, at=4:10, labels=2^(4:10))
alpha <- .1
curve(LRcrit2, add=TRUE, lwd=2)
alpha <- .025
curve(LRcrit2, add=TRUE, lwd=2)
text(x=7, y=c(11.7,6,3), labels=c(expression(alpha==.025),
  expression(alpha==.05), expression(alpha==.1)))

# relative likelihood
n <- 20
df <- n - 2
func <- function(x) dt(x, df=df, ncp=x)/dt(x, df=df, ncp=0)
curve(func, from=0, to=2.5)
abline(h=LRcrit1(n), v=qt(.975, df=df), lty=2)
func <- function(x){
  dt(x, df=df, ncp=x*4*df/(4*df+3))/dt(x, df=df, ncp=0)
}
curve(func, from=0, to=2.5)
abline(h=LRcrit1(n), v=qt(.975, df=df), lty=2)

n <- 200
df <- n - 2
func <- function(x) dt(x, df=df, ncp=x)/dt(x, df=df, ncp=0)
curve(func, from=1.9, to=2.1)
abline(h=LRcrit1(n), v=qt(.975, df=df), lty=2)
func <- function(x){
  dt(x, df=df, ncp=x*4*df/(4*df+3))/dt(x, df=df, ncp=0)
}
curve(func, from=1.9, to=2.1)
abline(h=LRcrit1(n), v=qt(.975, df=df), lty=2)

LRcrit3 <- function(n){
  df <- n-3
  dt(qt(.975, df=df), df=n-2, ncp=qt(.975, df=df))/
    dt(qt(.975, df=df), df=n-2, ncp=0)
}
LRcrit3(20)

LRcrit4 <- function(n){
  n <- 2^n
  df <- n-3
  dt(qt(.975, df=df), df=n-2, ncp=qt(.975, df=df))/
    dt(qt(.975, df=df), df=n-2, ncp=0)
}
curve(LRcrit4, xlim=c(3,10), xaxt="n", lwd=2,
      main="Critical value of likelihood ratio for\na t-test or simple regression model",
      xlab="Sample size", ylab="Likelihood ratio critical value")
axis(side=1, at=3:10, labels=2^(3:10))


# plot of multivariate t --------------------------------------------------


# pow(rho1_p=.3, rho2_p=.3, delta=.5, alpha1=.6, alpha2=.6, n=40)
# [1] 1.778164
#         [,1]    [,2]
# [1,]  1.0000 -0.2775
# [2,] -0.2775  1.0000
#        p1        p2         p      rho1      rho2     delta    rho1_p 
# 0.2737027 0.1363932 0.4100959 0.4610840 0.4610840 0.5000000 0.3000000 
#    rho2_p   delta_p    pcorr1    scorr1    pcorr2    scorr2    r_x1x2 
# 0.3000000 0.3650000 0.2805853 0.3571541 0.2805853 0.3571541 0.3000000 
#    r2true 
# 0.2834646 

rang <- seq(1, 5, length.out=200)
samp <- expand.grid(t1=rang, t2=rang)
mat <- cbind(c(1, -0.2775), c(-0.2775, 1))
z <- matrix(apply(samp, 1, dmvt, delta=c(1.778164, 1.778164), df=40-3,
                  sigma=mat, log=FALSE), ncol=length(rang), byrow=TRUE)

mat <- cbind(c(1, 0.3435837), c(0.3435837, 1))
z <- matrix(apply(samp, 1, dmvt, delta=c(3.807526, 3.807526), df=40-3,
                  sigma=mat, log=FALSE), ncol=length(rang), byrow=TRUE)


### build plot
png("/Users/Jake/Desktop/joint_t_dist.png", type="cairo",
    units="in", height=8, width=8, res=200, pointsize=18)
contour(y=rang, x=rang, z=z, drawlabels=FALSE,
        main="Joint distribution of t-statistics",
        xlab=expression(paste("t-statistic for ",beta[X1])),
        ylab=expression(paste("t-statistic for ",beta[X2])))
rect(xleft=c(qt(.975, 40-3), qt(.975, 40-3)), ytop=c(3.5, 3.5),
     xright=c(3.5, 3.5), ybottom=c(qt(.975, 40-3), -.5), border=NA,
     col=rgb(0,0,0,.15))
abline(v=qt(.975, 40-3), h=qt(.975, 40-3), lty=2, lwd=2)
text(x=rep(2.55, 4), y=c(3, 2.82, 1.9, 1.72),
     labels=c(expression(paste("Region ",p[2])), "(13.6%)",
              expression(paste("Region ",p[1])), "(27.4%)"))
text(x=c(.25, 1.9), y=c(2.15, .3), labels=expression(t[crit]))
dev.off()


# plots for true delta=1 --------------------------------------------------


xrange <- seq(0, .99, .01)
y1 <- sapply(xrange, function(x){
  pow(rho1=.1, rho2=.1, delta=1, alpha1=x, alpha2=x, n=100)["p"]
})
y2 <- sapply(xrange, function(x){
  pow(rho1=.3, rho2=.3, delta=1, alpha1=x, alpha2=x, n=100)["p"]
})
y3 <- sapply(xrange, function(x){
  pow(rho1=.5, rho2=.5, delta=1, alpha1=x, alpha2=x, n=100)["p"]
})
y <- cbind(y1, y2, y3)
matplot(x=xrange, y=y, type="l", col="black", lwd=2, lty=1,
        ylim=0:1, ylab="Probability", xlab="Reliabilities")


# plots for true rho1=0 ---------------------------------------------------


### delta=.1
xrange <- seq(0, 1, .01)
y1 <- sapply(xrange, function(x){
  pow(rho1_p=0, rho2=.1, delta=.1, alpha1=x, alpha2=x, n=100)["p"]
})
y2 <- sapply(xrange, function(x){
  pow(rho1_p=0, rho2=.3, delta=.1, alpha1=x, alpha2=x, n=100)["p"]
})
y3 <- sapply(xrange, function(x){
  pow(rho1_p=0, rho2=.5, delta=.1, alpha1=x, alpha2=x, n=100)["p"]
})
y <- cbind(y1, y2, y3)
matplot(x=xrange, y=y, type="l", col="black", lwd=2, lty=1, ylim=0:1,
        ylab="Probability", xlab="Reliabilities", main="Delta=.1")

### delta=.3
xrange <- seq(0, 1, .01)
y1 <- sapply(xrange, function(x){
  pow(rho1_p=0, rho2=.1, delta=.3, alpha1=x, alpha2=x, n=100)["p"]
})
y2 <- sapply(xrange, function(x){
  pow(rho1_p=0, rho2=.3, delta=.3, alpha1=x, alpha2=x, n=100)["p"]
})
y3 <- sapply(xrange, function(x){
  pow(rho1_p=0, rho2=.5, delta=.3, alpha1=x, alpha2=x, n=100)["p"]
})
y <- cbind(y1, y2, y3)
matplot(x=xrange, y=y, type="l", col="black", lwd=2, lty=1, ylim=0:1,
        ylab="Probability", xlab="Reliabilities", main="Delta=.3")

### delta=.5
xrange <- seq(0, 1, .01)
y1 <- sapply(xrange, function(x){
  pow(rho1_p=0, rho2=.1, delta=.5, alpha1=x, alpha2=x, n=100)["p"]
})
y2 <- sapply(xrange, function(x){
  pow(rho1_p=0, rho2=.3, delta=.5, alpha1=x, alpha2=x, n=100)["p"]
})
y3 <- sapply(xrange, function(x){
  pow(rho1_p=0, rho2=.5, delta=.5, alpha1=x, alpha2=x, n=100)["p"]
})
y <- cbind(y1, y2, y3)
matplot(x=xrange, y=y, type="l", col="black", lwd=2, lty=1, ylim=0:1,
        ylab="Probability", xlab="Reliabilities", main="Delta=.5")


# compare two likelihoods -------------------------------------------------

xrange <- seq(0, 1, .01)

layout(matrix(1:6, ncol=2, byrow=T))
par(mar=c(2,2,2,2)+.1)

# n=50
y1 <- sapply(xrange, function(x){
  pow(rho1_p=0, rho2=.5, delta=.5, alpha1=x, alpha2=x, n=50)["p"]
})
y2 <- sapply(xrange, function(x){
  pow(rho1_p=.21, rho2=.5, delta=.5, alpha1=x, alpha2=x, n=50)["p"]
})
y <- cbind(y1, y2)
matplot(x=xrange, y=y, type="l", col="black", lwd=2, lty=1, ylim=0:1,
        ylab="Probability", xlab="Reliabilities", main="")
func <- function(x){
  pow(rho1_p=.21, rho2=.5, delta=.5, alpha1=x, alpha2=x, n=50)["p"]/
    pow(rho1_p=0, rho2=.5, delta=.5, alpha1=x, alpha2=x, n=50)["p"]
}
plot(y=sapply(seq(0,1,.01), func), x=seq(0,1,.01), type="l", lwd=2,
     ylab="Likelihood ratio", xlab="Reliabilities", ylim=c(1,20))
abline(h=7, lty=2)

# n=100
y1 <- sapply(xrange, function(x){
  pow(rho1_p=0, rho2=.5, delta=.5, alpha1=x, alpha2=x, n=100)["p"]
})
y2 <- sapply(xrange, function(x){
  pow(rho1_p=.21, rho2=.5, delta=.5, alpha1=x, alpha2=x, n=100)["p"]
})
y <- cbind(y1, y2)
matplot(x=xrange, y=y, type="l", col="black", lwd=2, lty=1, ylim=0:1,
        ylab="Probability", xlab="Reliabilities", main="")
func <- function(x){
  pow(rho1_p=.21, rho2=.5, delta=.5, alpha1=x, alpha2=x, n=100)["p"]/
    pow(rho1_p=0, rho2=.5, delta=.5, alpha1=x, alpha2=x, n=100)["p"]
}
plot(y=sapply(seq(0,1,.01), func), x=seq(0,1,.01), type="l", lwd=2,
     ylab="Likelihood ratio", xlab="Reliabilities", ylim=c(1,20))
abline(h=7, lty=2)

# n=200
y1 <- sapply(xrange, function(x){
  pow(rho1_p=0, rho2=.5, delta=.5, alpha1=x, alpha2=x, n=200)["p"]
})
y2 <- sapply(xrange, function(x){
  pow(rho1_p=.21, rho2=.5, delta=.5, alpha1=x, alpha2=x, n=200)["p"]
})
y <- cbind(y1, y2)
matplot(x=xrange, y=y, type="l", col="black", lwd=2, lty=1, ylim=0:1,
        ylab="Probability", xlab="Reliabilities", main="")
func <- function(x){
  pow(rho1_p=.21, rho2=.5, delta=.5, alpha1=x, alpha2=x, n=200)["p"]/
    pow(rho1_p=0, rho2=.5, delta=.5, alpha1=x, alpha2=x, n=200)["p"]
}
plot(y=sapply(seq(0,1,.01), func), x=seq(0,1,.01), type="l", lwd=2,
     ylab="Likelihood ratio", xlab="Reliabilities", ylim=c(1,20))
abline(h=7, lty=2)

par(mar=c(5,4,4,2)+.1)
layout(1)


# 3 sample sizes ----------------------------------------------------------


# simple sizes leading to power = .5, .8, .95 
# when rho1 = .3 in simple regression where 2nd predictor is junk:
pow(rho1=.3, rho2=0, delta=0, alpha1=1, alpha2=1, n=44)["p"]
pow(rho1=.3, rho2=0, delta=0, alpha1=1, alpha2=1, n=84)["p"]
pow(rho1=.3, rho2=0, delta=0, alpha1=1, alpha2=1, n=136)["p"]
# n = 44, 84, 136

# in the case of simple correlation:
# n = 43, 84, 135

# check with hand calculations
n <- 30
df <- n - 2
rToT <- function(r) sqrt(r^2*df/(1-r^2))
pt(qt(.975, df=df), df=df, ncp=rToT(.3), lower.tail=F) +
  pt(qt(.025, df=df), df=df, ncp=rToT(.3))
y <- sapply(seq(0,1,.01), function(x){
  pt(qt(.975, df=df), df=df, ncp=rToT(x), lower.tail=F) +
    pt(qt(.025, df=df), df=df, ncp=rToT(x))
})
plot(y=y, x=seq(0,1,.01), type="l")

xrange <- 50:500
y <- sapply(xrange, function(x){
  pow(rho1=.21, rho2=0, delta=0, alpha1=1, alpha2=1, n=x)["p"]
})
plot(y=y, x=xrange, type="l")
xrange[which.min((y-.5)^2)] # 88
xrange[which.min((y-.8)^2)] # 175
xrange[which.min((y-.95)^2)] # 287



# explaining likelihood ratio ---------------------------------------------


y1 <- sapply(xrange, function(x){
  pow(rho1_p=0, rho2=.6, delta=.6, alpha1=x, alpha2=x, n=100)["p"]
})
y2 <- sapply(xrange, function(x){
  pow(rho1_p=.3, rho2=.6, delta=.6, alpha1=x, alpha2=x, n=100)["p"]
})
y <- cbind(y1, y2)
# critical LR as function of sample size (log2 scale)
LRcrit2 <- function(n){
  n <- 2^n
  df <- n-2
  (1 + qf(1-alpha, 1, df)/df)^(n/2)
}

png("/Users/Jake/Desktop/likelihoods.png", type="cairo",
    units="in", height=8, width=24, res=200, pointsize=30)
layout(rbind(1:3))
matplot(x=xrange, y=y, type="l", col="black", lwd=3, lty=1,
        ylim=0:1, ylab="Likelihoods",
        xlab=expression(paste("Reliabilities of ",X[1]," and ",X[2])),
        main=expression(paste("Probability of rejecting ",beta[X1]==0)))
text(x=c(.7, .4), y=c(.85, .38),
     labels=c(expression(rho[1.2]==.3), expression(rho[1.2]==0)))
func <- function(x){
  pow(rho1_p=.3, rho2=.6, delta=.6, alpha1=x, alpha2=x, n=100)["p"]/
    pow(rho1_p=0, rho2=.6, delta=.6, alpha1=x, alpha2=x, n=100)["p"]
}
plot(y=sapply(seq(0,1,.01), func), x=seq(0,1,.01), type="l", lwd=3,
     ylim=c(1,20), ylab="Likelihood ratio",
     main=expression(paste("Ratio of probabilities of rejecting ",beta[X1]==0)),
     xlab=expression(paste("Reliabilities of ",X[1]," and ",X[2])))
abline(h=7, lty=2)
alpha <- .05
curve(LRcrit2, xlim=c(4,10), ylim=c(1,20), xaxt="n", lwd=2,
      main=expression(paste("Critical value of likelihood ratio in a\nt-test or simple regression context")),
      #      main=expression(paste("     Critical value of likelihood ratio\nas a function of sample size and ",alpha,"-level")),
      xlab="Sample size", ylab="Likelihood ratio critical value")
axis(side=1, at=4:10, labels=2^(4:10))
alpha <- .1
curve(LRcrit2, add=TRUE, lwd=2)
alpha <- .025
curve(LRcrit2, add=TRUE, lwd=2)
text(x=7, y=c(11.7,6,3), labels=c(expression(alpha==.025),
  expression(alpha==.05), expression(alpha==.1)))
abline(h=7, lty=2)
layout(1)
dev.off()


# plots for true rho1_p = 0 -----------------------------------------------


params <- expand.grid(rho2=c(.3, .5, .7), delta=c(.3, .5, .7))
res <- 200
ps <- 20

png("/Users/Jake/Desktop/true_rho1p_a.png", type="cairo",
    units="in", height=8, width=8, res=res, pointsize=ps)
par(mar=numeric(4))
par(oma=c(4,4,4,3))
layout(matrix(1:9, nrow=3))

size <- 40
mapply(function(rho2, delta, n, rho1_p){
  func <- function(x){
    pow(rho1_p=rho1_p, rho2=rho2, delta=delta, alpha1=x, alpha2=x, n=n)["p"]/
      pow(rho1_p=0, rho2=rho2, delta=delta, alpha1=x, alpha2=x, n=n)["p"]
  }
  plot(y=sapply(seq(0,1,.01), func), x=seq(0,1,.01), type="l", lwd=3,
       ylab="Likelihood ratio", xlab="Reliabilities", ylim=c(1,20),
       yaxt="n", xaxt="n")
  abline(h=7, lty=2)
  if(delta==min(params$delta)) axis(side=2)
  if(rho2==max(params$rho2)){
    axis(side=1,at=seq(0,1,.2),labels=c("0",".2",".4",".6",".8","1"))
  } else axis(side=1, labels=FALSE)
  if(delta==max(params$delta)){
    mtext(bquote(rho[2] == .(rho2)), side=4, line=.5)
  }
  if(rho2==min(params$rho2)){
    mtext(bquote(delta == .(delta)), side=3, line=.5)
  }
}, rho2=params$rho2, delt=params$delta,
   MoreArgs=list(n=size, rho1_p=.3))
mtext("Reliabilities", side=1, line=2.5, outer=T, at=.5)
mtext("Likelihood Ratio", side=2, line=2.5, outer=T, at=.5)
mtext(expression(paste("Low Power (", rho[1.2]==.3,", ", n==40,")")),
      side=3, line=2, outer=T, at=.5, cex=1.3)
dev.off()

png("/Users/Jake/Desktop/true_rho1p_b.png", type="cairo",
    units="in", height=8, width=8, res=res, pointsize=ps)
par(mar=numeric(4))
par(oma=c(4,4,4,3))
layout(matrix(1:9, nrow=3))
size <- 80
mapply(function(rho2, delta, n, rho1_p){
  func <- function(x){
    pow(rho1_p=rho1_p, rho2=rho2, delta=delta, alpha1=x, alpha2=x, n=n)["p"]/
      pow(rho1_p=0, rho2=rho2, delta=delta, alpha1=x, alpha2=x, n=n)["p"]
  }
  plot(y=sapply(seq(0,1,.01), func), x=seq(0,1,.01), type="l", lwd=3,
       ylab="Likelihood ratio", xlab="Reliabilities", ylim=c(1,20),
       yaxt="n", xaxt="n")
  abline(h=7, lty=2)
  if(delta==min(params$delta)) axis(side=2)
  if(rho2==max(params$rho2)){
    axis(side=1,at=seq(0,1,.2),labels=c("0",".2",".4",".6",".8","1"))
  } else axis(side=1, labels=FALSE)
  if(delta==max(params$delta)){
    mtext(bquote(rho[2] == .(rho2)), side=4, line=.5)
  }
  if(rho2==min(params$rho2)){
    mtext(bquote(delta == .(delta)), side=3, line=.5)
  }
}, rho2=params$rho2, delta=params$delta,
MoreArgs=list(n=size, rho1_p=.3))
mtext("Reliabilities", side=1, line=2.5, outer=T, at=.5)
mtext("Likelihood Ratio", side=2, line=2.5, outer=T, at=.5)
mtext(expression(paste("Moderate Power (", rho[1.2]==.3,", ", n==80,")")),
      side=3, line=2, outer=T, at=.5, cex=1.3)
dev.off()

png("/Users/Jake/Desktop/true_rho1p_c.png", type="cairo",
    units="in", height=8, width=8, res=res, pointsize=ps)
par(mar=numeric(4))
par(oma=c(4,4,4,3))
layout(matrix(1:9, nrow=3))
size <- 160
mapply(function(rho2, delta, n, rho1_p){
  func <- function(x){
    pow(rho1_p=rho1_p, rho2=rho2, delta=delta, alpha1=x, alpha2=x, n=n)["p"]/
      pow(rho1_p=0, rho2=rho2, delta=delta, alpha1=x, alpha2=x, n=n)["p"]
  }
  plot(y=sapply(seq(0,1,.01), func), x=seq(0,1,.01), type="l", lwd=3,
       ylab="Likelihood ratio", xlab="Reliabilities", ylim=c(1,20),
       yaxt="n", xaxt="n")
  abline(h=7, lty=2)
  if(delta==min(params$delta)) axis(side=2)
  if(rho2==max(params$rho2)){
    axis(side=1,at=seq(0,1,.2),labels=c("0",".2",".4",".6",".8","1"))
  } else axis(side=1, labels=FALSE)
  if(delta==max(params$delta)){
    mtext(bquote(rho[2] == .(rho2)), side=4, line=.5)
  }
  if(rho2==min(params$rho2)){
    mtext(bquote(delta == .(delta)), side=3, line=.5)
  }
}, rho2=params$rho2, delta=params$delta,
MoreArgs=list(n=size, rho1_p=.3))
mtext("Reliabilities", side=1, line=2.5, outer=T, at=.5)
mtext("Likelihood Ratio", side=2, line=2.5, outer=T, at=.5)
mtext(expression(paste("High Power (", rho[1.2]==.3,", ", n==160,")")),
      side=3, line=2, outer=T, at=.5, cex=1.3)
layout(1)
dev.off()


# plots for true delta = 0 ------------------------------------------------


# varying reliabilities
png("/Users/Jake/Desktop/SIC.png", type="cairo",
    units="in", height=8, width=12, res=200, pointsize=20)
par(mar=c(4,4,1,1)+.1)
layout(rbind(1:3, 4:6))
xrange <- seq(0, .99999, length.out=100)
### plot likelihoods
sapply(c(80,160,320), function(n){
  y1 <- sapply(xrange, function(x){
    pow(rho1=.3, rho2=.3, delta=1, alpha1=x, alpha2=x, n=n)["p2"]
  })
  y2 <- sapply(xrange, function(x){
    pow(rho1=.3, rho2=.3, delta=.8, alpha1=x, alpha2=x, n=n)["p2"]
  })
  y3 <- sapply(xrange, function(x){
    pow(rho1=.3, rho2=.3, delta=.5, alpha1=x, alpha2=x, n=n)["p2"]
  })
  y4 <- sapply(xrange, function(x){
    pow(rho1=.3, rho2=.3, delta=.2, alpha1=x, alpha2=x, n=n)["p2"]
  })
  y <- cbind(y1,y2,y3,y4)
  matplot(x=xrange, y=y, type="l", col="black",
          lwd=2, lty=1, ylim=0:1, ylab="Likelihoods",
          xlab="Reliabilities",
          main="")
})
### plot likelihood ratios
sapply(c(80,160,320), function(n){
  y1 <- sapply(xrange, function(x){
    pow(rho1=.3, rho2=.3, delta=.8, alpha1=x, alpha2=x, n=n)["p2"]/
      pow(rho1=.3, rho2=.3, delta=1, alpha1=x, alpha2=x, n=n)["p2"]
  })
  y2 <- sapply(xrange, function(x){
    pow(rho1=.3, rho2=.3, delta=.5, alpha1=x, alpha2=x, n=n)["p2"]/
      pow(rho1=.3, rho2=.3, delta=1, alpha1=x, alpha2=x, n=n)["p2"]
  })
  y3 <- sapply(xrange, function(x){
    pow(rho1=.3, rho2=.3, delta=.2, alpha1=x, alpha2=x, n=n)["p2"]/
      pow(rho1=.3, rho2=.3, delta=1, alpha1=x, alpha2=x, n=n)["p2"]
  })
  y <- cbind(y1, y2, y3)
  matplot(x=xrange, y=log(y), type="l", col="black", yaxt="n",
          ylim=c(log(0.1071009), log(674.2675)), lwd=2, lty=1,
          ylab="Likelihood Ratios (log scale)", xlab="Reliabilities",
          main="")
  axis(side=2, at=log(c(1/7,1,7,20,100,500)),
       labels=c("1/7","1","7","20","100","500"))
  abline(h=log(c(1/7, 7)), lty=2)
})
dev.off()


# plots for true rho1.2 = rho2.1 ------------------------------------------


test <- expand.grid(x=c(.1,.3,.5), y=c(0,.2,.4))
test <- within(test, y <- x+y)
cbind(cos(acos(test$x)+acos(test$y)),cos(acos(test$x)-acos(test$y)))

xrange <- seq(0, 1, length.out=100)

### delta = .3
sapply(c(.1,.3,.5), function(r2){
  y1 <- sapply(xrange, function(x){
    pow(rho1_p=r2, rho2_p=r2, delta=.3, alpha1=x, alpha2=x, n=80)["p1"]
  })
  y2 <- sapply(xrange, function(x){
    pow(rho1_p=r2+.2, rho2_p=r2, delta=.3, alpha1=x, alpha2=x, n=80)["p1"]
  })
  y3 <- sapply(xrange, function(x){
    pow(rho1_p=r2+.4, rho2_p=r2, delta=.3, alpha1=x, alpha2=x, n=80)["p1"]
  })
  y <- cbind(y1,y2,y3)
  matplot(x=xrange, y=y, type="l", col="black",
          lwd=2, lty=1, ylim=0:1, ylab="Likelihoods",
          xlab="Reliabilities",
          main="")
})
### delta = .5
sapply(c(.1,.3,.5), function(r2){
  y1 <- sapply(xrange, function(x){
    pow(rho1_p=r2, rho2_p=r2, delta=.5, alpha1=x, alpha2=x, n=80)["p1"]
  })
  y2 <- sapply(xrange, function(x){
    pow(rho1_p=r2+.2, rho2_p=r2, delta=.5, alpha1=x, alpha2=x, n=80)["p1"]
  })
  y3 <- sapply(xrange, function(x){
    pow(rho1_p=r2+.4, rho2_p=r2, delta=.5, alpha1=x, alpha2=x, n=80)["p1"]
  })
  y <- cbind(y1,y2,y3)
  matplot(x=xrange, y=y, type="l", col="black",
          lwd=2, lty=1, ylim=0:1, ylab="Likelihoods",
          xlab="Reliabilities",
          main="")
})
### delta = .7
sapply(c(.1,.3,.5), function(r2){
  y1 <- sapply(xrange, function(x){
    pow(rho1_p=r2, rho2_p=r2, delta=.7, alpha1=x, alpha2=x, n=80)["p1"]
  })
  y2 <- sapply(xrange, function(x){
    pow(rho1_p=r2+.2, rho2_p=r2, delta=.7, alpha1=x, alpha2=x, n=80)["p1"]
  })
  y3 <- sapply(xrange, function(x){
    pow(rho1_p=r2+.4, rho2_p=r2, delta=.7, alpha1=x, alpha2=x, n=80)["p1"]
  })
  y <- cbind(y1,y2,y3)
  matplot(x=xrange, y=y, type="l", col="black",
          lwd=2, lty=1, ylim=0:1, ylab="Likelihoods",
          xlab="Reliabilities",
          main="")
})


# fig 1, take 2 -----------------------------------------------------------


png("/Users/Jake/Desktop/fig1.png", type="cairo",
    units="in", height=8, width=8, res=200, pointsize=16)
par(mar=c(3.5,3.5,2.5,.5)+.1)
layout(matrix(1:4, nrow=2))

# vary rho1_p, probabilities
xrange <- seq(0, 1, length.out=100)
y1 <- sapply(xrange, function(x){
  pow(rho1_p=0, rho2=.5, delta=.5, alpha1=x, alpha2=x, n=100)["p"]
})
y2 <- sapply(xrange, function(x){
  pow(rho1_p=.2, rho2=.5, delta=.5, alpha1=x, alpha2=x, n=100)["p"]
})
y3 <- sapply(xrange, function(x){
  pow(rho1_p=.5, rho2=.5, delta=.5, alpha1=x, alpha2=x, n=100)["p"]
})
y4 <- sapply(xrange, function(x){
  pow(rho1_p=.8, rho2=.5, delta=.5, alpha1=x, alpha2=x, n=100)["p"]
})
y <- cbind(y1,y2,y3,y4)
matplot(x=xrange, y=y, type="l", col="black", lwd=2, lty=1,
        ylim=0:1, ylab="", xlab="", main="")
mtext("Probability", side=2, line=2.25, cex=.85)
mtext(expression(paste("Reliabilities of ",X[1]," and ",X[2])),
      side=1, line=2.5, cex=.85)
mtext(expression(paste("Probability of rejecting ",beta[X1]==0)),
      side=3, line=.5, cex=.9)
# text(x=c(.7, .4), y=c(.85, .38),
#      labels=c(expression(rho[1.2]==.3), expression(rho[1.2]==0)))

# vary rho1_p, likelihood ratios
y1 <- sapply(xrange, function(x){
  pow(rho1_p=.2, rho2=.5, delta=.5, alpha1=x, alpha2=x, n=100)["p"]/
    pow(rho1_p=0, rho2=.5, delta=.5, alpha1=x, alpha2=x, n=100)["p"]
})
y2 <- sapply(xrange, function(x){
  pow(rho1_p=.5, rho2=.5, delta=.5, alpha1=x, alpha2=x, n=100)["p"]/
    pow(rho1_p=0, rho2=.5, delta=.5, alpha1=x, alpha2=x, n=100)["p"]
})
y3 <- sapply(xrange, function(x){
  pow(rho1_p=.8, rho2=.5, delta=.5, alpha1=x, alpha2=x, n=100)["p"]/
    pow(rho1_p=0, rho2=.5, delta=.5, alpha1=x, alpha2=x, n=100)["p"]
})
y <- cbind(y1,y2,y3)
matplot(x=xrange, y=y, type="l", col="black", lwd=2, lty=1,
        ylim=c(1,20), ylab="", xlab="", main="")
abline(h=c(7,8), lty=2)
mtext("Likelihood ratio", side=2, line=2.25, cex=.85)
mtext(expression(paste("Reliabilities of ",X[1]," and ",X[2])),
      side=1, line=2.5, cex=.85)
mtext(expression(paste("Likelihood ratio relative to ",rho[1.2]==0)),
      side=3, line=.5, cex=.9)

# vary delta, probabilities
xrange <- seq(0, .99999, length.out=100)
y1 <- sapply(xrange, function(x){
  pow(rho1=.5, rho2=.5, delta=1, alpha1=x, alpha2=x, n=100)["p2"]
})
y2 <- sapply(xrange, function(x){
  pow(rho1=.5, rho2=.5, delta=.8, alpha1=x, alpha2=x, n=100)["p2"]
})
y3 <- sapply(xrange, function(x){
  pow(rho1=.5, rho2=.5, delta=.5, alpha1=x, alpha2=x, n=100)["p2"]
})
y4 <- sapply(xrange, function(x){
  pow(rho1=.5, rho2=.5, delta=.2, alpha1=x, alpha2=x, n=100)["p2"]
})
y <- cbind(y1,y2,y3,y4)
matplot(x=xrange, y=y, type="l", col="black",
        lwd=2, lty=1, ylim=0:1, ylab="", xlab="", main="")
mtext("Probability", side=2, line=2.25, cex=.85)
mtext(expression(paste("Reliabilities of ",X[1]," and ",X[2])),
      side=1, line=2.5, cex=.85)
mtext(expression(paste("Probability of rejecting ",beta[X1]," = ",beta[X2]," = 0")),
      side=3, line=.5, cex=.9)

# vary delta, likelihood ratios
y1 <- sapply(xrange, function(x){
  pow(rho1=.5, rho2=.5, delta=.8, alpha1=x, alpha2=x, n=100)["p2"]/
    pow(rho1=.5, rho2=.5, delta=1, alpha1=x, alpha2=x, n=100)["p2"]
})
y2 <- sapply(xrange, function(x){
  pow(rho1=.5, rho2=.5, delta=.5, alpha1=x, alpha2=x, n=100)["p2"]/
    pow(rho1=.5, rho2=.5, delta=1, alpha1=x, alpha2=x, n=100)["p2"]
})
y3 <- sapply(xrange, function(x){
  pow(rho1=.5, rho2=.5, delta=.2, alpha1=x, alpha2=x, n=100)["p2"]/
    pow(rho1=.5, rho2=.5, delta=1, alpha1=x, alpha2=x, n=100)["p2"]
})
y <- cbind(y1, y2, y3)
matplot(x=xrange, y=log(y), type="l", col="black", yaxt="n", ylab="",
        ylim=c(log(1), log(max(y3))), lwd=2, lty=1, xlab="", main="")
axis(side=2, at=log(c(1,20,100,500)),
     labels=c("1","20","100","500"))
abline(h=log(c(7,8)), lty=2)
mtext("Likelihood ratio (log scale)", side=2, line=2.25, cex=.85)
mtext(expression(paste("Reliabilities of ",X[1]," and ",X[2])),
      side=1, line=2.5, cex=.85)
mtext(expression(paste("Likelihood ratio relative to ",delta==1)),
      side=3, line=.5, cex=.9)

dev.off()


# contour plots -----------------------------------------------------------


# round(exp(seq(log(20), log(300), length.out=7)))
# # [1]  20  31  49  77 122 191 300
# log(c(20,30,50,75,120,200,300))
# # [1] 2.995732 3.401197 3.912023 4.317488 4.787492 5.298317 5.703782

z <- list(z1,z2,z3,z4,z5,z6,z7,z8,z9,
          z10,z11,z12,z13,z14,z15,z16,z17,z18)
saveRDS(z, file="/Users/Jake/Desktop/Dropbox/Dropbox/SIC/z.rds")

z <- readRDS("/Users/Jake/Desktop/Dropbox/Dropbox/SIC/z.rds")
for(i in 1:18) assign(paste0("z",i), z[[i]])

# plot for incremental validity -------------------------------------------

nVals <- c(seq(20,100,10), 200, 300)
relRange <- seq(0, 1, length.out=50)
nRange <- seq(log(20), log(300), length.out=50)
samp <- expand.grid(n=nRange, rel=relRange)
# system.time({
#   z1 <- matrix(apply(samp, 1, function(x) pow(rho1_p=0, rho2=.3, delta=.3,
#     alpha1=x["rel"], alpha2=x["rel"], n=round(exp(x["n"])))["p"]),
#     ncol=50, byrow=TRUE)
#   z2 <- matrix(apply(samp, 1, function(x) pow(rho1_p=0, rho2=.5, delta=.5,
#     alpha1=x["rel"], alpha2=x["rel"], n=round(exp(x["n"])))["p"]),
#     ncol=50, byrow=TRUE)
#   z3 <- matrix(apply(samp, 1, function(x) pow(rho1_p=0, rho2=.7, delta=.7,
#     alpha1=x["rel"], alpha2=x["rel"], n=round(exp(x["n"])))["p"]),
#     ncol=50, byrow=TRUE)
# })
z1[45:50,] <- rbind(rep(.051, times=50), rep(.051, times=50),
                    rep(.051, times=50), rep(.051, times=50),
                    rep(.051, times=50), rep(.04999, times=50))
z1[1:2,] <- rbind(rep(.04999, times=50), rep(.051, times=50))
z2[49:50,] <- rbind(rep(.051, times=50), rep(.04999, times=50))
z2[1,] <- rep(.051, times=50)
z3[50,] <- rep(.04999, times=50)
z3[1,] <- rep(.051, times=50)

tiff(paste0(path,"Fig3.tif"),
    units="in", height=8, width=8*3, res=100, pointsize=35)
layout(rbind(1:3))
par(mar=c(4, 4, 2, 1)+.1)
par(oma=c(0,0,2,0))
contour(y=nRange, x=relRange, z=z1, yaxt="n", lwd=2, labcex=.55,
        levels=c(.05, .06, .07, .15, seq(.1, 1, .1)),
        xlab=expression(paste("Reliabilities ",(alpha[1]==alpha[2]))),
        ylab="Sample size (log scale)",
        main=expression(paste("A: ",delta," = ",rho[2]," = .3")))
axis(side=2, at=log(nVals), labels=FALSE, las=1)
mtext(c("20","50","100","300"), side=2, at=log(c(20,50,100,300)), line=1, cex=.75)
contour(y=nRange, x=relRange, z=z2, yaxt="n", lwd=2, labcex=.55,
        levels=c(.05, .06, .07, .15, seq(.1, 1, .1)),
        xlab=expression(paste("Reliabilities ",(alpha[1]==alpha[2]))),
        ylab="Sample size (log scale)",
        main=expression(paste("B: ",delta," = ",rho[2]," = .5")))
axis(side=2, at=log(nVals), labels=FALSE, las=1)
mtext(c("20","50","100","300"), side=2, at=log(c(20,50,100,300)), line=1, cex=.75)
contour(y=nRange, x=relRange, z=z3, yaxt="n", lwd=2, labcex=.55,
        levels=c(.05, .06, .15, seq(.1, 1, .1)),
        xlab=expression(paste("Reliabilities ",(alpha[1]==alpha[2]))),
        ylab="Sample size (log scale)",
        main=expression(paste("C: ",delta," = ",rho[2]," = .7")))
axis(side=2, at=log(nVals), labels=FALSE, las=1)
mtext(c("20","50","100","300"), side=2, at=log(c(20,50,100,300)), line=1, cex=.75)
mtext(expression(paste("Type 1 error rates for rejecting ",rho[1.2]," = 0 if ",beta[X1]," significant")),
      side=3, line=0, outer=TRUE, at=.5)
dev.off()

# plot for SIC ------------------------------------------------------------

nVals <- c(seq(20,100,10), 200, 300)
relRange <- seq(0, .999, length.out=50)
relRange2 <- c(relRange, c(.999999, .9999999))
nRange <- seq(log(20), log(300), length.out=50)
samp <- expand.grid(n=nRange, rel=relRange)
# system.time({
#   z4 <- matrix(apply(samp, 1, function(x) pow(rho1=.3, rho2=.3, delta=1,
#     alpha1=x["rel"], alpha2=x["rel"], n=round(exp(x["n"])))["p2"]),
#     ncol=50, byrow=TRUE)
#   z5 <- matrix(apply(samp, 1, function(x) pow(rho1=.5, rho2=.5, delta=1,
#     alpha1=x["rel"], alpha2=x["rel"], n=round(exp(x["n"])))["p2"]),
#     ncol=50, byrow=TRUE)
#   z6 <- matrix(apply(samp, 1, function(x) pow(rho1=.7, rho2=.7, delta=1,
#     alpha1=x["rel"], alpha2=x["rel"], n=round(exp(x["n"])))["p2"]),
#     ncol=50, byrow=TRUE)
# })
add <- rbind(rep(c(.05,.049), each=25), rep(c(.051,.05), each=25))
z4 <- rbind(z4, add)
z5 <- rbind(z5, add)
z6 <- rbind(z6, add)

tiff(paste0(path,"Fig4.tif"),
    units="in", height=8, width=8*3, res=100, pointsize=35)
layout(rbind(1:3))
par(mar=c(4, 4, 2, 1)+.1)
par(oma=c(0,0,2,0))
contour(y=nRange, x=relRange2, z=z4, yaxt="n", lwd=2, labcex=.55,
        levels=c(.01, .05, seq(.1, 1, .1)), las=1,
        xlab=expression(paste("Reliabilities ",(alpha[1]==alpha[2]))),
        ylab="Sample size (log scale)",
        main=expression(paste("A: ",rho[1]," = ",rho[2]," = .3")))
axis(side=2, at=log(nVals), labels=FALSE, las=1)
mtext(c("20","50","100","300"), side=2, at=log(c(20,50,100,300)), line=1, cex=.75)
contour(y=nRange, x=relRange2, z=z5, yaxt="n", lwd=2, labcex=.55,
        levels=c(.01, .05, seq(.1, 1, .1)), las=1,
        xlab=expression(paste("Reliabilities ",(alpha[1]==alpha[2]))),
        ylab="Sample size (log scale)",
        main=expression(paste("B: ",rho[1]," = ",rho[2]," = .5")))
axis(side=2, at=log(nVals), labels=FALSE, las=1)
mtext(c("20","50","100","300"), side=2, at=log(c(20,50,100,300)), line=1, cex=.75)
contour(y=nRange, x=relRange2, z=z6, yaxt="n", lwd=2, labcex=.55,
        levels=c(.01, .05, seq(.1, 1, .1)), las=1,
        xlab=expression(paste("Reliabilities ",(alpha[1]==alpha[2]))),
        ylab="Sample size (log scale)",
        main=expression(paste("C: ",rho[1]," = ",rho[2]," = .7")))
axis(side=2, at=log(nVals), labels=FALSE, las=1)
mtext(c("20","50","100","300"), side=2, at=log(c(20,50,100,300)), line=1, cex=.75)
mtext(expression(paste("Type 1 error rates for rejecting ",rho[1.2]," = 0 OR ",rho[2.1],"= 0 if both predictors significant")),
      side=3, line=0, outer=TRUE, at=.5)
dev.off()

# plot for improved measurement -------------------------------------------

nVals <- c(20,30,50,75,120,200,300)
relRange <- seq(0, .999, length.out=60)
nRange <- seq(log(20), log(300), length.out=60)
samp <- expand.grid(n=nRange, rel=relRange)
# system.time({
#   z7 <- matrix(apply(samp, 1, function(x) pow(rho1_p=.1, rho2_p=.1, delta_p=0,
#     alpha1=x["rel"], alpha2=x["rel"], n=round(exp(x["n"])))["p1"]),
#     ncol=60, byrow=TRUE)
#   z8 <- matrix(apply(samp, 1, function(x) pow(rho1_p=.3, rho2_p=.3, delta_p=0,
#     alpha1=x["rel"], alpha2=x["rel"], n=round(exp(x["n"])))["p1"]),
#     ncol=60, byrow=TRUE)
#   z9 <- matrix(apply(samp, 1, function(x) pow(rho1_p=.5, rho2_p=.5, delta_p=0,
#     alpha1=x["rel"], alpha2=x["rel"], n=round(exp(x["n"])))["p1"]),
#     ncol=60, byrow=TRUE)
# })

png("/Users/Jake/Desktop/Dropbox/Dropbox/SIC/IM_contour.png", type="cairo",
    units="in", height=8, width=8*3, res=200, pointsize=31.5)
layout(rbind(1:3))
par(mar=c(4, 4, 2, 1)+.1)
par(oma=c(0,0,2,0))
contour(y=nRange, x=relRange, z=z7, yaxt="n", lwd=2, labcex=.65,
        levels=c(.01, .05, seq(.1, 1, .1)),
        xlab=expression(paste("Reliabilities of ",X[1]," and ",X[2])),
        ylab="Sample size (log scale)",
        main="Y,T partial correlations = .1")
axis(side=2, at=log(nVals), labels=nVals)
contour(y=nRange, x=relRange, z=z8, yaxt="n", lwd=2, labcex=.65,
        levels=c(.01, .05, seq(.1, 1, .1)),
        xlab=expression(paste("Reliabilities of ",X[1]," and ",X[2])),
        ylab="Sample size (log scale)",
        main="Y,T partial correlations = .3")
axis(side=2, at=log(nVals), labels=nVals)
contour(y=nRange, x=relRange, z=z9, yaxt="n", lwd=2, labcex=.65,
        levels=c(.01, .05, seq(.1, 1, .1)),
        xlab=expression(paste("Reliabilities of ",X[1]," and ",X[2])),
        ylab="Sample size (log scale)",
        main="Y,T partial correlations = .5")
axis(side=2, at=log(nVals), labels=nVals)
mtext(expression(paste("Type 1 error rates for testing Y,",T[1]," partial correlation = Y,",T[2]," partial correlation")),
      side=3, line=0, outer=TRUE, at=.5)
dev.off()

# 3x3 plot for IM ---------------------------------------------------------

nVals <- c(20,30,50,75,120,200,300)
relRange <- seq(0, .999, length.out=50)
nRange <- seq(log(20), log(300), length.out=50)
samp <- expand.grid(n=nRange, rel=relRange)
# system.time({ # elapsed time = 7.5 minutes
#   z10 <- matrix(apply(samp, 1, function(x) pow(rho1_p=.1, rho2_p=.1, delta=.3,
#     alpha1=x["rel"], alpha2=x["rel"], n=round(exp(x["n"])))["p1"]),
#     ncol=50, byrow=TRUE)
#   z11 <- matrix(apply(samp, 1, function(x) pow(rho1_p=.1, rho2_p=.1, delta=.5,
#     alpha1=x["rel"], alpha2=x["rel"], n=round(exp(x["n"])))["p1"]),
#     ncol=50, byrow=TRUE)
#   z12 <- matrix(apply(samp, 1, function(x) pow(rho1_p=.1, rho2_p=.1, delta=.7,
#     alpha1=x["rel"], alpha2=x["rel"], n=round(exp(x["n"])))["p1"]),
#     ncol=50, byrow=TRUE)
#   z13 <- matrix(apply(samp, 1, function(x) pow(rho1_p=.3, rho2_p=.3, delta=.3,
#     alpha1=x["rel"], alpha2=x["rel"], n=round(exp(x["n"])))["p1"]),
#     ncol=50, byrow=TRUE)
#   z14 <- matrix(apply(samp, 1, function(x) pow(rho1_p=.3, rho2_p=.3, delta=.5,
#     alpha1=x["rel"], alpha2=x["rel"], n=round(exp(x["n"])))["p1"]),
#     ncol=50, byrow=TRUE)
#   z15 <- matrix(apply(samp, 1, function(x) pow(rho1_p=.3, rho2_p=.3, delta=.7,
#     alpha1=x["rel"], alpha2=x["rel"], n=round(exp(x["n"])))["p1"]),
#     ncol=50, byrow=TRUE)
#   z16 <- matrix(apply(samp, 1, function(x) pow(rho1_p=.5, rho2_p=.5, delta=.3,
#     alpha1=x["rel"], alpha2=x["rel"], n=round(exp(x["n"])))["p1"]),
#     ncol=50, byrow=TRUE)
#   z17 <- matrix(apply(samp, 1, function(x) pow(rho1_p=.5, rho2_p=.5, delta=.5,
#     alpha1=x["rel"], alpha2=x["rel"], n=round(exp(x["n"])))["p1"]),
#     ncol=50, byrow=TRUE)
#   z18 <- matrix(apply(samp, 1, function(x) pow(rho1_p=.5, rho2_p=.5, delta=.7,
#     alpha1=x["rel"], alpha2=x["rel"], n=round(exp(x["n"])))["p1"]),
#     ncol=50, byrow=TRUE)
# })

png("/Users/Jake/Desktop/Dropbox/Dropbox/SIC/IM_contour_3x3.png", type="cairo",
    units="in", height=8*3, width=8*3, res=200, pointsize=31.5)
layout(matrix(1:9, nrow=3, ncol=3))
par(mar=c(4, 4, 2, 1)+.1)
par(oma=c(0,0,2,0))
contour(y=nRange, x=relRange, z=z10, yaxt="n", lwd=2, labcex=.65,
        levels=c(.01, .05, seq(.1, 1, .1)),
        xlab=expression(paste("Reliabilities of ",X[1]," and ",X[2])),
        ylab="Sample size (log scale)",
        main="Y,T partial correlations = .1")
axis(side=2, at=log(nVals), labels=nVals)
contour(y=nRange, x=relRange, z=z11, yaxt="n", lwd=2, labcex=.65,
        levels=c(.01, .05, seq(.1, 1, .1)),
        xlab=expression(paste("Reliabilities of ",X[1]," and ",X[2])),
        ylab="Sample size (log scale)",
        main="Y,T partial correlations = .3")
axis(side=2, at=log(nVals), labels=nVals)
contour(y=nRange, x=relRange, z=z12, yaxt="n", lwd=2, labcex=.65,
        levels=c(.01, .05, seq(.1, 1, .1)),
        xlab=expression(paste("Reliabilities of ",X[1]," and ",X[2])),
        ylab="Sample size (log scale)",
        main="Y,T partial correlations = .5")
axis(side=2, at=log(nVals), labels=nVals)
contour(y=nRange, x=relRange, z=z13, yaxt="n", lwd=2, labcex=.65,
        levels=c(.01, .05, seq(.1, 1, .1)),
        xlab=expression(paste("Reliabilities of ",X[1]," and ",X[2])),
        ylab="Sample size (log scale)",
        main="Y,T partial correlations = .1")
axis(side=2, at=log(nVals), labels=nVals)
contour(y=nRange, x=relRange, z=z14, yaxt="n", lwd=2, labcex=.65,
        levels=c(.01, .05, seq(.1, 1, .1)),
        xlab=expression(paste("Reliabilities of ",X[1]," and ",X[2])),
        ylab="Sample size (log scale)",
        main="Y,T partial correlations = .3")
axis(side=2, at=log(nVals), labels=nVals)
contour(y=nRange, x=relRange, z=z15, yaxt="n", lwd=2, labcex=.65,
        levels=c(.01, .05, seq(.1, 1, .1)),
        xlab=expression(paste("Reliabilities of ",X[1]," and ",X[2])),
        ylab="Sample size (log scale)",
        main="Y,T partial correlations = .5")
axis(side=2, at=log(nVals), labels=nVals)
contour(y=nRange, x=relRange, z=z16, yaxt="n", lwd=2, labcex=.65,
        levels=c(.01, .05, seq(.1, 1, .1)),
        xlab=expression(paste("Reliabilities of ",X[1]," and ",X[2])),
        ylab="Sample size (log scale)",
        main="Y,T partial correlations = .1")
axis(side=2, at=log(nVals), labels=nVals)
contour(y=nRange, x=relRange, z=z17, yaxt="n", lwd=2, labcex=.65,
        levels=c(.01, .05, seq(.1, 1, .1)),
        xlab=expression(paste("Reliabilities of ",X[1]," and ",X[2])),
        ylab="Sample size (log scale)",
        main="Y,T partial correlations = .3")
axis(side=2, at=log(nVals), labels=nVals)
contour(y=nRange, x=relRange, z=z18, yaxt="n", lwd=2, labcex=.65,
        levels=c(.01, .05, seq(.1, 1, .1)),
        xlab=expression(paste("Reliabilities of ",X[1]," and ",X[2])),
        ylab="Sample size (log scale)",
        main="Y,T partial correlations = .5")
axis(side=2, at=log(nVals), labels=nVals)
mtext(expression(paste("Type 1 error rates for testing Y,",T[1]," partial correlation = Y,",T[2]," partial correlation")),
      side=3, line=0, outer=TRUE, at=.5)
dev.off()

# 1x3 plot for IM ---------------------------------------------------------

nVals <- c(seq(20,100,10), 200, 300)

tiff(paste0(path,"Fig5.tif"), units="in",
    height=8, width=8*3, res=100, pointsize=35)
layout(rbind(1:3))
par(mar=c(4, 4, 2, 1)+.1)
par(oma=c(0,0,2,0))
contour(y=nRange, x=relRange, z=z11, yaxt="n", lwd=2, labcex=.65,
        levels=c(.01, .05, seq(.1, 1, .1)),
        xlab=expression(paste("Reliabilities ",(alpha[1]==alpha[2]))),
        ylab="Sample size (log scale)",
        main=expression(paste("A: ", rho[1.2]," = ",rho[2.1]," = ",.1)))
axis(side=2, at=log(nVals), labels=FALSE)
mtext(c("20","50","100","300"), side=2, at=log(c(20,50,100,300)), line=1, cex=.75)
contour(y=nRange, x=relRange, z=z14, yaxt="n", lwd=2, labcex=.65,
        levels=c(.01, .05, seq(.1, 1, .1)),
        xlab=expression(paste("Reliabilities ",(alpha[1]==alpha[2]))),
        ylab="Sample size (log scale)",
        main=expression(paste("B: ", rho[1.2]," = ",rho[2.1]," = ",.3)))
axis(side=2, at=log(nVals), labels=FALSE)
mtext(c("20","50","100","300"), side=2, at=log(c(20,50,100,300)), line=1, cex=.75)
contour(y=nRange, x=relRange, z=z17, yaxt="n", lwd=2, labcex=.65,
        levels=c(.01, .05, seq(.1, 1, .1)),
        xlab=expression(paste("Reliabilities ",(alpha[1]==alpha[2]))),
        ylab="Sample size (log scale)",
        main=expression(paste("C: ", rho[1.2]," = ",rho[2.1]," = ",.5)))
axis(side=2, at=log(nVals), labels=FALSE)
mtext(c("20","50","100","300"), side=2, at=log(c(20,50,100,300)), line=1, cex=.75)
mtext(expression(paste("Type 1 error rates for rejecting ",rho[1.2]," = ",rho[2.1]," if only one predictor significant")),
      side=3, line=0, outer=TRUE, at=.5)
dev.off()

# ice cream example -------------------------------------------------------

fToC <- function(f) (f-32)/1.8
n <- 120

# set.seed(4624635)
set.seed(56246)
tempF <- rnorm(n, mean=80, sd=10)
plot(density(tempF))
heat <- sqrt(.5)*scale(tempF) + sqrt(.5)*rnorm(n)
heat <- round(6*(heat - min(heat))/diff(range(heat))+1)

tiff(paste0(path,"Fig1.tif"), type="cairo",
    units="in", height=4, width=4, res=400, pointsize=25/2)
par(mar=c(4, 4, 1, 1)+.1)
plot(heat ~ tempF, pch=19, xaxt="n",
  xlab="Recorded temperature", ylab="Subjective heat")
text(y=6.5, x=60, labels="Reliability = 0.40")
axis(side=1, labels=FALSE)
mtext(paste(seq(50,100,10),"F\n",round(fToC(seq(50,100,10))),"C"),
  side=1, at=seq(50,100,10), line=1.6, cex=.9)
dev.off()

cor(heat, tempF)^2

deaths <- rpois(n, lambda=exp(1 + .05*(tempF - mean(tempF))))
plot(table(deaths))
plot(deaths ~ tempF)
cor(deaths, tempF)

sales <- round(1000 + 20*(tempF - mean(tempF)) + rnorm(n, sd=200))
plot(density(sales))
plot(sales ~ tempF)
cor(sales, tempF)

cor(cbind(deaths, sales, tempF))
cor(cbind(deaths, sales, heat))

pcor(cbind(deaths, sales, tempF))
pcor(cbind(deaths, sales, heat))

tiff(paste0(path,"Fig2.tif"),
    units="in", height=8, width=8*3, res=100, pointsize=35)
layout(rbind(1:3))
par(mar=c(4, 4, 3, 2.25)+.1)
plot(deaths ~ sales, pch=19, mgp=c(2,1,0), 
  ylab="Swimming pool deaths per day",
  xlab="Ice cream cones sold per day",
  main="A: Simple correlation\nr = 0.49, p < .001")
abline(lm(deaths ~ sales), lwd=4)
y <- resid(lm(deaths ~ heat))
x <- resid(lm(sales ~ heat))
plot(y ~ x, pch=19, mgp=c(2,1,0),
  ylab="Swimming pool deaths per day\n(adjusted for heat)",
  xlab="Ice cream cones sold (adjusted for heat)",
  main="B: Controlling for subjective heat\nPartial r = 0.33, p < .001")
abline(lm(y ~ x), lwd=4)
y <- resid(lm(deaths ~ tempF))
x <- resid(lm(sales ~ tempF))
plot(y ~ x, pch=19, mgp=c(2,1,0),
  ylab="Swimming pool deaths per day\n(adjusted for temperature)",
  xlab="Ice cream cones sold (adjusted for temperature)",
  main="C: Controlling for recorded temperature\nPartial r = -0.02, p = .81")
abline(lm(y ~ x), lwd=4)
dev.off()

# two measures of temperature ---------------------------------------------

set.seed(567453)
tempEst <- round(tempF + rnorm(n, sd=15))
plot(tempEst ~ tempF)
summary(lm(tempEst ~ tempF))

summary(lm(deaths ~ heat + tempEst))
pcor(cbind(deaths, heat, tempEst))
avPlots(lm(deaths ~ heat + tempEst))



# toy sim, IV -------------------------------------------------------------

# controlling for heat
sim1 <- function(n, alpha){
  tempF <- rnorm(n, mean=80, sd=10)
  heat <- sqrt(alpha)*scale(tempF) + sqrt(1-alpha)*rnorm(n)
  heat <- round(6*(heat - min(heat))/diff(range(heat))+1)
  deaths <- rpois(n, lambda=exp(1 + .05*(tempF - mean(tempF))))
  sales <- round(1000 + 20*(tempF - mean(tempF)) + rnorm(n, sd=200))
  coef(summary(lm(deaths ~ sales + heat)))["sales","Pr(>|t|)"] < .05
}
mean(replicate(10000, sim1(n=120, alpha=.5)))
# 0.9237

# controlling for temperature
sim2 <- function(n, alpha){
  tempF <- rnorm(n, mean=80, sd=10)
  deaths <- rpois(n, lambda=exp(1 + .05*(tempF - mean(tempF))))
  sales <- round(1000 + 20*(tempF - mean(tempF)) + rnorm(n, sd=200))
  coef(summary(lm(deaths ~ sales + tempF)))["sales","Pr(>|t|)"] < .05
}
mean(replicate(10000, sim2(n=120, alpha=.5)))


# toy sim, IC -------------------------------------------------------------

sim3 <- function(n, alpha){
  tempF <- rnorm(n, mean=80, sd=10)
  heat <- sqrt(alpha)*scale(tempF) + sqrt(1-alpha)*rnorm(n)
  heat <- round(6*(heat - min(heat))/diff(range(heat))+1)
  tempEst <- round(tempF + rnorm(n, sd=15))
  deaths <- rpois(n, lambda=exp(1 + .05*(tempF - mean(tempF))))
  all(coef(summary(lm(deaths ~ heat + tempEst)))[c("heat","tempEst"),"Pr(>|t|)"] < .05)
}
mean(replicate(10000, sim3(n=120, alpha=.5)))
# 0.682
