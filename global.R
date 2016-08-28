library(mvtnorm) # for pmvt
library(corpcor) # for cor2pcor

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
  # x1-x2 partial correlation
  pcorr3 <- (r["x1","x2"] - r["y","x1"]*r["y","x2"])/
    sqrt(1 - r["y","x1"]^2)/sqrt(1 - r["y","x2"]^2)
  # ncp for testing partial correlation between y and x1
  t1 <- sqrt(nu*pcorr1^2/(1 - pcorr1^2))
  # ncp for testing partial correlation between y and x1
  t2 <- sqrt(nu*pcorr2^2/(1 - pcorr2^2))
  
  ### get elements of sigma
  # residual variance
  err <- 1 - (r["y","x1"]^2 + r["y","x2"]^2 - 2*r["y","x1"]*r["y","x2"]*r["x1","x2"])/
    (1 - r["t1","t2"]^2*alpha1*alpha2)
  # variances of b1 and b2
  vb <- err/(n-1)/(1 - r["x1","x2"]^2)
  # cov(b1,b2)
  covb <- err*r["x1","x2"]/(n-1)/(r["x1","x2"]^2 - 1)
  
  # return probability of observing significant x1 slope
  # p1: probability that x1 significant & x2 n.s.
  # p2: probability that x1 significant & x2 significant
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
           r_x1x2=r["x1","x2"], pcorr3=pcorr3, r2true=r2true, err=err))
}