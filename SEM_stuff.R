library("OpenMx") # for SEM model fitting
library("umx") # for umxExpCov
library("semPlot") # for path diagrams
source("/Users/Jake/Desktop/Google Drive/SEM_class/rewriteSemplotFunctions.R")


set.seed(56246)
n <- 120
tempF <- rnorm(n, mean=80, sd=10)
heat <- sqrt(.5)*scale(tempF) + sqrt(.5)*rnorm(n)
heat <- round(6*(heat - min(heat))/diff(range(heat))+1)
deaths <- rpois(n, lambda=exp(1 + .05*(tempF - mean(tempF))))
sales <- round(1000 + 20*(tempF - mean(tempF)) + rnorm(n, sd=200))
summary(lm(deaths ~ heat + sales))
# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -1.4569014  0.6608626  -2.205 0.029442 *  
# heat         0.5335614  0.1660621   3.213 0.001697 ** 
# sales        0.0024622  0.0006507   3.784 0.000245 ***
summary(lm(deaths ~ tempF + sales))


# reproduce multiple regression -------------------------------------------

mod0 <- mxTryHard(mxModel(name="mod0",
  type="RAM",
  manifestVars=c("deaths","heat","sales"),
  latentVars=c("L_heat","L_sales"),
  # indicators
  mxPath(from=c("L_heat","L_sales"), to=c("heat","sales"),
    connect="single", free=FALSE, value=1),
  # residual variances
  mxPath(from=c("deaths","heat","sales"), arrows=2, free=c(T,F,F),
         value=c(1, 0, 0)),
  # fix scale of latents
  mxPath(from=c("L_heat","L_sales"), arrows=2, free=FALSE, value=1),
  # allow latents to covary
  mxPath(from="L_heat", to="L_sales", arrows=2),
  # regress Y on Ts
  mxPath(from=c("L_heat","L_sales"), to=c("deaths")),
  # data source
  mxData(observed=cov(data.frame(deaths,heat,sales)), type="cov", numObs=100)))
summary(mod0)

mod0b <- mxTryHard(mxModel(name="mod0b",
  type="RAM",
  manifestVars=c("deaths","heat","sales"),
  latentVars=c("L_heat","L_sales"),
  # indicators
  mxPath(from=c("L_heat","L_sales"), to=c("heat","sales"),
         connect="single", free=FALSE, value=1),
  # residual variances
  mxPath(from=c("deaths","heat","sales"), arrows=2, free=c(T,F,F),
         value=c(1, 0, 0)),
  # fix scale of latents
  mxPath(from=c("L_heat","L_sales"), arrows=2, free=FALSE, value=1),
  # allow latents to covary
  mxPath(from="L_heat", to="L_sales", arrows=2),
  # regress Y on Ts
  mxPath(from=c("L_heat","L_sales"), to=c("deaths"), free=c(T,F), value=c(1,0)),
  # data source
  mxData(observed=cov(data.frame(deaths,heat,sales)), type="cov", numObs=100)))
umxCompare(mod0, mod0b)
#   Model EP &Delta; -2LL &Delta; df       p     AIC Compare with Model
# 1  mod0  4                                 9637618                   
# 2 mod0b  3    11.430585          1 < 0.001 9637627               mod0
sqrt(11.430585) # 3.380915

# add measurement error ---------------------------------------------------

sales2 <- sales/200
mod1 <- mxTryHard(mxModel(name="mod1",
  type="RAM",
  manifestVars=c("deaths","heat","sales2"),
  latentVars=c("L_heat","L_sales"),
  # indicators
  mxPath(from=c("L_heat","L_sales"), to=c("heat","sales2"),
         connect="single", free=FALSE, value=1),
  # residual variances
  mxPath(from=c("deaths","heat","sales2"), arrows=2, free=c(T,F,F),
         value=c(1, (1-.4)*var(heat), (1-.4)*var(sales2))),
  # fix scale of latents
  mxPath(from=c("L_heat","L_sales"), arrows=2, free=FALSE, value=1),
  # allow latents to covary
  mxPath(from="L_heat", to="L_sales", arrows=2),
  # regress Y on Ts
  mxPath(from=c("L_heat","L_sales"), to=c("deaths")),
  # data source
  mxData(observed=cov(data.frame(deaths,heat,sales2)), type="cov", numObs=100)))
summary(mod1)

mod1b <- mxTryHard(mxModel(name="mod1b",
  type="RAM",
  manifestVars=c("deaths","heat","sales2"),
  latentVars=c("L_heat","L_sales"),
  # indicators
  mxPath(from=c("L_heat","L_sales"), to=c("heat","sales2"),
        connect="single", free=FALSE, value=1),
  # residual variances
  mxPath(from=c("deaths","heat","sales2"), arrows=2, free=c(T,F,F),
         value=c(1, (1-.4)*var(heat), (1-1)*var(sales2))),
  # fix scale of latents
  mxPath(from=c("L_heat","L_sales"), arrows=2, free=FALSE, value=1),
  # allow latents to covary
  mxPath(from="L_heat", to="L_sales", arrows=2),
  # regress Y on Ts
  mxPath(from=c("L_heat","L_sales"), to=c("deaths"), free=c(T,F), value=c(1,0)),
  # data source
  mxData(observed=cov(data.frame(deaths,heat,sales2)), type="cov", numObs=100)))
umxCompare(mod1, mod1b)
#   Model EP &Delta; -2LL &Delta; df     p       AIC Compare with Model
# 1  mod1  4                               -1.198244                   
# 2 mod1b  3    0.0015238          1 0.969 -3.196720               mod1
sqrt(0.5492853) # 0.7411378

# plot model --------------------------------------------------------------

tiff(paste0(path,"Fig9.tif"), units="in",
    height=5, width=6, res=200, pointsize=18)
semPaths(semPlotModel_MxRAMModel(mod1), whatLabels=c("name"),
         sizeMan=10, sizeLat=10, mar=c(10,6,10,6), layout="tree2",
         edge.color="black", edge.label.cex=1.5, intercepts=FALSE,
         color=list(man="pink", lat="lightblue"), style="ram",
         nodeLabels=c("NEO\nscore","Drug\nuse","HEXACO\nscore","NEO\nconstruct","HEXACO\nconstruct"),
         edgeLabels=list(expression(beta[1]), expression(beta[2]),"1.0","1.0",
           expression(sigma^2), expression((1-alpha[1])*s[NEO]^2),
           expression((1-alpha[2])*s[HEXACO]^2), "1.0",
           expression(delta), "1.0"))
dev.off()


# partial correlation model -----------------------------------------------

sales2 <- sales/200
mod2 <- mxTryHard(mxModel(name="mod2",
  type="RAM",
  manifestVars=c("deaths","heat","sales2"),
  latentVars=c("L_heat","L_sales"),
  # indicators
  mxPath(from=c("L_heat","L_sales"), to=c("heat","sales2"),
         connect="single", free=FALSE, value=1),
  # residual variances
  mxPath(from=c("deaths","heat","sales2"), arrows=2, free=c(T,F,F),
         value=c(1, .6*var(heat), .1*var(sales2))),
  # fix scale of latents
  mxPath(from=c("L_heat","L_sales"), arrows=2, free=c(F,T), value=c(1,1)),
  # regress sales on heat
  mxPath(from="L_heat", to="L_sales"),
  # regress Y on Ts
  mxPath(from=c("L_heat","L_sales"), to=c("deaths")),
  # data source
  mxData(observed=cov(data.frame(deaths,heat,sales2)), type="cov", numObs=100)))
summary(mod2)

png("/Users/Jake/Desktop/SEM.png",
    type="cairo", units="in", height=5, width=6, res=200, pointsize=18)
semPaths(semPlotModel_MxRAMModel(mod2), whatLabels=c("std"),
         sizeMan=10, sizeLat=10, mar=c(10,6,10,6), layout="tree",
         edge.color="black", edge.label.cex=1.5, intercepts=FALSE,
         color=list(man="pink", lat="lightblue"), style="lisrel")
dev.off()

# z as f(reliability) -----------------------------------------------------

getZ <- function(alpha){
  mod <- mxTryHard(mxModel(name="mod",
     type="RAM",
     manifestVars=c("deaths","heat","sales2"),
     latentVars=c("L_heat","L_sales"),
     mxPath(from=c("L_heat","L_sales"), to=c("heat","sales2"),
            connect="single", free=FALSE, value=1),
     mxPath(from=c("deaths","heat","sales2"), arrows=2, free=c(T,F,F),
            value=c(1, (1-alpha)*var(heat), (1-1)*var(sales2))),
     mxPath(from=c("L_heat","L_sales"), arrows=2, free=FALSE, value=1),
     mxPath(from="L_heat", to="L_sales", arrows=2),
     mxPath(from=c("L_heat","L_sales"), to=c("deaths"),
            labels=c("heatCoef","salesCoef")),
     mxData(observed=cov(data.frame(deaths,heat,sales2)),
            type="cov", numObs=length(deaths))))
  return(mod$output$estimate["salesCoef"]/
           mod$output$standardErrors["salesCoef",])
}

z <- sapply(seq(0,1,.1), getZ)

png("/Users/Jake/Desktop/z_as_fAlpha.png", type="cairo",
    units="in", height=8, width=8, res=200, pointsize=25)
par(mar=c(4, 4, 2.5, 1)+.1)
plot(y=z, x=seq(0,1,.1), type="o", pch=20, ylab="Wald z-statistic",
     main=expression(paste("Partial effect of ice cream sales ",(beta[2]))),
     xlab=expression(paste("Assumed reliability of heat measure ",(alpha[1]))))
abline(h=qnorm(.975), lty=2)
dev.off()

# chi-square as f(alpha) --------------------------------------------------

set.seed(567453)
tempEst <- round(tempF + rnorm(n, sd=15))

summary(lm(deaths ~ heat + tempEst))

tempEst2 <- tempEst/20

getX2 <- function(alpha){
  modA <- mxRun(mxModel(name="modA",
     type="RAM",
     manifestVars=c("deaths","heat","tempEst2"),
     latentVars=c("L_heat","L_temp"),
     mxPath(from=c("L_heat","L_temp"), to=c("heat","tempEst2"),
            connect="single", free=FALSE, value=1),
     mxPath(from=c("deaths","heat","tempEst2"), arrows=2, free=c(T,F,F),
            value=c(1, (1-alpha)*var(heat), (1-alpha)*var(tempEst2))),
     mxPath(from=c("L_heat","L_temp"), arrows=2, free=FALSE, value=1),
     mxPath(from="L_heat", to="L_temp", arrows=2),
     mxPath(from=c("L_heat","L_temp"), to=c("deaths"),
            labels=c("heatCoef","tempEstCoef")),
     mxData(observed=cov(data.frame(deaths,heat,tempEst2)),
            type="cov", numObs=length(deaths))))
  modC <- mxRun(mxModel(name="modC",
    type="RAM",
    manifestVars=c("deaths","heat","tempEst2"),
    latentVars=c("L"),
    mxPath(from=c("L"), to=c("heat","tempEst2"), free=FALSE, value=1),
    mxPath(from=c("deaths","heat","tempEst2"), arrows=2, free=c(T,F,F),
           value=c(1, (1-alpha)*var(heat), (1-alpha)*var(tempEst2))),
    mxPath(from=c("L"), arrows=2, free=FALSE, value=1),
    mxPath(from=c("L"), to=c("deaths")),
    mxData(observed=cov(data.frame(deaths,heat,tempEst2)),
           type="cov", numObs=length(deaths))))
  umxCompare(modA,modC)[2,3]
}

LR <- sapply(c(seq(0,.9,.1), .99), getX2)
LR <- sapply(seq(0,.99,.01), getX2)

png("/Users/Jake/Desktop/LR_as_fAlpha.png", type="cairo",
    units="in", height=8, width=8, res=200, pointsize=25)
par(mar=c(4, 4, 3, 1.5)+.1)
plot(y=log(LR), x=seq(0,.99,.01), type="l", pch=20, yaxt="n", lwd=2,
     ylab=expression(paste(chi^2," difference (log scale)")),
     main="", mgp=c(2.5,1,0),
     xlab=expression(paste("Assumed reliabilities ",(alpha[1]==alpha[2]))))
title("Subjective heat and estimated temperature:", line=2, cex.main=.9)
title(expression(paste(chi^2," difference test of ",delta==1," and ",beta[1]==beta[2])),
      line=1, cex.main=1.1)
axis(side=2, at=log(c(1,7,50,400,3000)), labels=c(1,7,50,400,3000))
abline(h=log(qchisq(.975, df=2)), lty=2)
dev.off()