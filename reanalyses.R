library("Hmisc") # for spss.get()
# library("OpenMx") # for SEM model fitting
# library("umx") # for umxExpCov and other cool stuff
library("semTools") # for reliability()
library("semPlot") # for path diagrams
library("lavaan") # for SEMs

# run this to make semPlot play nice with OpenMx
source("/Users/Jake/Desktop/Google Drive/SEM_class/rewriteSemplotFunctions.R")

# read in panel dataset ---------------------------------------------------

path <- "/Users/Jake/Desktop/Dropbox/Dropbox/SIC/"
panel <- read.csv(paste0(path,"anes2008_2009panel/anes2008_2009panel_rawdata.txt"))

# examine vars ------------------------------------------------------------

table(panel$derw10a1) # for whom will R vote?
# 1 = mccain, 2 = obama, 3 = someone else
# -6  -5  -1   1   2   3 
# 80   3  38 518 536  58 

table(panel$der16) # For whom did R vote for president in 2008 election
# 13 = mccain, 18 = obama, 3 = someone else
# -6  -2  -1   1   2   3 
# 27  28 113 497 542  26 

# good agreement between prospective and retrospective
with(panel, table(der16, derw10a1))
#      derw10a1
# der16  -6  -5  -1   1   2   3
#    -6  27   0   0   0   0   0
#    -2   6   0   0  12   6   4
#    -1   4   0   0  57  43   9
#    1   22   1  19 425  13  17
#    2   21   2  18  16 470  15
#    3    0   0   1   8   4  13

# grab names of AMP response variables
n <- nchar(names(panel))
name <- names(panel)[substr(names(panel), n-6, n) == "_choice"]

# verify that respondents complete only 1 of the AMPs, w9 or w10
apply(panel[,name], 1, function(x) sum(!x %in% 1:2))

### explicit attitude measures

# wave 2 warm/cold
table(panel$w2d11) # do you feel warm, cold, or neither toward Blacks?
# 1 = warm, 2 = cold, 3 = neither
table(panel$w2d12) # How warm does R feel to blacks
# 1 = extremely warm, 2 = moderately warm, 3 = a little warm
table(panel$w2d13) # How cold does R feel to blacks
# 1 = extremely cold, 2 = moderately cold, 3 = a little cold

# wave 10 warm/cold
table(panel$w10d11) # do you feel warm, cold, or neither toward Blacks?
table(panel$w10d11_warm) # How warm does R feel to blacks
table(panel$w10d11_cold) # How cold does R feel to blacks

# wave 20 warm/cold
table(panel$w20d11) # do you feel warm, cold, or neither toward Blacks?
table(panel$w20d12) # How warm does R feel to blacks
table(panel$w20d13) # How cold does R feel to blacks

# sympathy for blacks
# 1 = always, 2 = most time, 3 = half time, 4 = once in a while, 5 = never
table(panel$w9zb24)
table(panel$w11zb24)
table(panel$w17x24)

# admiration for blacks
# 1 = always, 2 = most time, 3 = half time, 4 = once in a while, 5 = never
table(panel$w9zb25)
table(panel$w11zb25)
table(panel$w17x25)

# Do blacks have too much or too little pol influence
# 1 = too much, 2 = just about right, 3 = too little
table(panel$w9zb23) # factor, not numeric!
table(panel$w11zb23)
table(panel$w17x26)

### control variables:

# party identification
table(panel$der08w1)
table(panel$der08w9)
table(panel$der08w10)
table(panel$der08w11)
table(panel$der08w17)
table(panel$der08w19)

# liberalism/conservatism
table(panel$der09w1)
table(panel$der09w2)
table(panel$der09w6)
table(panel$der09w10)
table(panel$der09w11)

# gender
table(panel$rgenderr)
# age
hist(panel$cpyourself_age, col="gray")
# race
table(panel$der04)
# education
table(panel$der05)
# income
table(panel$der06)

# build dataset -----------------------------------------------------------

# grab names of AMP response variables
n <- nchar(names(panel))
name <- names(panel)[substr(names(panel), n-6, n) == "_choice"]

# construct dataset with only variables relevant to AMP analysis
clean <- panel[,c(
  "caseid","der16", # president vote
  "w2d11","w2d12","w2d13", # wave 2 warm/cold
  "w10d11","w10d11_warm","w10d11_cold", # wave 10 warm/cold
  "w20d11","w20d12","w20d13", # wave 20 warm/cold
  "w9zb24","w11zb24","w17x24", # sympathy
  "w9zb25","w11zb25","w17x25", # admiration
  "w9zb23","w11zb23","w17x26", # influence
  # party identification
  "der08w1","der08w9","der08w10","der08w11","der08w17","der08w19",
  # liberal/conservative
  "der09w1","der09w2","der09w6","der09w10","der09w11",
  "rgenderr","cpyourself_age","der04","der05","der06")] # demographics
names(clean) <- c(
  "caseid","vote",
  "warmCold_w2","howWarm_w2","howCold_w2",
  "warmCold_w10","howWarm_w10","howCold_w10",
  "warmCold_w20","howWarm_w20","howCold_w20",
  "sympathy_w9","sympathy_w11","sympathy_w17",
  "admiration_w9","admiration_w11","admiration_w17",
  "influence_w9","influence_w11","influence_w17",
  "party_w1","party_w9","party_w10","party_w11","party_w17","party_w19",
  "libCon_w1","libCon_w2","libCon_w6","libCon_w10","libCon_w11",
  "gender","age","race","education","income")
# not sure why this one is a factor, but cast to numeric
clean$influence_w9 <- as.numeric(as.character(clean$influence_w9))

# replace missing codes (negative values) with NA
clean[clean < 0] <- NA
# and remove 2.5% of Rs who voted for other than mccain/obama
clean$vote[clean$vote==3] <- NA

# grab the AMP responses
choices <- panel[,name]
choices <- data.frame(lapply(choices, function(x) as.numeric(as.character(x))))
choices[choices < 1] <- NA

# compute new variables
key <- c("1"=1,"2"=-1,"3"=0)
clean <- within(clean, {
  vote <- vote - 1
  black <- rowMeans(choices[c(1:24, 49:72)]==1, na.rm=TRUE)
  white <- rowMeans(choices[c(25:48, 73:96)]==1, na.rm=TRUE)
  diff <- black - white
  sympathy <- rowMeans(cbind(sympathy_w9, sympathy_w11, sympathy_w17), na.rm=TRUE)
  sympathy <- 6 - sympathy
  admiration <- rowMeans(cbind(admiration_w9, admiration_w11, admiration_w17), na.rm=TRUE)
  admiration <- 6 - admiration
  influence <- rowMeans(cbind(influence_w9, influence_w11, influence_w17), na.rm=TRUE)
  warm_w2 <- key[as.character(warmCold_w2)]
  warm_w2 <- warm_w2*cbind(howCold_w2,0,howWarm_w2)[cbind(seq(nrow(clean)),warm_w2+2)]
  warm_w10 <- key[as.character(warmCold_w10)]
  warm_w10 <- warm_w10*cbind(howCold_w10,0,howWarm_w10)[cbind(seq(nrow(clean)),warm_w10+2)]
  warm_w20 <- key[as.character(warmCold_w20)]
  warm_w20 <- warm_w20*cbind(howCold_w20,0,howWarm_w20)[cbind(seq(nrow(clean)),warm_w20+2)]
  warm <- rowMeans(cbind(warm_w2, warm_w10, warm_w20), na.rm=TRUE)
  explicit <- rowMeans(scale(cbind(warm,sympathy,admiration,influence)),
                       na.rm=TRUE)
  party <- rowMeans(cbind(party_w1,party_w9,party_w10,party_w11,party_w17,party_w19), na.rm=TRUE)
  libCon <- rowMeans(cbind(libCon_w1,libCon_w2,libCon_w6,libCon_w10,libCon_w11), na.rm=TRUE)
})

# merge in pre-computed D-scores
iat <- spss.get(paste0(path,"IATscores/ANES0809Panel_IAT.por"))
names(iat)[match("CASEID", names(iat))] <- "caseid"
clean <- merge(clean, iat[,c("caseid","IAT.D")], all.x=TRUE)

# mean center most predictors
index <- match(c("vote","diff","race","IAT.D"), names(clean))
clean[,-index] <- scale(clean[,-index], scale=FALSE)

# dummy code race
dummy <- contr.treatment(4)
colnames(dummy) <- c("race1","race2","race3")
clean <- data.frame(clean, dummy[clean$race,])

# check number of missing obs by variable
apply(clean, 2, function(x) sum(is.na(x)))

# check simple correlations
round(cor(clean[,c("diff","IAT.D","explicit","party","libCon","gender",
  "age","education","income")], use="pairwise.complete.obs"), 2)

# AMP regression models ---------------------------------------------------

# note two simplifications compared to Payne et al.:
# single implicit predictor (% diff) rather than black % + white %
# single vote outcome (mccain vs. obama) rather than 2 vote dummies

# model using only AMP scores as predictors
summary(glm(vote ~ diff, family=binomial, data=clean))

# control for explicit measures
summary(glm(vote ~ diff + explicit, family=binomial, data=clean))

# w/ demographics, without controlling for explicit
summary(glm(vote ~ diff +
              party + libCon + gender + age + race + education + income,
            family=binomial, data=clean))

# w/ demographics, controlling for explicit
summary(glm(vote ~ diff + explicit +
              party + libCon + gender + age + race + education + income,
            family=binomial, data=clean))

# AMP SEM -----------------------------------------------------------------


ind <- c("warm","sympathy","admiration","influence")
demo <- c("party","libCon","gender","age","education","income",
          "race1","race2","race3")
# the model has a hard time with missing data!
include <- !is.na(clean$diff) & !is.na(clean$vote)
ampdat <- within(na.omit(clean[include, c("vote","diff",ind,demo)]), {
  # scale variables to have similar SDs
  diff <- diff*3
  influence <- influence*2
  party <- party/2
  libCon <- libCon/2
  age <- age/10
  income <- income/4
  race1 <- race1*3
  race2 <- race2*3
  race3 <- race3*3
  vote <- mxFactor(vote, levels=0:1)
})

# subsample
ampdat <- ampdat[sample(seq(nrow(ampdat)), 400),]

# check SDs
sqrt(diag(var(ampdat, na.rm=TRUE)))
# check number of missing obs by variable (should be none)
apply(ampdat, 2, function(x) sum(is.na(x)))

# test out with reliability = .4 (no demographic controls)
alpha <- .55
modA <- mxTryHard(mxOption(model=mxModel(name="modA",
  type="RAM",
  manifestVars=c("vote","diff",ind),
  latentVars=c("implicit","explicit"),
  # indicators
  mxPath(from="explicit", to=ind),
  mxPath(from="implicit", to="diff", free=FALSE, value=1),
  # residual variances
  mxPath(from=ind, arrows=2),
  mxPath(from="diff", arrows=2, free=FALSE,
         value=(1-alpha)*var(ampdat$diff, na.rm=TRUE)),
  mxPath(from="vote", arrows=2, free=FALSE, values=1),
  # fix scale of latents
  mxPath(from=c("implicit","explicit"), arrows=2, free=FALSE, value=1),
  # allow latents to covary
  mxPath(from="implicit", to="explicit", arrows=2),
  # regress Y on Ts
  mxPath(from=c("implicit","explicit"), to="vote"),
  # means/intercepts and thresholds
  mxPath(from="one", to=c("implicit","explicit","vote"), free=FALSE, values=0),
  mxPath(from="one", to=c("diff",ind)),
  mxThreshold(vars="vote", nThresh=1, values=0),
  # data source
  mxData(observed=ampdat, type="raw")),
  # enable parallel (arguments to mxOption)
  key="Number of Threads", value=(omxDetectCores()-1)),
  # number of times to try new starting values (argument to mxTryHard)
  extraTries=99)

summary(modA)
umxSummary(modA)

modC <- mxTryHard(mxOption(model=mxModel(name="modC",
  type="RAM",
  manifestVars=c("vote","diff",ind),
  latentVars=c("attitude"),
  # indicators
  mxPath(from="attitude", to=ind),
  mxPath(from="attitude", to="diff", free=FALSE, value=1),
  # residual variances
  mxPath(from=ind, arrows=2),
  mxPath(from="diff", arrows=2, free=FALSE,
  value=(1-alpha)*var(ampdat$diff, na.rm=TRUE)),
  mxPath(from="vote", arrows=2, free=FALSE, values=1),
  # fix scale of latent
  mxPath(from=c("attitude"), arrows=2, free=FALSE, value=1),
  # regress Y on latent
  mxPath(from=c("attitude"), to="vote"),
  # means/intercepts and thresholds
  mxPath(from="one", to=c("attitude","vote"), free=FALSE, values=0),
  mxPath(from="one", to=c("diff",ind)),
  mxThreshold(vars="vote", nThresh=1, values=0),
  # data source
  mxData(observed=ampdat, type="raw")),
  # enable parallel (arguments to mxOption)
  key="Number of Threads", value=(omxDetectCores()-1)),
  # number of times to try new starting values (argument to mxTryHard)
  extraTries = 99)

summary(modC)

umxCompare(modA,modC)

loads <- modA@matrices$A$values[ind,"explicit"]
errors <- diag(modA@matrices$S$values[ind,ind])
round(sum(loads)^2/(sum(loads)^2 + sum(errors)), 2)
# reliability = 0.86

# define function to get chi-square differences as f(alpha)
getX2 <- function(alpha){
  modA <- mxTryHard(mxOption(model=mxModel(name="modA",
    type="RAM",
    manifestVars=c("vote","diff",ind),
    latentVars=c("implicit","explicit"),
    mxPath(from="explicit", to=ind),
    mxPath(from="implicit", to="diff", free=FALSE, value=1),
    mxPath(from=ind, arrows=2),
    mxPath(from="diff", arrows=2, free=FALSE,
          value=(1-alpha)*var(ampdat$diff, na.rm=TRUE)),
    mxPath(from="vote", arrows=2, free=FALSE, values=1),
    mxPath(from=c("implicit","explicit"), arrows=2, free=FALSE, value=1),
    mxPath(from="implicit", to="explicit", arrows=2),
    mxPath(from=c("implicit","explicit"), to="vote"),
    mxPath(from="one", to=c("implicit","explicit","vote"), free=FALSE, values=0),
    mxPath(from="one", to=c("diff",ind)),
    mxThreshold(vars="vote", nThresh=1, values=0),
    mxData(observed=ampdat, type="raw")),
    key="Number of Threads", value=(omxDetectCores()-1)),
    extraTries=99)
  modC <- mxTryHard(mxOption(model=mxModel(name="modC",
    type="RAM",
    manifestVars=c("vote","diff",ind),
    latentVars=c("attitude"),
    mxPath(from="attitude", to=ind),
    mxPath(from="attitude", to="diff", free=FALSE, value=1),
    mxPath(from=ind, arrows=2),
    mxPath(from="diff", arrows=2, free=FALSE,
      value=(1-alpha)*var(ampdat$diff, na.rm=TRUE)),
    mxPath(from="vote", arrows=2, free=FALSE, values=1),
    mxPath(from=c("attitude"), arrows=2, free=FALSE, value=1),
    mxPath(from=c("attitude"), to="vote"),
    mxPath(from="one", to=c("attitude","vote"), free=FALSE, values=0),
    mxPath(from="one", to=c("diff",ind)),
    mxThreshold(vars="vote", nThresh=1, values=0),
    mxData(observed=ampdat, type="raw")),
    key="Number of Threads", value=(omxDetectCores()-1)),
    extraTries = 99)
  umxCompare(modA,modC)[2,3]
}

system.time({
  X2diffs <- sapply(seq(from=.01, to=.99, length.out=11), getX2)
}) # elapsed = 1.5 minutes

png("/Users/Jake/Desktop/AMP_LR_as_fAlpha_subsample.png",
    units="in", height=8, width=8, res=200, pointsize=25)
plot(y=log(X2diffs), x=seq(from=.01, to=.99, length.out=11),
     type="o", pch=19, yaxt="n", lwd=2, ylim=c(0,max(log(X2diffs))),
     ylab=expression(paste(chi^2," difference (log scale)")),
     main="AMP", mgp=c(2.5,1,0),
     xlab="Assumed reliability of implicit measure")
axis(side=2, at=log(c(1,7,50,400,3000)), labels=c(1,7,50,400,3000))
abline(h=c(0, log(qchisq(.975, df=2))), lty=2)
dev.off()

### this time with demographic controls

alpha <- .4
modA_demo <- mxTryHard(mxOption(model=mxModel(name="modA_demo",
  type="RAM",
  manifestVars=c("vote","diff",ind,demo),
  latentVars=c("implicit","explicit"),
  # indicators
  mxPath(from="explicit", to=ind),
  mxPath(from="implicit", to="diff", free=FALSE, value=1),
  # residual variances
  mxPath(from=ind, arrows=2),
  mxPath(from="diff", arrows=2, free=FALSE,
    value=(1-alpha)*var(ampdat$diff, na.rm=TRUE)),
  mxPath(from="vote", arrows=2, free=FALSE, values=1),
  # fix scale of latents
  mxPath(from=c("implicit","explicit"), arrows=2, free=FALSE, value=1),
  # allow latents to covary
  #mxPath(from="implicit", to="explicit", arrows=2),
  # regress Y on Ts
  mxPath(from=c("implicit","explicit"), to="vote"),
  # means/intercepts and thresholds
  mxPath(from="one", to=c("implicit","explicit","vote"), free=FALSE, values=0),
  mxPath(from="one", to=c("diff",ind)),
  mxThreshold(vars="vote", nThresh=1, values=0),
  # demographic controls
  mxPath(from=demo, to="vote"),
  mxPath(from=demo, arrows=2),
  mxPath(from=c(demo,"implicit","explicit"), arrows=2, connect="unique.bivariate"),
  # data source
  mxData(observed=ampdat, type="raw")),
  # enable parallel (arguments to mxOption)
  key="Number of Threads", value=(omxDetectCores()-1)),
  # number of times to try new starting values (argument to mxTryHard)
  extraTries=99)

# define function to get chi-square differences as f(alpha)
getX2_demo <- function(alpha){
  modA <- mxTryHard(mxOption(model=mxModel(name="modA",
    type="RAM",
    manifestVars=c("vote","diff",ind,demo),
    latentVars=c("implicit","explicit"),
    mxPath(from="explicit", to=ind),
    mxPath(from="implicit", to="diff", free=FALSE, value=1),
    mxPath(from=ind, arrows=2),
    mxPath(from="diff", arrows=2, free=FALSE,
    value=(1-alpha)*var(ampdat$diff, na.rm=TRUE)),
    mxPath(from="vote", arrows=2, free=FALSE, values=1),
    mxPath(from=c("implicit","explicit"), arrows=2, free=FALSE, value=1),
    mxPath(from="implicit", to="explicit", arrows=2),
    mxPath(from=c("implicit","explicit"), to="vote"),
    mxPath(from="one", to=c("implicit","explicit","vote"), free=FALSE, values=0),
    mxPath(from="one", to=c("diff",ind)),
    mxThreshold(vars="vote", nThresh=1, values=0),
    mxPath(from=demo, to="vote"),
    mxPath(from=demo, arrows=2),
    mxData(observed=ampdat, type="raw")),
    key="Number of Threads", value=(omxDetectCores()-1)),
    extraTries=99)
  modC <- mxTryHard(mxOption(model=mxModel(name="modC",
    type="RAM",
    manifestVars=c("vote","diff",ind,demo),
    latentVars=c("attitude"),
    mxPath(from="attitude", to=ind),
    mxPath(from="attitude", to="diff", free=FALSE, value=1),
    mxPath(from=ind, arrows=2),
    mxPath(from="diff", arrows=2, free=FALSE,
      value=(1-alpha)*var(ampdat$diff, na.rm=TRUE)),
    mxPath(from="vote", arrows=2, free=FALSE, values=1),
    mxPath(from=c("attitude"), arrows=2, free=FALSE, value=1),
    mxPath(from=c("attitude"), to="vote"),
    mxPath(from="one", to=c("attitude","vote"), free=FALSE, values=0),
    mxPath(from="one", to=c("diff",ind)),
    mxThreshold(vars="vote", nThresh=1, values=0),
    mxPath(from=demo, to="vote"),
    mxPath(from=demo, arrows=2),
    mxData(observed=ampdat, type="raw")),
    key="Number of Threads", value=(omxDetectCores()-1)),
    extraTries = 99)
  umxCompare(modA,modC)[2,3]
}

system.time({
  X2diffs_demo <- sapply(seq(from=.01, to=.99, length.out=11), getX2_demo)
}) # elapsed = 1.5 minutes

# path diagrams -----------------------------------------------------------

# for plotting purposes, must change ordinal outcome to continuous

modA2 <- mxTryHard(mxOption(model=mxModel(name="modA2",
  type="RAM",
  manifestVars=c("vote2","diff",ind),
  latentVars=c("implicit","explicit"),
  # indicators
  mxPath(from="explicit", to=ind),
  mxPath(from="implicit", to="diff", free=FALSE, value=1),
  # residual variances
  mxPath(from=ind, arrows=2),
  mxPath(from="diff", arrows=2, free=FALSE,
         value=(1-alpha)*var(ampdat$diff, na.rm=TRUE)),
  mxPath(from="vote2", arrows=2),
  # fix scale of latents
  mxPath(from=c("implicit","explicit"), arrows=2, free=FALSE, value=1),
  # allow latents to covary
  mxPath(from="implicit", to="explicit", arrows=2),
  # regress Y on Ts
  mxPath(from=c("implicit","explicit"), to="vote2"),
  # means/intercepts and thresholds
  mxPath(from="one", to=c("implicit","explicit"), free=FALSE, values=0),
  mxPath(from="one", to=c("diff",ind,"vote2")),
  # data source
  mxData(observed=ampdat[,-match("vote",names(ampdat))], type="raw")),
  # enable parallel (arguments to mxOption)
  key="Number of Threads", value=(omxDetectCores()-1)))

modC2 <- mxTryHard(mxOption(model=mxModel(name="modC2",
  type="RAM",
  manifestVars=c("vote2","diff",ind),
  latentVars=c("attitude"),
  # indicators
  mxPath(from="attitude", to=ind),
  mxPath(from="attitude", to="diff", free=FALSE, value=1),
  # residual variances
  mxPath(from=ind, arrows=2),
  mxPath(from="diff", arrows=2, free=FALSE,
         value=(1-alpha)*var(ampdat$diff, na.rm=TRUE)),
  mxPath(from="vote2", arrows=2),
  # fix scale of latent
  mxPath(from=c("attitude"), arrows=2, free=FALSE, value=1),
  # regress Y on latent
  mxPath(from=c("attitude"), to="vote2"),
  # means/intercepts and thresholds
  mxPath(from="one", to=c("attitude"), free=FALSE, values=0),
  mxPath(from="one", to=c("diff",ind,"vote2")),
  # data source
  mxData(observed=ampdat[,-match("vote",names(ampdat))], type="raw")),
  # enable parallel (arguments to mxOption)
  key="Number of Threads", value=(omxDetectCores()-1)))

# plot path diagram
semPaths(semPlotModel_MxRAMModel(modA2), whatLabels=c("par"),
         sizeMan=6, sizeLat=9, intercepts=FALSE, style="lisrel",
         edge.color="black", edge.label.cex=1)

semPaths(semPlotModel_MxRAMModel(modC2), whatLabels=c("par"),
         sizeMan=6, sizeLat=9, intercepts=FALSE, style="lisrel",
         edge.color="black", edge.label.cex=1)

# IAT regression models ---------------------------------------------------

# model using only D scores as predictors
summary(glm(vote ~ IAT.D, family=binomial, data=clean))

# control for explicit measures
summary(glm(vote ~ IAT.D + explicit, family=binomial, data=clean))

# w/ demographics, without controlling for explicit
summary(glm(vote ~ IAT.D +
              party + libCon + gender + age + race + education + income,
            family=binomial, data=clean))

# w/ demographics, controlling for explicit
summary(glm(vote ~ IAT.D + explicit +
              party + libCon + gender + age + race + education + income,
            family=binomial, data=clean))

# IAT SEM -----------------------------------------------------------------

ind <- c("warm","sympathy","admiration","influence")
# the model has a hard time with missing data!
include <- !is.na(clean$IAT.D) & !is.na(clean$vote)
iatdat <- within(na.omit(clean[include, c("vote","IAT.D",ind)]), {
  # scale variables to have similar SDs
  D <- IAT.D
  IAT.D <- NULL
  influence <- influence
  vote <- mxFactor(vote, levels=0:1)
})
# check SDs
sqrt(diag(var(iatdat, na.rm=TRUE)))
# check number of missing obs by variable (should be none)
apply(iatdat, 2, function(x) sum(is.na(x)))

# test out with reliability = .5
alpha <- .4
modA <- mxTryHard(mxOption(model=mxModel(name="modA",
  type="RAM",
  manifestVars=c("vote","D",ind),
  latentVars=c("implicit","explicit"),
  # indicators
  mxPath(from="explicit", to=ind),
  mxPath(from="implicit", to="D", free=FALSE, value=1),
  # residual variances
  mxPath(from=ind, arrows=2),
  mxPath(from="D", arrows=2, free=FALSE,
  value=(1-alpha)*var(iatdat$D, na.rm=TRUE)),
  mxPath(from="vote", arrows=2, free=FALSE, values=1),
  # fix scale of latents
  mxPath(from=c("implicit","explicit"), arrows=2, free=FALSE, value=1),
  # allow latents to covary
  mxPath(from="implicit", to="explicit", arrows=2),
  # regress Y on Ts
  mxPath(from=c("implicit","explicit"), to="vote"),
  # means/intercepts and thresholds
  mxPath(from="one", to=c("implicit","explicit","vote"), free=FALSE, values=0),
  mxPath(from="one", to=c("D",ind)),
  mxThreshold(vars="vote", nThresh=1, values=0),
  # data source
  mxData(observed=iatdat, type="raw")),
  # enable parallel (arguments to mxOption)
  key="Number of Threads", value=(omxDetectCores()-1)),
  # number of times to try new starting values (argument to mxTryHard)
  extraTries=99)

summary(modA)
umxSummary(modA)

modC <- mxTryHard(mxOption(model=mxModel(name="modC",
  type="RAM",
  manifestVars=c("vote","D",ind),
  latentVars=c("attitude"),
  # indicators
  mxPath(from="attitude", to=ind),
  mxPath(from="attitude", to="D", free=FALSE, value=1),
  # residual variances
  mxPath(from=ind, arrows=2),
  mxPath(from="D", arrows=2, free=FALSE,
    value=(1-alpha)*var(iatdat$D, na.rm=TRUE)),
  mxPath(from="vote", arrows=2, free=FALSE, values=1),
  # fix scale of latent
  mxPath(from=c("attitude"), arrows=2, free=FALSE, value=1),
  # regress Y on latent
  mxPath(from=c("attitude"), to="vote"),
  # means/intercepts and thresholds
  mxPath(from="one", to=c("attitude","vote"), free=FALSE, values=0),
  mxPath(from="one", to=c("D",ind)),
  mxThreshold(vars="vote", nThresh=1, values=0),
  # data source
  mxData(observed=iatdat, type="raw")),
  # enable parallel (arguments to mxOption)
  key="Number of Threads", value=(omxDetectCores()-1)),
  # number of times to try new starting values (argument to mxTryHard)
  extraTries = 99)

summary(modC)

umxCompare(modA,modC)

loads <- modA@matrices$A$values[ind,"explicit"]
errors <- diag(modA@matrices$S$values[ind,ind])
round(sum(loads)^2/(sum(loads)^2 + sum(errors)), 2)
# reliability = 0.77

# define function to get chi-square differences as f(alpha)
getX2 <- function(alpha){
  modA <- mxTryHard(mxOption(model=mxModel(name="modA",
    type="RAM",
    manifestVars=c("vote","D",ind),
    latentVars=c("implicit","explicit"),
    mxPath(from="explicit", to=ind),
    mxPath(from="implicit", to="D", free=FALSE, value=1),
    mxPath(from=ind, arrows=2),
    mxPath(from="D", arrows=2, free=FALSE,
      value=(1-alpha)*var(iatdat$D, na.rm=TRUE)),
    mxPath(from="vote", arrows=2, free=FALSE, values=1),
    mxPath(from=c("implicit","explicit"), arrows=2, free=FALSE, value=1),
    mxPath(from="implicit", to="explicit", arrows=2),
    mxPath(from=c("implicit","explicit"), to="vote"),
    mxPath(from="one", to=c("implicit","explicit","vote"), free=FALSE, values=0),
    mxPath(from="one", to=c("D",ind)),
    mxThreshold(vars="vote", nThresh=1, values=0),
    mxData(observed=iatdat, type="raw")),
    key="Number of Threads", value=(omxDetectCores()-1)),
    extraTries=99)
  modC <- mxTryHard(mxOption(model=mxModel(name="modC",
    type="RAM",
    manifestVars=c("vote","D",ind),
    latentVars=c("attitude"),
    mxPath(from="attitude", to=ind),
    mxPath(from="attitude", to="D", free=FALSE, value=1),
    mxPath(from=ind, arrows=2),
    mxPath(from="D", arrows=2, free=FALSE,
      value=(1-alpha)*var(iatdat$D, na.rm=TRUE)),
    mxPath(from="vote", arrows=2, free=FALSE, values=1),
    mxPath(from=c("attitude"), arrows=2, free=FALSE, value=1),
    mxPath(from=c("attitude"), to="vote"),
    mxPath(from="one", to=c("attitude","vote"), free=FALSE, values=0),
    mxPath(from="one", to=c("D",ind)),
    mxThreshold(vars="vote", nThresh=1, values=0),
    mxData(observed=iatdat, type="raw")),
    key="Number of Threads", value=(omxDetectCores()-1)),
    extraTries = 99)
  umxCompare(modA,modC)[2,3]
}

system.time({
  X2diffs2 <- sapply(seq(from=.01, to=.99, length.out=11), getX2)
}) # elapsed = 1 minute

png("/Users/Jake/Desktop/IAT_LR_as_fAlpha.png",
    units="in", height=8, width=8, res=200, pointsize=25)
plot(y=log(X2diffs2), x=seq(from=.01, to=.99, length.out=11),
     type="o", pch=19, yaxt="n", lwd=2, ylim=c(0,max(log(X2diffs2))),
     ylab=expression(paste(chi^2," difference (log scale)")),
     main="IAT", mgp=c(2.5,1,0),
     xlab="Assumed reliability of implicit measure")
axis(side=2, at=log(c(1,7,50,400,3000)), labels=c(1,7,50,400,3000))
abline(h=c(0, log(qchisq(.975, df=2))), lty=2)
dev.off()

# IPIP, all 11 factors ----------------------------------------------------

list.files(paste0(path,"ipip_data/"))

bri <- read.table(paste0(path,"ipip_data/bri.dat"), sep="\t", header=TRUE)
bri_clus <- read.table(paste0(path,"ipip_data/bri_clus.dat"), sep="\t", header=TRUE)
hexaco <- read.table(paste0(path,"ipip_data/hexaco.dat"), sep="\t", header=TRUE)
neo <- read.table(paste0(path,"ipip_data/neo.dat"), sep="\t", header=TRUE)

ipip <- Reduce(merge, list(bri_clus, neo[,1:36], hexaco))

# scale NEO variables so all manifests have similar SDs
index <- apply(expand.grid(c("n","e","o","a","c"), 1:6), 1, paste, collapse="")
ipip[,index] <- ipip[,index]/6

# multiple regression of BRI drugs on neo + hexaco
modA <- lm(drugs ~ n + e + o + a + c +
             hones + emoti + extra + agree + consc + openn, data=ipip)
modC <- lm(drugs ~ n + e + o + a + c, data=ipip)
anova(modC, modA)
#   Res.Df    RSS Df Sum of Sq      F    Pr(>F)    
# 1    598 209.14                                  
# 2    592 190.25  6    18.892 9.7974 2.611e-10 ***

# again w/ SEM
h <- names(hexaco)[2:25]
modA <- OpenMx::mxTryHard(mxOption(model=mxModel(name="modA",
  type="RAM",
  manifestVars=c("drugs", names(neo)[7:36], h),
  latentVars=c("n","e","o","a","c",
               "hones","emoti","extra","agree","consc","openn"),
  # indicators
  mxPath(from="n", to=paste0("n",1:6)),
  mxPath(from="e", to=paste0("e",1:6)),
  mxPath(from="o", to=paste0("o",1:6)),
  mxPath(from="a", to=paste0("a",1:6)),
  mxPath(from="c", to=paste0("c",1:6)),
  mxPath(from="hones", to=h[substr(h,1,1)=="h"]),
  mxPath(from="emoti", to=h[substr(h,1,1)=="e"]),
  mxPath(from="extra", to=h[substr(h,1,1)=="x"]),
  mxPath(from="agree", to=h[substr(h,1,1)=="a"]),
  mxPath(from="consc", to=h[substr(h,1,1)=="c"]),
  mxPath(from="openn", to=h[substr(h,1,1)=="o"]),
  # residual variances
  mxPath(from=c("drugs", names(neo)[7:36], h), arrows=2),
  # fix scale of latents
  mxPath(from=c("n","e","o","a","c","hones","emoti","extra","agree",
                "consc","openn"), arrows=2, free=FALSE, value=1),
  # allow latents to covary
  mxPath(from=c("n","e","o","a","c","hones","emoti","extra","agree",
                "consc","openn"), arrows=2, connect="unique.bivariate"),
  # regress Y on Ts
  mxPath(from=c("n","e","o","a","c","hones","emoti","extra","agree",
                "consc","openn"), to="drugs"),
  # data source
  mxData(observed=cov(ipip[,c("drugs", names(neo)[7:36], h)]),
         type="cov", numObs=nrow(ipip))),
  # enable parallel (arguments to mxOption)
  key="Number of Threads", value=(omxDetectCores()-1)),
  # number of times to try new starting values (argument to mxTryHard)
  extraTries=99)
# WON'T CONVERGE

summary(modA)

# check SDs
sqrt(diag(cov(ipip[,c("drugs", names(neo)[7:36], h)])))

semPaths(semPlotModel_MxRAMModel(modA), whatLabels=c("par"),
         sizeMan=3, sizeLat=9, sizeInt=1,  mar=c(5,5,5,5),
         edge.color="black", edge.label.cex=1, structural=TRUE,
         color=list(man="pink", lat="lightblue"))

omxGraphviz(modA, "/Users/Jake/Desktop/test.dot")

# IPIP, subsets of factors ------------------------------------------------

# multiple regression of BRI drugs on neo + hexaco
summary(lm(drugs ~ n + e + o + a + c + hones, data=ipip))
summary(lm(drugs ~ hones + emoti + extra + agree + consc + openn, data=ipip))

# NEO variables + HEXACO honesty
h <- names(hexaco)[2:5]
modA <- mxTryHard(mxOption(model=mxModel(name="modA",
  type="RAM",
  manifestVars=c("drugs", names(neo)[7:36], h),
  latentVars=c("n","e","o","a","c","hones"),
  # indicators
  mxPath(from="n", to=paste0("n",1:6)),
  mxPath(from="e", to=paste0("e",1:6)),
  mxPath(from="o", to=paste0("o",1:6)),
  mxPath(from="a", to=paste0("a",1:6)),
  mxPath(from="c", to=paste0("c",1:6)),
  mxPath(from="hones", to=h),
  # residual variances
  mxPath(from=c("drugs", names(neo)[7:36], h), arrows=2),
  # fix scale of latents
  mxPath(from=c("n","e","o","a","c","hones"),
         arrows=2, free=FALSE, value=1),
  # allow latents to covary
  mxPath(from=c("n","e","o","a","c","hones"),
         arrows=2, connect="unique.bivariate"),
  # regress Y on Ts
  mxPath(from=c("n","e","o","a","c","hones"), to="drugs"),
  # data source
  mxData(observed=cov(ipip[,c("drugs", names(neo)[7:36], h)]),
    type="cov", numObs=nrow(ipip))),
  # enable parallel (arguments to mxOption)
  key="Number of Threads", value=(omxDetectCores()-1)),
  # number of times to try new starting values (argument to mxTryHard)
  extraTries=99)
# WON'T CONVERGE

summary(modA)

# HEXACO variables only
h <- names(hexaco)[2:25]
modA <- mxTryHard(mxOption(model=mxModel(name="modA",
  type="RAM",
  manifestVars=c("drugs", h),
  latentVars=c("hones","emoti","extra","agree","consc","openn"),
  # indicators
  mxPath(from="hones", to=h[substr(h,1,1)=="h"]),
  mxPath(from="emoti", to=h[substr(h,1,1)=="e"]),
  mxPath(from="extra", to=h[substr(h,1,1)=="x"]),
  mxPath(from="agree", to=h[substr(h,1,1)=="a"]),
  mxPath(from="consc", to=h[substr(h,1,1)=="c"]),
  mxPath(from="openn", to=h[substr(h,1,1)=="o"]),
  # residual variances
  mxPath(from=c("drugs", h), arrows=2),
  # fix scale of latents
  mxPath(from=c("hones","emoti","extra","agree","consc","openn"),
         arrows=2, free=FALSE, value=1),
  # allow latents to covary
  mxPath(from=c("hones","emoti","extra","agree","consc","openn"),
         arrows=2, connect="unique.bivariate"),
  # regress Y on Ts
  mxPath(from=c("hones","emoti","extra","agree","consc","openn"), to="drugs"),
  # data source
  mxData(observed=cov(ipip[,c("drugs", h)]), type="cov", numObs=nrow(ipip))),
  # enable parallel (arguments to mxOption)
  key="Number of Threads", value=(omxDetectCores()-1)),
  # number of times to try new starting values (argument to mxTryHard)
  extraTries=99)
# WON'T CONVERGE
# still won't converge if "drugs" outcome changed to "religion"

summary(modA)

# NEO-PI variables only
h <- names(hexaco)[2:25]
modA <- mxTryHard(mxOption(model=mxModel(name="modA",
  type="RAM",
  manifestVars=c("drugs", names(neo)[7:36]),
  latentVars=c("n","e","o","a","c"),
  # indicators
  mxPath(from="n", to=paste0("n",1:6)),
  mxPath(from="e", to=paste0("e",1:6)),
  mxPath(from="o", to=paste0("o",1:6)),
  mxPath(from="a", to=paste0("a",1:6)),
  mxPath(from="c", to=paste0("c",1:6)),
  # residual variances
  mxPath(from=c("drugs", names(neo)[7:36]), arrows=2),
  # fix scale of latents
  mxPath(from=c("n","e","o","a","c"), arrows=2, free=FALSE, value=1),
  # allow latents to covary
  mxPath(from=c("n","e","o","a","c"), arrows=2, connect="unique.bivariate"),
  # regress Y on Ts
  mxPath(from=c("n","e","o","a","c"), to="drugs"),
  # data source
  mxData(observed=cov(ipip[,c("drugs", names(neo)[7:36])]),
  type="cov", numObs=nrow(ipip))),
  # enable parallel (arguments to mxOption)
  key="Number of Threads", value=(omxDetectCores()-1)),
  # number of times to try new starting values (argument to mxTryHard)
  extraTries=99)
# WON'T CONVERGE

summary(modA)

# NEO extraversion + HEXACO extraversion
h <- names(hexaco)[2:25]
h <- h[substr(h,1,1)=="x"]
modA <- mxTryHard(mxOption(model=mxModel(name="modA",
  type="RAM",
  manifestVars=c("drugs", paste0("e",1:6), h),
  latentVars=c("e","extra"),
  # indicators
  mxPath(from="e", to=paste0("e",1:6)),
  mxPath(from="extra", to=h),
  # residual variances
  mxPath(from=c("drugs", paste0("e",1:6), h), arrows=2),
  # fix scale of latents
  mxPath(from=c("e","extra"), arrows=2, free=FALSE, value=1),
  # allow latents to covary
  mxPath(from=c("e","extra"), arrows=2, connect="unique.bivariate"),
  # regress Y on Ts
  mxPath(from=c("e","extra"), to="drugs"),
  # data source
  mxData(observed=cov(ipip[,c("drugs", paste0("e",1:6), h)]),
    type="cov", numObs=nrow(ipip))),
  # enable parallel (arguments to mxOption)
  key="Number of Threads", value=(omxDetectCores()-1)),
  # number of times to try new starting values (argument to mxTryHard)
  extraTries=99)
# WON'T CONVERGE

summary(modA)

# NEO extraversion + HEXACO extraversion, measurement model only
h <- names(hexaco)[2:25]
h <- h[substr(h,1,1)=="x"]
modA <- mxTryHard(mxOption(model=mxModel(name="modA",
  type="RAM",
  manifestVars=c(paste0("e",1:6), h),
  latentVars=c("e","extra"),
  # indicators
  mxPath(from="e", to=paste0("e",1:6)),
  mxPath(from="extra", to=h),
  # residual variances
  mxPath(from=c(paste0("e",1:6), h), arrows=2),
  # fix scale of latents
  mxPath(from=c("e","extra"), arrows=2, free=FALSE, value=1),
  # allow latents to covary
  mxPath(from=c("e","extra"), arrows=2, connect="unique.bivariate"),
  # data source
  mxData(observed=cov(ipip[,c(paste0("e",1:6), h)]),
    type="cov", numObs=nrow(ipip))),
  # enable parallel (arguments to mxOption)
  key="Number of Threads", value=(omxDetectCores()-1)),
  # number of times to try new starting values (argument to mxTryHard)
  extraTries=99)
# WON'T CONVERGE

summary(modA)

# NEO variables only (no AC), measurement model only
modA <- mxTryHard(mxOption(model=mxModel(name="modA",
  type="RAM",
  manifestVars=c(names(neo)[7:26]),
  latentVars=c("n","e","o"),
  # indicators
  mxPath(from="n", to=paste0("n",1:6)),
  mxPath(from="e", to=paste0("e",1:6)),
  mxPath(from="o", to=paste0("o",1:6)),
  # residual variances
  mxPath(from=c(names(neo)[7:26]), arrows=2),
  # fix scale of latents
  mxPath(from=c("n","e","o"), arrows=2, free=FALSE, value=1),
  # allow latents to covary
  mxPath(from=c("n","e","o"), arrows=2, connect="unique.bivariate"),
  # data source
#   mxData(observed=cov(ipip[,c(names(neo)[7:26])]),
#     type="cov", numObs=nrow(ipip))),
  mxData(observed=ipip[,c(names(neo)[7:26])],
         type="raw"),
  # means/intercepts and thresholds
  mxPath(from="one", to=c("n","e","o"), free=FALSE, values=0),
  mxPath(from="one", to=names(neo)[7:26])),
  # enable parallel (arguments to mxOption)
  key="Number of Threads", value=(omxDetectCores()-1)),
  # number of times to try new starting values (argument to mxTryHard)
  extraTries=99)
# WON'T CONVERGE

### try with lavaan

model <- "n =~ n1 + n2 + n3 + n4 + n5 + n6
          e =~ e1 + e2 + e3 + e4 + e5 + e6
          o =~ o1 + o2 + o3 + o4 + o5 + o6"
fit <- cfa(model, std.lv=TRUE, data=ipip[,c(names(neo)[7:26])])

summary(fit)

# IPIP in lavaan ----------------------------------------------------------

# NEO variables + HEXACO honesty

summary(lm(drugs ~ n + e + o + a + c + hones, data=ipip))

mod <- sem("
           # measurement models
           n =~ n1 + n2 + n3 + n4 + n5 + n6
           e =~ e1 + e2 + e3 + e4 + e5 + e6
           o =~ o1 + o2 + o3 + o4 + o5 + o6
           a =~ a1 + a2 + a3 + a4 + a5 + a6
           c =~ c1 + c2 + c3 + c4 + c5 + c6
           honesty =~ hsinc + hfair + hgree + hmode
           # regressions
           drugs ~ n + e + o + a + c + honesty
           ", data=ipip)
summary(mod)

semPaths(mod, whatLabels=c("par"), layout="tree",
         sizeMan=3, sizeLat=5, sizeInt=1, mar=c(3,3,3,3),
         edge.color="black", edge.label.cex=1, # structural=TRUE,
         color=list(man="pink", lat="lightblue"))

# HEXACO variables only

summary(lm(drugs ~ emoti + extra + agree + consc + openn + hones, data=ipip))

mod <- sem("
           # measurement models
           x =~ xexpr + xsocb + xsoci + xlive
           e =~ efear + eanxi + edepe + esent
           o =~ oaesa + oinqu + ocrea + ounco
           a =~ aforg + agent + aflex + apati
           c =~ corga + cdili + cperf + cprud
           honesty =~ hsinc + hfair + hgree + hmode
           # regressions
           drugs ~ x + e + o + a + c + honesty
           ", data=ipip)
summary(mod)

semPaths(mod, whatLabels=c("par"), layout="tree",
         sizeMan=3, sizeLat=5, sizeInt=1, mar=c(5,5,5,5),
         edge.color="black", edge.label.cex=1, # structural=TRUE,
         color=list(man="pink", lat="lightblue"))

# E: emotional stability / neuroticism

summary(lm(drugs ~ n + emoti, data=ipip)) # both significant

eMod <- sem("
            # measurement models
            neo_emo =~ n1 + n2 + n3 + n4 + n5 + n6
            hex_emo =~ efear + eanxi + edepe + esent
            # regressions
            drugs ~ neo_emo + hex_emo
            ", data=ipip, std.lv=TRUE)
cov2cor(lavInspect(eMod, "cov.lv")) # r = .723
loads <- colSums(lavInspect(eMod, "est")$lambda[,-3])^2
loads[1]/(loads[1] + sum(diag(lavInspect(eMod, "est")$theta)[1:6])) # .87
loads[2]/(loads[2] + sum(diag(lavInspect(eMod, "est")$theta)[7:10])) # .64
summary(eMod) # both significant

semPaths(eMod, whatLabels=c("par"), layout="tree",
         sizeMan=5, sizeLat=5, sizeInt=1, mar=c(5,5,5,5),
         edge.color="black", edge.label.cex=1, # structural=TRUE,
         color=list(man="pink", lat="lightblue"))


# X: extraversion

summary(lm(drugs ~ e + extra, data=ipip)) # only NEO sig

xMod <- sem("
            # measurement models
            neo_extra =~ e1 + e2 + e3 + e4 + e5 + e6
            hex_extra =~ xexpr + xsocb + xsoci + xlive
            # regressions
            drugs ~ neo_extra + hex_extra
            ", data=ipip, std.lv=TRUE)
cov2cor(lavInspect(xMod, "cov.lv")) # r > 1, non-converged
loads <- colSums(lavInspect(xMod, "est")$lambda[,-3])^2
loads[1]/(loads[1] + sum(diag(lavInspect(xMod, "est")$theta)[1:6]))
loads[2]/(loads[2] + sum(diag(lavInspect(xMod, "est")$theta)[7:10]))
summary(xMod)

# A: agreeableness

summary(lm(drugs ~ a + agree, data=ipip)) # both sig

aMod <- sem("
            # measurement models
            neo_agree =~ a1 + a2 + a3 + a4 + a5 + a6
            hex_agree =~ aforg + agent + aflex + apati
            # regressions
            drugs ~ neo_agree + hex_agree
            ", data=ipip, std.lv=TRUE)
cov2cor(lavInspect(aMod, "cov.lv")) # r = .767
loads <- colSums(lavInspect(aMod, "est")$lambda[,-3])^2
loads[1]/(loads[1] + sum(diag(lavInspect(aMod, "est")$theta)[1:6])) # .73
loads[2]/(loads[2] + sum(diag(lavInspect(aMod, "est")$theta)[7:10])) # .79
summary(aMod) # only NEO sig

# C: conscientousness

summary(lm(drugs ~ c + consc, data=ipip)) # only neo sig

cMod <- sem("
            # measurement models
            neo_consc =~ c1 + c2 + c3 + c4 + c5 + c6
            hex_consc =~ corga + cdili + cperf + cprud
            # regressions
            drugs ~ neo_consc + hex_consc
            ", data=ipip, std.lv=TRUE)
cov2cor(lavInspect(cMod, "cov.lv")) # r = .951
loads <- colSums(lavInspect(cMod, "est")$lambda[,-3])^2
loads[1]/(loads[1] + sum(diag(lavInspect(cMod, "est")$theta)[1:6])) # .83
loads[2]/(loads[2] + sum(diag(lavInspect(cMod, "est")$theta)[7:10])) # .67
summary(cMod) # neither sig

# O: openness

summary(lm(drugs ~ o + openn, data=ipip)) # only NEO sig

oMod <- sem("
            # measurement models
            neo_open =~ o1 + o2 + o3 + o4 + o5 + o6
            hex_open =~ oaesa + oinqu + ocrea + ounco
            # regressions
            drugs ~ neo_open + hex_open
            ", data=ipip, std.lv=TRUE)
cov2cor(lavInspect(oMod, "cov.lv")) # r = .977
loads <- colSums(lavInspect(oMod, "est")$lambda[,-3])^2
loads[1]/(loads[1] + sum(diag(lavInspect(oMod, "est")$theta)[1:6])) # .79
loads[2]/(loads[2] + sum(diag(lavInspect(oMod, "est")$theta)[7:10])) # .76
summary(oMod) # neither sig

# IPIP, predicting all BRI clusters ---------------------------------------

getSEM <- function(x){
  cat(paste(x,"\n"))
  
  # NEO variables + HEXACO honesty
  form <- as.formula(paste(x,"~ n + e + o + a + c + hones"))
  hReg <- coef(summary(lm(form, data=ipip)))["hones","t value"]
  hMod <- sem(paste0("
           # measurement models
           n =~ n1 + n2 + n3 + n4 + n5 + n6
           e =~ e1 + e2 + e3 + e4 + e5 + e6
           o =~ o1 + o2 + o3 + o4 + o5 + o6
           a =~ a1 + a2 + a3 + a4 + a5 + a6
           c =~ c1 + c2 + c3 + c4 + c5 + c6
           honesty =~ hsinc + hfair + hgree + hmode
           # regressions
           ",x," ~ n + e + o + a + c + honesty
           "), data=ipip, std.lv=TRUE)
  hSEM <- unname(lavInspect(hMod, "est")$beta[x,"honesty"] /
    sqrt(diag(vcov(hMod)))[paste0(x,"~honesty")])
  
  # E: emotional stability / neuroticism
  form <- as.formula(paste(x,"~ n + emoti"))
  eReg <- unname(coef(summary(lm(form, data=ipip)))[2:3,"t value"])
  eMod <- sem(paste0("
              # measurement models
              neo_emo =~ n1 + n2 + n3 + n4 + n5 + n6
              hex_emo =~ efear + eanxi + edepe + esent
              # regressions
              ",x," ~ neo_emo + hex_emo
              "), data=ipip, std.lv=TRUE)
  eSEM <- unname(lavInspect(eMod, "est")$beta[x,c("neo_emo","hex_emo")] /
    sqrt(diag(vcov(eMod)))[paste0(x,c("~neo_emo","~hex_emo"))])
  
  # X: extraversion
  form <- as.formula(paste(x,"~ e + extra"))
  xReg <- unname(coef(summary(lm(form, data=ipip)))[2:3,"t value"])
  xMod <- sem(paste0("
              # measurement models
              neo_extra =~ e1 + e2 + e3 + e4 + e5 + e6
              hex_extra =~ xexpr + xsocb + xsoci + xlive
              # regressions
              ",x," ~ neo_extra + hex_extra
              "), data=ipip, std.lv=TRUE)
  xSEM <- unname(lavInspect(xMod, "est")$beta[x,c("neo_extra","hex_extra")] /
    sqrt(diag(vcov(xMod)))[paste0(x,c("~neo_extra","~hex_extra"))])
  
  # A: agreeableness
  form <- as.formula(paste(x,"~ a + agree"))
  aReg <- unname(coef(summary(lm(form, data=ipip)))[2:3,"t value"])
  aMod <- sem(paste0("
              # measurement models
              neo_agree =~ a1 + a2 + a3 + a4 + a5 + a6
              hex_agree =~ aforg + agent + aflex + apati
              # regressions
              ",x," ~ neo_agree + hex_agree
              "), data=ipip, std.lv=TRUE)
  aSEM <- unname(lavInspect(aMod, "est")$beta[x,c("neo_agree","hex_agree")] /
    sqrt(diag(vcov(aMod)))[paste0(x,c("~neo_agree","~hex_agree"))])
  
  # C: conscientousness
  form <- as.formula(paste(x,"~ c + consc"))
  cReg <- unname(coef(summary(lm(form, data=ipip)))[2:3,"t value"])
  cMod <- sem(paste0("
              # measurement models
              neo_consc =~ c1 + c2 + c3 + c4 + c5 + c6
              hex_consc =~ corga + cdili + cperf + cprud
              # regressions
              ",x," ~ neo_consc + hex_consc
              "), data=ipip, std.lv=TRUE)
  cSEM <- unname(lavInspect(cMod, "est")$beta[x,c("neo_consc","hex_consc")] /
    sqrt(diag(vcov(cMod)))[paste0(x,c("~neo_consc","~hex_consc"))])
  
  # O: openness
  form <- as.formula(paste(x,"~ o + openn"))
  oReg <- unname(coef(summary(lm(form, data=ipip)))[2:3,"t value"])
  oMod <- sem(paste0("
            # measurement models
            neo_open =~ o1 + o2 + o3 + o4 + o5 + o6
            hex_open =~ oaesa + oinqu + ocrea + ounco
            # regressions
            ",x," ~ neo_open + hex_open
            "), data=ipip, std.lv=TRUE)
  oSEM <- unname(lavInspect(oMod, "est")$beta[x,c("neo_open","hex_open")] /
    sqrt(diag(vcov(oMod)))[paste0(x,c("~neo_open","~hex_open"))])

  ans <- c(hConverge=lavInspect(hMod, "converged"),
           hReg=hReg, hSEM=hSEM,
           eConverge=lavInspect(eMod, "converged"),
           eNeoReg=eReg[1], eNeoSEM=eSEM[1],
           eHexReg=eReg[2], eHexSEM=eSEM[2],
           xConverge=lavInspect(xMod, "converged"),
           xNeoReg=xReg[1], xNeoSEM=xSEM[1],
           xHexReg=xReg[2], xHexSEM=xSEM[2],
           aConverge=lavInspect(aMod, "converged"),
           aNeoReg=aReg[1], aNeoSEM=aSEM[1],
           aHexReg=aReg[2], aHexSEM=aSEM[2],
           cConverge=lavInspect(cMod, "converged"),
           cNeoReg=cReg[1], cNeoSEM=cSEM[1],
           cHexReg=cReg[2], cHexSEM=cSEM[2],
           oConverge=lavInspect(oMod, "converged"),
           oNeoReg=oReg[1], oNeoSEM=oSEM[1],
           oHexReg=oReg[2], oHexSEM=oSEM[2])
  return(ans)
}

# test on 2 outcomes
sapply(c("drugs","tv"), getSEM)

# do it on all outcomes!
system.time({
  results <- sapply(names(bri_clus)[-1], getSEM)
}) # 90 seconds

# check convergence
rowSums(results[paste0(c("h","e","x","a","c","o"),"Converge"),])
# all converged

### plot results

# HEXACO honesty
tiff(paste0(path,"Fig8.tif"),
    units="in", height=5, width=6, res=200, pointsize=15)
plot(y=abs(results["hSEM",]), x=abs(results["hReg",]), pch=20,
     ylim=range(abs(results["hReg",])),
     main="Effect of HEXACO 'Honesty/Humility'\ncontrolling for 5 NEO PI-R factors",
     ylab="SEM absolute test statistics",
     xlab="Regression absolute test statistics")
legend("topleft", lty=1:3, lwd=3:1, bty="n",
       legend=c("Identity line","Regression line","Critical value"))
abline(0,1,lwd=3)
abline(lm(abs(results["hSEM",]) ~ abs(results["hReg",])), lwd=2, lty=2)
abline(h=qnorm(.975), v=qt(.975, df=597), lty=3)
dev.off()

# NEO over HEXACO

tiff(paste0(path,"Fig7.tif"),
    units="in", height=5, width=7.5, res=200, pointsize=15)
layout(matrix(1:6, nrow=2, byrow=TRUE))
par(mar=c(4,3,3,1)+.1)
par(oma=c(0,0,2,0))
main <- c("A: Emotionality \n/ Neuroticism","B: Extraversion",
          "C: Agreeableness","D: Conscientiousness","E: Openness\nto experience")
names(main) <- c("e","x","a","c","o")
for(x in c("e","x","a","c","o")){
  plot(y=abs(results[paste0(x,"NeoSEM"),]), mgp=2:0,
       x=abs(results[paste0(x,"NeoReg"),]), pch=20,
       ylim=range(abs(results[paste0(x,"NeoReg"),])),
       main=main[x],
       ylab="SEM absolute test statistics",
       xlab="Regression absolute test statistics")
  abline(0,1,lwd=3)
  abline(lm(abs(results[paste0(x,"NeoSEM"),]) ~ abs(results[paste0(x,"NeoReg"),])),
         lwd=2, lty=2)
  abline(h=qnorm(.975), v=qt(.975, df=601), lty=3)
}
par(mar=c(4,3,3,1)+.1)
plot(y=0:1, x=0:1, cex=0, bty="n", ylab="", xlab="", yaxt="n", xaxt="n")
legend("topleft", lty=1:3, lwd=3:1, bty="n",
       legend=c("Identity line","Regression line","Critical value"))
mtext("Effects of NEO factors, controlling for HEXACO factors", outer=TRUE)
dev.off()

# HEXACO over NEO

tiff(paste0(path,"Fig6.tif"),
    units="in", height=5, width=7.5, res=200, pointsize=15)
layout(matrix(1:6, nrow=2, byrow=TRUE))
par(mar=c(4,3,3,1)+.1)
par(oma=c(0,0,2,0))
main <- c("A: Emotionality \n/ Neuroticism","B: Extraversion",
          "C: Agreeableness","D: Conscientiousness","E: Openness\nto experience")
names(main) <- c("e","x","a","c","o")
for(x in c("e","x","a","c","o")){
  plot(y=abs(results[paste0(x,"HexSEM"),]), mgp=2:0,
       x=abs(results[paste0(x,"HexReg"),]), pch=20,
       ylim=range(abs(results[paste0(x,"HexReg"),])),
       main=main[x],
       ylab="SEM absolute test statistics",
       xlab="Regression absolute test statistics")
  abline(0,1,lwd=3)
  abline(lm(abs(results[paste0(x,"HexSEM"),]) ~ abs(results[paste0(x,"HexReg"),])),
         lwd=2, lty=2)
  abline(h=qnorm(.975), v=qt(.975, df=601), lty=3)
}
par(mar=c(4,3,3,1)+.1)
plot(y=0:1, x=0:1, cex=0, bty="n", ylab="", xlab="", yaxt="n", xaxt="n")
legend("topleft", lty=1:3, lwd=3:1, bty="n",
       legend=c("Identity line","Regression line","Critical value"))
mtext("Effects of HEXACO factors, controlling for NEO factors", outer=TRUE)
dev.off()

# single indicator analysis for a single behavior -------------------------

# scale down sum-scores
ipip[,c("n","e","o","a","c")] <- ipip[,c("n","e","o","a","c")]/6/6
sqrt(diag(var(ipip[,c("n","e","o","a","c","hones","drugs")])))
#         n         e         o         a         c     hones     drugs 
# 0.6447893 0.5522065 0.5865386 0.4631684 0.5190752 0.4516519 0.6487788 

alpha <- .8
dv <- "drugs"

# function to compute z-statistics as f(alpha) for given DV
getZ <- function(dv, alpha, allowInadmissible){
  # NEO variables + HEXACO honesty
  mod <- sem(paste("
                   # measurement models
                   N  =~ 1*n
                   E  =~ 1*e
                   O  =~ 1*o
                   A  =~ 1*a
                   C  =~ 1*c
                   H =~ 1*hones
                   n ~~ (1 -",alpha,")*",var(ipip$n),"*n
                   e ~~ (1 -",alpha,")*",var(ipip$e),"*e
                   o ~~ (1 -",alpha,")*",var(ipip$o),"*o
                   a ~~ (1 -",alpha,")*",var(ipip$a),"*a
                   c ~~ (1 -",alpha,")*",var(ipip$c),"*c
                   hones ~~ (1 -",alpha,")*",var(ipip$hones),"*hones
                   # regressions
                   ",dv," ~ N + E + O + A + C + H
                   "), std.lv=TRUE, data=ipip)
  if(lavInspect(mod, "converged")){
    zH <- try(lavInspect(mod, "est")$beta[dv,"H"]/
                    sqrt(diag(vcov(mod)))[paste0(dv,"~H")],
                  silent=TRUE)
    if(class(zH) == "try-error") zH <- NA
  } else zH <- NA
  if(!allowInadmissible & det(inspect(mod,"cov.lv")) < 0){
    zH <- NA
  }
  
  # E: emotional stability / neuroticism
  eMod <- sem(paste("
                    ### loadings
                    neo_emo  =~ 1*n
                    hex_emo =~ 1*emoti
                    ### residual variances
                    n ~~ (1 -",alpha,")*",var(ipip$n),"*n
                    emoti ~~ (1 -",alpha,")*",var(ipip$emoti),"*emoti
                    ### regressions
                    ",dv," ~ neo_emo + hex_emo
                    "), data=ipip, std.lv=TRUE)
  if(lavInspect(eMod, "converged")){
    zE_neo <- try(lavInspect(eMod, "est")$beta[dv,"neo_emo"]/
                sqrt(diag(vcov(eMod)))[paste0(dv,"~neo_emo")],
              silent=TRUE)
    if(class(zE_neo) == "try-error") zE_neo <- NA
    zE_hex <- try(lavInspect(eMod, "est")$beta[dv,"hex_emo"]/
                    sqrt(diag(vcov(eMod)))[paste0(dv,"~hex_emo")],
                  silent=TRUE)
    if(class(zE_hex) == "try-error") zE_hex <- NA
  } else {
    zE_neo <- NA
    zE_hex <- NA
  }
  if(!allowInadmissible & det(inspect(eMod,"cov.lv")) < 0){
    zE_neo <- NA
    zE_hex <- NA
  }

  # X: extraversion
  xMod <- sem(paste("
                    ### loadings
                    neo_extra =~ 1*e
                    hex_extra =~ 1*extra
                    ### residual variances
                    e ~~ (1 -",alpha,")*",var(ipip$e),"*e
                    extra ~~ (1 -",alpha,")*",var(ipip$extra),"*extra
                    ### regressions
                    ",dv," ~ neo_extra + hex_extra
                    "), data=ipip, std.lv=TRUE)
  if(lavInspect(xMod, "converged")){
    zX_neo <- try(lavInspect(xMod, "est")$beta[dv,"neo_extra"]/
                    sqrt(diag(vcov(xMod)))[paste0(dv,"~neo_extra")],
                  silent=TRUE)
    if(class(zX_neo) == "try-error") zX_neo <- NA
    zX_hex <- try(lavInspect(xMod, "est")$beta[dv,"hex_extra"]/
                    sqrt(diag(vcov(xMod)))[paste0(dv,"~hex_extra")],
                  silent=TRUE)
    if(class(zX_hex) == "try-error") zX_hex <- NA
  } else {
    zX_neo <- NA
    zX_hex <- NA
  }
  if(!allowInadmissible & det(inspect(xMod,"cov.lv")) < 0){
    zX_neo <- NA
    zX_hex <- NA
  }

  # A: agreeableness
  aMod <- sem(paste("
                    ### loadings
                    neo_agree =~ 1*a
                    hex_agree =~ 1*agree
                    ### residual variances
                    a ~~ (1 -",alpha,")*",var(ipip$a),"*a
                    agree ~~ (1 -",alpha,")*",var(ipip$agree),"*agree
                    ### regressions
                    ",dv," ~ neo_agree + hex_agree
                    "), data=ipip, std.lv=TRUE)
  if(lavInspect(aMod, "converged")){
    zA_neo <- try(lavInspect(aMod, "est")$beta[dv,"neo_agree"]/
                    sqrt(diag(vcov(aMod)))[paste0(dv,"~neo_agree")],
                  silent=TRUE)
    if(class(zA_neo) == "try-error") zA_neo <- NA
    zA_hex <- try(lavInspect(aMod, "est")$beta[dv,"hex_agree"]/
                    sqrt(diag(vcov(aMod)))[paste0(dv,"~hex_agree")],
                  silent=TRUE)
    if(class(zA_hex) == "try-error") zA_hex <- NA
  } else {
    zA_neo <- NA
    zA_hex <- NA
  }
  if(!allowInadmissible & det(inspect(aMod,"cov.lv")) < 0){
    zA_neo <- NA
    zA_hex <- NA
  }

  # C: conscientousness
  cMod <- sem(paste("
                    ### loadings
                    neo_consc =~ 1*c
                    hex_consc =~ 1*consc
                    ### residual variances
                    c ~~ (1 -",alpha,") *",var(ipip$c),"*c
                    consc ~~ (1 -",alpha,")*",var(ipip$consc),"*consc
                    ### regressions
                    ",dv," ~ neo_consc + hex_consc
                    "), data=ipip, std.lv=TRUE)
  if(lavInspect(cMod, "converged")){
    zC_neo <- try(lavInspect(cMod, "est")$beta[dv,"neo_consc"]/
                    sqrt(diag(vcov(cMod)))[paste0(dv,"~neo_consc")],
                  silent=TRUE)
    if(class(zC_neo) == "try-error") zC_neo <- NA
    zC_hex <- try(lavInspect(cMod, "est")$beta[dv,"hex_consc"]/
                    sqrt(diag(vcov(cMod)))[paste0(dv,"~hex_consc")],
                  silent=TRUE)
    if(class(zC_hex) == "try-error") zC_hex <- NA
  } else {
    zC_neo <- NA
    zC_hex <- NA
  }
  if(!allowInadmissible & det(inspect(cMod,"cov.lv")) < 0){
    zC_neo <- NA
    zC_hex <- NA
  }

  # O: openness
  oMod <- sem(paste("
                    ### loadings
                    neo_open =~ 1*o
                    hex_open =~ 1*openn
                    ### residual variances
                    o ~~ (1 -",alpha,")*",var(ipip$o),"*o
                    openn ~~ (1 -",alpha,")*",var(ipip$openn),"*openn
                    ### regressions
                    ",dv," ~ neo_open + hex_open
                    "), data=ipip, std.lv=TRUE)
  if(lavInspect(oMod, "converged")){
    zO_neo <- try(lavInspect(oMod, "est")$beta[dv,"neo_open"]/
                    sqrt(diag(vcov(oMod)))[paste0(dv,"~neo_open")],
                  silent=TRUE)
    if(class(zO_neo) == "try-error") zO_neo <- NA
    zO_hex <- try(lavInspect(oMod, "est")$beta[dv,"hex_open"]/
                    sqrt(diag(vcov(oMod)))[paste0(dv,"~hex_open")],
                  silent=TRUE)
    if(class(zO_hex) == "try-error") zO_hex <- NA
  } else {
    zO_neo <- NA
    zO_hex <- NA
  }
  if(!allowInadmissible & det(inspect(oMod,"cov.lv")) < 0){
    zO_neo <- NA
    zO_hex <- NA
  }
  
  return(c(zH=unname(zH),
           zE_neo=unname(zE_neo), zE_hex=unname(zE_hex),
           zX_neo=unname(zX_neo), zX_hex=unname(zX_hex),
           zA_neo=unname(zA_neo), zA_hex=unname(zA_hex),
           zC_neo=unname(zC_neo), zC_hex=unname(zC_hex),
           zO_neo=unname(zO_neo), zO_hex=unname(zO_hex)))
  }

# "drugs" DV, allowing inadmissible solutions
siResults <- sapply(seq(.1,1,.1), function(x){
  getZ(dv="drugs", alpha=x, allowInadmissible=TRUE)
})
round(siResults, 2)
#         [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10]
# zH     -3.43 -2.91 -1.57 -0.17 -1.25  0.03 -1.15 -2.19 -3.13 -3.99
# zE_neo -5.07 -4.19 -2.83    NA  2.53  4.21  4.54  4.39  4.10  3.79
# zE_hex  2.69  2.79  2.30    NA -2.91 -5.71 -7.19 -8.05 -8.65 -9.12
# zX_neo  1.65  1.15  0.64  0.12 -0.41 -0.92 -1.34  1.69  2.46  2.90
# zX_hex  5.03  4.63  4.19  3.71  3.20  2.65  1.98 -1.39 -1.09 -0.56
# zA_neo -2.96 -1.88 -0.81  0.22  0.81 -2.00 -3.34 -4.36 -5.31 -6.20
# zA_hex -6.12 -5.08 -4.02 -2.89 -1.29  0.92  0.10 -0.86 -1.79 -2.68
# zC_neo -1.14 -0.76 -0.37  0.02  0.40  0.76  0.42 -1.49 -1.83 -2.15
# zC_hex -3.90 -3.50 -3.07 -2.62 -2.15 -1.65 -0.45  0.73  0.26 -0.19
# zO_neo  4.49  3.67  2.82  1.95  1.07  0.20 -0.59  1.46  2.26  2.98
# zO_hex  7.11  6.27  5.36  4.39  3.37  2.31  1.14 -0.21  0.81  1.78

# "drugs" DV, not allowing inadmissible solutions
siResults <- sapply(seq(.1,1,.1), function(x){
  getZ(dv="drugs", alpha=x, allowInadmissible=FALSE)
})
round(siResults, 2)
#        [,1] [,2] [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10]
# zH       NA   NA   NA -0.17    NA  0.03 -1.15 -2.19 -3.13 -3.99
# zE_neo   NA   NA   NA    NA  2.53  4.21  4.54  4.39  4.10  3.79
# zE_hex   NA   NA   NA    NA -2.91 -5.71 -7.19 -8.05 -8.65 -9.12
# zX_neo   NA   NA   NA    NA    NA    NA    NA  1.69  2.46  2.90
# zX_hex   NA   NA   NA    NA    NA    NA    NA -1.39 -1.09 -0.56
# zA_neo   NA   NA   NA    NA    NA -2.00 -3.34 -4.36 -5.31 -6.20
# zA_hex   NA   NA   NA    NA    NA  0.92  0.10 -0.86 -1.79 -2.68
# zC_neo   NA   NA   NA    NA    NA    NA    NA -1.49 -1.83 -2.15
# zC_hex   NA   NA   NA    NA    NA    NA    NA  0.73  0.26 -0.19
# zO_neo   NA   NA   NA    NA    NA    NA    NA  1.46  2.26  2.98
# zO_hex   NA   NA   NA    NA    NA    NA    NA -0.21  0.81  1.78

# interpolate more values and plot
# elapsed time = 50 seconds
xaxis <- seq(.1,1,.01)
system.time({
  siResults <- sapply(xaxis, function(x){
    getZ(dv="drugs", alpha=x, allowInadmissible=FALSE)
  })
})

# a few low alpha values weirdly converged for H factor
# set them to NA for less confusion
h2 <- siResults["zH",]
h2[1:45] <- NA

# plot!
tiff(paste0(path,"Fig10.tif"),
    units="in", height=5, width=7.5, res=200, pointsize=14)
layout(matrix(1:6, nrow=2, byrow=TRUE))
par(mar=c(3,4,3,1)+.1)
par(oma=c(0,0,0,2))
x <- c(NA,.7,.8, .8,.85,.7,.8,.8,.85,.8,.85)
y <- c(NA, 3,-5,3.5, -3,-5, 3,-3,  3, 3, -3)
sapply(c(2,4,10,6,8), function(i){
  plot(x=c(.2,1), y=c(-9,4.5), cex=0, mgp=2:0,
       ylab="z-statistics from regressing\ndrug use on both predictors",
       xlab="Assumed reliability of predictors")
  polygon(x=c(.1,1.1,1.1,.1,.1), y=c(1.96,1.96,-1.96,-1.96,1.96),
          col=rgb(0,0,0,.2), border=NA)
  abline(h=c(0, 1.96, -1.96), lty=c(1,2,2))
  matlines(y=t(siResults[c(i,i+1),]), x=xaxis,
          type="l", col=c("red","blue"), lty=1, lwd=3)
  mtext(switch(as.character(i),
               "2"="A: Emotionality/Neuroticism",
               "4"="B: Extraversion",
               "10"="C: Openness to experience",
               "6"="D: Agreeableness",
               "8"="E: Conscientiousness"), cex=.8, line=.5)
  text(x=x[c(i,i+1)], y=y[c(i,i+1)], col=c("red","blue"),
       labels=c("NEO","HEXACO"))
})
plot(x=c(.2,1), y=c(-9,4.5), cex=0, mgp=2:0,
     ylab="z-statistic for\nHonesty/Humility predictor",
     xlab="Assumed reliability of predictors")
polygon(x=c(.1,1.1,1.1,.1,.1), y=c(1.96,1.96,-1.96,-1.96,1.96),
        col=rgb(0,0,0,.2), border=NA)
abline(h=c(0, 1.96, -1.96), lty=c(1,2,2))
lines(y=h2, x=xaxis, lwd=3, col="blue")
mtext("F: HEXACO 'Honesty/Humility'\ncontrolling for 5 NEO factors", cex=.8)
dev.off()

# do again with other outcomes
# achievement, religion, work, anger
system.time({
  siResults2 <- sapply(xaxis, function(x){
    getZ(dv="anger", alpha=x, allowInadmissible=FALSE)
  })
})

layout(matrix(1:6, nrow=2, byrow=TRUE))
par(mar=c(3,4,3,1)+.1)
par(oma=c(0,0,0,1))
sapply(c(2,4,10,6,8), function(i){
  matplot(y=abs(t(siResults2[c(i,i+1),])), x=xaxis, mgp=2:0,
          type="l", col=c("red","blue"), lty=1, lwd=3,
          ylab="Absolute z-statistic\nfor predicting 'drugs' outcome",
          xlab="Assumed reliability of all measures")
  abline(h=c(0, 1.96), lty=1:2)
  mtext(switch(as.character(i),
               "2"="Emotionality/Neuroticism",
               "4"="Extraversion",
               "10"="Openness to experience",
               "6"="Agreeableness",
               "8"="Conscientiousness"), cex=.8, line=.5)
})
plot(y=abs(siResults2["zH",]), x=xaxis, type="l", lwd=3, mgp=2:0, col="blue",
     ylab="Absolute z-statistic\nfor predicting 'drugs' outcome",
     xlab="Assumed reliability of all measures")
abline(h=c(0, 1.96), lty=1:2)
mtext("HEXACO 'Honesty/Humility'\ncontrolling for 5 NEO factors", cex=.8)

