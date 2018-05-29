## WHITE SHARK DETERRENT TRIAL ANALYSIS
## Corey Bradshaw March-May 2018 (data from Charlie Huveneers)
## Flinders University, College of Science and Engineering
## accompanies online report: https://www.dpi.nsw.gov.au/__data/assets/pdf_file/0007/815164/Shark-response-to-personal-deterrents_Flinders.pdf

## Remove everything
rm(list = ls())

## source functions & libraries
source("/.../new_lmer_AIC_tables3.R") # update source-file folder
source("/.../r.squared.R") # update source-file folder
library(lme4)
library(boot)
library(Hmisc)
library(ggplot2)

## functions
delta.IC <- function(x) x - min(x) ## where x is a vector of an IC
weight.IC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dIC


## DISTANCE OF APPROACH
## set working directory
setwd("~/.../")

# import distance data
dist.dat <- read.table("dist_data.csv", header=T, sep=",")
dist.dat$trialset <- ordered(dist.dat$Trial_set)
dist.dat$trial <- ordered(dist.dat$Trial)

# remove zero distances
dist.no0dist <- subset(dist.dat, dist > 0)
dist.no0dist$ldist <- log10(dist.no0dist$dist)

## GLMM
# remove unknown sharks
dist.nounk <- subset(dist.no0dist, ID != "unknown")
dist.nounk$ID <- factor(dist.nounk$ID)

# model set
m1 <- "ldist ~ det + trialset + (1|ID)"
m2 <- "ldist ~ det + (1|ID)"
m3 <- "ldist ~ trialset + (1|ID)"
m4 <- "ldist ~ 1 + (1|ID)"

## Make model vector
mod.vec <- c(m1,m2,m3,m4)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]), data=dist.nounk, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- r.squared(fit)$AIC
  Rm[i] <- 100*r.squared(fit)$Marginal # marginal R-squared
  Rc[i] <- 100*r.squared(fit)$Conditional # conditional R-squared
  print(i)
}

dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)

sumtable <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,10),round(Rm,4),round(Rc,4))
colnames(sumtable) <- c("model","k","LL","AICc","dAICc","wAICc","Rc","Rm")
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,5],decreasing=F),1:8]
summary.table

fit.sat <- lmer(as.formula(mod.vec[1]), data=dist.nounk, na.action=na.omit)
summary(fit.sat)

# redo treating trialset as covariate (ignoring 'time' per se, but using trialset as temporal marker)
dist.nounk.tscov <- dist.nounk
dist.nounk.tscov$trialset <- as.integer(dist.nounk.tscov$trialset)

# model set
m1 <- "ldist ~ det + trialset + det*trialset + (1|ID)"
m2 <- "ldist ~ det + trialset + (1|ID)"
m3 <- "ldist ~ det + (1|ID)"
m4 <- "ldist ~ trialset + (1|ID)"
m5 <- "ldist ~ 1 + (1|ID)"

## Make model vector
mod.vec <- c(m1,m2,m3,m4,m5)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]), data=dist.nounk.tscov, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- r.squared(fit)$AIC
  Rm[i] <- 100*r.squared(fit)$Marginal # marginal R-squared
  Rc[i] <- 100*r.squared(fit)$Conditional # conditional R-squared
  print(i)
}

dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)

sumtable <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,10),round(Rm,4),round(Rc,4))
colnames(sumtable) <- c("model","k","LL","AICc","dAICc","wAICc","Rc","Rm")
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,5],decreasing=F),1:8]
summary.table

fit.sat <- lmer(as.formula(mod.vec[1]), data=dist.nounk.tscov, na.action=na.omit)
summary(fit.sat)


## DISTANCE OF APPROACH (SS & controls only)
setwd("~/Documents/Papers/Fish/Sharks/White sharks/deterrant trials/")
dist2.dat <- subset(dist.dat, det == "SS" | det == "Control")
dist2.dat$trialset <- as.integer(dist2.dat$trialset)

# remove zero distances
dist2.no0dist <- subset(dist2.dat, dist > 0)
dist2.no0dist$ldist <- log10(dist2.no0dist$dist)


## GLMM
# remove unknown sharks
dist2.nounk <- subset(dist2.no0dist, ID != "unknown")
dist2.nounk$ID <- factor(dist2.nounk$ID)

# model set
m1 <- "ldist ~ det + trialset + (1|ID)"
m2 <- "ldist ~ det + (1|ID)"
m3 <- "ldist ~ trialset + (1|ID)"
m4 <- "ldist ~ 1 + (1|ID)"

## Make model vector
mod.vec <- c(m1,m2,m3,m4)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]), data=dist2.nounk, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- r.squared(fit)$AIC
  Rm[i] <- 100*r.squared(fit)$Marginal # marginal R-squared
  Rc[i] <- 100*r.squared(fit)$Conditional # conditional R-squared
  print(i)
}

dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)

sumtable <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,10),round(Rm,4),round(Rc,4))
colnames(sumtable) <- c("model","k","LL","AICc","dAICc","wAICc","Rc","Rm")
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,5],decreasing=F),1:8]
summary.table

fit.sat <- lmer(as.formula(mod.vec[1]), data=dist2.nounk, na.action=na.omit)
summary(fit.sat)


## TIME TO FEEDING
##########################
## set working directory
setwd("~/.../")

# import distance data
TF.dat <- read.table("TF_data.csv", header=T, sep=",")
TF.dat$trialset <- ordered(TF.dat$Trial_set)
TF.dat$trial <- factor(TF.dat$Trial)

# remove zero distances
TF.dat$ltime <- log10(TF.dat$time)

## GLMM
# remove unknown sharks
TF.nounk <- subset(TF.dat, ID != "unknown")
TF.nounk$ID <- factor(TF.nounk$ID)
TF.nounk$stime <- scale(log10(TF.nounk$time), center=T, scale=T)

# model set
m1 <- "stime ~ det + trialset + (1|ID)"
m2 <- "stime ~ det + (1|ID)"
m3 <- "stime ~ trialset + (1|ID)"
m4 <- "stime ~ 1 + (1|ID)"

## Make model vector
mod.vec <- c(m1,m2,m3,m4)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]), data=TF.nounk, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- r.squared(fit)$AIC
  Rm[i] <- 100*r.squared(fit)$Marginal # marginal R-squared
  Rc[i] <- 100*r.squared(fit)$Conditional # conditional R-squared
  print(i)
}

dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)

sumtable <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,10),round(Rm,4),round(Rc,4))
colnames(sumtable) <- c("model","k","LL","AICc","dAICc","wAICc","Rc","Rm")
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,5],decreasing=F),1:8]
summary.table

fit.sat <- lmer(as.formula(mod.vec[1]), data=TF.nounk, na.action=na.omit)
summary(fit.sat)


# redo with trialset as integer
TF.nounk$trialset <- as.integer(TF.nounk$trialset)

# model set
m1 <- "stime ~ det + trialset + det*trialset + (1|ID)"
m2 <- "stime ~ det + trialset + (1|ID)"
m3 <- "stime ~ det + (1|ID)"
m4 <- "stime ~ trialset + (1|ID)"
m5 <- "stime ~ 1 + (1|ID)"

## Make model vector
mod.vec <- c(m1,m2,m3,m4,m5)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]), data=TF.nounk, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- r.squared(fit)$AIC
  Rm[i] <- 100*r.squared(fit)$Marginal # marginal R-squared
  Rc[i] <- 100*r.squared(fit)$Conditional # conditional R-squared
  print(i)
}

dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)

sumtable <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,10),round(Rm,4),round(Rc,4))
colnames(sumtable) <- c("model","k","LL","AICc","dAICc","wAICc","Rc","Rm")
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,5],decreasing=F),1:8]
summary.table

fit.sat <- lmer(as.formula(mod.vec[1]), data=TF.nounk, na.action=na.omit)
summary(fit.sat)




# eat or not (binomial model)
# import data
setwd("~/.../deterrant trials/")
E.dat <- read.table("Binomial_data.csv", header=T, sep=",")
E.dat$trialset <- as.integer(E.dat$Trial_set)
E.dat$trial <- factor(E.dat$Trial)

## GLMM
# remove unknown sharks
E.nounk <- subset(E.dat, ID != "unknown")
E.nounk$ID <- factor(E.nounk$ID)

# model set
m1 <- "eat ~ det + trialset + (1|ID)"
m2 <- "eat ~ det + (1|ID)"
m3 <- "eat ~ trialset + (1|ID)"
m4 <- "eat ~ 1 + (1|ID)"

## Make model vector
mod.vec <- c(m1,m2,m3,m4)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glmer(as.formula(mod.vec[i]), family=binomial(link="logit"),data=E.nounk, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- r.squared(fit)$AIC
  Rm[i] <- 100*r.squared(fit)$Marginal # marginal R-squared
  Rc[i] <- 100*r.squared(fit)$Conditional # conditional R-squared
  print(i)
}

dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)

sumtable <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,10),round(Rm,4),round(Rc,4))
colnames(sumtable) <- c("model","k","LL","AICc","dAICc","wAICc","Rc","Rm")
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,5],decreasing=F),1:8]
summary.table

fit.sat <- glmer(as.formula(mod.vec[2]), family=binomial(link="logit"), data=E.nounk, na.action=na.omit)
summary(fit.sat)


## SS & controls only
E2.nounk <- subset(E.nounk, det == "SS" | det == "Control")

# model set
m1 <- "eat ~ det + trialset + (1|ID)"
m2 <- "eat ~ det + (1|ID)"
m3 <- "eat ~ trialset + (1|ID)"
m4 <- "eat ~ 1 + (1|ID)"

## Make model vector
mod.vec <- c(m1,m2,m3,m4)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glmer(as.formula(mod.vec[i]), family=binomial(link="logit"),data=E2.nounk, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- r.squared(fit)$AIC
  Rm[i] <- 100*r.squared(fit)$Marginal # marginal R-squared
  Rc[i] <- 100*r.squared(fit)$Conditional # conditional R-squared
  print(i)
}

dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)

sumtable <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,10),round(Rm,4),round(Rc,4))
colnames(sumtable) <- c("model","k","LL","AICc","dAICc","wAICc","Rc","Rm")
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,5],decreasing=F),1:8]
summary.table

fit.sat <- glmer(as.formula(mod.vec[2]), family=binomial(link="logit"), data=E2.nounk, na.action=na.omit)
summary(fit.sat)



## Number of approaches
setwd("~/Documents/Papers/Fish/Sharks/White sharks/deterrant trials/")
A.dat <- read.table("NoApp_data.csv", header=T, sep=",")
A.dat$trialset <- as.integer(A.dat$Trial_set)
A.dat$trial <- factor(A.dat$Trial)

hist(A.dat$appr)
range(A.dat$appr, na.rm=T)
A.dat$lappr <- log10(scale((A.dat$appr), center=F, scale=T))

## GLMM
# remove unknown sharks
A.nounk <- subset(A.dat, ID != "unknown")
A.nounk$ID <- factor(A.nounk$ID)
hist(A.nounk$appr)
A.nounk$lappr <- log10(scale((A.nounk$appr), center=F, scale=T))
hist((A.nounk$lappr))

# model set
m1 <- "lappr ~ det + trialset + det*trialset + (1|ID)"
m2 <- "lappr ~ det + trialset + (1|ID)"
m3 <- "lappr ~ det + (1|ID)"
m4 <- "lappr ~ trialset + (1|ID)"
m5 <- "lappr ~ 1 + (1|ID)"

## Make model vector
mod.vec <- c(m1,m2,m3,m4,m5)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]),data=A.nounk, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- r.squared(fit)$AIC
  Rm[i] <- 100*r.squared(fit)$Marginal # marginal R-squared
  Rc[i] <- 100*r.squared(fit)$Conditional # conditional R-squared
  print(i)
}

dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)

sumtable <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,10),round(Rm,4),round(Rc,4))
colnames(sumtable) <- c("model","k","LL","AICc","dAICc","wAICc","Rc","Rm")
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,5],decreasing=F),1:8]
summary.table

fit.sat <- lmer(as.formula(mod.vec[2]), data=A.nounk, na.action=na.omit)
summary(fit.sat)


## Number of approaches (SS & control only)
A2.dat <- subset(A.dat, det == "SS" | det == "Control")

## GLMM
# remove unknown sharks
A2.nounk <- subset(A2.dat, ID != "unknown")
A2.nounk$ID <- factor(A2.nounk$ID)
hist(A2.nounk$appr)
A2.nounk$lappr <- log10(scale((A2.nounk$appr), center=F, scale=T))
hist((A2.nounk$lappr))

# model set
#m1 <- "lappr ~ det + trialset + det*trialset + (1|ID)"
m1 <- "lappr ~ det + trialset + (1|ID)"
m2 <- "lappr ~ det + (1|ID)"
m3 <- "lappr ~ trialset + (1|ID)"
m4 <- "lappr ~ 1 + (1|ID)"

## Make model vector
mod.vec <- c(m1,m2,m3,m4)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]),data=A2.nounk, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- r.squared(fit)$AIC
  Rm[i] <- 100*r.squared(fit)$Marginal # marginal R-squared
  Rc[i] <- 100*r.squared(fit)$Conditional # conditional R-squared
  print(i)
}

dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)

sumtable <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,10),round(Rm,4),round(Rc,4))
colnames(sumtable) <- c("model","k","LL","AICc","dAICc","wAICc","Rc","Rm")
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,5],decreasing=F),1:8]
summary.table

fit.sat <- lmer(as.formula(mod.vec[1]), data=A2.nounk, na.action=na.omit)
summary(fit.sat)



###########################################################
# POWER ANALYSES
###########################################################

## DISTANCE
# MBAND VERSUS CONTROL
distMB.dat <- subset(dist.dat, det == "Mband" | det == "Control")
distMB.dat$trialset <- as.integer(distMB.dat$trialset)
distMB.dat$det <- factor(distMB.dat$det)

# remove zero distances
distMB.no0dist <- subset(distMB.dat, dist > 0)
distMB.no0dist$ldist <- log10(distMB.no0dist$dist)

# remove unknown sharks
distMB.nounk <- subset(distMB.no0dist, ID != "unknown")
distMB.nounk$ID <- factor(distMB.nounk$ID)
distMB.nounk$det2 <- ifelse(distMB.nounk$det == "Mband", "Treat", distMB.nounk$det)
n.size <- dim(distMB.nounk)[1]
barplot(xtabs(distMB.nounk$dist ~ distMB.nounk$det)/table(distMB.nounk$det))

distMB.C <- subset(distMB.nounk, det == "Control")
distMB.Cse <- sd(distMB.C$dist)/sqrt(dim(distMB.C)[1])
distMB.T <- subset(distMB.nounk, det == "Mband")
distMB.Tse <- sd(distMB.T$dist)/sqrt(dim(distMB.T)[1])
distMB.se <- c(distMB.Cse, distMB.Tse)

MB.plotdat <- data.frame(treatment=attr(table(distMB.nounk$det), 'names'),
                         distance=as.numeric(xtabs(distMB.nounk$dist ~ distMB.nounk$det)/table(distMB.nounk$det)),
                         se=distMB.se)
ggplot(MB.plotdat) +
  geom_bar( aes(x=treatment, y=distance), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=treatment, ymin=distance-se, ymax=distance+se), width=0.4, colour="orange", alpha=0.9, size=1.3)


# MLEASH VERSUS CONTROL
distML.dat <- subset(dist.dat, det == "Mleash" | det == "Control")
distML.dat$trialset <- as.integer(distML.dat$trialset)
distML.dat$det <- factor(distML.dat$det)

# remove zero distances
distML.no0dist <- subset(distML.dat, dist > 0)
distML.no0dist$ldist <- log10(distML.no0dist$dist)

# remove unknown sharks
distML.nounk <- subset(distML.no0dist, ID != "unknown")
distML.nounk$ID <- factor(distML.nounk$ID)
distML.nounk$det2 <- ifelse(distML.nounk$det == "Mleash", "Treat", distML.nounk$det)
n.size <- dim(distML.nounk)[1]
barplot(xtabs(distML.nounk$dist ~ distML.nounk$det)/table(distML.nounk$det))

distML.C <- subset(distML.nounk, det == "Control")
distML.Cse <- sd(distML.C$dist)/sqrt(dim(distML.C)[1])
distML.T <- subset(distML.nounk, det == "Mleash")
distML.Tse <- sd(distML.T$dist)/sqrt(dim(distML.T)[1])
distML.se <- c(distML.Cse, distML.Tse)

ML.plotdat <- data.frame(treatment=attr(table(distML.nounk$det), 'names'),
                         distance=as.numeric(xtabs(distML.nounk$dist ~ distML.nounk$det)/table(distML.nounk$det)),
                         se=distML.se)
ggplot(ML.plotdat) +
  geom_bar( aes(x=treatment, y=distance), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=treatment, ymin=distance-se, ymax=distance+se), width=0.4, colour="orange", alpha=0.9, size=1.3)


# RPELA VERSUS CONTROL
distRP.dat <- subset(dist.dat, det == "Rpela" | det == "Control")
distRP.dat$trialset <- as.integer(distRP.dat$trialset)
distRP.dat$det <- factor(distRP.dat$det)

# remove zero distances
distRP.no0dist <- subset(distRP.dat, dist > 0)
distRP.no0dist$ldist <- log10(distRP.no0dist$dist)

# remove unknown sharks
distRP.nounk <- subset(distRP.no0dist, ID != "unknown")
distRP.nounk$ID <- factor(distRP.nounk$ID)
distRP.nounk$det2 <- ifelse(distRP.nounk$det == "Rpela", "Treat", distRP.nounk$det)
n.size <- dim(distRP.nounk)[1]
barplot(xtabs(distRP.nounk$dist ~ distRP.nounk$det)/table(distRP.nounk$det))

distRP.C <- subset(distRP.nounk, det == "Control")
distRP.Cse <- sd(distRP.C$dist)/sqrt(dim(distRP.C)[1])
distRP.T <- subset(distRP.nounk, det == "Rpela")
distRP.Tse <- sd(distRP.T$dist)/sqrt(dim(distRP.T)[1])
distRP.se <- c(distRP.Cse, distRP.Tse)

RP.plotdat <- data.frame(treatment=attr(table(distRP.nounk$det), 'names'),
                         distance=as.numeric(xtabs(distRP.nounk$dist ~ distRP.nounk$det)/table(distRP.nounk$det)),
                         se=distRP.se)
ggplot(RP.plotdat) +
  geom_bar( aes(x=treatment, y=distance), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=treatment, ymin=distance-se, ymax=distance+se), width=0.4, colour="orange", alpha=0.9, size=1.3)


# SS VERSUS CONTROL
distSS.dat <- subset(dist.dat, det == "SS" | det == "Control")
distSS.dat$trialset <- as.integer(distSS.dat$trialset)
distSS.dat$det <- factor(distSS.dat$det)

# remove zero distances
distSS.no0dist <- subset(distSS.dat, dist > 0)
distSS.no0dist$ldist <- log10(distSS.no0dist$dist)

# remove unknown sharks
distSS.nounk <- subset(distSS.no0dist, ID != "unknown")
distSS.nounk$ID <- factor(distSS.nounk$ID)
distSS.nounk$det2 <- ifelse(distSS.nounk$det == "SS", "Treat", distSS.nounk$det)
n.size <- dim(distSS.nounk)[1]
barplot(xtabs(distSS.nounk$dist ~ distSS.nounk$det)/table(distSS.nounk$det))

distSS.C <- subset(distSS.nounk, det == "Control")
distSS.Cse <- sd(distSS.C$dist)/sqrt(dim(distSS.C)[1])
distSS.T <- subset(distSS.nounk, det == "SS")
distSS.Tse <- sd(distSS.T$dist)/sqrt(dim(distSS.T)[1])
distSS.se <- c(distSS.Cse, distSS.Tse)

SS.plotdat <- data.frame(treatment=attr(table(distSS.nounk$det), 'names'),
                         distance=as.numeric(xtabs(distSS.nounk$dist ~ distSS.nounk$det)/table(distSS.nounk$det)),
                         se=distSS.se)
ggplot(SS.plotdat) +
  geom_bar( aes(x=treatment, y=distance), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=treatment, ymin=distance-se, ymax=distance+se), width=0.4, colour="orange", alpha=0.9, size=1.3)


# WAX VERSUS CONTROL
distWX.dat <- subset(dist.dat, det == "Wax" | det == "Control")
distWX.dat$trialset <- as.integer(distWX.dat$trialset)
distWX.dat$det <- factor(distWX.dat$det)

# remove zero distances
distWX.no0dist <- subset(distWX.dat, dist > 0)
distWX.no0dist$ldist <- log10(distWX.no0dist$dist)

# remove unknown sharks
distWX.nounk <- subset(distWX.no0dist, ID != "unknown")
distWX.nounk$ID <- factor(distWX.nounk$ID)
distWX.nounk$det2 <- ifelse(distWX.nounk$det == "Wax", "Treat", distWX.nounk$det)
n.size <- dim(distWX.nounk)[1]
barplot(xtabs(distWX.nounk$dist ~ distWX.nounk$det)/table(distWX.nounk$det))

distWX.C <- subset(distWX.nounk, det == "Control")
distWX.Cse <- sd(distWX.C$dist)/sqrt(dim(distWX.C)[1])
distWX.T <- subset(distWX.nounk, det == "Wax")
distWX.Tse <- sd(distWX.T$dist)/sqrt(dim(distWX.T)[1])
distWX.se <- c(distWX.Cse, distWX.Tse)

WX.plotdat <- data.frame(treatment=attr(table(distWX.nounk$det), 'names'),
                         distance=as.numeric(xtabs(distWX.nounk$dist ~ distWX.nounk$det)/table(distWX.nounk$det)),
                         se=distWX.se)
ggplot(WX.plotdat) +
  geom_bar( aes(x=treatment, y=distance), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=treatment, ymin=distance-se, ymax=distance+se), width=0.4, colour="orange", alpha=0.9, size=1.3)


## choose which deterrant-control pair
distPAIR.nounk <- distMB.nounk
#distPAIR.nounk <- distML.nounk
#distPAIR.nounk <- distRP.nounk
#distPAIR.nounk <- distSS.nounk
#distPAIR.nounk <- distWX.nounk

n.size <- dim(distPAIR.nounk)[1]

# create percentage increase vector for treatment
inc.vec <- seq(1,1.5,0.05)
iter <- 1000
itdiv <- iter/100

# model set
m1 <- "ldistnew ~ det + trialset + det*trialset + (1|ID)"
m2 <- "ldistnew ~ det + trialset + (1|ID)"
m3 <- "ldistnew ~ det + (1|ID)"
m4 <- "ldistnew ~ trialset + (1|ID)"
m5 <- "ldistnew ~ 1 + (1|ID)"

## Make model vector
mod.vec <- c(m1,m2,m3,m4,m5)

## Define n.mod
n.mod <- length(mod.vec)

# storage for each increment
topmod.most <- mod3AICw.med <- mod3AICw.lo <- mod3AICw.up <- mod3Rc.med <- mod3Rc.lo <- mod3Rc.up <- mod3ER.med <- mod3ER.lo <- mod3ER.up <- rep(0,length(inc.vec))

# increment increase in response for treatment
for (n in 1:length(inc.vec)) {

  distPAIR.nounk$distnew <- ifelse(distPAIR.nounk$det2 == "Treat", distPAIR.nounk$dist * inc.vec[n], distPAIR.nounk$dist)
  distPAIR.nounk$ldistnew <- log10(distPAIR.nounk$distnew)
  
  # model outcome storage vectors
  topmod.vec <- topmodwAICc.vec <- topmodRc.vec <- mod3wAICc.vec <- mod3ER.vec <- mod3Rc.vec <- rep(0,iter)
  
  # resample dataset for iter iterations
  for (i in 1:iter) {
    resamp.sub <- sample(1:n.size, n.size, replace=TRUE)
    distPAIR.resamp <- distPAIR.nounk[resamp.sub,]
    
    # Model fitting and logLik output loop
    Modnum <- length(mod.vec)
    LL.vec <- SaveCount <- AICc.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
    mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
    mod.num <- seq(1,Modnum,1)
    
    for(m in 1:Modnum) {
      fit <- lmer(as.formula(mod.vec[m]), data=distPAIR.resamp, na.action=na.omit)
      assign(paste("fit",m,sep=""), fit)
      mod.list[[m]] <- fit
      LL.vec[m] <- as.numeric(logLik(fit))
      k.vec[m] <- attr(logLik(fit),"df")
      AICc.vec[m] <- r.squared(fit)$AIC
      Rm[m] <- 100*r.squared(fit)$Marginal # marginal R-squared
      Rc[m] <- 100*r.squared(fit)$Conditional # conditional R-squared
    }
    
    dAICc <- delta.IC(AICc.vec)
    wAICc <- weight.IC(dAICc)
    
    sumtable <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,10),round(Rm,4),round(Rc,4))
    colnames(sumtable) <- c("model","k","LL","AICc","dAICc","wAICc","Rc","Rm")
    row.names(sumtable) <- mod.vec
    summary.table <- sumtable[order(sumtable[,5],decreasing=F),1:8]
    summary.table
    
    topmod.vec[i] <- summary.table[1,1]
    topmodwAICc.vec[i] <- summary.table[1,6]
    topmodRc.vec[i] <- summary.table[1,7]
    
    mod3wAICc.vec[i] <- summary.table[which(summary.table[,1] == 3), 6]
    mod3Rc.vec[i] <- summary.table[which(summary.table[,1] == 3), 7]
    mod3ER.vec[i] <- mod3wAICc.vec[i] / summary.table[which(summary.table[,1] == 5), 6]
    
    if (i %% itdiv==0) print(i)
  }

  topmod.most.sub <- as.numeric(attr(rev(sort(table(topmod.vec))), 'names')[1])
  topmod.most[n] <- mod.vec[topmod.most.sub]
  
  mod3AICw.med[n] <- median(mod3wAICc.vec, na.rm=T) 
  mod3AICw.lo[n] <- quantile(mod3wAICc.vec, probs=0.025, na.rm=T) 
  mod3AICw.up[n] <- quantile(mod3wAICc.vec, probs=0.975, na.rm=T) 
  mod3ER.med[n] <- median(mod3ER.vec, na.rm=T) 
  mod3ER.lo[n] <- quantile(mod3ER.vec, probs=0.025, na.rm=T) 
  mod3ER.up[n] <- quantile(mod3ER.vec, probs=0.975, na.rm=T) 
  mod3Rc.med[n] <- median(mod3Rc.vec, na.rm=T) 
  mod3Rc.lo[n] <- quantile(mod3Rc.vec, probs=0.025, na.rm=T) 
  mod3Rc.up[n] <- quantile(mod3Rc.vec, probs=0.975, na.rm=T)
  
  print("-----------------------")
  print(paste("increment number =", n))
  print("-----------------------")
}

par(mfrow=c(1,3))
plot(100*(inc.vec-1), mod3AICw.med, type="l", lty=2, lwd=2, xlab="% incr in treatment", ylab="model 3 wAICc", ylim=c(min(mod3AICw.lo), max(mod3AICw.up)))
lines(100*(inc.vec-1), mod3AICw.lo, lty=1, lwd=1, col="red")
lines(100*(inc.vec-1), mod3AICw.up, lty=1, lwd=1, col="red")

plot(100*(inc.vec-1), log10(mod3ER.med), type="l", lty=2, lwd=2, xlab="% incr in treatment", ylab="model 3 log10(ER)", ylim=c(log10(min(mod3ER.lo)), log10(max(mod3ER.up))))
lines(100*(inc.vec-1), log10(mod3ER.lo), lty=1, lwd=1, col="red")
lines(100*(inc.vec-1), log10(mod3ER.up), lty=1, lwd=1, col="red")

plot(100*(inc.vec-1), mod3Rc.med, type="l", lty=2, lwd=2, xlab="% incr in treatment", ylab="model 3 Rc", ylim=c(min(mod3Rc.lo), max(mod3Rc.up)))
lines(100*(inc.vec-1), mod3Rc.lo, lty=1, lwd=1, col="red")
lines(100*(inc.vec-1), mod3Rc.up, lty=1, lwd=1, col="red")
par(mfrow=c(1,1))

print(data.frame(inc.vec, topmod.most))

results.out <- data.frame(inc.vec, mod3AICw.med, mod3AICw.up, mod3AICw.lo, mod3ER.med, mod3ER.up, mod3ER.lo, mod3Rc.med, mod3Rc.up, mod3Rc.lo)
common.model.out <- data.frame(inc.vec, topmod.most)


##########################################################################
## TIME TO FEEDING
# MBAND VERSUS CONTROL
tfMB.dat <- subset(TF.dat, det == "Mband" | det == "Control")
tfMB.dat$trialset <- as.integer(tfMB.dat$trialset)
tfMB.dat$det <- factor(tfMB.dat$det)

# remove unknown sharks
tfMB.nounk <- subset(tfMB.dat, ID != "unknown")
tfMB.nounk$ID <- factor(tfMB.nounk$ID)
tfMB.nounk$stime <- as.numeric(scale(log10(tfMB.nounk$time), center=T, scale=T))
tfMB.nounk$det2 <- ifelse(tfMB.nounk$det == "Mband", "Treat", "Control")
n.size <- dim(tfMB.nounk)[1]
barplot(xtabs(tfMB.nounk$time ~ tfMB.nounk$det)/table(tfMB.nounk$det))

tfMB.C <- subset(tfMB.nounk, det == "Control")
tfMB.Cse <- sd(tfMB.C$time)/sqrt(dim(tfMB.C)[1])
tfMB.T <- subset(tfMB.nounk, det == "Mband")
tfMB.Tse <- sd(tfMB.T$time)/sqrt(dim(tfMB.T)[1])
tfMB.se <- c(tfMB.Cse, tfMB.Tse)

MB.plotdat <- data.frame(treatment=attr(table(tfMB.nounk$det), 'names'),
                         TimeToFeeding=as.numeric(xtabs(tfMB.nounk$time ~ tfMB.nounk$det)/table(tfMB.nounk$det)),
                         se=tfMB.se)
ggplot(MB.plotdat) +
  geom_bar( aes(x=treatment, y=TimeToFeeding), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=treatment, ymin=TimeToFeeding-se, ymax=TimeToFeeding+se), width=0.4, colour="orange", alpha=0.9, size=1.3)


# MLEASH VERSUS CONTROL
tfML.dat <- subset(TF.dat, det == "Mleash" | det == "Control")
tfML.dat$trialset <- as.integer(tfML.dat$trialset)
tfML.dat$det <- factor(tfML.dat$det)

# remove unknown sharks
tfML.nounk <- na.omit(subset(tfML.dat, ID != "unknown"))
tfML.nounk$ID <- factor(tfML.nounk$ID)
tfML.nounk$stime <- as.numeric(scale(log10(tfML.nounk$time), center=T, scale=T))
tfML.nounk$det2 <- ifelse(tfML.nounk$det == "Mleash", "Treat", "Control")
n.size <- dim(tfML.nounk)[1]
barplot(xtabs(tfML.nounk$time ~ tfML.nounk$det)/table(tfML.nounk$det))

tfML.C <- subset(tfML.nounk, det == "Control")
tfML.Cse <- sd(tfML.C$time)/sqrt(dim(tfML.C)[1])
tfML.T <- subset(tfML.nounk, det == "Mleash")
tfML.Tse <- sd(tfML.T$time)/sqrt(dim(tfML.T)[1])
tfML.se <- c(tfML.Cse, tfML.Tse)

ML.plotdat <- data.frame(treatment=attr(table(tfML.nounk$det), 'names'),
                         TimeToFeeding=as.numeric(xtabs(tfML.nounk$time ~ tfML.nounk$det)/table(tfML.nounk$det)),
                         se=tfML.se)
ggplot(ML.plotdat) +
  geom_bar( aes(x=treatment, y=TimeToFeeding), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=treatment, ymin=TimeToFeeding-se, ymax=TimeToFeeding+se), width=0.4, colour="orange", alpha=0.9, size=1.3)


# RPELA VERSUS CONTROL
tfRP.dat <- subset(TF.dat, det == "Rpela" | det == "Control")
tfRP.dat$trialset <- as.integer(tfRP.dat$trialset)
tfRP.dat$det <- factor(tfRP.dat$det)

# remove unknown sharks
tfRP.nounk <- na.omit(subset(tfRP.dat, ID != "unknown"))
tfRP.nounk$ID <- factor(tfRP.nounk$ID)
tfRP.nounk$stime <- as.numeric(scale(log10(tfRP.nounk$time), center=T, scale=T))
tfRP.nounk$det2 <- ifelse(tfRP.nounk$det == "Rpela", "Treat", "Control")
n.size <- dim(tfRP.nounk)[1]
barplot(xtabs(tfRP.nounk$time ~ tfRP.nounk$det)/table(tfRP.nounk$det))

tfRP.C <- subset(tfRP.nounk, det == "Control")
tfRP.Cse <- sd(tfRP.C$time)/sqrt(dim(tfRP.C)[1])
tfRP.T <- subset(tfRP.nounk, det == "Rpela")
tfRP.Tse <- sd(tfRP.T$time)/sqrt(dim(tfRP.T)[1])
tfRP.se <- c(tfRP.Cse, tfRP.Tse)

RP.plotdat <- data.frame(treatment=attr(table(tfRP.nounk$det), 'names'),
                         TimeToFeeding=as.numeric(xtabs(tfRP.nounk$time ~ tfRP.nounk$det)/table(tfRP.nounk$det)),
                         se=tfRP.se)
ggplot(RP.plotdat) +
  geom_bar( aes(x=treatment, y=TimeToFeeding), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=treatment, ymin=TimeToFeeding-se, ymax=TimeToFeeding+se), width=0.4, colour="orange", alpha=0.9, size=1.3)


# SS VERSUS CONTROL
tfSS.dat <- subset(TF.dat, det == "SS" | det == "Control")
tfSS.dat$trialset <- as.integer(tfSS.dat$trialset)
tfSS.dat$det <- factor(tfSS.dat$det)

# remove unknown sharks
tfSS.nounk <- na.omit(subset(tfSS.dat, ID != "unknown"))
tfSS.nounk$ID <- factor(tfSS.nounk$ID)
tfSS.nounk$stime <- as.numeric(scale(log10(tfSS.nounk$time), center=T, scale=T))
tfSS.nounk$det2 <- ifelse(tfSS.nounk$det == "SS", "Treat", "Control")
n.size <- dim(tfSS.nounk)[1]
barplot(xtabs(tfSS.nounk$time ~ tfSS.nounk$det)/table(tfSS.nounk$det))

tfSS.C <- subset(tfSS.nounk, det == "Control")
tfSS.Cse <- sd(tfSS.C$time)/sqrt(dim(tfSS.C)[1])
tfSS.T <- subset(tfSS.nounk, det == "SS")
tfSS.Tse <- sd(tfSS.T$time)/sqrt(dim(tfSS.T)[1])
tfSS.se <- c(tfSS.Cse, tfSS.Tse)

SS.plotdat <- data.frame(treatment=attr(table(tfSS.nounk$det), 'names'),
                         TimeToFeeding=as.numeric(xtabs(tfSS.nounk$time ~ tfSS.nounk$det)/table(tfSS.nounk$det)),
                         se=tfSS.se)
ggplot(SS.plotdat) +
  geom_bar( aes(x=treatment, y=TimeToFeeding), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=treatment, ymin=TimeToFeeding-se, ymax=TimeToFeeding+se), width=0.4, colour="orange", alpha=0.9, size=1.3)


# WAX VERSUS CONTROL
tfWX.dat <- subset(TF.dat, det == "Wax" | det == "Control")
tfWX.dat$trialset <- as.integer(tfWX.dat$trialset)
tfWX.dat$det <- factor(tfWX.dat$det)

# remove unknown sharks
tfWX.nounk <- na.omit(subset(tfWX.dat, ID != "unknown"))
tfWX.nounk$ID <- factor(tfWX.nounk$ID)
tfWX.nounk$stime <- as.numeric(scale(log10(tfWX.nounk$time), center=T, scale=T))
tfWX.nounk$det2 <- ifelse(tfWX.nounk$det == "Wax", "Treat", "Control")
n.size <- dim(tfWX.nounk)[1]
barplot(xtabs(tfWX.nounk$time ~ tfWX.nounk$det)/table(tfWX.nounk$det))

tfWX.C <- subset(tfWX.nounk, det == "Control")
tfWX.Cse <- sd(tfWX.C$time)/sqrt(dim(tfWX.C)[1])
tfWX.T <- subset(tfWX.nounk, det == "Wax")
tfWX.Tse <- sd(tfWX.T$time)/sqrt(dim(tfWX.T)[1])
tfWX.se <- c(tfWX.Cse, tfWX.Tse)

WX.plotdat <- data.frame(treatment=attr(table(tfWX.nounk$det), 'names'),
                         TimeToFeeding=as.numeric(xtabs(tfWX.nounk$time ~ tfWX.nounk$det)/table(tfWX.nounk$det)),
                         se=tfWX.se)
ggplot(WX.plotdat) +
  geom_bar( aes(x=treatment, y=TimeToFeeding), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=treatment, ymin=TimeToFeeding-se, ymax=TimeToFeeding+se), width=0.4, colour="orange", alpha=0.9, size=1.3)


## choose which deterrant-control pair
tfPAIR.nounk <- tfMB.nounk
#tfPAIR.nounk <- tfML.nounk
#tfPAIR.nounk <- tfRP.nounk
#tfPAIR.nounk <- tfSS.nounk
#tfPAIR.nounk <- tfWX.nounk

n.size <- dim(tfPAIR.nounk)[1]

# create percentage increase vector for treatment
inc.vec <- seq(1,1.5,0.05)
iter <- 1000
itdiv <- iter/100

# model set
m1 <- "stimenew ~ det + trialset + det*trialset + (1|ID)"
m2 <- "stimenew ~ det + trialset + (1|ID)"
m3 <- "stimenew ~ det + (1|ID)"
m4 <- "stimenew ~ trialset + (1|ID)"
m5 <- "stimenew ~ 1 + (1|ID)"

## Make model vector
mod.vec <- c(m1,m2,m3,m4,m5)

## Define n.mod
n.mod <- length(mod.vec)

# storage for each increment
topmod.most <- mod3AICw.med <- mod3AICw.lo <- mod3AICw.up <- mod3Rc.med <- mod3Rc.lo <- mod3Rc.up <- mod3ER.med <- mod3ER.lo <- mod3ER.up <- rep(0,length(inc.vec))

# increment increase in response for treatment
for (n in 1:length(inc.vec)) {
  
  tfPAIR.nounk$timenew <- ifelse(tfPAIR.nounk$det2 == "Treat", tfPAIR.nounk$time * inc.vec[n], tfPAIR.nounk$time)
  tfPAIR.nounk$stimenew <- as.numeric(scale(log10(tfPAIR.nounk$timenew), center=T, scale=T))
  
  # model outcome storage vectors
  topmod.vec <- topmodwAICc.vec <- topmodRc.vec <- mod3wAICc.vec <- mod3ER.vec <- mod3Rc.vec <- rep(0,iter)
  
  # resample dataset for iter iterations
  for (i in 1:iter) {
    resamp.sub <- sample(1:n.size, n.size, replace=TRUE)
    tfPAIR.resamp <- tfPAIR.nounk[resamp.sub,]
    
    # Model fitting and logLik output loop
    Modnum <- length(mod.vec)
    LL.vec <- SaveCount <- AICc.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
    mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
    mod.num <- seq(1,Modnum,1)
    
    for(m in 1:Modnum) {
      fit <- lmer(as.formula(mod.vec[m]), data=tfPAIR.resamp, na.action=na.omit)
      assign(paste("fit",m,sep=""), fit)
      mod.list[[m]] <- fit
      LL.vec[m] <- as.numeric(logLik(fit))
      k.vec[m] <- attr(logLik(fit),"df")
      AICc.vec[m] <- r.squared(fit)$AIC
      Rm[m] <- 100*r.squared(fit)$Marginal # marginal R-squared
      Rc[m] <- 100*r.squared(fit)$Conditional # conditional R-squared
    }
    
    dAICc <- delta.IC(AICc.vec)
    wAICc <- weight.IC(dAICc)
    
    sumtable <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,10),round(Rm,4),round(Rc,4))
    colnames(sumtable) <- c("model","k","LL","AICc","dAICc","wAICc","Rc","Rm")
    row.names(sumtable) <- mod.vec
    summary.table <- sumtable[order(sumtable[,5],decreasing=F),1:8]
    summary.table
    
    topmod.vec[i] <- summary.table[1,1]
    topmodwAICc.vec[i] <- summary.table[1,6]
    topmodRc.vec[i] <- summary.table[1,7]
    
    mod3wAICc.vec[i] <- summary.table[which(summary.table[,1] == 3), 6]
    mod3Rc.vec[i] <- summary.table[which(summary.table[,1] == 3), 7]
    mod3ER.vec[i] <- mod3wAICc.vec[i] / summary.table[which(summary.table[,1] == 5), 6]
    
    if (i %% itdiv==0) print(i)
  }
  
  topmod.most.sub <- as.numeric(attr(rev(sort(table(topmod.vec))), 'names')[1])
  topmod.most[n] <- mod.vec[topmod.most.sub]
  
  mod3AICw.med[n] <- median(mod3wAICc.vec, na.rm=T) 
  mod3AICw.lo[n] <- quantile(mod3wAICc.vec, probs=0.025, na.rm=T) 
  mod3AICw.up[n] <- quantile(mod3wAICc.vec, probs=0.975, na.rm=T) 
  mod3ER.med[n] <- median(mod3ER.vec, na.rm=T) 
  mod3ER.lo[n] <- quantile(mod3ER.vec, probs=0.025, na.rm=T) 
  mod3ER.up[n] <- quantile(mod3ER.vec, probs=0.975, na.rm=T) 
  mod3Rc.med[n] <- median(mod3Rc.vec, na.rm=T) 
  mod3Rc.lo[n] <- quantile(mod3Rc.vec, probs=0.025, na.rm=T) 
  mod3Rc.up[n] <- quantile(mod3Rc.vec, probs=0.975, na.rm=T)
  
  print("-----------------------")
  print(paste("increment number =", n))
  print("-----------------------")
}

par(mfrow=c(1,3))
plot(100*(inc.vec-1), mod3AICw.med, type="l", lty=2, lwd=2, xlab="% incr in treatment", ylab="model 3 wAICc", ylim=c(min(mod3AICw.lo), max(mod3AICw.up)))
lines(100*(inc.vec-1), mod3AICw.lo, lty=1, lwd=1, col="red")
lines(100*(inc.vec-1), mod3AICw.up, lty=1, lwd=1, col="red")

plot(100*(inc.vec-1), log10(mod3ER.med), type="l", lty=2, lwd=2, xlab="% incr in treatment", ylab="model 3 log10(ER)", ylim=c(log10(min(mod3ER.lo)), log10(max(mod3ER.up))))
lines(100*(inc.vec-1), log10(mod3ER.lo), lty=1, lwd=1, col="red")
lines(100*(inc.vec-1), log10(mod3ER.up), lty=1, lwd=1, col="red")

plot(100*(inc.vec-1), mod3Rc.med, type="l", lty=2, lwd=2, xlab="% incr in treatment", ylab="model 3 Rc", ylim=c(min(mod3Rc.lo), max(mod3Rc.up)))
lines(100*(inc.vec-1), mod3Rc.lo, lty=1, lwd=1, col="red")
lines(100*(inc.vec-1), mod3Rc.up, lty=1, lwd=1, col="red")
par(mfrow=c(1,1))

print(data.frame(inc.vec, topmod.most))

results.out <- data.frame(inc.vec, mod3AICw.med, mod3AICw.up, mod3AICw.lo, mod3ER.med, mod3ER.up, mod3ER.lo, mod3Rc.med, mod3Rc.up, mod3Rc.lo)
common.model.out <- data.frame(inc.vec, topmod.most)



##########################################################################
## NUMBER OF APPROACHES
# MBAND VERSUS CONTROL
appMB.dat <- subset(A.dat, det == "Mband" | det == "Control")
appMB.dat$trialset <- as.integer(appMB.dat$trialset)
appMB.dat$det <- factor(appMB.dat$det)

# remove unknown sharks
appMB.nounk <- subset(appMB.dat, ID != "unknown")
appMB.nounk$ID <- factor(appMB.nounk$ID)
appMB.nounk$lappr <- log10(scale((appMB.nounk$appr), center=F, scale=T))
appMB.nounk$det2 <- ifelse(appMB.nounk$det == "Mband", "Treat", "Control")
n.size <- dim(appMB.nounk)[1]
barplot(xtabs(appMB.nounk$app ~ appMB.nounk$det)/table(appMB.nounk$det))

appMB.C <- subset(appMB.nounk, det == "Control")
appMB.Cse <- sd(appMB.C$appr)/sqrt(dim(appMB.C)[1])
appMB.T <- subset(appMB.nounk, det == "Mband")
appMB.Tse <- sd(appMB.T$appr)/sqrt(dim(appMB.T)[1])
appMB.se <- c(appMB.Cse, appMB.Tse)

MB.plotdat <- data.frame(treatment=attr(table(appMB.nounk$det), 'names'),
                         NoApproaches=as.numeric(xtabs(appMB.nounk$appr ~ appMB.nounk$det)/table(appMB.nounk$det)),
                         se=appMB.se)
ggplot(MB.plotdat) +
  geom_bar( aes(x=treatment, y=NoApproaches), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=treatment, ymin=NoApproaches-se, ymax=NoApproaches+se), width=0.4, colour="orange", alpha=0.9, size=1.3)


# MLEASH VERSUS CONTROL
appML.dat <- subset(A.dat, det == "Mleash" | det == "Control")
appML.dat$trialset <- as.integer(appML.dat$trialset)
appML.dat$det <- factor(appML.dat$det)

# remove unknown sharks
appML.nounk <- subset(appML.dat, ID != "unknown")
appML.nounk$ID <- factor(appML.nounk$ID)
appML.nounk$lappr <- log10(scale((appML.nounk$appr), center=F, scale=T))
appML.nounk$det2 <- ifelse(appML.nounk$det == "Mleash", "Treat", "Control")
n.size <- dim(appML.nounk)[1]
barplot(xtabs(appML.nounk$app ~ appML.nounk$det)/table(appML.nounk$det))

appML.C <- subset(appML.nounk, det == "Control")
appML.Cse <- sd(appML.C$appr)/sqrt(dim(appML.C)[1])
appML.T <- subset(appML.nounk, det == "Mleash")
appML.Tse <- sd(appML.T$appr)/sqrt(dim(appML.T)[1])
appML.se <- c(appML.Cse, appML.Tse)

ML.plotdat <- data.frame(treatment=attr(table(appML.nounk$det), 'names'),
                         NoApproaches=as.numeric(xtabs(appML.nounk$appr ~ appML.nounk$det)/table(appML.nounk$det)),
                         se=appML.se)
ggplot(ML.plotdat) +
  geom_bar( aes(x=treatment, y=NoApproaches), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=treatment, ymin=NoApproaches-se, ymax=NoApproaches+se), width=0.4, colour="orange", alpha=0.9, size=1.3)


# RPELA VERSUS CONTROL
appRP.dat <- subset(A.dat, det == "Rpela" | det == "Control")
appRP.dat$trialset <- as.integer(appRP.dat$trialset)
appRP.dat$det <- factor(appRP.dat$det)

# remove unknown sharks
appRP.nounk <- subset(appRP.dat, ID != "unknown")
appRP.nounk$ID <- factor(appRP.nounk$ID)
appRP.nounk$lappr <- log10(scale((appRP.nounk$appr), center=F, scale=T))
appRP.nounk$det2 <- ifelse(appRP.nounk$det == "Rpela", "Treat", "Control")
n.size <- dim(appRP.nounk)[1]
barplot(xtabs(appRP.nounk$app ~ appRP.nounk$det)/table(appRP.nounk$det))

appRP.C <- subset(appRP.nounk, det == "Control")
appRP.Cse <- sd(appRP.C$appr)/sqrt(dim(appRP.C)[1])
appRP.T <- subset(appRP.nounk, det == "Rpela")
appRP.Tse <- sd(appRP.T$appr)/sqrt(dim(appRP.T)[1])
appRP.se <- c(appRP.Cse, appRP.Tse)

RP.plotdat <- data.frame(treatment=attr(table(appRP.nounk$det), 'names'),
                         NoApproaches=as.numeric(xtabs(appRP.nounk$appr ~ appRP.nounk$det)/table(appRP.nounk$det)),
                         se=appRP.se)
ggplot(RP.plotdat) +
  geom_bar( aes(x=treatment, y=NoApproaches), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=treatment, ymin=NoApproaches-se, ymax=NoApproaches+se), width=0.4, colour="orange", alpha=0.9, size=1.3)


# SS VERSUS CONTROL
appSS.dat <- subset(A.dat, det == "SS" | det == "Control")
appSS.dat$trialset <- as.integer(appSS.dat$trialset)
appSS.dat$det <- factor(appSS.dat$det)

# remove unknown sharks
appSS.nounk <- subset(appSS.dat, ID != "unknown")
appSS.nounk$ID <- factor(appSS.nounk$ID)
appSS.nounk$lappr <- log10(scale((appSS.nounk$appr), center=F, scale=T))
appSS.nounk$det2 <- ifelse(appSS.nounk$det == "SS", "Treat", "Control")
n.size <- dim(appSS.nounk)[1]
barplot(xtabs(appSS.nounk$app ~ appSS.nounk$det)/table(appSS.nounk$det))

appSS.C <- subset(appSS.nounk, det == "Control")
appSS.Cse <- sd(appSS.C$appr)/sqrt(dim(appSS.C)[1])
appSS.T <- subset(appSS.nounk, det == "SS")
appSS.Tse <- sd(appSS.T$appr)/sqrt(dim(appSS.T)[1])
appSS.se <- c(appSS.Cse, appSS.Tse)

SS.plotdat <- data.frame(treatment=attr(table(appSS.nounk$det), 'names'),
                         NoApproaches=as.numeric(xtabs(appSS.nounk$appr ~ appSS.nounk$det)/table(appSS.nounk$det)),
                         se=appSS.se)
ggplot(SS.plotdat) +
  geom_bar( aes(x=treatment, y=NoApproaches), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=treatment, ymin=NoApproaches-se, ymax=NoApproaches+se), width=0.4, colour="orange", alpha=0.9, size=1.3)


# WAX VERSUS CONTROL
appWX.dat <- subset(A.dat, det == "Wax" | det == "Control")
appWX.dat$trialset <- as.integer(appWX.dat$trialset)
appWX.dat$det <- factor(appWX.dat$det)

# remove unknown sharks
appWX.nounk <- subset(appWX.dat, ID != "unknown")
appWX.nounk$ID <- factor(appWX.nounk$ID)
appWX.nounk$lappr <- log10(scale((appWX.nounk$appr), center=F, scale=T))
appWX.nounk$det2 <- ifelse(appWX.nounk$det == "Wax", "Treat", "Control")
n.size <- dim(appWX.nounk)[1]
barplot(xtabs(appWX.nounk$app ~ appWX.nounk$det)/table(appWX.nounk$det))

appWX.C <- subset(appWX.nounk, det == "Control")
appWX.Cse <- sd(appWX.C$appr)/sqrt(dim(appWX.C)[1])
appWX.T <- subset(appWX.nounk, det == "Wax")
appWX.Tse <- sd(appWX.T$appr)/sqrt(dim(appWX.T)[1])
appWX.se <- c(appWX.Cse, appWX.Tse)

WX.plotdat <- data.frame(treatment=attr(table(appWX.nounk$det), 'names'),
                         NoApproaches=as.numeric(xtabs(appWX.nounk$appr ~ appWX.nounk$det)/table(appWX.nounk$det)),
                         se=appWX.se)
ggplot(WX.plotdat) +
  geom_bar( aes(x=treatment, y=NoApproaches), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=treatment, ymin=NoApproaches-se, ymax=NoApproaches+se), width=0.4, colour="orange", alpha=0.9, size=1.3)


## choose which deterrant-control pair
appPAIR.nounk <- appMB.nounk
#appPAIR.nounk <- appML.nounk
#appPAIR.nounk <- appRP.nounk
#appPAIR.nounk <- appSS.nounk
#appPAIR.nounk <- appWX.nounk

n.size <- dim(appPAIR.nounk)[1]

# create percentage increase vector for treatment
inc.vec <- seq(1,1.5,0.05)
iter <- 1000
itdiv <- iter/100

# model set
m1 <- "lapprnew ~ det2 + trialset + det2*trialset + (1|ID)"
m2 <- "lapprnew ~ det2 + trialset + (1|ID)"
m3 <- "lapprnew ~ det2 + (1|ID)"
m4 <- "lapprnew ~ trialset + (1|ID)"
m5 <- "lapprnew ~ 1 + (1|ID)"

## Make model vector
mod.vec <- c(m1,m2,m3,m4,m5)

## Define n.mod
n.mod <- length(mod.vec)

# storage for each increment
topmod.most <- mod3AICw.med <- mod3AICw.lo <- mod3AICw.up <- mod3Rc.med <- mod3Rc.lo <- mod3Rc.up <- mod3ER.med <- mod3ER.lo <- mod3ER.up <- rep(0,length(inc.vec))

# increment increase in response for treatment
for (n in 1:length(inc.vec)) {
  
  appPAIR.nounk$apprnew <- ifelse(appPAIR.nounk$det2 == "Treat", round((appPAIR.nounk$appr * inc.vec[n]), 0), appPAIR.nounk$appr)
  appPAIR.nounk$lapprnew <- as.numeric(log10(scale(appPAIR.nounk$apprnew, center=F, scale=T)))
  
  # model outcome storage vectors
  topmod.vec <- topmodwAICc.vec <- topmodRc.vec <- mod3wAICc.vec <- mod3ER.vec <- mod3Rc.vec <- rep(0,iter)
  
  # resample dataset for iter iterations
  for (i in 1:iter) {
    resamp.sub <- sample(1:n.size, n.size, replace=TRUE)
    appPAIR.resamp <- appPAIR.nounk[resamp.sub,]
    
    # Model fitting and logLik output loop
    Modnum <- length(mod.vec)
    LL.vec <- SaveCount <- AICc.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
    mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
    mod.num <- seq(1,Modnum,1)
    
    for(m in 1:Modnum) {
      fit <- lmer(as.formula(mod.vec[m]), data=appPAIR.resamp, na.action=na.omit)
      assign(paste("fit",m,sep=""), fit)
      mod.list[[m]] <- fit
      LL.vec[m] <- as.numeric(logLik(fit))
      k.vec[m] <- attr(logLik(fit),"df")
      AICc.vec[m] <- r.squared(fit)$AIC
      Rm[m] <- 100*r.squared(fit)$Marginal # marginal R-squared
      Rc[m] <- 100*r.squared(fit)$Conditional # conditional R-squared
    }
    
    dAICc <- delta.IC(AICc.vec)
    wAICc <- weight.IC(dAICc)
    
    sumtable <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,10),round(Rm,4),round(Rc,4))
    colnames(sumtable) <- c("model","k","LL","AICc","dAICc","wAICc","Rc","Rm")
    row.names(sumtable) <- mod.vec
    summary.table <- sumtable[order(sumtable[,5],decreasing=F),1:8]
    summary.table
    
    topmod.vec[i] <- summary.table[1,1]
    topmodwAICc.vec[i] <- summary.table[1,6]
    topmodRc.vec[i] <- summary.table[1,7]
    
    mod3wAICc.vec[i] <- summary.table[which(summary.table[,1] == 3), 6]
    mod3Rc.vec[i] <- summary.table[which(summary.table[,1] == 3), 7]
    mod3ER.vec[i] <- mod3wAICc.vec[i] / summary.table[which(summary.table[,1] == 5), 6]
    
    if (i %% itdiv==0) print(i)
  }
  
  topmod.most.sub <- as.numeric(attr(rev(sort(table(topmod.vec))), 'names')[1])
  topmod.most[n] <- mod.vec[topmod.most.sub]
  
  mod3AICw.med[n] <- median(mod3wAICc.vec, na.rm=T) 
  mod3AICw.lo[n] <- quantile(mod3wAICc.vec, probs=0.025, na.rm=T) 
  mod3AICw.up[n] <- quantile(mod3wAICc.vec, probs=0.975, na.rm=T) 
  mod3ER.med[n] <- median(mod3ER.vec, na.rm=T) 
  mod3ER.lo[n] <- quantile(mod3ER.vec, probs=0.025, na.rm=T) 
  mod3ER.up[n] <- quantile(mod3ER.vec, probs=0.975, na.rm=T) 
  mod3Rc.med[n] <- median(mod3Rc.vec, na.rm=T) 
  mod3Rc.lo[n] <- quantile(mod3Rc.vec, probs=0.025, na.rm=T) 
  mod3Rc.up[n] <- quantile(mod3Rc.vec, probs=0.975, na.rm=T)
  
  print("-----------------------")
  print(paste("increment number =", n))
  print("-----------------------")
}

par(mfrow=c(1,3))
plot(100*(inc.vec-1), mod3AICw.med, type="l", lty=2, lwd=2, xlab="% incr in treatment", ylab="model 3 wAICc", ylim=c(min(mod3AICw.lo), max(mod3AICw.up)))
lines(100*(inc.vec-1), mod3AICw.lo, lty=1, lwd=1, col="red")
lines(100*(inc.vec-1), mod3AICw.up, lty=1, lwd=1, col="red")

plot(100*(inc.vec-1), log10(mod3ER.med), type="l", lty=2, lwd=2, xlab="% incr in treatment", ylab="model 3 log10(ER)", ylim=c(log10(min(mod3ER.lo)), log10(max(mod3ER.up))))
lines(100*(inc.vec-1), log10(mod3ER.lo), lty=1, lwd=1, col="red")
lines(100*(inc.vec-1), log10(mod3ER.up), lty=1, lwd=1, col="red")

plot(100*(inc.vec-1), mod3Rc.med, type="l", lty=2, lwd=2, xlab="% incr in treatment", ylab="model 3 Rc", ylim=c(min(mod3Rc.lo), max(mod3Rc.up)))
lines(100*(inc.vec-1), mod3Rc.lo, lty=1, lwd=1, col="red")
lines(100*(inc.vec-1), mod3Rc.up, lty=1, lwd=1, col="red")
par(mfrow=c(1,1))

print(data.frame(inc.vec, topmod.most))

results.out <- data.frame(inc.vec, mod3AICw.med, mod3AICw.up, mod3AICw.lo, mod3ER.med, mod3ER.up, mod3ER.lo, mod3Rc.med, mod3Rc.up, mod3Rc.lo)
common.model.out <- data.frame(inc.vec, topmod.most)



##########################################################################
## PROBABILITY OF TAKING BAIT (EAT)
library(tcltk) # status bars

# MBAND VERSUS CONTROL
eatMB.dat <- subset(E.dat, det == "Mband" | det == "Control")
eatMB.dat$trialset <- as.integer(eatMB.dat$trialset)
eatMB.dat$det <- factor(eatMB.dat$det)

# remove unknown sharks
eatMB.nounk <- subset(eatMB.dat, ID != "unknown")
eatMB.nounk$ID <- factor(eatMB.nounk$ID)
eatMB.nounk$det2 <- ifelse(eatMB.nounk$det == "Mband", "Treat", "Control")
n.size <- dim(eatMB.nounk)[1]
barplot(xtabs(eatMB.nounk$eat ~ eatMB.nounk$det)/table(eatMB.nounk$det2))


# MLEASH VERSUS CONTROL
eatML.dat <- subset(E.dat, det == "Mleash" | det == "Control")
eatML.dat$trialset <- as.integer(eatML.dat$trialset)
eatML.dat$det <- factor(eatML.dat$det)

# remove unknown sharks
eatML.nounk <- subset(eatML.dat, ID != "unknown")
eatML.nounk$ID <- factor(eatML.nounk$ID)
eatML.nounk$det2 <- ifelse(eatML.nounk$det == "Mleash", "Treat", "Control")
n.size <- dim(eatML.nounk)[1]
barplot(xtabs(eatML.nounk$eat ~ eatML.nounk$det)/table(eatML.nounk$det2))


# RPELA VERSUS CONTROL
eatRP.dat <- subset(E.dat, det == "Rpela" | det == "Control")
eatRP.dat$trialset <- as.integer(eatRP.dat$trialset)
eatRP.dat$det <- factor(eatRP.dat$det)

# remove unknown sharks
eatRP.nounk <- subset(eatRP.dat, ID != "unknown")
eatRP.nounk$ID <- factor(eatRP.nounk$ID)
eatRP.nounk$det2 <- ifelse(eatRP.nounk$det == "Rpela", "Treat", "Control")
n.size <- dim(eatRP.nounk)[1]
barplot(xtabs(eatRP.nounk$eat ~ eatRP.nounk$det)/table(eatRP.nounk$det2))


# SS VERSUS CONTROL
eatSS.dat <- subset(E.dat, det == "SS" | det == "Control")
eatSS.dat$trialset <- as.integer(eatSS.dat$trialset)
eatSS.dat$det <- factor(eatSS.dat$det)

# remove unknown sharks
eatSS.nounk <- subset(eatSS.dat, ID != "unknown")
eatSS.nounk$ID <- factor(eatSS.nounk$ID)
eatSS.nounk$det2 <- ifelse(eatSS.nounk$det == "SS", "Treat", "Control")
n.size <- dim(eatSS.nounk)[1]
barplot(xtabs(eatSS.nounk$eat ~ eatSS.nounk$det)/table(eatSS.nounk$det2))


# WAX VERSUS CONTROL
eatWX.dat <- subset(E.dat, det == "Wax" | det == "Control")
eatWX.dat$trialset <- as.integer(eatWX.dat$trialset)
eatWX.dat$det <- factor(eatWX.dat$det)

# remove unknown sharks
eatWX.nounk <- subset(eatWX.dat, ID != "unknown")
eatWX.nounk$ID <- factor(eatWX.nounk$ID)
eatWX.nounk$det2 <- ifelse(eatWX.nounk$det == "Wax", "Treat", "Control")
n.size <- dim(eatWX.nounk)[1]
barplot(xtabs(eatWX.nounk$eat ~ eatWX.nounk$det)/table(eatWX.nounk$det2))

# ignore errors
library(plyr)
safe.glmer <- failwith(NULL, glmer, quiet=TRUE)
safe.logLik <- failwith(NULL, logLik, quiet=TRUE)

## choose which deterrant-control pair
eatPAIR.nounk <- eatMB.nounk
#eatPAIR.nounk <- eatML.nounk
#eatPAIR.nounk <- eatRP.nounk
#eatPAIR.nounk <- eatSS.nounk
#eatPAIR.nounk <- eatWX.nounk

n.size <- dim(eatPAIR.nounk)[1]

# estimate binomial probability
prob <- sum(eatPAIR.nounk$eat)/n.size

# set up increment vector
prob.inc.vec <- seq(from=prob,to=0.50*prob,by=-((prob - 0.50*prob)/10))
inc.decr <- seq(0,50,5)

# iterations
iter <- 1000
itdiv <- iter/100

# model set
m1 <- "eatnew ~ det2 + trialset + det2*trialset + (1|ID)"
m2 <- "eatnew ~ det2 + trialset + (1|ID)"
m3 <- "eatnew ~ det2 + (1|ID)"
m4 <- "eatnew ~ trialset + (1|ID)"
m5 <- "eatnew ~ 1 + (1|ID)"

## Make model vector
mod.vec <- c(m1,m2,m3,m4,m5)

## Define n.mod
n.mod <- length(mod.vec)

# storage for each increment
topmod.most <- mod3AICw.med <- mod3AICw.lo <- mod3AICw.up <- mod3Rc.med <- mod3Rc.lo <- mod3Rc.up <- mod3ER.med <- mod3ER.lo <- mod3ER.up <- rep(0,length(prob.inc.vec))

pb <- txtProgressBar(min=1, max=length(prob.inc.vec), style=3)

# increment increase in response for treatment
for (n in 1:length(prob.inc.vec)) {
 
  if (prob.inc.vec[n] == prob.inc.vec[1]) {
    eatPAIR.nounk$eatnew <- eatPAIR.nounk$eat
  }

  if (prob.inc.vec[n] < prob.inc.vec[1]) {
    eatPAIR.nounk$eatnew <- eatPAIR.nounk$eat
    treat.set <- subset(eatPAIR.nounk, det2 == "Treat")
    cntrl.set <- subset(eatPAIR.nounk, det2 == "Control")
    ltreatset <- dim(treat.set)[1]
    treat.set$eatnew <- rbinom(ltreatset,1,prob.inc.vec[n])
    eatPAIR.nounk2 <- rbind(treat.set, cntrl.set)
  }
  
  # model outcome storage vectors
  topmod.vec <- topmodwAICc.vec <- topmodRc.vec <- mod3wAICc.vec <- mod3ER.vec <- mod3Rc.vec <- rep(0,iter)
  pb2 <- txtProgressBar(min=1, max=iter, style=3)
  
  # resample dataset for iter iterations
  for (i in 1:iter) {
    resamp.sub <- sample(1:n.size, n.size, replace=TRUE)
    eatPAIR.resamp <- eatPAIR.nounk2[resamp.sub,]
    
    # Model fitting and logLik output loop
    Modnum <- length(mod.vec)
    LL.vec <- SaveCount <- AICc.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
    mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
    mod.num <- seq(1,Modnum,1)
    
    for(m in 1:Modnum) {
      fit <- safe.glmer(as.formula(mod.vec[m]), family=binomial(link="logit"), data=eatPAIR.resamp, na.action=na.omit)
      assign(paste("fit",m,sep=""), fit)
      mod.list[[m]] <- fit
      LL.vec[m] <- ifelse(is.null(safe.logLik(fit)), NA, as.numeric(safe.logLik(fit)))
      k.vec[m] <- ifelse(is.null(safe.logLik(fit)), NA, attr(logLik(fit),"df"))
      AICc.vec[m] <- ifelse(is.null(safe.logLik(fit)), NA, r.squared(fit)$AIC)
      Rm[m] <- ifelse(is.null(safe.logLik(fit)), NA, 100*r.squared(fit)$Marginal) # marginal R-squared
      Rc[m] <- ifelse(is.null(safe.logLik(fit)), NA, 100*r.squared(fit)$Conditional) # conditional R-squared
    }
    
    dAICc <- delta.IC(AICc.vec)
    wAICc <- weight.IC(dAICc)
    
    sumtable <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,10),round(Rm,4),round(Rc,4))
    colnames(sumtable) <- c("model","k","LL","AICc","dAICc","wAICc","Rc","Rm")
    row.names(sumtable) <- mod.vec
    summary.table <- sumtable[order(sumtable[,5],decreasing=F),1:8]
    summary.table
    
    topmod.vec[i] <- ifelse(is.na(summary.table[1,3]), NA, summary.table[1,1])
    topmodwAICc.vec[i] <- summary.table[1,6]
    topmodRc.vec[i] <- summary.table[1,7]
    
    mod3wAICc.vec[i] <- summary.table[which(summary.table[,1] == 3), 6]
    mod3Rc.vec[i] <- summary.table[which(summary.table[,1] == 3), 7]
    mod3ER.vec[i] <- mod3wAICc.vec[i] / summary.table[which(summary.table[,1] == 5), 6]
    
    #if (i %% itdiv==0) print(i)
    setTxtProgressBar(pb2, i) # iterations status bar
  }
  close(pb2)
  
  topmod.most.sub <- as.numeric(attr(rev(sort(table(topmod.vec))), 'names')[1])
  topmod.most[n] <- mod.vec[topmod.most.sub]
  
  mod3AICw.med[n] <- median(mod3wAICc.vec, na.rm=T) 
  mod3AICw.lo[n] <- quantile(mod3wAICc.vec, probs=0.025, na.rm=T) 
  mod3AICw.up[n] <- quantile(mod3wAICc.vec, probs=0.975, na.rm=T) 
  mod3ER.med[n] <- median(mod3ER.vec, na.rm=T) 
  mod3ER.lo[n] <- quantile(mod3ER.vec, probs=0.025, na.rm=T) 
  mod3ER.up[n] <- quantile(mod3ER.vec, probs=0.975, na.rm=T) 
  mod3Rc.med[n] <- median(mod3Rc.vec, na.rm=T) 
  mod3Rc.lo[n] <- quantile(mod3Rc.vec, probs=0.025, na.rm=T) 
  mod3Rc.up[n] <- quantile(mod3Rc.vec, probs=0.975, na.rm=T)
  
  #print("-----------------------")
  #print(paste("increment number =", n))
  #print("-----------------------")
  setTxtProgressBar(pb, n) # increments status bar
}
close(pb)

par(mfrow=c(1,3))
plot(prob-prob.inc.vec, mod3AICw.med, type="l", lty=2, lwd=2, xlab="binomial prob decrease", ylab="model 3 wAICc", ylim=c(min(mod3AICw.lo), max(mod3AICw.up)))
lines(prob-prob.inc.vec, mod3AICw.lo, lty=1, lwd=1, col="red")
lines(prob-prob.inc.vec, mod3AICw.up, lty=1, lwd=1, col="red")

plot(prob-prob.inc.vec, log10(mod3ER.med), type="l", lty=2, lwd=2, xlab="binomial prob decrease", ylab="model 3 log10(ER)", ylim=c(log10(min(mod3ER.lo)), log10(max(mod3ER.up))))
lines(prob-prob.inc.vec, log10(mod3ER.lo), lty=1, lwd=1, col="red")
lines(prob-prob.inc.vec, log10(mod3ER.up), lty=1, lwd=1, col="red")

plot(prob-prob.inc.vec, mod3Rc.med, type="l", lty=2, lwd=2, xlab="binomial prob decrease", ylab="model 3 Rc", ylim=c(min(mod3Rc.lo), max(mod3Rc.up)))
lines(prob-prob.inc.vec, mod3Rc.lo, lty=1, lwd=1, col="red")
lines(prob-prob.inc.vec, mod3Rc.up, lty=1, lwd=1, col="red")
par(mfrow=c(1,1))

print(data.frame(prob-prob.inc.vec, topmod.most))

results.out <- data.frame(inc.decr, prob-prob.inc.vec, mod3AICw.med, mod3AICw.up, mod3AICw.lo, mod3ER.med, mod3ER.up, mod3ER.lo, mod3Rc.med, mod3Rc.up, mod3Rc.lo)
common.model.out <- data.frame(inc.decr, prob-prob.inc.vec, topmod.most)
