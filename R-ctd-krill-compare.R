#compares krill profiles to ctd data using density calculated from 10x50m extracted integration intervals
#author: Lisa-Marie Harrison
#date: 24/08/2015

setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/BROKE-West")
density <- read.csv("brokewest_krill_ctd.csv", header = T)
source("C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/R code/R-mixed-models/calc_conditional_marginal_Rsquared.R")
library(car)
library(caret)
library(nlme)
library(lme4)
library(flux)

#get glm.spl
source("C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/R code/R-broke-west/R-set-up-fluoro.R")

#find stations where krill data is available and subset glm.spl
stn <- unique(density$stn)
fluoro <- glm.spl[glm.spl$stn %in% stn, ]

#transform krill to be on the same scale as ctd data
p <- rep(NA, nrow(fluoro))
for (i in 1:nrow(density)) {
  w <- which(fluoro$stn == density$stn[i] & (abs(fluoro$z - density$depth[i]) <= 1))[1]
  p[w] <- density$p[i]
}

#number of stations included
length(unique(fluoro$stn))

#plot log krill against environmental variables
plot(cbind(fluoro[c(1:2, 6:9)], log(p)))

#plot krill presence/absence against environmental variables
pa <- rep(NA, length(p))
pa[p > 0] <- 1
pa[p == 0] <- 0
plot(cbind(fluoro[c(1:2, 6:9)], pa))

#boxplots for presence/absence
par(mfrow = c(1, 6))
boxplot(fluoro$oxy ~ pa)
boxplot(fluoro$sal ~ pa)
boxplot(fluoro$z ~ pa)
boxplot(fluoro$par ~ pa)
boxplot(fluoro$temp ~ pa)
boxplot(fluoro$l.obs ~ pa)

d <- data.frame(cbind(pa, fluoro$oxy, fluoro$sal, fluoro$z, fluoro$par, fluoro$temp, p, fluoro$stn))
colnames(d) <- c("pa", "oxy", "sal", "z", "par", "temp", "p", "stn")
d <- na.omit(d)

#-------------------- BINOMIAL GLM FOR PRESENCE/ABSENCE -----------------------#

pa.lm <- glm(pa ~ sal + z + oxy + par + stn - 1, dat = d, family = "binomial")
summary(pa.lm)

#mixed model with station random effect
pa.lm <- glmer(pa ~ sal + z + par + (1|stn), dat = d, family = "binomial")


#table of false and true 0 and 1
table(na.omit(d)$pa, round(fitted(pa.lm)))

#calculate variance inflation factors (<5 = good)
vif(pa.lm)

#calculate sensitivity and specificity
sensitivity(as.factor(round(fitted(pa.lm))), as.factor(na.omit(d)$pa))
specificity(as.factor(round(fitted(pa.lm))), as.factor(na.omit(d)$pa))

#plot a ROC curve for the binomial glm
roc.curve <- function(s, print = FALSE) {
  Ps <- (S > s)*1
  FP <- sum((Ps == 1)*(Y == 0))/sum(Y == 0)
  TP <- sum((Ps == 1)*(Y == 1))/sum(Y == 1)
  if (print) {
    print(table(Observed = Y, Predicted = Ps))
  }
  vect <- c(FP, TP)
  names(vect) <- c("FPR", "TPR")
  return(vect)
}

S <- predict(pa.lm, type = "response")
threshold <- 0.5
Y <- na.omit(d)$pa
roc.curve(threshold, print = TRUE)
ROC.curve <- Vectorize(roc.curve)
M.ROC <- ROC.curve(seq(0, 1, by = 0.01))

par(mfrow = c(1, 1))
plot(M.ROC[1, ], M.ROC[2, ], lwd = 2, type = "l", xlab = "False Positive Rate", ylab = "True Positive Rate")
title("ROC curve")

#calculate the area under the ROC curve (0.5 = bad, 1 = perfect)
auc(M.ROC[1,], M.ROC[2,])


#------------------ LINEAR MODELS FOR DENSITY GIVEN PRESENCE ------------------#

#plot log(density) and depth, the best predictor
par(mfrow = c(1, 1))
plot(fluoro$z, log(p), col = fluoro$stn, pch = 19)

#subset data frame to get only presence
dat <- d[d$pa == 1, ]
dat$stn <- as.factor(dat$stn)

#simple linear model
p.lm <- lm(log(p) ~ z, dat = dat)
summary(p.lm)


#model with station random intercept on log scale
p.lm <- lme(log(p) ~ z, dat = dat, random =~ 1 | stn, na.action = "na.omit")
summary(p.lm)

par(mfrow = c(1, 1))
plot(dat$z, log(dat$p), xlab = "depth", ylab = "log(krill density)", pch = 19, col = dat$stn, ylim = c(-3, 10))
y <- p.lm$coefficients$fixed[1] + p.lm$coefficients$fixed[2]*dat$z
x <- dat$z
xy <- cbind(x, y)
xy <- xy[order(xy[, 1]), ]
points(xy[, 1], xy[, 2], col = "black", type = "l", lwd = 4)
title("log(krill) vs depth")

for(i in 1:nlevels(dat$stn)) {
  y <- p.lm$coefficients$fixed[1] + p.lm$coefficients$random$stn[i] + p.lm$coefficients$fixed[2]*dat$z[dat$stn == levels(dat$stn)[i]]
  x <- dat$z[dat$stn == levels(dat$stn)[i]]
  xy <- cbind(x, y)
  xy <- xy[order(xy[, 1]), ]
  points(xy[, 1], xy[, 2], type = "l", col = i)
}

#plot of model with random intercept and slope on log scale
p.lm <- lme(log(p) ~ z, dat = dat, random =~ 1 + z | stn, na.action = "na.omit")
r.squared.lme(p.lm)

plot(dat$z, log(dat$p), xlab = "depth", ylab = "log(krill density)", pch = 19, col = dat$stn, ylim = c(-2, 10))
y <- p.lm$coefficients$fixed[1] + p.lm$coefficients$fixed[2]*dat$z
x <- dat$z
xy <- cbind(x, y)
xy <- xy[order(xy[, 1]), ]
points(xy[, 1], xy[, 2], col = "black", type = "l", lwd = 4)
title("log(krill) vs depth")

for(i in 1:nlevels(dat$stn)) {
  y <- p.lm$coefficients$fixed[1] + p.lm$coefficients$random$stn[i, 1] + (p.lm$coefficients$fixed[2] + p.lm$coefficients$random$stn[i, 2])*dat$z[dat$stn == levels(dat$stn)[i]]
  x <- dat$z[dat$stn == levels(dat$stn)[i]]
  xy <- cbind(x, y)
  xy <- xy[order(xy[, 1]), ]
  points(xy[, 1], xy[, 2], type = "l", col = i)
}




plot(fluoro$z, log(p), col = "white")
for (i in unique(fluoro$stn)) {
  log_p <- na.omit(log(p)[fluoro$stn == i])
  depth <- fluoro$z[fluoro$stn == i][which(!is.na(log(p)[fluoro$stn == i]))]
  points(depth, log_p, type = "l")
}








