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
boxplot(fluoro$oxy ~ pa, main = "Oxygen")
boxplot(fluoro$sal ~ pa, main = "Salinity")
boxplot(fluoro$z ~ pa, main = "Depth")
boxplot(fluoro$par ~ pa, main = "PAR")
boxplot(fluoro$temp ~ pa, main = "temp")
boxplot(fluoro$l.obs ~ pa, main = "l.obs")

d <- data.frame(cbind(pa, fluoro$oxy, fluoro$sal, fluoro$z, fluoro$par, fluoro$temp, p, fluoro$stn))
colnames(d) <- c("pa", "oxy", "sal", "z", "par", "temp", "p", "stn")
d <- na.omit(d)

#-------------------- BINOMIAL GLM FOR PRESENCE/ABSENCE -----------------------#

pa.lm <- glm(pa ~ z + temp + sal, dat = d, family = "binomial")
summary(pa.lm)

#mixed model with station random effect

d$sal <- scale(d$sal)
d$z <- scale(d$z)
d$par <- scale(d$par)
d$temp <- scale(d$temp)
d$oxy <- scale(d$oxy)

pa.lm <- glmer(pa ~ z + temp + sal + (1|stn), data = d, family = "binomial")
summary(pa.lm)

#table of false and true 0 and 1
table(na.omit(d)$pa, round(fitted(pa.lm)))

#calculate variance inflation factors (<5 = good)
vif(pa.lm)

#calculate sensitivity and specificity
sensitivity(as.factor(round(fitted(pa.lm))), as.factor(na.omit(d)$pa), positive = "1", negative = "0")
specificity(as.factor(round(fitted(pa.lm))), as.factor(na.omit(d)$pa), positive = "1", negative = "0")

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

png(file = "roc.png", width = 1000, height = 750, res = 100)
par(mfrow = c(1, 1), mar = c(5, 5, 1, 1))
plot(M.ROC[1, ], M.ROC[2, ], lwd = 2, type = "l", xlab = "False Positive Rate", ylab = "True Positive Rate", cex.lab = 2, cex.axis = 2)
#title("ROC curve")
lines(c(0, 1), c(0, 1), col = "red")
dev.off()

#calculate the area under the ROC curve (0.5 = bad, 0.8 = good, 0.9 = excellent, 1 = perfect)
auc(M.ROC[1,], M.ROC[2,])


#------------------ LINEAR MODELS FOR DENSITY GIVEN PRESENCE ------------------#

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

for (i in 1:nlevels(dat$stn)) {
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

for (i in 1:nlevels(dat$stn)) {
  y <- p.lm$coefficients$fixed[1] + p.lm$coefficients$random$stn[i, 1] + (p.lm$coefficients$fixed[2] + p.lm$coefficients$random$stn[i, 2])*dat$z[dat$stn == levels(dat$stn)[i]]
  x <- dat$z[dat$stn == levels(dat$stn)[i]]
  xy <- cbind(x, y)
  xy <- xy[order(xy[, 1]), ]
  points(xy[, 1], xy[, 2], type = "l", col = i)
}


#------------------------ PLOTS ON THE NATURAL SCALE --------------------------#

#plot of model with random intercept and slope on log scale
p.lm <- lme(log(p) ~ z, dat = dat, random =~ 1 + z | stn, na.action = "na.omit")
r.squared.lme(p.lm)

plot(dat$z, (dat$p), xlab = "depth", ylab = "log(krill density)", pch = 19, col = dat$stn)
for (i in 1:nlevels(dat$stn)) {
  y <- exp(p.lm$coefficients$fixed[1] + p.lm$coefficients$random$stn[i, 1] + (p.lm$coefficients$fixed[2] + p.lm$coefficients$random$stn[i, 2])*dat$z[dat$stn == levels(dat$stn)[i]])
  x <- dat$z[dat$stn == levels(dat$stn)[i]]
  xy <- cbind(x, y)
  xy <- xy[order(xy[, 1]), ]
  points(xy[, 1], xy[, 2], type = "l", col = i)
}

y <- exp(p.lm$coefficients$fixed[1] + p.lm$coefficients$fixed[2]*dat$z)
x <- dat$z
xy <- cbind(x, y)
xy <- xy[order(xy[, 1]), ]
points(xy[, 1], xy[, 2], col = "black", type = "l", lwd = 4)
title("log(krill) vs depth")



#-------------------------- KRILL VS PHYTOPLANKTON ----------------------------#


#subset data frame to get only stations with 8 or more data points
d <- data.frame(cbind(pa, fluoro$oxy, fluoro$sal, fluoro$z, fluoro$par, fluoro$temp, p, fluoro$stn, fluoro$l.obs, fluoro$obs))
colnames(d) <- c("pa", "oxy", "sal", "z", "par", "temp", "p", "stn", "l.obs", "obs")
d <- na.omit(d)
d <- d[d$stn %in% sort(unique(d$stn))[which(table(d$stn) >= 10)], ]
dat <- d[d$pa == 1, ]
dat$stn <- as.factor(dat$stn)
dat <- dat[dat$stn %in% sort(unique(dat$stn))[which(table(dat$stn) >= 10)], ]



#plot of model with random slope on log scale using l.obs
p.lm <- lme(log(p) ~ exp(l.obs), random =~ exp(l.obs) | stn, data = dat, na.action = na.omit)
summary(p.lm)
r.squared.lme(p.lm)

plot(dat$l.obs, log(dat$p), xlab = "l.obs", ylab = "log(krill density)", pch = 19, col = dat$stn)
y <- (p.lm$coefficients$fixed[1] + p.lm$coefficients$fixed[2]*exp(dat$l.obs))
x <- dat$l.obs
xy <- cbind(x, y)
xy <- xy[order(xy[, 1]), ]
points(xy[, 1], xy[, 2], col = "black", type = "l", lwd = 4)
title("log(krill) vs l.obs")

dat$stn <- as.factor(dat$stn)
for (i in 1:length(unique(dat$stn))) {
  y <- p.lm$coefficients$fixed[1] + p.lm$coefficients$random$stn[i, 1] + (p.lm$coefficients$fixed[2] + p.lm$coefficients$random$stn[i, 2])*exp(dat$l.obs)[dat$stn == sort(unique(dat$stn))[i]]
  x <- dat$l.obs[dat$stn == sort(unique(dat$stn))[i]]
  xy <- cbind(x, y)
  xy <- xy[order(xy[, 1]), ]
  points(xy[, 1], xy[, 2], type = "l", col = i)
}



#plot of model with random slope on natural scale using l.obs
p.lm <- lme(log(p) ~ obs, random =~ obs - 1| stn, data = dat, na.action = na.omit)
summary(p.lm)
r.squared.lme(p.lm)

plot(dat$l.obs, (dat$p), xlab = "l.obs", ylab = "krill density", pch = 19, col = dat$stn)
y <- exp(p.lm$coefficients$fixed[1] + p.lm$coefficients$fixed[2]*dat$obs)
x <- dat$l.obs
xy <- cbind(x, y)
xy <- xy[order(xy[, 1]), ]
points(xy[, 1], xy[, 2], col = "black", type = "l", lwd = 4)
title("log(krill) vs l.obs")

dat$stn <- as.factor(dat$stn)
for (i in 1:length(unique(dat$stn))) {
  y <- exp(p.lm$coefficients$fixed[1] + (p.lm$coefficients$fixed[2] + p.lm$coefficients$random$stn[i, 1])*dat$obs[dat$stn == sort(unique(dat$stn))[i]])
  x <- dat$l.obs[dat$stn == sort(unique(dat$stn))[i]]
  xy <- cbind(x, y)
  xy <- xy[order(xy[, 1]), ]
  points(xy[, 1], xy[, 2], type = "l", col = i)
}







#plot of model with random slope on log scale using l.obs
p.lm <- lme(log(p) ~ (obs), random =~ obs*2 - 1 | stn, data = dat, na.action = na.omit)
summary(p.lm)
r.squared.lme(p.lm)

plot(log(dat$obs), log(dat$p), xlab = "l.obs", ylab = "krill density", pch = 19, col = dat$stn)
y <- (p.lm$coefficients$fixed[1] + p.lm$coefficients$fixed[2]*dat$obs)
x <- log(dat$obs)
xy <- cbind(x, y)
xy <- xy[order(xy[, 1]), ]
points(xy[, 1], xy[, 2], col = "black", type = "l", lwd = 4)
title("log(krill) vs l.obs")

dat$stn <- as.factor(dat$stn)
for (i in 1:length(unique(dat$stn))) {
  y <- p.lm$coefficients$fixed[1] + (p.lm$coefficients$fixed[2] + p.lm$coefficients$random$stn[i, 1])*dat$obs[dat$stn == sort(unique(dat$stn))[i]]
  x <- log(dat$obs)[dat$stn == sort(unique(dat$stn))[i]]
  xy <- cbind(x, y)
  xy <- xy[order(xy[, 1]), ]
  points(xy[, 1], xy[, 2], type = "l", col = i)
}



#plot of model with random slope on log scale using l.obs for presentation
p.lm <- lme(log(p) ~ exp(0.5*l.obs), random =~ exp(0.5*l.obs) - 1| stn, data = dat, na.action = na.omit)
summary(p.lm)
r.squared.lme(p.lm)

png(file = "lobs.png", width = 1000, height = 750, res = 100)
par(mar = c(5, 5, 1, 1))
plot(0.5*dat$l.obs, log(dat$p), xlab = "log(phytoplankton density)", ylab = "log(krill density)", pch = 19, col = "grey", cex.lab = 2, cex.axis = 2)
y <- (p.lm$coefficients$fixed[1] + p.lm$coefficients$fixed[2]*exp(0.5*dat$l.obs))
x <- dat$l.obs
xy <- cbind(x, y)
xy1 <- xy[order(xy[, 1]), ]

dat$stn <- as.factor(dat$stn)
for (i in 1:length(unique(dat$stn))) {
  y <- p.lm$coefficients$fixed[1] + (p.lm$coefficients$fixed[2] + p.lm$coefficients$random$stn[i, 1])*exp(0.5*dat$l.obs)[dat$stn == sort(unique(dat$stn))[i]]
  x <- dat$l.obs[dat$stn == sort(unique(dat$stn))[i]]
  xy <- cbind(x, y)
  xy <- xy[order(xy[, 1]), ]
  points(xy[, 1], xy[, 2], type = "l")
}
points(xy1[, 1], xy1[, 2], col = "red", type = "l", lwd = 4)
dev.off()




#subset data frame to get only stations with 8 or more data points
d <- data.frame(cbind(pa, fluoro$oxy, fluoro$sal, fluoro$z, fluoro$par, fluoro$temp, p, fluoro$stn, fluoro$l.obs, fluoro$obs))
colnames(d) <- c("pa", "oxy", "sal", "z", "par", "temp", "p", "stn", "l.obs", "obs")
d <- na.omit(d)
d <- d[d$stn %in% sort(unique(d$stn))[which(table(d$stn) >= 7)], ]
dat <- d[d$pa == 1, ]
dat$stn <- as.factor(dat$stn)
dat <- dat[dat$stn %in% sort(unique(dat$stn))[which(table(dat$stn) >= 7)], ]



p.lm <- lme(log(p) ~ exp(0.5*l.obs), random =~ exp(0.5*l.obs) - 1 | stn, data = dat, na.action = na.omit)
summary(p.lm)
r.squared.lme(p.lm)

png(file = "lobs.png", width = 1000, height = 750, res = 100)
par(mar = c(5, 5, 1, 1))
plot(log(dat$obs), log(dat$p), xlab = "log(phytoplankton density)", ylab = "log(krill density)", pch = 19, col = "grey", cex.lab = 2, cex.axis = 2)
y <- (p.lm$coefficients$fixed[1] + p.lm$coefficients$fixed[2]*exp(0.5*dat$l.obs))
x <- (dat$l.obs)

xy <- cbind(x, y)
xy1 <- xy[order(xy[, 1]), ]

dat$stn <- as.factor(dat$stn)
for (i in 1:length(unique(dat$stn))) {
  y <- p.lm$coefficients$fixed[1] + (p.lm$coefficients$fixed[2] + p.lm$coefficients$random$stn[i, 1])*exp(0.5*dat$l.obs)[dat$stn == sort(unique(dat$stn))[i]]
  x <- (dat$l.obs)[dat$stn == sort(unique(dat$stn))[i]]
  xy <- cbind(x, y)
  xy <- xy[order(xy[, 1]), ]
  points(xy[, 1], xy[, 2], type = "l")
}
points(xy1[, 1], xy1[, 2], col = "red", type = "l", lwd = 4)
dev.off()


#using gamm
library(mgcv)
p.lm <- gamm(log(p) ~ s(l.obs), random = list(stn =~ 1), dat = dat)
summary(p.lm$gam)


plot(dat$l.obs, log(dat$p), xlab = "log(phytoplankton density)", ylab = "log(krill density)", pch = 19, col = "grey", cex.lab = 2, cex.axis = 2)
y <- (p.lm$lme$coefficients$fixed[1] + p.lm$coefficients$fixed[2]*(dat$l.obs))
x <- dat$l.obs
xy <- cbind(x, y)
xy1 <- xy[order(xy[, 1]), ]

dat$stn <- as.factor(dat$stn)
for (i in 1:length(unique(dat$stn))) {
  y <- p.lm$lme$coefficients$fixed[1] + p.lm$lme$coefficients$random$stn[i, 1] + p.lm$lme$coefficients$fixed[2]*
  x <- dat$l.obs[dat$stn == sort(unique(dat$stn))[i]]
  xy <- cbind(x, y)
  xy <- xy[order(xy[, 1]), ]
  points(xy[, 1], xy[, 2], type = "l")
}
points(xy1[, 1], xy1[, 2], col = "red", type = "l", lwd = 4)


a <- predict(p.lm, dat = dat, type="lpmatrix")



