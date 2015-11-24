#compares krill profiles to ctd data using density calculated from 10x50m extracted integration intervals
#author: Lisa-Marie Harrison
#date: 24/08/2015

if (Sys.info()[4] == "SCI-6246") {
  setwd(dir = "C:/Users/43439535/Documents/Lisa/phd/Mixed models")
} else {
  setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models")
}

density <- read.csv("Data/brokewest_krill_ctd.csv", header = T)
library(car)
library(caret)
library(nlme)
library(lme4)
library(flux)
library(itsadug)
library(mgcv)

#source required functions
function_list <- c("setUpFluoro.R",
                   "calc_conditional_marginal_Rsquared.R",
                   "calc_asreml_conditional_marginal_Rsquared.R",
                   "rocCurve.R",
                   "krillPresenceAbsence")

for (f in function_list) {
  source(paste("R code/R-functions-southern-ocean/", f, sep = ""))
}

#find stations where krill data is available and subset glm.spl
dat <- read.csv(file = "Data/procCTD.csv", header= T)
names(dat) <- c("survey", "stn", "lat", "long", "start.time", "end.time", "depth", "transmittance", "cond", "temp", "sal", "par", "oxygen", "fluoro", "x2", "ice", "wm")

glm.spl <- setUpFluoro(dat)
stn <- unique(density$stn)
fluoro <- glm.spl[glm.spl$stn %in% stn, ]

#transform krill to be on the same scale as ctd data

transformKrill <- function (krill, fluoro) {
  
  p <- rep(NA, nrow(fluoro))
  for (i in 1:nrow(krill)) {
    w <- which(fluoro$stn == krill$stn[i] & (abs(fluoro$z - krill$depth[i]) <= 1))[1]
    p[w] <- krill$p[i]
  }
  return(p)
  
}

p <- transformKrill(density, fluoro)

#plot log krill against environmental variables
plot(cbind(fluoro[c(1:2, 6:9)], log(p)))

#plot krill presence/absence against environmental variables
pa <- krillPresenceAbsence(p)

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

d$sal  <- scale(d$sal)
d$z    <- scale(d$z)
d$par  <- scale(d$par)
d$temp <- scale(d$temp)
d$oxy  <- scale(d$oxy)

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
ROC.curve <- Vectorize(rocCurve(model = pa.lm, s = 0.5, data = pa, print = TRUE))
M.ROC     <- ROC.curve(seq(0, 1, by = 0.01))

par(mfrow = c(1, 1), mar = c(5, 5, 1, 1))
plot(M.ROC[1, ], M.ROC[2, ], lwd = 2, type = "l", xlab = "False Positive Rate", ylab = "True Positive Rate", cex.lab = 2, cex.axis = 2)
title("ROC curve")
lines(c(0, 1), c(0, 1), col = "red")

#calculate the area under the ROC curve (0.5 = bad, 0.8 = good, 0.9 = excellent, 1 = perfect)
auc(M.ROC[1,], M.ROC[2,])


#-------------------------- KRILL VS PHYTOPLANKTON ----------------------------#

#subset data frame to get only stations with 5 or more data points
d <- data.frame(cbind(pa, fluoro$oxy, fluoro$sal, fluoro$z, fluoro$par, fluoro$temp, p, fluoro$stn, fluoro$l.obs, fluoro$obs))
colnames(d) <- c("pa", "oxy", "sal", "z", "par", "temp", "p", "stn", "l.obs", "obs")
d <- na.omit(d)
dat <- d[d$pa == 1, ]
dat$stn <- as.factor(dat$stn)
dat <- dat[dat$stn %in% sort(unique(dat$stn))[which(table(dat$stn) >= 5)], ]

dat$pwr <- 2^dat$l.obs

p.lm <- lme(log(p) ~ pwr, random =~ pwr -1 | stn, data = dat, na.action = na.omit)
summary(p.lm)
r.squared.lme(p.lm)

par(mar = c(5, 5, 1, 1))
plot(log(dat$obs), log(dat$p), xlab = "log(phytoplankton density)", ylab = "log(krill density)", pch = 19, col = "grey", cex.lab = 2, cex.axis = 2)
y <- (p.lm$coefficients$fixed[1] + p.lm$coefficients$fixed[2]*dat$pwr)
x <- (dat$l.obs)
xy <- cbind(x, y)
xy1 <- xy[order(xy[, 1]), ]

dat$stn <- as.factor(dat$stn)
for (i in 1:length(unique(dat$stn))) {
  y <- p.lm$coefficients$fixed[1] + (p.lm$coefficients$fixed[2] + p.lm$coefficients$random$stn[i, 1])*dat$pwr[dat$stn == sort(unique(dat$stn))[i]]
  x <- (dat$l.obs)[dat$stn == sort(unique(dat$stn))[i]]
  xy <- cbind(x, y)
  xy <- xy[order(xy[, 1]), ]
  points(xy[, 1], xy[, 2], type = "l")
}
points(xy1[, 1], xy1[, 2], col = "red", type = "l", lwd = 4)


#using gamm
p.lm <- gam(log(p) ~ s(l.obs) + s(stn, bs = "re"), dat = dat, gamma = 2)
summary(p.lm)
plot(p.lm)

plot_smooth(p.lm, view = "l.obs", cond = list(stn = unique(dat$stn)[1]), se = FALSE)
for (i in sort(unique(dat$stn))) {
  plot_smooth(p.lm, view = "l.obs", cond = list(stn = i), se = FALSE, add = TRUE)
}

p.lm <- gamm(log(p) ~ s(l.obs), random = list(stn =~ -1 + l.obs), dat = dat, control = list(niter = 100000))
pd <- plot(p.lm$gam)

plot(log(dat$obs), log(dat$p), pch = 19, col = "grey")
points(pd[[1]]$x, pd[[1]]$fit, type = "l")
for (i in sort(unique(dat$stn))) {
  points(pd[[1]]$x, (pd[[1]]$fit-0.471522814)*p.lm$lme$coefficients$random$stn[grep(rownames(p.lm$lme$coefficients$random$stn), pattern = i), 2] + p.lm$lme$coefficients$random$stn[grep(rownames(p.lm$lme$coefficients$random$stn), pattern = i), 1], type = "l")
}


#using asreml
dat_full <- as.data.frame(matrix(NA, nrow = length(unique(dat$stn))*length(unique(dat$z)), ncol = ncol(dat)))
names(dat_full) <- names(dat)
dat_full$z <- rep(unique(dat$z), length(unique(dat$stn)))
dat_full$stn <- sort(rep(unique(dat$stn), length(unique(dat$z))))
for (i in 1:nrow(dat)) {
  w <- which(dat_full$stn == dat$stn[i] & dat_full$z == dat$z[i])
  dat_full[w, ] <- dat[i, ]
}
dat_full$z.fact <- as.factor(dat_full$z)
dat_full <- dat_full[order(as.numeric(as.character(dat_full$stn))), ]



p.lm <- asreml(log(p) ~ obs, random =~ obs:stn, data = dat_full, na.method.X = 'include')
summary(p.lm)

plot(dat_full$l.obs, log(dat_full$p), xlab = "log(phytoplankton density)", ylab = "log(krill density)", pch = 19, col = "grey", cex.lab = 2, cex.axis = 2)
y <- p.lm$coefficients$fixed[2] + p.lm$coefficients$fixed[1]*dat_full$obs2
x <- dat_full$l.obs
xy <- cbind(x, y)
xy1 <- xy[order(xy[, 1]), ]

for (i in 1:length(unique(dat_full$stn))) {
  y <- fitted(p.lm)[dat_full$stn == sort(unique(dat_full$stn))[i]]
  x <- dat_full$l.obs[dat_full$stn == sort(unique(dat_full$stn))[i]]
  xy <- cbind(x, y)
  xy <- xy[order(xy[, 1]), ]
  points(xy[, 1], xy[, 2], type = "l")
}
points(xy1[, 1], xy1[, 2], col = "red", type = "l", lwd = 4)

calcRsquared(p.lm, "obs:stn")
