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
library(AICcmodavg)

#source required functions
function_list <- c("setUpFluoro.R",
                   "calc_conditional_marginal_Rsquared.R",
                   "calc_asreml_conditional_marginal_Rsquared.R",
                   "rocCurve.R",
                   "krillPresenceAbsence.R")

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

#bottom depth

depths <- read.csv("C:/Users/Lisa/Documents/phd/southern ocean/depths/brokeWestDepths.csv")
stn_coords <- read.csv("~/phd/southern ocean/Mixed models/Data/dat.stn.csv", header = T)
stn_coords <- stn_coords[, 2:4]
stn_coords <- stn_coords[!duplicated(stn_coords[, 1]), ]

depths$Latitude[depths$Latitude == 999] <- NA
depths$Longitude[depths$Longitude == 999] <- NA

locV = NA
for(i in 1:nrow(stn_coords))
{
  dists=distHaversine(p1=c(stn_coords$Longitude[i],stn_coords$Latitude[i]),p2=cbind(depths$Longitude, depths$Latitude))
  minDistLOC=which.min(dists)
  message('Minimum distance = ',dists[minDistLOC],' m')  
  locV[i]=minDistLOC
}

stn_coords$depth <- depths$Depth[locV]

d <- data.frame(cbind(pa, fluoro$oxy, fluoro$sal, fluoro$z, fluoro$par, fluoro$temp, p, fluoro$stn, fluoro$obs))
colnames(d) <- c("pa", "oxy", "sal", "z", "par", "temp", "p", "stn", "obs")
d <- na.omit(d)
d$depth <- stn_coords$depth[match(d$stn, stn_coords$Cast.Number)]

#bottomm depth plots
boxplot(-d$depth ~ d$pa)

unscaled <- d
#-------------------- BINOMIAL GLM FOR PRESENCE/ABSENCE -----------------------#
 
#glm
pa.lm <- glm(pa ~ z + temp + sal + par, dat = d, family = "binomial")
summary(pa.lm)

#calculate variance inflation factors (<5 = good)
vif(pa.lm)

#scale or model doesn't converge
d <- cbind(d[, c(1, 7:8)], apply(d[, c(2:6, 9)], 2, scale))

#mixed model with station random effect
pa.lm <- glmer(pa ~ z + temp + sal + par -1 +(1|stn), data = d, family = "binomial")
summary(pa.lm)

#calculate sensitivity and specificity
sensitivity(as.factor(round(fitted(pa.lm))), as.factor(na.omit(d)$pa), positive = "1", negative = "0")
specificity(as.factor(round(fitted(pa.lm))), as.factor(na.omit(d)$pa), positive = "1", negative = "0")

#plot a ROC curve for the binomial glm
M.ROC <- rocCurve(model = pa.lm, threshold = 0.5, data = pa, print = TRUE)

par(mfrow = c(1, 1), mar = c(5, 5, 1, 1))
plot(M.ROC[1, ], M.ROC[2, ], lwd = 2, type = "l", xlab = "False Positive Rate", ylab = "True Positive Rate", cex.lab = 2, cex.axis = 2)
title("ROC curve")
lines(c(0, 1), c(0, 1), col = "red")

#calculate the area under the ROC curve (0.5 = bad, 0.8 = good, 0.9 = excellent, 1 = perfect)
auc(M.ROC[1,], M.ROC[2,])

#partial plots
par(mfrow = c(2, 2))

predict_pa_re <- expand.grid(seq(min(d$z), max(d$z), length.out = 100), unique(d$stn))
predict_pa <- data.frame("z" = predict_pa_re$Var1, "stn" = predict_pa_re$Var2, "temp" = 0, "sal" = 0, "par" = 0)
pred_z <- predict(pa.lm, newdata = predict_pa, allow.new.level = T, type = "response")
pred_z <- aggregate(pred_z, list(predict_pa$z), FUN = mean)
plot(pred_z$Group.1 * sd(unscaled$z) + mean(unscaled$z), pred_z$x, main = "Depth (m)", type = "l", xlab = "Depth", ylab = "Probability of krill", ylim = c(0, 1))

predict_pa_re <- expand.grid(seq(min(d$temp), max(d$temp), length.out = 100), unique(d$stn))
predict_pa <- data.frame("z" = 0, "stn" = predict_pa_re$Var2, "temp" = predict_pa_re$Var1, "sal" = 0, "par" = 0)
pred_temp <- predict(pa.lm, newdata = predict_pa, allow.new.level = T, type = "response")
pred_temp <- aggregate(pred_temp, list(predict_pa$temp), FUN = mean)
plot(pred_z$Group.1 * sd(unscaled$temp) + mean(unscaled$temp), pred_temp$x, main = "Temperature (Degrees celcius)", ylim = c(0, 1), type = "l", xlab = "Temperature", ylab = "Probability of krill")

predict_pa_re <- expand.grid(seq(min(d$sal), max(d$sal), length.out = 100), unique(d$stn))
predict_pa <- data.frame("z" = 0, "stn" = predict_pa_re$Var2, "temp" = 0, "sal" = predict_pa_re$Var1, "par" = 0)
pred_sal <- predict(pa.lm, newdata = predict_pa, allow.new.level = T, type = "response")
pred_sal <- aggregate(pred_sal, list(predict_pa$sal), FUN = mean)
plot(pred_z$Group.1 * sd(unscaled$sal) + mean(unscaled$sal), pred_sal$x, main = "Salinity", ylim = c(0, 1), type = "l", xlab = "Salinity", ylab = "Probability of krill")

predict_pa_re <- expand.grid(seq(min(d$par), max(d$par), length.out = 100), unique(d$stn))
predict_pa <- data.frame("z" = 0, "stn" = predict_pa_re$Var2, "temp" = 0, "sal" = 0, "par" = predict_pa_re$Var1)
pred_par <- predict(pa.lm, newdata = predict_pa, allow.new.level = T, type = "response")
pred_par <- aggregate(pred_par, list(predict_pa$par), FUN = mean)
plot(pred_z$Group.1 * sd(unscaled$par) + mean(unscaled$par), pred_par$x, main = "PAR", ylim = c(0, 1), type = "l", xlab = "PAR", ylab = "Probability of krill")






#------------------------ cross validation drop 1 station ---------------------------------#

truth <- NULL
pred  <- NULL
for (i in unique(d$stn)) {
  
  pa.lm <- glmer(pa ~ z + temp + sal + par - 1 + (1|stn), data = d[d$stn != i, ], family = "binomial")
  pred <- c(pred, predict(pa.lm, newdata = d[d$stn == i, ], allow.new.levels = T, type = "response"))
  truth <- c(truth, d$pa[d$stn == i])
  
}

pred <- round(pred)

table(pred, truth)

sensitivity(data = as.factor(pred), reference = as.factor(truth), positive = "1", negative = "0") #true positive rate
specificity(data = as.factor(pred), reference = as.factor(truth), positive = "1", negative = "0") #true negative rate


#--------------------- simulating random effects for cross validation ------------------------#

#Need to add random effect back in. Can't leave it out because data transformation in binomial glm will skew results
#Using example from Simon Wotherspoon

ilogit <- function(x) 1/(1+exp(-x))

truth <- NULL
pred  <- NULL
for (i in unique(d$stn)) {
  
  pa.lm <- glmer(pa ~ z + temp + sal + par - 1 + (1|stn), data = d[d$stn != i, ], family = "binomial")
  stn_sd <- unlist(lapply(VarCorr(pa.lm), function(m) sqrt(diag(m))))
  fixed_effect <- rowSums(sweep(cbind(d$z[d$stn == i], d$temp[d$stn == i], d$sal[d$stn == i], d$par[d$stn == i]),MARGIN=2,fixef(pa.lm),`*`))

  #resample to add in station random effect using extracted random effect sd
  pred_sample <- NULL
  for (j in 1:1000) {
    pred_sample <- cbind(pred_sample, ilogit(fixed_effect + rnorm(1, mean = 0, sd = stn_sd)))
  }
  
  pred <- c(pred, rowMeans(pred_sample)) #use average of resampled predictions
  truth <- c(truth, d$pa[d$stn == i])
  
}

pred <- round(pred)

table(pred, truth)

sensitivity(data = as.factor(pred), reference = as.factor(truth), positive = "1", negative = "0") #true positive rate
specificity(data = as.factor(pred), reference = as.factor(truth), positive = "1", negative = "0") #true negative rate


#------------------------ cross validation drop 1 station without random effect -------------------------------#

truth <- NULL
pred  <- NULL
for (i in unique(d$stn)) {
  
  pa.lm <- glm(pa ~ z + temp + sal + par - 1, data = d[d$stn != i, ], family = "binomial")
  pred <- c(pred, predict(pa.lm, newdata = d[d$stn == i, ], allow.new.levels = T, type = "response"))
  truth <- c(truth, d$pa[d$stn == i])
  
}

pred <- round(pred)

table(pred, truth)

sensitivity(data = as.factor(pred), reference = as.factor(truth), positive = "1", negative = "0") #true positive rate
specificity(data = as.factor(pred), reference = as.factor(truth), positive = "1", negative = "0") #true negative rate



#-------------------------- KRILL VS PHYTOPLANKTON ----------------------------#

d <- data.frame(cbind(pa, fluoro$oxy, fluoro$sal, fluoro$z, fluoro$par, fluoro$temp, p, fluoro$stn, fluoro$l.obs, fluoro$obs))
colnames(d) <- c("pa", "oxy", "sal", "z", "par", "temp", "p", "stn", "l.obs", "obs")
d <- na.omit(d)
dat <- d[d$pa == 1, ]
dat$stn <- as.factor(dat$stn)
dat <- dat[dat$stn %in% sort(unique(dat$stn))[which(table(dat$stn) >= 5)], ]

dat$lat <- stn_coords$Latitude[match(dat$stn, stn_coords$Cast.Number)]
dat$long <- stn_coords$Longitude[match(dat$stn, stn_coords$Cast.Number)]

locV = NA
for(i in 1:nrow(dat))
{
  dists=distHaversine(p1=c(dat$long[i], dat$lat[i]),p2=cbind(depths$Longitude, depths$Latitude))
  minDistLOC=which.min(dists)
  message('Minimum distance = ',dists[minDistLOC],' m')  
  locV[i]=minDistLOC
}

dat$depth <- depths$Depth[locV]

plot(dat$depth[!duplicated(dat$stn)], aggregate(dat$p, list(dat$stn), mean)$x, xlab = "bottom depth (m)",
     ylab = "mean stn krill density", pch = 19)



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
p.lm <- gam(log(p) ~ s(l.obs) + s(l.obs, stn, bs = "re"), dat = dat, gamma = 2)
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

#oxygen with linear relationship

p.lm <- lme(log(p) ~ oxy, random =~ 1 + oxy | stn, data = dat, na.action = na.omit, control = list(opt='optim'))
summary(p.lm)
r.squared.lme(p.lm)


par(mar = c(5, 5, 1, 1))
plot(dat$oxy, log(dat$p), xlab = "oxygen", ylab = "log(krill density)", pch = 19, col = "grey", cex.lab = 2, cex.axis = 2)
y <- (p.lm$coefficients$fixed[1] + p.lm$coefficients$fixed[2]*dat$oxy)
x <- (dat$oxy)
xy <- cbind(x, y)
xy1 <- xy[order(xy[, 1]), ]

dat$stn <- as.factor(dat$stn)
for (i in 1:length(unique(dat$stn))) {
  y <- p.lm$coefficients$fixed[1] + p.lm$coefficients$random$stn[i, 1] + (p.lm$coefficients$fixed[2] + p.lm$coefficients$random$stn[i, 2])*dat$oxy[dat$stn == sort(unique(dat$stn))[i]]
  x <- (dat$oxy)[dat$stn == sort(unique(dat$stn))[i]]
  xy <- cbind(x, y)
  xy <- xy[order(xy[, 1]), ]
  points(xy[, 1], xy[, 2], type = "l")
}
points(xy1[, 1], xy1[, 2], col = "red", type = "l", lwd = 4)


#fluoro with linear relationship

p.lm <- lme(log(p) ~ l.obs, random =~ 1 + l.obs | stn, data = dat, na.action = na.omit, control = list(opt='optim'))
summary(p.lm)
r.squared.lme(p.lm)


par(mar = c(5, 5, 1, 1))
plot(dat$l.obs, log(dat$p), xlab = "l.fluoro", ylab = "log(krill density)", pch = 19, col = "grey", cex.lab = 2, cex.axis = 2)
y <- (p.lm$coefficients$fixed[1] + p.lm$coefficients$fixed[2]*dat$l.obs)
x <- (dat$l.obs)
xy <- cbind(x, y)
xy1 <- xy[order(xy[, 1]), ]

dat$stn <- as.factor(dat$stn)
for (i in 1:length(unique(dat$stn))) {
  y <- p.lm$coefficients$fixed[1] + p.lm$coefficients$random$stn[i, 1] + (p.lm$coefficients$fixed[2] + p.lm$coefficients$random$stn[i, 2])*dat$l.obs[dat$stn == sort(unique(dat$stn))[i]]
  x <- (dat$l.obs)[dat$stn == sort(unique(dat$stn))[i]]
  xy <- cbind(x, y)
  xy <- xy[order(xy[, 1]), ]
  points(xy[, 1], xy[, 2], type = "l")
}
points(xy1[, 1], xy1[, 2], col = "red", type = "l", lwd = 4)


#--------------- fluoro and oxy with linear relationship and interaction term ---------------#

dat$oxy <- scale(dat$oxy)
dat$obs <- scale(dat$obs)

p.lm <- lme(log(p) ~ obs * oxy, random =~ 1 + oxy + obs | stn, data = dat, na.action = na.omit, 
            control = list(opt='optim'))
summary(p.lm)
r.squared.lme(p.lm)


#3D plot of obs*oxy interaction
interaction_data <- expand.grid(seq(min(dat$oxy), max(dat$oxy), length.out = 20), seq(min(dat$obs), max(dat$obs), length.out = 20), unique(dat$stn))
colnames(interaction_data) <- c("oxy", "obs", "stn")
pred_interaction <- predict(p.lm, newdata = interaction_data)

pred_fixed <- aggregate(pred_interaction, list(interaction_data$obs, interaction_data$oxy), FUN = mean)
pred_fixed$x[exp(pred_fixed$x) > 500] <- NA

scatter3d(pred_fixed$Group.1, pred_fixed$Group.2, exp(pred_fixed$x), xlab = "Phyto", ylab = "Oxygen", zlab = "Krill density")

#extract fitted including only fixed effects
y <- p.lm$fitted[, 1]

#plot residuals against fitted and covariates
par(mfrow = c(1, 5))
plot(fitted(p.lm), resid(p.lm, type = "normalized"))
plot(dat$obs, resid(p.lm, type = "normalized"))
plot(dat$oxy, resid(p.lm, type = "normalized"))
plot(log(dat$p), y, main = "observed vs fitted with only fixed effects")
plot(log(dat$p), fitted(p.lm), main = "observed vs fitted")

#3d scatter plot of residuals against covariates
scatter3d(dat$obs, resid(p.lm, type = "normalized"), dat$oxy)

#plot of only oxygen
par(mar = c(5, 5, 1, 1))
plot(dat$oxy, log(dat$p), xlab = "oxygen", ylab = "log(krill density)", pch = 19, col = "grey", cex.lab = 2, cex.axis = 2)
y <- (p.lm$coefficients$fixed[1] + p.lm$coefficients$fixed[3]*dat$oxy)
x <- (dat$oxy)
xy <- cbind(x, y)
xy1 <- xy[order(xy[, 1]), ]

dat$stn <- as.factor(dat$stn)
for (i in 1:length(unique(dat$stn))) {
  y <- p.lm$coefficients$fixed[1] + p.lm$coefficients$random$stn[i, 1] + 
    (p.lm$coefficients$fixed[3] + p.lm$coefficients$random$stn[i, 2])*dat$oxy[dat$stn == sort(unique(dat$stn))[i]] 
  x <- (dat$oxy)[dat$stn == sort(unique(dat$stn))[i]]
  xy <- cbind(x, y)
  xy <- xy[order(xy[, 1]), ]
  points(xy[, 1], xy[, 2], type = "l")
}
points(xy1[, 1], xy1[, 2], col = "red", type = "l", lwd = 4)

#plot of only l.obs
par(mar = c(5, 5, 1, 1))
plot(dat$l.obs, log(dat$p), xlab = "l.obs", ylab = "log(krill density)", pch = 19, col = "grey", cex.lab = 2, cex.axis = 2)
y <- (p.lm$coefficients$fixed[1] + p.lm$coefficients$fixed[2]*dat$l.obs)
x <- (dat$l.obs)
xy <- cbind(x, y)
xy1 <- xy[order(xy[, 1]), ]

dat$stn <- as.factor(dat$stn)
for (i in 1:length(unique(dat$stn))) {
  y <- p.lm$coefficients$fixed[1] + p.lm$coefficients$random$stn[i, 1] + 
    (p.lm$coefficients$fixed[2] + p.lm$coefficients$random$stn[i, 3])*dat$l.obs[dat$stn == sort(unique(dat$stn))[i]] 
  x <- (dat$l.obs)[dat$stn == sort(unique(dat$stn))[i]]
  xy <- cbind(x, y)
  xy <- xy[order(xy[, 1]), ]
  points(xy[, 1], xy[, 2], type = "l")
}
points(xy1[, 1], xy1[, 2], col = "red", type = "l", lwd = 4)

#----------------------------- cross validation drop 1 station --------------------------#

pred <- NULL
pred_se <- NULL
truth <- NULL
for (i in unique(dat$stn)) {

  p.lm <- lme(log(p) ~ obs * oxy, random =~ 1 + oxy + obs | stn, data = dat[dat$stn != i, ], na.action = na.exclude, 
              control = list(opt='optim'), weights = varExp(form =~ oxy))
  new_predictions <- predictSE.lme(p.lm, newdata = dat[dat$stn == i, ], na.action = na.omit)
  pred <- c(pred, exp(new_predictions$fit))
  pred_se <- c(pred_se, new_predictions$se.fit)
  truth <- c(truth, dat$p[dat$stn == i])
  
}

rmse  <- sqrt(sum(na.omit((pred - truth)^2))/length(truth)) #calculate root mean square error

rmse/(max(truth) - min(truth)) #rmse % of true range

plot(log(truth), log(pred))


#----------------------------- cross validation drop 1 point --------------------------#

pred <- NULL
truth <- NULL
for (i in 1:nrow(dat)) {
  
  p.lm <- lme(log(p) ~ obs * oxy, random =~ 1 + oxy + obs | stn, data = dat[-i, ], na.action = na.omit, 
              control = list(opt='optim'), weights = varExp(form =~ oxy))
  new_predictions <- predict(p.lm, newdata = dat[i, ], na.action = na.omit)
  pred <- c(pred, new_predictions)
  truth <- c(truth, log(dat$p[i]))
  
}

rmse  <- sqrt(sum(na.omit((pred - truth)^2))/length(truth)) #calculate root mean square error

rmse/(max(truth) - min(truth)) #rmse % of true range

plot(truth, pred)


#------------------------ cross validation drop 1 station with simulated random effects --------------------------#

pred <- NULL
truth <- NULL
for (i in unique(dat$stn)) {
  
  p.lm <- lme(log(p) ~ obs * oxy, random =~ 1 + oxy + obs | stn, data = dat[dat$stn != i, ], na.action = na.exclude, 
              control = list(opt='optim'), weights = varExp(form =~ oxy))
  
  #resample using extracted random effect sds to simulate random effects
  #manually create predictions
  total <- NULL
  for (j in 1:1000) {
    
    #intercept
    intercept <- p.lm$coefficients$fixed[1] + rnorm(1, 0, sd = as.numeric(VarCorr(p.lm)[1, 2]))
    
    #oxy
    oxygen <- (p.lm$coefficients$fixed[3] + rnorm(1, 0, sd = as.numeric(VarCorr(p.lm)[2, 2])))*dat$oxy[dat$stn == i]
    
    #obs
    obs <- (p.lm$coefficients$fixed[2] + rnorm(1, 0, sd = as.numeric(VarCorr(p.lm)[3, 2])))*dat$obs[dat$stn == i]
    
    #interaction
    interaction <- p.lm$coefficients$fixed[4]*(dat$oxy*dat$obs)[dat$stn == i]
    
    combined <- exp(intercept + oxygen + obs + interaction)
    
    combined[combined > 500] <- NA #catch unreasonably large values so they don't skew resample estimates
    
    total <- cbind(total, combined)
    
  }
  
  pred  <- c(pred, rowMeans(total, na.rm = TRUE))
  truth <- c(truth, dat$p[dat$stn == i])
  
}

rmse  <- sqrt(sum(na.omit((pred - truth)^2))/length(truth)) #calculate root mean square error

rmse/(max(truth) - min(truth)) #rmse % of true range

plot(log(truth), log(pred))

