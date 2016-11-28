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
library(rgl)
library(colorRamps)  
library(chron)
library(boot) #inv.logit

#source required functions
function_list <- c("setUpFluoro.R",
                   "calc_conditional_marginal_Rsquared.R",
                   "calc_asreml_conditional_marginal_Rsquared.R",
                   "rocCurve.R",
                   "krillPresenceAbsence.R",
                   "crossValROC.R")

for (f in function_list) {
  source(paste("R code/R-functions-southern-ocean/", f, sep = ""))
}

#find stations where krill data is available and subset glm.spl
dat <- read.csv(file = "Data/procCTD.csv", header= T)
names(dat) <- c("survey", "stn", "lat", "long", "start.time", "end.time", "depth", "transmittance", "cond", "temp", "sal", "par", "oxygen", "fluoro", "x2", "ice", "wm")

ctd_time <- read.csv("Data/ctd_times.csv", header = T)

glm.spl <- setUpFluoro(dat, scale = FALSE)
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

d$time <- NA
d$day <- NA
for (i in 1:nrow(d)) {
  d$time[i] <- chron(times. = ctd_time$start.time[ctd_time$stn == d$stn[i]], format = "h:m:s")
  d$day[i]  <- ctd_time$julian_day[ctd_time$stn == d$stn[i]]
}

d$time <- chron(times. = d$time , format = "h:m:s")
d$hour <- hours(d$time)

#d$depth <- stn_coords$depth[match(d$stn, stn_coords$Cast.Number)]

d$log_p <- log(d$p)
d$log_p[is.infinite(d$log_p)] <- NA

d$obs[d$obs < 0] <- NA

d$l.obs <- log(d$obs)
d$l.obs[is.infinite(d$l.obs)] <- NA

unscaled <- d
#-------------------- BINOMIAL GLM FOR PRESENCE/ABSENCE -----------------------#
 
#glm
pa.lm <- glm(pa ~ z + temp + sal + par, dat = d, family = "binomial")
summary(pa.lm)

#calculate variance inflation factors (<5 = good)
vif(pa.lm)

#scale or model doesn't converge
d <- cbind(d[, c(1, 7:8)], apply(d[, c(2:6, 9, 10, 14)], 2, scale))

#mixed model with station random effect
pa.lm <- glmer(pa ~ z + temp + sal + par +(1|stn), data = d, family = "binomial")
summary(pa.lm)

#calculate sensitivity and specificity
sensitivity(as.factor(round(fitted(pa.lm))), as.factor(d$pa), positive = "1", negative = "0")
specificity(as.factor(round(fitted(pa.lm))), as.factor(d$pa), positive = "1", negative = "0")

#plot a ROC curve for the binomial glm
M.ROC <- rocCurve(model = pa.lm, threshold = 0.5, data = pa, print = TRUE)

par(mfrow = c(1, 1), mar = c(5, 5, 1, 1))
plot(M.ROC[1, ], M.ROC[2, ], lwd = 2, type = "l", xlab = "False Positive Rate", ylab = "True Positive Rate", cex.lab = 2, cex.axis = 2)
title("ROC curve")
lines(c(0, 1), c(0, 1), col = "red")

#calculate the area under the ROC curve (0.5 = bad, 0.8 = good, 0.9 = excellent, 1 = perfect)
auc(M.ROC[1,], M.ROC[2,])

#partial plots
pdf("C:/Users/43439535/Dropbox/uni/hurdle paper/figures/fig_1.pdf", width = 10, height = 9)
par(mar = c(4.1,4.1,3.1,2.1), mfrow = c(2, 2), lwd = 2)
par(oma = c(3, 6, 0, 0))

#depth
predict_pa_re <- expand.grid(seq(min(d$z), max(d$z), length.out = 100), unique(d$stn))
predict_pa <- data.frame("z" = predict_pa_re$Var1, "stn" = predict_pa_re$Var2, "temp" = 0, "sal" = 0, "par" = 0, "day" = 0)
pred_z <- predict(pa.lm, newdata = predict_pa, allow.new.level = T, type = "response")
pred_z <- aggregate(pred_z, list(predict_pa$z), FUN = mean)
plot(pred_z$Group.1 * sd(unscaled$z) + mean(unscaled$z), pred_z$x, type = "l", xlab = "Depth (m)", ylab = "", ylim = c(0, 1), bty = "l", cex.lab = 1.5, cex.axis = 1.5)
legend("topleft", "(a)", bty = "n", cex = 1.5, x.intersp = 0, y.intersp = 0)
mm <- model.matrix(~z, predict_pa)
y <- mm%*%fixef(pa.lm)[1:2]
pvar1 <- diag(mm %*% tcrossprod(vcov(pa.lm)[1:2, 1:2], mm))
tlo = inv.logit(y - 1.96*sqrt(pvar1))
thi = inv.logit(y + 1.96*sqrt(pvar1))
tlo <- aggregate(tlo, list(predict_pa$z), FUN = mean)
thi <- aggregate(thi, list(predict_pa$z), FUN = mean)
points(tlo$Group.1  * sd(unscaled$z) + mean(unscaled$z), tlo$V1, lty = 2, type = "l")
points(thi$Group.1  * sd(unscaled$z) + mean(unscaled$z), thi$V1, lty = 2, type = "l")


#temperature
predict_pa_re <- expand.grid(seq(min(d$temp), max(d$temp), length.out = 100), unique(d$stn))
predict_pa <- data.frame("z" = 0, "stn" = predict_pa_re$Var2, "temp" = predict_pa_re$Var1, "sal" = 0, "par" = 0, "day" = 0)
pred_temp <- predict(pa.lm, newdata = predict_pa, allow.new.level = T, type = "response")
pred_temp <- aggregate(pred_temp, list(predict_pa$temp), FUN = mean)
plot(pred_temp$Group.1 * sd(unscaled$temp) + mean(unscaled$temp), pred_temp$x, ylim = c(0, 1),type = "l", cex.lab = 1.5, xlab = expression(Temperature~(~degree~C)), ylab = "", bty = "l", cex.axis = 1.5)
legend("topleft", "(b)", bty = "n", cex = 1.5, x.intersp = 0, y.intersp = 0)
mm <- model.matrix(~temp, predict_pa)
y <- mm%*%fixef(pa.lm)[c(1, 3)]
pvar1 <- diag(mm %*% tcrossprod(vcov(pa.lm)[c(1, 3), c(1, 3)], mm))
tlo = inv.logit(y - 1.96*sqrt(pvar1))
thi = inv.logit(y + 1.96*sqrt(pvar1))
tlo <- aggregate(tlo, list(predict_pa$temp), FUN = mean)
thi <- aggregate(thi, list(predict_pa$temp), FUN = mean)
points(tlo$Group.1  * sd(unscaled$temp) + mean(unscaled$temp), tlo$V1, lty = 2, type = "l")
points(thi$Group.1  * sd(unscaled$temp) + mean(unscaled$temp), thi$V1, lty = 2, type = "l")


#salinity
predict_pa_re <- expand.grid(seq(min(d$sal), max(d$sal), length.out = 100), unique(d$stn))
predict_pa <- data.frame("z" = 0, "stn" = predict_pa_re$Var2, "temp" = 0, "sal" = predict_pa_re$Var1, "par" = 0, "day" = 0)
pred_sal <- predict(pa.lm, newdata = predict_pa, allow.new.level = T, type = "response")
pred_sal <- aggregate(pred_sal, list(predict_pa$sal), FUN = mean)
plot(pred_sal$Group.1 * sd(unscaled$sal) + mean(unscaled$sal), pred_sal$x, ylim = c(0, 1), type = "l", xlab = "Salinity (ppm)", ylab = "", bty = "l", cex.lab = 1.5, cex.axis = 1.5)
legend("topleft", "(c)", bty = "n", cex = 1.5, x.intersp = 0, y.intersp = 0)
mm <- model.matrix(~sal, predict_pa)
y <- mm%*%fixef(pa.lm)[c(1, 4)]
pvar1 <- diag(mm %*% tcrossprod(vcov(pa.lm)[c(1, 4), c(1, 4)], mm))
tlo = inv.logit(y - 1.96*sqrt(pvar1))
thi = inv.logit(y + 1.96*sqrt(pvar1))
tlo <- aggregate(tlo, list(predict_pa$sal), FUN = mean)
thi <- aggregate(thi, list(predict_pa$sal), FUN = mean)
points(tlo$Group.1  * sd(unscaled$sal) + mean(unscaled$sal), tlo$V1, lty = 2, type = "l")
points(thi$Group.1  * sd(unscaled$sal) + mean(unscaled$sal), thi$V1, lty = 2, type = "l")

#par
predict_pa_re <- expand.grid(seq(min(d$par), max(d$par), length.out = 100), unique(d$stn))
predict_pa <- data.frame("z" = 0, "stn" = predict_pa_re$Var2, "temp" = 0, "sal" = 0, "par" = predict_pa_re$Var1, "day" = 0)
pred_par <- predict(pa.lm, newdata = predict_pa, allow.new.level = T, type = "response")
pred_par <- aggregate(pred_par, list(predict_pa$par), FUN = mean)
plot(pred_par$Group.1 * sd(unscaled$par) + mean(unscaled$par), pred_par$x, ylim = c(0, 1), type = "l",xlab = expression("PAR" ~ (mu~E ~ m^{-2} ~ s^{-1})), ylab = "", bty = "l", cex.lab = 1.5, cex.axis = 1.5)
legend("topleft", "(d)", bty = "n", cex = 1.5, x.intersp = 0, y.intersp = 0)
mm <- model.matrix(~par, predict_pa)
y <- mm%*%fixef(pa.lm)[c(1, 5)]
pvar1 <- diag(mm %*% tcrossprod(vcov(pa.lm)[c(1, 5), c(1, 5)], mm))
tlo = inv.logit(y - 1.96*sqrt(pvar1))
thi = inv.logit(y + 1.96*sqrt(pvar1))
tlo <- aggregate(tlo, list(predict_pa$par), FUN = mean)
thi <- aggregate(thi, list(predict_pa$par), FUN = mean)
points(tlo$Group.1  * sd(unscaled$par) + mean(unscaled$par), tlo$V1, lty = 2, type = "l")
points(thi$Group.1  * sd(unscaled$par) + mean(unscaled$par), thi$V1, lty = 2, type = "l")


mtext("Probability of krill presence", side = 2, outer = TRUE, line = 2, cex = 2)

dev.off()


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

#cross validation ROC
crossValROC(pred, truth)

pred <- round(pred)

table(pred, truth)

sensitivity(data = as.factor(pred), reference = as.factor(truth), positive = "1", negative = "0") #true positive rate
specificity(data = as.factor(pred), reference = as.factor(truth), positive = "1", negative = "0") #true negative rate


#-------------------------- KRILL VS PHYTOPLANKTON ----------------------------#

#--------------- fluoro and oxy with linear relationship and interaction term ---------------#

d <- data.frame(cbind(pa, fluoro$oxy, fluoro$sal, fluoro$z, fluoro$par, fluoro$temp, p, fluoro$stn, fluoro$obs))
colnames(d) <- c("pa", "oxy", "sal", "z", "par", "temp", "p", "stn", "obs")
d <- na.omit(d)


d$time <- NA
d$day <- NA
for (i in 1:nrow(d)) {
  d$time[i] <- chron(times. = ctd_time$start.time[ctd_time$stn == d$stn[i]], format = "h:m:s")
  d$day[i] <- ctd_time$julian_day[ctd_time$stn == d$stn[i]]
}

d$time <- chron(times. = d$time , format = "h:m:s")
d$hour <- hours(d$time)

dat <- d[d$pa == 1 & round(fitted(pa.lm)) == 1, ]
dat <- dat[dat$stn %in% sort(unique(dat$stn))[which(table(dat$stn) >= 5)], ]
dat$stn <- as.factor(dat$stn)

dat$obs[dat$obs < 0] <- NA
dat$l.obs <- log(dat$obs)
dat$l.obs[is.infinite(dat$l.obs)] <- NA
dat_unscaled <- dat
dat_unscaled$obs[dat_unscaled$obs < 0] <- NA

dat$oxy <- scale(dat$oxy)
dat$obs <- scale(dat$obs)
dat$l.obs <- scale(dat$l.obs)

p.lm <- lme(log(p) ~ l.obs * oxy, random =~ 1 | stn, data = dat, na.action = na.omit,
            control = list(opt='optim'))
summary(p.lm)
r.squared.lme(p.lm)


#deviance explained
model_deviance <- -2*p.lm$logLik

p.null <- lm(log(p) ~ 1, data = dat, na.action = na.omit)
null_deviance <- -2*logLik(p.null)[1]

(null_deviance - model_deviance)/null_deviance


#3D plot of obs*oxy interaction
interaction_data <- expand.grid(seq(min(dat$oxy), max(dat$oxy), length.out = 200), seq(min(na.omit(dat$l.obs)), max(na.omit(dat$l.obs)), length.out = 200), unique(dat$stn))
colnames(interaction_data) <- c("oxy", "l.obs", "stn")
pred_interaction <- predict(p.lm, newdata = interaction_data)

pred_fixed <- aggregate(exp(pred_interaction), list(interaction_data$l.obs, interaction_data$oxy), FUN = mean)
pred_fixed <- na.omit(pred_fixed)

plot_dat <- data.frame("x" = exp(pred_fixed$Group.1*sd(na.omit(dat_unscaled$l.obs)) + mean(na.omit(dat_unscaled$l.obs))), "y" = pred_fixed$Group.2*sd(dat_unscaled$oxy) + mean(dat_unscaled$oxy), "z" = pred_fixed$x)

#interactive dot plot
plot3d(plot_dat$x, plot_dat$y, plot_dat$z, xlab = "Phytoplankton Fluoresence", ylab = "Dissolved Oxygen", zlab = "Krill density (g/m2)")

#static plot for paper

wireframe(z ~ x * y, data = plot_dat, xlab = expression("Phytoplankton" ~ (mu~g ~ L^{-1})), ylab = expression("Dissolved oxygen" ~ (mu~mol ~ L^{-1})), zlab = expression("Krill density"~(gm^-2)),
          perspective = FALSE, colorkey = FALSE, scales = list(arrows=FALSE,tick.number = 10, x = list(distance = 1.5), y = list(distance = 1.5), col = "black"),
          drape = T,  col.regions = colorRampPalette( c("lightblue", "darkblue"))(100), col = "transparent", par.settings = list(axis.line = list(col = 'transparent')))


#interactive surface plot
jet.colors <- colorRampPalette(matlab.like(50))
colorjet <- jet.colors(5)
open3d()
rgl.surface(x=unique(pred_fixed$Group.1)*10, z=unique(pred_fixed$Group.2)*10, y=pred_fixed$x, 
            color=colorjet[ findInterval(exp(pred_fixed$x), seq(min(na.omit(exp(pred_fixed$x))), max(na.omit(exp(pred_fixed$x))), length=100))])
axes3d()
title3d(xlab = "Phytoplankton Fluoresence", zlab = "Dissolved Oxygen", ylab = "Krill density (g/m2)")


#extract fitted including only fixed effects
y <- p.lm$fitted[, 1]

#plot residuals against fitted and covariates
par(mfrow = c(1, 5))
plot(fitted(p.lm), resid(p.lm, type = "normalized"))
plot(na.omit(dat)$l.obs, resid(p.lm, type = "normalized"))
plot(na.omit(dat)$oxy, resid(p.lm, type = "normalized"))
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
  points(exp(xy[, 1]), exp(xy[, 2]), type = "l")
}
points(exp(xy1[, 1]), exp(xy1[, 2]), col = "red", type = "l", lwd = 4)



#------------------------ cross validation drop 1 station with simulated random effects --------------------------#

pred <- NULL
truth <- NULL
for (i in unique(dat$stn[dat$stn != 47])) {
  
  p.lm <- lme(log(p) ~ l.obs * oxy, random =~ 1 + oxy + l.obs | stn, data = dat, na.action = na.omit, 
              control = list(opt='optim'))
  
  #resample using extracted random effect sds to simulate random effects
  #manually create predictions
  total <- NULL
  for (j in 1:1000) {
    
    #intercept
    intercept <- p.lm$coefficients$fixed[1] + rnorm(1, 0, sd = as.numeric(VarCorr(p.lm)[1, 2]))
    
    #oxy
    oxygen <- (p.lm$coefficients$fixed[3] + rnorm(1, 0, sd = as.numeric(VarCorr(p.lm)[2, 2])))*dat$oxy[dat$stn == i]
    
    #obs
    l.obs <- (p.lm$coefficients$fixed[2] + rnorm(1, 0, sd = as.numeric(VarCorr(p.lm)[3, 2])))*dat$l.obs[dat$stn == i]
    
    #interaction
    interaction <- p.lm$coefficients$fixed[4]*(dat$oxy*dat$l.obs)[dat$stn == i]
    
    combined <- exp(intercept + oxygen + l.obs + interaction)
    
    combined[combined > 500] <- NA #catch unreasonably large values so they don't skew resample estimates
    
    total <- cbind(total, combined)
    
  }
  
  pred  <- c(pred, rowMeans(total, na.rm = TRUE))
  truth <- c(truth, dat$p[dat$stn == i])
  
}

rmse  <- sqrt(sum(na.omit((pred - truth)^2))/length(truth)) #calculate root mean square error

rmse/(max(truth) - min(truth)) #rmse % of true range

plot(log(truth), log(pred))

