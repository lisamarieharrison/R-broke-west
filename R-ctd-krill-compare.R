#compares krill profiles to ctd data using density calculated from 2x25m extracted integration intervals
#author: Lisa-Marie Harrison
#date: 13/08/2015

setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/BROKE-West")
density <- read.csv("brokewest_krill_ctd.csv", header = T)
library(lattice)

#get glm.spl
source("C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/R code/R-broke-west/R-set-up-fluoro.R")

#find stations where krill data is available and subset glm.spl
stn <- unique(density$stn)
fluoro <- glm.spl[glm.spl$stn %in% stn, ]

#latplot of krill at each station
krill_stn <- rep(stn, each = 125)
lat.plot <- xyplot(density$p ~ density$depth | density$stn, xlab = "depth (m)",
                   ylab = "krill density (g/m2)", outer = FALSE, type = "l")
update(lat.plot)

#histogram of krill density
hist(density$p, main = "Histogram of krill density", xlab = "krill density (g/m2)")

#transform krill to be on the same scale as ctd data
p <- 0
for (i in 1:nrow(fluoro)) {
  p[i] <- density$p[density$stn == fluoro$stn[i]][which.min(abs(density$depth[density$stn == fluoro$stn[i]] - fluoro$z[i]))]
}

#plot krill against fluoro
plot(exp(fluoro$l.obs), p, ylab = "krill density (g/m2)", pch = 19,
     xlab = "phytoplankton fluorescence", main = "Krill vs Phytoplankton")

#plot log krill against environmental variables
plot(cbind(fluoro[c(2, 6:9)], log(p)))

#plot krill presence/absence against environmental variables
pa <- rep(0, length(p))
pa[p > 0] <- 1
plot(cbind(fluoro[c(2, 6:9)], pa))


par(mfrow = c(1, 3))
boxplot(fluoro$oxy ~ pa)
boxplot(fluoro$sal ~ pa)
boxplot(fluoro$z ~ pa)

d <- data.frame(cbind(pa, fluoro$oxy, fluoro$sal, fluoro$z))
colnames(d) <- c("pa", "oxy", "sal", "z")

pa.lm <- glm(pa ~ oxy + sal + z, dat = d, family = "binomial")
summary(pa.lm)






