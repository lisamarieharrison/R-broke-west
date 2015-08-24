#compares krill profiles to phytoplankton fluorescence at each ctd station
#author: Lisa-Marie Harrison
#date: 30/03/2015

setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/BROKE-West/Echoview/Extracted data/ctd/density")
krill_files <- list.files()
density = unlist(lapply(krill_files, read.csv))
density[density > 5000] <- NA #remove noise
library(lattice)

#get glm.spl
source("C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/R code/R-broke-west/R-set-up-fluoro.R")

#find stations where krill data is available and subset glm.spl
stn <- as.numeric(gsub(".*stn_(.*)\\..*","\\1", krill_files, perl=T))
fluoro <- glm.spl[glm.spl$stn %in% stn, ]
fluoro <- fluoro[order(fluoro$stn), ]
density <- density[-c(1:125)] #remove station 1 because no fluoro data

#latplot of krill at each station
krill_stn <- rep(stn, each = 125)
lat.plot <- xyplot(density ~ rep(seq(2, 250, by = 2), 104) | krill_stn, xlab = "depth (m)",
                   ylab = "krill density (g/m2)", outer = FALSE, type = "l")
update(lat.plot)

#histogram of krill density
hist(density, main = "Histogram of krill density", xlab = "krill density (g/m2)")




#plot log krill against environmental variables
plot(cbind(fluoro[c(1:2, 6:9)], log(density)))

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

d <- data.frame(cbind(pa, fluoro$oxy, fluoro$sal, fluoro$z, fluoro$par, fluoro$temp, p, fluoro$stn))
colnames(d) <- c("pa", "oxy", "sal", "z", "par", "temp", "p", "stn")

pa.lm <- glm(pa ~ sal + z + temp, dat = d, family = "binomial")
summary(pa.lm)

#table of false and true 0 and 1
table(na.omit(d)$pa, round(fitted(pa.lm)))

#calculate variance inflation factors (<5 = good)
vif(pa.lm)

#calculate sensitivity and specificity
sensitivity(as.factor(round(fitted(pa.lm))), as.factor(na.omit(d)$pa))
specificity(as.factor(round(fitted(pa.lm))), as.factor(na.omit(d)$pa))

