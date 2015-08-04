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

#latplot of krill at each station
krill_stn <- rep(stn, each = 125)
lat.plot <- xyplot(density ~ rep(seq(2, 250, by = 2), 104) | krill_stn, xlab = "depth (m)",
                   ylab = "krill density (g/m2)", outer = FALSE, type = "l")
update(lat.plot)

#histogram of krill density
hist(density, main = "Histogram of krill density", xlab = "krill density (g/m2)")

#plot krill against fluoro
plot(exp(fluoro$l.obs), density, ylab = "krill density (g/m2)", pch = 19,
     xlab = "phytoplankton fluorescence")
title("Krill vs Phytoplankton")

#plot log krill against log fluoro
plot(log(density), fluoro$l.obs, xlab = "log krill density (g/m2)", pch = 19,
     ylab = "log phytoplankton fluorescence")
title("log Krill vs log Phytoplankton")

#optional kernal smoothing of krill density
density <- ksmooth(c(1:length(density)), density, bandwidth = 11)$y


