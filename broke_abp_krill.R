# BROKE air-breathing predator preliminary analysis
#integration interval used for acoustic data is 10x50m
#Lisa-Marie Harrison
#02/09/2014

setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/")
dat   <- read.csv(file = "BROKE-West/Air breathing predator/broke_baleen.csv", header = T, fill = T)
track <- read.csv(file  = "BROKE-West/Echoview/integrated data/broke_cruise_track.csv", header = T)
library(chron)

#remove null rows from acoustic files
files <- list.files("C:/Users/Lisa/Documents/phd/southern ocean/BROKE-West/Echoview/Extracted data/250x2000 integration", full.names = T)
for (i in files) {
  
  dat <- read.csv(i, header = T)
  w <- !dat$ï..Region_ID == -9999
  dat <- dat[w, ]
  
  write.csv(dat, i, row.names = F)
  
}

#read all acoustic data files and combine into one
acoustic_38 <- rep(0, 85)
files_38 <- list.files("C:/Users/Lisa/Documents/phd/southern ocean/BROKE-West/Echoview/Extracted data/250x2000 integration", full.names = T, pattern = "38kHz")
for (i in files_38) {
  dat <- read.csv(i, header = T)
  acoustic_38 <- rbind(acoustic_38, dat)  
}
acoustic_38 <- acoustic_38[-1, ]
colnames(acoustic_38) <- colnames(dat)

acoustic_120 <- rep(0, 85)
files_120 <- list.files("C:/Users/Lisa/Documents/phd/southern ocean/BROKE-West/Echoview/Extracted data/250x2000 integration", full.names = T, pattern = "120kHz")
for (i in files_120) {
  dat <- read.csv(i, header = T)
  acoustic_120 <- rbind(acoustic_120, dat)  
}
acoustic_120 <- acoustic_120[-1, ]
colnames(acoustic_120) <- colnames(dat)


#calculate 120kHz - 38kHz for each 10x50 window
sv_38 <- acoustic_38$Sv_mean
sv_120 <- acoustic_120$Sv_mean
sv_38[sv_38 > 500 | sv_38 < -500] <- NA
sv_120[sv_120 > 500 | sv_120 < -100] <- NA
noise <- is.na(sv_120)
sv_diff <- sv_120 - sv_38


#remove 120 - 38 kHz values outside of [1.02, 14.75] because these are unlikely to be krill
#dB difference window is from Potts AAD report for KAOS data
sv_diff[sv_diff < 2.5 | sv_diff > 14.7] <- NA
sv_120[is.na(sv_diff)] <- NA

sv <- 10^(sv_120/10)

mvbs <- 10*log10(sv)
mvbs[mvbs == -Inf] <- NA

#convert to density using target strength (kg/m2 per interval)
p <- 250*10 ^((mvbs - -42.22)/10)*1000

#add zero krill intervals back in
p[is.na(p)] <- 0

#remove noise values
p[p > 5000] <- NA

#calculate interval length (m)
interval_length <- 0
for (k in 1:nrow(acoustic_38)) {
  interval_length[k] <- (acoustic_38$Dist_E[k] - acoustic_38$Dist_S[k])
}

#calculate interval weighting
interval_weight <- interval_length/sum(na.omit(interval_length))

#calculate mean weighted density
survey_mean <- sum(na.omit(p*interval_weight))




#plot locations - full data set
plot(dat$longitude, dat$latitude, xlab = "Longitude", ylab = "Latitude")
title("Locations of air-breathing predator sightings - full data set")

#check time range of sightings
dat$time <- chron(times. = dat$time, format = "h:m:s")
acoustic_38$Time_S <- chron(times. = acoustic_38$Time_S, format = "h:m:s")
acoustic_38$Time_E <- chron(times. = acoustic_38$Time_E, format = "h:m:s")
hist(dat$time)

#subset data to only include times that overlap with the acoustic survey
dat$date <- chron(dates. = as.character(dat$date), format = "d/m/y")
acoustic_38$Date_S <- chron(dates. = as.character(acoustic_38$Date_S), format = "ymd", out.format = "d/m/y")
dat_sub <- dat[dat$date %in% acoustic_38$Date_S, ]


#table of species seen during acoustic survey
sp.seen <- table(dat_sub$species, dat_sub$best.group.size)
sightings <- unname(apply(sp.seen, 1, function(x)sum(x != 0)))
ind <- unname(rowSums(sp.seen))
sp <- cbind(sightings, ind)
rownames(sp) <- levels(dat_sub$species)
sp <- sp[sightings != 0, ]
colnames(sp) <- c("Sightings", "Individuals")


#plot of track with sighting locations superimposed
plot(acoustic_38$Lon_S[acoustic_38$Lat_S < 900], acoustic_38$Lat_S[acoustic_38$Lat_S < 900], 
     xlab = "Longitude", ylab = "Latitude")
points(dat_sub$longitude, dat_sub$latitude, col = "red", pch = 19)
title("Cruise track (black) with sighting locations superimposed (red)")



#-------------------------- AIR BREATHING PREDATOR ----------------------------#

#assign predator sightings to intervals
#pred is a 0 1 vector of whether there any sightings during an interval
#multiple sightings are not reported, only presence/absence of sightings
#because intervals only contain 1 ping, the predator is assigned to the closest ping
pred <- rep(NA, nrow(acoustic_38))
for (i in 1:length(dat_sub$date)) {

  w <- which(acoustic_38$Date_S == dat_sub$date[i] & acoustic_38$Time_S <= dat_sub$time[i] & acoustic_38$Time_E >= dat_sub$time[i])
  if (is.integer(w)) {
    pred[w] <- 1
  }
}



#plot density (summed through all depths) for each time interval
#times of predator observations are shown with a red vertical line
plot(p, ylab = "krill density (g/m2)")
abline(v = which(pred > 0), col = "red")
title("Krill interval density with predator sightings in red")


#assign predator sightings to intervals
#pred_num is a vector of the number of sightings during an interval
individuals <- rep(0, nrow(acoustic_38))
for (i in 1:length(dat_sub$date)) {
  
  w <- which(acoustic_38$Date_S == dat_sub$date[i] & acoustic_38$Time_S <= dat_sub$time[i] & acoustic_38$Time_E >= dat_sub$time[i])
  if (is.integer(w)) {
    individuals[w] <- individuals[w] + dat_sub$best.group.size[i]
  }
}

#plot of krill biomass against number of sightings at that time
plot(p, individuals, xlab = "krill density", ylab = "number of sightings")
title("krill interval density against the number of predator in that interval")


#--------------------------- SCHOOLS DETECTION --------------------------------#

#using EV schools detection
start.date <- sc$Date_S
start.time <- sc$Time_S
end.date   <- sc$Date_E
end.time   <- sc$Time_E
height     <- sc$Height_mean
depth      <- sc$Depth_mean

#each integration interval is marked as 0 if no krill present or 1 if krill present
#inaccurate because the exported schools detection only gives average location for the school
krill <- matrix(0, ncol = length(int_td)/25, nrow = 25) 
s <- chron(times. = as.character(start.time), format = "h:m:s")
e <- chron(times. = as.character(end.time), format = "h:m:s")
for (i in 1:length(start.date)) {
  w <- which(int_date == start.date[i] & as.character(int_time) >= s[i] & as.character(int_time) <= e[i] & 
               int_depth >= (depth[i] - height[i]/2) & int_depth <= (depth[i] + height[i]/2))
  krill[w] <- 1
}

#plot krill biomass with predator locations in black
p_ext_pred <- matrix(rep(pred, 25), nrow = 25, byrow = T) #add 1s to first row at intervals where there are predators
p_ext_pred[p_ext_pred == 0 ] <- NA
image(t(krill)[,nrow(krill):1])
image(t(p_ext_pred)[, nrow(p_ext_pred):1], add = T, col = "black")


#---------------------- MODELS FOR PREDATOR PRESENCE --------------------------#


#glm for biomass ~ predator presence/absence
biomass.glm <- glm(pred ~ p, family = binomial)
summary(biomass.glm)





