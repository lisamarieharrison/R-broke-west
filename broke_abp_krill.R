# BROKE air-breathing predator preliminary analysis
#integration interval used for acoustic data is 10x50m
#Lisa-Marie Harrison
#02/09/2014

setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/")
dat   <- read.csv(file = "BROKE-West raw data/Air breathing predator/broke_seabirds.csv", header = T, fill = T)
track <- read.csv(file  = "BROKE-West raw data/Echoview/integrated data/broke_cruise_track.csv", header = T)
acoustic_38  <- read.csv(file = "BROKE-West raw data/Echoview/integrated data/broke_38khz_integration_hrp.csv", header = T)
acoustic_120 <- read.csv(file = "BROKE-West raw data/Echoview/integrated data/broke_120khz_integration_hrp.csv", header = T)
library(chron)

#convert latitude/longitude from dms to decimal
dat$latitude <- as.numeric(substr(dat$latitude, 1, 2)) + as.numeric(substr(dat$latitude, 5, 8))/60
dat$longitude <- c(as.numeric(substr(dat$longitude[1:54], 1, 3)) + as.numeric(substr(dat$longitude[1:54], 6, 9))/60, as.numeric(substr(dat$longitude[55:length(dat$longitude)], 1, 2)) + as.numeric(substr(dat$longitude[55:length(dat$longitude)], 5, 8))/60)

#plot locations - full data set
plot(long, -lat, xlab = "Longitude", ylab = "Latitude")
title("Locations of air-breathing predator sightings - full data set")

#check time range of sightings
hist(dat$time, xlab = "Time (hhmm)", main = "Histogram of observation times (UTC)")

#subset data to only include times that overlap with the acoustic survey
obs  <- c(1127:1168)
num.obs <- length(obs)
dat_sub <- dat[obs, ]

#table of species seen during acoustic survey
sp.seen <- table(dat_sub$species, dat_sub$count)
sightings <- unname(apply(sp.seen, 1, function(x)sum(x != 0)))
ind <- unname(rowSums(sp.seen))
sp <- cbind(sightings, ind)
rownames(sp) <- levels(dat_sub$species)
sp <- sp[sightings != 0, ]
colnames(sp) <- c("Sightings", "Individuals")


#plot of track with sighting locations superimposed
gps <- unique(track[,c('Latitude','Longitude')]) #filter unique track values
plot(gps$Longitude, gps$Latitude, xlab = "Longitude", ylab = "Latitude")
points(dat_sub$longitude, -dat_sub$latitude, col = "red", pch = 19)
title("Cruise track (black) with sighting locations superimposed (red)")


#calculate 120kHz - 38kHz for each 10x50 window
sv_diff <- acoustic_120$Sv_mean - acoustic_38$Sv_mean
sv_diff[sv_diff < -500 | sv_diff > 500] <- NA #remove no data areas

int_depth <- acoustic_38$Depth_mean
int_time  <- acoustic_38$Time_M
int_date  <- acoustic_38$Date_M
int_td <- paste(int_date, int_time) #date time for each integration interval

hist(sv_diff)
abline(v = c(2, 12), col = "red")
#write.csv(sv_diff, "C:/Users/Lisa/Documents/phd/southern ocean/BROKE-West raw data/Echoview/integrated data/sv_diff.csv", row.names = F)

#remove 120 - 38 kHz values outside of [2, 12] because these are unlikely to be krill
sv_diff[sv_diff < 2 | sv_diff > 12] <- NA
nasc_120 <- acoustic_120$PRC_NASC
nasc_120[sv_diff == NA] <- 0
nasc_120[int_depth < 15] <- NA #remove the first interval because of noise
nasc_120[acoustic_120$Sv_noise == -999] <- NA #remove possible error values

#convert to density using target strength
#use 0.028*sv_120 (Demer & Hewitt 1995) as an approximation of krill lengths
p <- nasc_120*0.028 #average volumetric density/interval
b <- (p*10*50)/1000 #biomass in kg/interval

#sum density through depths for each time point to find biomass in water column at each time
pt <- colSums(matrix(b, nrow = 25), na.rm = T)
plot(pt)


#------------------------------------------------------------------------------#

#find date and middle time during an integration interval for unique intervals
int_d <- unlist(strsplit(unique(int_td), "  "))[seq(1, length(unlist(strsplit(unique(int_td), "  "))), by = 2)]
int_t <- unlist(strsplit(unique(int_td), "  "))[seq(2, length(unlist(strsplit(unique(int_td), "  "))), by = 2)]

#find start and end times for each unique interval (rather than down all depths)
int.start.time <- acoustic_38$Time_S[seq(1, length(acoustic_38$Time_S), by = 25)]
int.end.time   <- acoustic_38$Time_E[seq(1, length(acoustic_38$Time_E), by = 25)]
t_s <- chron(times. = as.character(int.start.time), format = "h:m:s")
t_e <- chron(times. = as.character(int.end.time), format = "h:m:s")

#find time at middle of integration interval and remove : and .
#use for x-axis of plots
t_m <- acoustic_38$Time_M[seq(1, length(acoustic_38$Time_E), by = 25)]
t_m <- gsub("[: -]", "" , t_m, perl=TRUE)
t_m <- gsub("[. -]", "" , t_m, perl=TRUE)
t_axt <- as.numeric(substr(t_m, start = 1, stop = 4))

#-------------------------- AIR BREATHING PREDATOR ----------------------------#

#assign predator sightings to intervals
#pred is a 0 1 vector of whether there any sightings during an interval
#multiple sightings are not reported, only presence/absence of sightings
#because intervals only contain 1 ping, the predator is assigned to the closest ping
pred <- rep(0, length(unique(int_td)))
for (i in 1:length(dat_sub$date)) {

  w <- which(int_d == as.character(dat_sub$date[i]) & t_s >= (chron(times. = dat_sub$t[i], format = "hh:mm:ss", )  - 1/24/60)
             & t_s <= (chron(times. = dat_sub$t[i], format = "hh:mm:ss", )  + 1/24/60))[1]
  pred[w] <- 1
}



#plot density (summed through all depths) for each time interval
#times of predator observations are shown with a red vertical line
plot(pt, xlab = "time of day", ylab = "krill biomass", xaxt = "n")
axis(1, at = seq(1, length(pt), by = 100), labels = t_axt[seq(1, length(pt), by = 100)])
title("Krill biomass summed in top 250m and timing of predator sightings in red")
abline(v = which(pred == 1), col = "red")


#assign predator sightings to intervals
#pred_num is a vector of the number of sightings during an interval
pred_num <- rep(0, length(unique(int_td)))
for (i in 1:length(dat_sub$d)) {
  t_s <- chron(times. = as.character(t_s), format = "h:m:s")
  t_e <- chron(times. = as.character(t_e), format = "h:m:s")
  
  w <- which(int_d == dat_sub$d[i] & t_e >= as.character(dat_sub$t[i])
             & t_s <= as.character(dat_sub$t[i]))[1]
  pred_num[w] <- pred_num[w] + 1
}

#plot of krill biomass against number of sightings at that time
plot(pt, pred_num, xlab = "krill biomass", ylab = "number of sightings")
title("krill biomass (top 250m) against the number of predator sightings at that time")


#assign predator sightings to intervals
#pred_count is a vector of the count during sightings in an interval
pred_count <- rep(0, length(unique(int_td)))
for (i in 1:length(dat_sub$d)) {
  t_s <- chron(times. = as.character(t_s), format = "h:m:s")
  t_e <- chron(times. = as.character(t_e), format = "h:m:s")
  
  w <- which(int_d == dat_sub$d[i] & t_e >= as.character(dat_sub$t[i])
             & t_s <= as.character(dat_sub$t[i]))[1]
  pred_count[w] <- pred_count[w] + dat_sub$Count[i]
}



#plot of krill biomass against number of individuals seen at that time
plot(pt, pred_count, xlab = "krill biomass", ylab = "number of individuals")
title("krill biomass (top 250m) against the number of predator individuals seen at that time")

#same plot as above but with only presences included
plot(pt[pred_count > 0], pred_count[pred_count> 0], xlab = "krill biomass", ylab = "number of individuals")
title("krill biomass (top 250m) against the number of predator individuals seen at that time")


#image of krill biomass along transect
p_mat <- matrix(p, nrow = 25)
image(t(p_mat)[,nrow(p_mat):1], main = "krill biomass in each 10x50m integration interval")


#plot krill biomass with predator locations in black
p_ext_pred <- matrix(rep(pred, 25), nrow = 25, byrow = T) #add 1s to first row at intervals where there are predators
p_ext_pred[p_ext_pred == 0 ] <- NA
image(t(p_mat)[, nrow(p_mat):1], main = "krill biomass with predator locations overlayed")
image(t(p_ext_pred)[, nrow(p_ext_pred):1], add = T, col = "black")

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
biomass.glm <- glm(pred ~ pt, family = binomial)
summary(biomass.glm)

plot(pt, pred, xlab = "krill biomass", ylab = "predator presence/absence")
title("Observed data with fitted values (red line) for biomass GLM pred ~ biomass")
points(pt, fitted(biomass.glm), col = "red", type = "l")

#--------------------- ADDING UNDERWAY FLUORESCENCE DATA ----------------------#

#get chlorophyll a values during survey
chla <- phyto$chl_a
chla.time <- phyto$LocalTime
chla.date <- phyto$DayNum
t <- chron(times. = as.character(chla.time), format = "h:m:s")


chl.int <- rep(0, length(unique(int_td)))
for (i in 1:length(chla)) {
  
  w <- which(dat_sub$d == chla.date[i] & t_s <= t[i] & t_e >= t[i])
  pred_num[w] <- chla[i]
}






