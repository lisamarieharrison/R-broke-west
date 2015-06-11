# BROKE air-breathing predator preliminary analysis
#integration interval used for acoustic data is 10x50m
#Lisa-Marie Harrison
#02/09/2014

setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/")
dat   <- read.csv(file = "BROKE-West/Air breathing predator/broke_baleen.csv", header = T, fill = T)
track <- read.csv(file  = "BROKE-West/Echoview/integrated data/broke_cruise_track.csv", header = T)
library(chron)

#remove null rows from acoustic files
files <- list.files("C:/Users/Lisa/Documents/phd/southern ocean/BROKE-West/Echoview/Extracted data/2x25 integration", full.names = T)
for (i in files) {
  
  dat <- read.csv(i, header = T)
  w <- !dat$Sv_mean == 9999
  dat <- dat[w, ]
  
  write.csv(dat, i, row.names = F)
  
}

transect <- "01" #specify transect number as a character

#read all acoustic data files and combine into one
acoustic_38 <- matrix(0, ncol = 85)
acoustic_120 <- matrix(0, ncol = 85)
files_38 <- list.files("C:/Users/Lisa/Documents/phd/southern ocean/BROKE-West/Echoview/Extracted data/2x25 integration", full.names = T, pattern = paste("Transect", transect, ".*38kHz", sep = ""))
files_120 <- list.files("C:/Users/Lisa/Documents/phd/southern ocean/BROKE-West/Echoview/Extracted data/2x25 integration", full.names = T, pattern = paste("Transect", transect, ".*120kHz", sep = ""))

for (i in 1:length(files_38)) {
  dat_38 <- read.csv(files_38[i], header = T)
  dat_120 <- read.csv(files_120[i], header = T)
  
  dat_38 <- dat_38[dat_38$Interval %in% unique(dat_120$Interval), ]
  dat_120 <- dat_120[dat_120$Interval %in% unique(dat_38$Interval), ]
  
  if(nrow(dat_120) != nrow(dat_38)) {
    
    compare <- table(dat_120$Interval, dat_120$Layer) == table(dat_38$Interval, dat_38$Layer)
    inds <- which(compare == FALSE, arr.ind=TRUE)
    for (j in 1:nrow(inds)) {
      dat_120 <- dat_120[!(dat_120$Interval == unique(rbind(dat_38, dat_120)$Interval)[inds[j, 1]] & dat_120$Layer == unique(dat_38$Layer)[inds[j, 2]]), ]
      dat_38 <- dat_38[!(dat_38$Interval == unique(rbind(dat_38, dat_120)$Interval)[inds[j, 1]] & dat_38$Layer == unique(dat_38$Layer)[inds[j, 2]]), ]
      
    }
  }
  
  colnames(acoustic_38) <- colnames(dat_38)
  acoustic_38 <- rbind(acoustic_38, dat_38)
  colnames(acoustic_120) <- colnames(dat_120)
  acoustic_120 <- rbind(acoustic_120, dat_120) 
}
acoustic_38 <- acoustic_38[-1, ]
acoustic_120 <- acoustic_120[-1, ]

max_layer <- max(acoustic_38$Layer)
change_loc <- which(acoustic_38$Layer == max_layer)
acoustic_38$unique_interval <- rep(1, nrow(acoustic_38))
for (i in 1:(length(change_loc) - 1)) {
  acoustic_38$unique_interval[(change_loc[i] + 1):change_loc[i+1]] <- i + 1
}

#calculate 120kHz - 38kHz for each 10x50 window
sv_38 <- acoustic_38$Sv_mean
sv_120 <- acoustic_120$Sv_mean
sv_38[sv_38 > 500 | sv_38 < -500] <- NA
sv_120[sv_120 > 500 | sv_120 < -80] <- NA
sv_diff <- sv_120 - sv_38


#remove 120 - 38 kHz values outside of [2, 16] because these are unlikely to be krill
sv_diff[sv_diff <= 2 | sv_diff >= 16] <- NA
sv_120[is.na(sv_diff)] <- NA

sv <- 10^(sv_120/10)

mvbs <- 0
for (i in 1:length(unique(acoustic_38$unique_interval))) {
  mvbs[i] = 10*log10(sum(na.omit(sv[acoustic_38$unique_interval == unique(acoustic_38$unique_interval)[i]])))
}
mvbs[mvbs == -Inf] <- NA


abc <- 10 ^((mvbs)/10)*2

deg2rad <- function(deg) {
  #converts degrees to radians
  #input: degree coordinate
  #returns: radian coordinate 
  
  return(deg*pi/180)
}

gcd.hf <- function(lat1, long1, lat2, long2) {
  #calculates distance between two coordinates using the Haversine formula (hf)
  #input: radian latitude and longitude coordinates
  #returns: distance between coordinates in m
  
  R <- 6371 # Earth mean radius [km]
  delta.long <- (deg2rad(long2) - deg2rad(long1))
  delta.lat  <- (deg2rad(lat2) - deg2rad(lat1))
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(min(1, sqrt(a)))
  d = R * c * 1000
  return(d) 
  
}


#calculate interval length (m) and time
interval_length <- 0
interval_time <- 0
int_matrix <- acoustic_38[acoustic_38$Layer == 5, ]
for (k in 1:nrow(int_matrix)) {
  interval_length[k] <- gcd.hf(int_matrix$Lat_S[k], int_matrix$Lon_S[k], int_matrix$Lat_E[k], int_matrix$Lon_E[k])
  time_start <- chron(times. = int_matrix$Time_S[k], format = "h:m:s")
  time_end <- chron(times. = int_matrix$Time_E[k], format = "h:m:s")
  interval_time[k] <- as.character(time_end - time_start)
  if(time_start > time_end) {
    interval_time[k] <- "00:00:00"
  }
}
interval_length[is.nan(interval_length)] <- 0

abc_nm <- 0
j <- 1
n_int <- 0
int_time <- rep("00:00:00", round((sum(interval_length)/2000)))
for (i in 1:(sum(interval_length)/2000)) {
  abc_nm[i]  <- 0
  cumulative_length <- 0
  n_int[i] <- 0
  while (cumulative_length < 2000) {
    abc_nm[i] <- abc_nm[i] + abc[j]
    j <- j + 1
    cumulative_length <- cumulative_length + interval_length[j]
    n_int[i] <- n_int[i] + 1
    int_time[i] <- as.character(chron(times.= int_time[i], format = "h:m:s") + chron(times. = interval_time[j], format = "h:m:s"))
    if (j >= length(abc)) {
      stop()
    }
  }
}

nasc <- (abc_nm/n_int)*4*pi*1852^2
nasc <- nasc[-!chron(times. = int_time, format = "h:m:s") > chron(times. = "01:00:00", format = "h:m:s")]

#method 1
p <- nasc*0.155
p[p > 5000] <- NA
p[is.na(p)] <- 0
mean(p)

#method 2
p <- nasc*0.6697
p[p > 5000] <- NA
p[is.na(p)] <- 0
mean(p)


#------------------------------------------------------------------------------#

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





