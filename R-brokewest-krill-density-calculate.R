# BROKE-West krill density calculation using methods from Jarvis et al 2010
#integration interval used for acoustic data is 2x25ping
#Lisa-Marie Harrison
#12/06/2015

cluster <- read.csv("C:/Users/Lisa/Documents/phd/southern ocean/BROKE-West/Krill target strength/EvToCluster.csv", header = T)
library(chron)
library(plyr)



transect <- "10" #specify transect number as a character

#read all acoustic data files and combine into one
acoustic_38 <- matrix(0, ncol = 86)
acoustic_120 <- matrix(0, ncol = 85)
files_38 <- list.files("C:/Users/Lisa/Documents/phd/southern ocean/BROKE-West/Echoview/Extracted data/2x25 integration", full.names = T, pattern = paste("Transect", transect, ".*38kHz", sep = ""))
files_120 <- list.files("C:/Users/Lisa/Documents/phd/southern ocean/BROKE-West/Echoview/Extracted data/2x25 integration", full.names = T, pattern = paste("Transect", transect, ".*120kHz", sep = ""))

for (i in 1:length(files_38)) {
  dat_38 <- read.csv(files_38[i], header = T)
  dat_120 <- read.csv(files_120[i], header = T)
  
  ev_files <- list.files("C:/Users/Lisa/Documents/phd/southern ocean/BROKE-West/Echoview/Extracted data/2x25 integration", full.names = T, pattern = paste("38kHz", sep = ""))
  
  dat_38$cluster_alloc <- rep(cluster$cluster[which(ev_files == files_38[i])], nrow(dat_38))
  
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

sv <- 10^(sv_120/10)*2*4*pi*1852^2

sv_mat <- as.data.frame(cbind(sv, acoustic_38$unique_interval, acoustic_38$cluster_alloc))
colnames(sv_mat) <- c("sv", "unique_interval", "cluster")
abc <- ddply(sv_mat, "unique_interval", numcolwise(sum), na.rm = TRUE)$sv
abc[abc == -Inf] <- NA
abc_cluster <- ddply(sv_mat, "unique_interval", numcolwise(mean), na.rm = TRUE)$cluster

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
acoustic_38$Lat_S[acoustic_38$Lat_S == 999] <- NA
acoustic_38$Lon_S[acoustic_38$Lon_S == 999] <- NA
acoustic_38$Lat_E[acoustic_38$Lat_E == 999] <- NA
acoustic_38$Lon_E[acoustic_38$Lon_E == 999] <- NA
int_matrix <- acoustic_38[acoustic_38$Layer == min(acoustic_38$Layer), ]
for (k in 1:nrow(int_matrix)) {
  interval_length[k] <- gcd.hf(int_matrix$Lat_S[k], int_matrix$Lon_S[k], int_matrix$Lat_E[k], int_matrix$Lon_E[k])
  time_start <- chron(times. = int_matrix$Time_S[k], format = "h:m:s")
  time_end <- chron(times. = int_matrix$Time_E[k], format = "h:m:s")
  interval_time[k] <- as.character(time_end - time_start)
  if(time_start > time_end) {
    interval_time[k] <- "00:00:00"
  }
}
interval_length[is.nan(interval_length) | is.na(interval_length)] <- 0
abc[interval_length == 0] <- NA
abc[is.na(abc)] <- 0

set_edsu_length <- 2000 #choose the edsu length in m
interval_length[interval_length > set_edsu_length] <- 0


j <- 1
n_int <- 0
int_time <- rep("00:00:00", 5000)
nasc_length <- 0
nasc_cluster <- 0
bin_include <- 0
for (i in 1:5000) {
  cumulative_length <- 0
  n_int[i] <- 0
  nasc_cluster[i] <- NA
  while ((cumulative_length + interval_length[j]) < set_edsu_length) {
    bin_include[j] <- i
    j <- j + 1
    cumulative_length <- cumulative_length + interval_length[j]
    nasc_length[i] <- cumulative_length
    nasc_cluster[i] <- abc_cluster[j]
    n_int[i] <- n_int[i] + 1
    int_time[i] <- as.character(chron(times.= int_time[i], format = "h:m:s") + chron(times. = interval_time[j], format = "h:m:s"))
    if (j > length(abc)) {
      stop()
    }
  }
}
int_time <- int_time[1:i]
int_time[int_time == "NA"] <- "00:00:00"


bin <- as.data.frame(cbind(interval_length, bin_include))
colnames(bin) <- c("interval_length", "bin_include")
nasc_length <- ddply(bin, "bin_include", numcolwise(sum), na.rm = TRUE)$interval_length

nasc_weight <- 0
for (i in unique(bin_include)) {
  nasc_weight[bin_include == i] <- interval_length[bin_include == i]/nasc_length[i]
}

abc_nm <- 0
for (i in 1:length(nasc_length)) {
  abc_nm[i] <- sum(abc[bin_include == i]*nasc_weight[bin_include == i])
}

conversion_factor_1 <- c(0.1587, 0.1548, 0.1516)
conversion_factor_2 <- c(0.6373, 0.6101, 0.7617)

#method 1
nasc <- (abc_nm)*(conversion_factor_1[nasc_cluster])
p <- nasc[-!chron(times. = int_time, format = "h:m:s") > chron(times. = "01:00:00", format = "h:m:s")]
p[is.na(p)] <- 0
p[p > 5000] <- NA
mean(na.omit(p))

#method 2
nasc <- (abc_nm)*(conversion_factor_2[nasc_cluster])
p <- nasc[-!chron(times. = int_time, format = "h:m:s") > chron(times. = "01:00:00", format = "h:m:s")]
p[is.na(p)] <- 0
p[p > 5000] <- NA
mean(na.omit(p))


