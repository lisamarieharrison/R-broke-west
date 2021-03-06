#extract krill acoustic data around ctd station for kaos data
#author: Lisa-Marie Harrison
#date: 18/06/2015


setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/BROKE-West/")
ctd <- read.csv("C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/Data/ctd_information.csv", header = T)
ev_times <- read.csv("C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/Data/ev_file_times.csv", header = T)
library(chron)

file.create("brokewest_krill_ctd.csv")

for (i in 1:length(unique(ctd$stn))) {
  
  stn <- unique(ctd$stn)[i]
  day <- ctd$day[ctd$stn == stn]
  if (nchar(day) == 1) {
    day <- paste(0, day, sep = "")
  }
  date <- paste(ctd$year[ctd$stn == stn], 0 , ctd$month[ctd$stn == stn], day, sep = "")
  start_time <- ctd$start_time[i]
  
  while (nchar(start_time) < 6) {
    start_time <- paste(0, start_time, sep = "")
  }
  
  time <- chron(times. = as.character(start_time), format = "hms")
  
  f <- which(ev_times$start_date == date & chron(times. = ev_times$start_time) < time & chron(times. = ev_times$end_time) > (time - 0.05))
    
  krill_file <- list.files("C:/Users/Lisa/Documents/phd/southern ocean/BROKE-West/Echoview/Extracted data/10x50 integration/", full.names = T)[c(f*2, f*2 - 1)]
  if (length(krill_file) != 0) {
    krill_120 <- read.csv(krill_file[2], header = T)
    krill_38 <- read.csv(krill_file[1], header = T)
  } else {
    next()
  }
  
  start_time <- as.character(unique(ctd$start_time[ctd$stn == stn]))
  if (nchar(start_time) != 6) {
    start_time <- paste(0, start_time, sep = "")
  }
  
  start_time <- chron(times. = start_time, format = "hms")
  
  bottom_time <- as.character(unique(ctd$bottom_time[ctd$stn == stn]))
  if (nchar(bottom_time) != 6) {
    bottom_time <- paste(0, bottom_time, sep = "")
  }
  
  bottom_time <- chron(times. = bottom_time, format = "hms")
  
  
  krill_38$Time_S <- chron(times. = krill_38$Time_S, format = "h:m:s")
  krill_38$Time_E <- chron(times. = krill_38$Time_E, format = "h:m:s")
  krill_120$Time_S <- chron(times. = krill_120$Time_S, format = "h:m:s")
  krill_120$Time_E <- chron(times. = krill_120$Time_E, format = "h:m:s")
  
  krill_38 <- krill_38[krill_38$Time_S > (start_time - 0.025) & krill_38$Time_E < start_time & krill_38$Date_E == date, ]  
  krill_120 <- krill_120[krill_120$Time_S > (start_time - 0.025) & krill_120$Time_E < start_time & krill_120$Date_E == date, ]  
  
  #just get the closest interval
  krill_38 <- krill_38[abs(krill_38$Sv_mean) != 9999 & is.na(krill_38$Sv_mean) == FALSE, ]
  krill_120 <- krill_120[abs(krill_120$Sv_mean) != 9999 & is.na(krill_38$Sv_mean) == FALSE, ]
  
  
  if (nrow(krill_38) == 0) {
    next()
  }
  
  krill_38 <- krill_38[krill_38$Interval %in% c((max(krill_38$Interval) - 100): max(krill_38$Interval)), ]
  krill_120 <- krill_120[krill_120$Interval %in% c((max(krill_38$Interval) - 100): max(krill_120$Interval)), ]

  
  #remove layers of -1 
  krill_38 <- krill_38[krill_38$Layer > 0, ]
  krill_120 <- krill_120[krill_120$Layer > 0, ]  
  
  if (length(unique(krill_38$Layer)) <= 1) {
    next()
  }
  
  #calculate 120kHz - 38kHz for each 10x50 window
  sv_38 <- krill_38$Sv_mean
  sv_120 <- krill_120$Sv_mean
  sv_38[sv_38 > 500 | sv_38 < -500] <- NA
  sv_120[sv_120 > 500 | sv_120 < -90] <- NA  
  
  sv <- 10^(sv_120/10)
  mvbs_120 <- 10*log10(aggregate(sv, by = list(krill_38$Layer), FUN = "mean", na.rm = TRUE)$x)
  
  sv1_38 <- 10^(sv_38/10)
  mvbs_38 <- 10*log10(aggregate(sv1_38, by = list(krill_38$Layer), FUN = "mean", na.rm = TRUE)$x)
  
  sv_diff <- mvbs_120 - mvbs_38
  
  mvbs_120[sv_diff < 2.5 | sv_diff > 14.7] <- NA
  p <- 10 ^((mvbs_120 - -42.22)/10)*1000*10
  p[is.na(mvbs_120)] <- 0
  
  
  #   #remove 120 - 38 kHz values outside of [1.02, 14.75] because these are unlikely to be krill
  #   #dB difference window is from Potts AAD report for KAOS data
  #   sv_diff[sv_diff < 1.02 | sv_diff > 14.75] <- NA
  #   sv_120[is.na(sv_diff)] <- NA
  #   
  #   #convert to density using target strength (kg/m2 per interval) and average across intervals
  #   p <- 2*10 ^((sv_120 - -42.22)/10)*1000
  #   p[is.na(p)] <- 0
  #   p <- aggregate(p, by = list(krill_38$Layer), FUN = "mean")$x  
  
  depth <- round(aggregate(krill_38$Depth_mean, by = list(krill_38$Layer), FUN = "mean"))$x
  
  p[p == -Inf] <- 0
  
  out <- cbind(rep(stn, length(p)), p, depth)
  
  write.table(out, "brokewest_krill_ctd.csv", sep = ",", row.names = F, col.names = F, append = T)  
  
  print(i)
  
}

dat <- read.csv("brokewest_krill_ctd.csv", header = F)
names(dat) <- c("stn", "p", "depth")
write.csv(dat, "brokewest_krill_ctd.csv", row.names = F)
