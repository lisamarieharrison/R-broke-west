#calculates krill density at each interval for CTD depths
#author: Lisa-Marie Harrison
#date: 26/03/2015

setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/BROKE-West/Echoview/Extracted data/ctd")
library(fields)

for (i in 1:118) {
  
  while(paste("stn_", i, "_extracted_120khz.csv", sep = "") %in% list.files("120kHz/") == FALSE) {
    i <- i + 1
  }
  
  k120 <- read.csv(paste("120kHz/stn_", i, "_extracted_120khz.csv", sep = ""), header = F, skip = 1)
  k38 <- read.csv(paste("38kHz/stn_", i, "_extracted_38khz.csv", sep = ""), header = F, skip = 1)
  
  
  #rename first 12 columns
  names(k120)[1:12] <- c("Ping_index", "Distance_gps", "Distance_vl", "Ping_date", "Ping_time", 
                         "Ping_milliseconds", "Latitude", "Longitude", "Depth_start", "Depth_stop", "Range_start",
                         "Range_stop")
  
  names(k38)[1:12] <- c("Ping_index", "Distance_gps", "Distance_vl", "Ping_date", "Ping_time", 
                        "Ping_milliseconds", "Latitude", "Longitude", "Depth_start", "Depth_stop", "Range_start",
                        "Range_stop")
  
  #for each row, find the krill density
  sv_120 <- k120[, 14:ncol(k120)]
  sv_120[sv_120 > 500 | sv_120 < -80] <- NA
  sv_38 <- k38[, 14:ncol(k38)]
  sv_38[sv_38 > 500 | sv_38 < -500] <- NA
  
  
  
  #bin data into 125 evenly spaced bins
  meanBins <- function(x) {
    
    y <- stats.bin(x = 1:484, y = x, breaks = seq(1, 485, length.out = 126))$stats[2, ]
    
    return(y)
    
  }
  
  
  sv_38 <- replace(sv_38, sv_38 == -9.9e+37, NA)
  sv <- 10^(sv_38/10)
  mvbs_38 <- 10*log10(apply(sv, 2, mean, na.rm = TRUE))
  
  
  sv_120 <- replace(sv_120, sv_120 == -9.9e+37, NA)
  sv <- 10^(sv_120/10)
  mvbs_120 <- 10*log10(apply(sv, 2, mean, na.rm = TRUE))
  
  
  sv_38_bin <- meanBins(mvbs_38)
  sv_120_bin <- meanBins(mvbs_120)
  
  #calculate difference window
  sv_diff <- sv_120_bin - sv_38_bin
  sv_diff[sv_diff <= 2.5 | sv_diff >= 14.7] <- NA
  sv_120_bin[is.na(sv_diff)] <- NA
  
  #convert to density using target strength (g/m2 per interval)
  p <- 10 ^((sv_120_bin - -42.22)/10)*1000*2
  p[is.na(p)] <- 0
  
  write.csv(p, paste("density/krill_density_gm2_stn_", i, ".csv", sep = ""), row.names = F)
  
  message(paste("Finished calculating density for station", i))
  
}








