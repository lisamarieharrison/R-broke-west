#calculates krill density at each interval for CTD depths
#author: Lisa-Marie Harrison
#date: 26/03/2015

setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/BROKE-West/Echoview/Extracted data/ctd")
k120 <- read.csv("120kHz/stn_2_extracted_120khz.csv", header = F, skip = 1)
k38 <- read.csv("38kHz/stn_2_extracted_38khz.csv", header = F, skip = 1)


#rename first 12 columns
names(k120)[1:12] <- c("Ping_index", "Distance_gps", "Distance_vl", "Ping_date", "Ping_time", 
                      "Ping_milliseconds", "Latitude", "Longitude", "Depth_start", "Depth_stop", "Range_start",
                      "Range_stop")

names(k38)[1:12] <- c("Ping_index", "Distance_gps", "Distance_vl", "Ping_date", "Ping_time", 
                       "Ping_milliseconds", "Latitude", "Longitude", "Depth_start", "Depth_stop", "Range_start",
                       "Range_stop")

#for each row, find the krill density
sv_120 <- k120[, 14:ncol(k120)]
sv_120[sv_120 == -9.9e+37] <- NA
sv_38 <- k38[, 14:ncol(k38)]
sv_38[sv_38 == -9.9e+37] <- NA


#calculate difference window
sv_diff <- k120[, 14:ncol(k120)] - k38[, 14:ncol(k38)]
sv_diff[sv_diff < 1.02 | sv_diff > 14.75] <- NA
sv_120[is.na(sv_diff)] <- NA

sv <- 10^(sv_120/10)

#bin data into 125 evenly spaced bins
x <- apply(sv, 1, meanBins)

meanBins <- function(x) {
  
  y <- stats.bin(x = 1:484, y = x, breaks = seq(1, 485, length.out = 126))$stats[2, ]
  
  return(y)
  
}

mean_120 <- stats.bin(x = 1:484, y = k120[1, 14:ncol(k120)], breaks = seq(1, 485, length.out = 126))$stats[2, ]



#convert to density using target strength (kg/m2 per interval)
p <- 2*10 ^((sv_120 - -42.22)/10)*1000


mvbs  <- 10*log10(aggregate(matrix(sv, ncol = 1), by = list(rep(c(1:(length(sv)/max(acoustic_38$Layer))), each = max(acoustic_38$Layer))), sum, na.rm = T)$V1/max(acoustic_38$Layer))
mvbs[mvbs == -Inf] <- NA




