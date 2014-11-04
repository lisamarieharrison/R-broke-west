mat.dat <- readMat("C:/Users/Lisa/Documents/phd/southern ocean/BROKE-West raw data/a0603dop.mat")

adcp <- cbind(rep(mat.dat$time, each= 60), rep(mat.dat$date, each = 60), rep(mat.dat$lat, each = 60), rep(mat.dat$lon, each = 60), rep(mat.dat$bindep, 2603), rep(mat.dat$shipu, each = 60), rep(mat.dat$shipv, each = 60), matrix(mat.dat$u, ncol = 1, byrow = F), matrix(mat.dat$v, ncol = 1, byrow = F), matrix(mat.dat$speed, ncol = 1, byrow = F))
colnames(adcp) <- c("time", "date", "lat", "long", "depth", "shipu", "shipv", "u", "v", "speed")

#convert dates and times to date and time format

date <- 0
time <- 0
for (i in 1:length(adcp$date)) {
  if (nchar(adcp$date[i]) == 5) {
    date[i] <- as.character(chron(dates. = as.character(paste("0", adcp$date[i], sep = "")), format = "dmy", out.format = "d/m/y"))
  } else {
    date[i] <- as.character(chron(dates. = as.character(adcp$date[i]), format = "dmy", out.format = "d/m/y"))
  } 
  if (nchar(adcp$time[i]) == 5) {
    time[i] <- as.character(chron(times. = as.character(paste("0", adcp$time[i], sep = "")), format = "hms", out.format = "h:m:s"))
  } 
  if (nchar(adcp$time[i]) == 6) {
    time[i] <- as.character(chron(times. = as.character(adcp$time[i]), format = "hms", out.format = "h:m:s"))
  }
  if (nchar(adcp$time[i]) == 1) {
    time[i] <- as.character(chron(times. = as.character(paste("00000", adcp$time[i]), sep = ""), format = "hms", out.format = "h:m:s"))
  }
  if (nchar(adcp$time[i]) == 4) {
    time[i] <- as.character(chron(times. = as.character(paste("00", adcp$time[i]), sep = ""), format = "hms", out.format = "h:m:s"))
  }
}

