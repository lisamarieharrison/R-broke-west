#gets the adcp current at each BROKE-West station
#author: Lisa-Marie Harrison

setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/Data/")
adcp <- read.csv("adcp_data.csv", header = T)
dat.cut <- read.csv("rstnCTD.csv", header= T)
date <- read.csv("stn_coordinates.csv", header = T)

dat.cut$date <- rep(date$date[date$station %in% unique(dat.cut$stn)], each = 125)

stn.current <- 0
for (i in 1:nrow(dat.cut)) {
  #narrow find similar date and time
  w <- which(adcp$date == unique(dat.cut$date[dat.cut$stn == dat.cut$stn[i]]) & (hours(chron(times. = adcp$time, format = "h:m:s") - unique(chron(times. = dat.cut$start.time, format = "h:m:s")[dat.cut$stn == dat.cut$stn[i]])) == 0) & (minutes(chron(times. = adcp$time, format = "h:m:s") - unique(chron(times. = dat.cut$start.time, format = "h:m:s")[dat.cut$stn == dat.cut$stn[i]])) < 20) & ((abs(adcp$depth - dat.cut$profile.depth[i])) < 2))
  if (length(w) > 0) {
  stn.current[i] <- adcp$speed[w][1]
  } else {
    stn.current[i] <- NA
  }
  if (i %% 100 == 0) print(i)
}
stn.current[is.nan(stn.current)] <- NA

dat.cut$current <- stn.current




