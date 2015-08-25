#collate start and end date and times of all BROKE-West files into a single csv
#author: Lisa-Marie Harrison
#date: 25/08/2015

setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/BROKE-West/Echoview/Extracted data/2x25 integration")
files <- list.files(pattern = "38kHz")
library(R.utils)

start_date <- 0
end_date   <- 0
start_time <- 0
end_time   <- 0

for (i in 1:length(files)) {
  
  n_lines <- countLines(files[i])
  line_1 <- read.csv(files[i], nrows = 1, header = T)
  last_line <- read.csv(files[i], skip = n_lines - 1, header = F)
  
  start_date[i] <- line_1$Date_S
  end_date[i] <- last_line$V22
  
  start_time[i] <- as.character(line_1$Time_S)
  end_time[i] <- as.character(last_line$V23)
  
}

dat <- data.frame(start_date, end_date, start_time, end_time)

write.csv(dat, "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/Data/ev_file_times_2m.csv", row.names = F)







