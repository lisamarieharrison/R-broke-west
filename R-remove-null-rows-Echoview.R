#removes null rows (-9999) from exported csv files from Echoview
#author: Lisa-Marie Harrison
#date: 10/06/2015

files <- list.files("C:/Users/Lisa/Documents/phd/southern ocean/BROKE-West/Echoview/Extracted data/250x2000 integration", full.names = T)

for (i in files) {
  
  dat <- read.csv(i, header = T)
  w <- !dat$ï..Region_ID == -9999
  dat <- dat[w, ]
  
  write.csv(dat, i, row.names = F)
  
}