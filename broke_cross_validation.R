#fits asreml model developed in simulation to BROKE-West data
#uses the distance between stations rather than latitude and longitude
#runs cross-validation by dropping one station at a time
#author: Lisa-Marie Harrison
#date: 11/11/2014

setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/Data")
dat <- read.csv(file = "procCTD.csv", header= T)
library(asreml)
library(nlme)
library(lattice)

names(dat) <- c("survey", "stn", "lat", "long", "start.time", "end.time", "depth", "transmittance", "cond", "temp", "sal", "par", "oxygen", "fluoro", "x2", "ice", "wm")

#remove null values
dat$sal[dat$sal == -9] <- NA
dat$temp[dat$temp == -9] <- NA
dat$par[dat$par == -9] <- NA
dat$fluoro[dat$fluoro == -9] <- NA

#compute log transformed fluoro values
dat$l.fluoro <- log(dat$fluoro)
dat$l.fluoro[is.nan(dat$l.fluoro)] <- NA

#plot latitude and longitude of each station
plot(dat$long[seq(1, nrow(dat), 125)], dat$lat[seq(1, nrow(dat), 125)], col = "white")
text(dat$long[seq(1, nrow(dat), 125)], dat$lat[seq(1, nrow(dat), 125)], dat$stn[seq(1, nrow(dat), 125)])

#get included stations for each run (dropping one arm at a time)
run1 <- setdiff(2:120, 27:44)
run2 <- setdiff(2:120, 45:59)
run3 <- setdiff(2:120, 60:71)
run4 <- setdiff(2:120, 72:85)
run5 <- setdiff(2:120, 86:102)
run6 <- setdiff(2:120, 103:120)

runs <- list(run1, run2, run3, run4, run5, run6)

pval_run <- list()
se_run <- list()

#-------------run cross-validation, dropping one arm at a time ----------------#

for (s in 1:6) {
  
  dat.cut <- dat[dat$stn %in% runs[s][[1]], ]
  
  #get latitude and longitude for each station
  n.station <- length(unique(dat.cut$stn))
  lat  <- dat.cut$lat[duplicated(dat.cut$stn) == FALSE]
  long <- dat.cut$long[duplicated(dat.cut$stn) == FALSE]
  
  #function to convert degrees to radians
  deg2rad <- function(deg) {
    return(deg*pi/180)
  }
  
  #Calculates the distance between two points with radian latitude/longitude using Haversine formula (hf)
  gcd.hf <- function(lat1, long1, lat2, long2) {
    R <- 6371 # Earth mean radius [km]
    delta.long <- (long2 - long1)
    delta.lat <- (lat2 - lat1)
    a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
    c <- 2 * asin(min(1,sqrt(a)))
    d = R * c
    return(d) # Distance in km
  }
  
  dist_x <- matrix(0, ncol = n.station, nrow = n.station)
  for (i in 1:n.station) {
    for (k in 1:n.station) {
      dist_x[i, k] <- gcd.hf(deg2rad(lat[i]), deg2rad(long[i]), deg2rad(lat[k]), deg2rad(long[i]))/100
    }
  }
  
  dist_y <- matrix(0, ncol = n.station, nrow = n.station)
  for (i in 1:n.station) {
    for (k in 1:n.station) {
      dist_y[i, k] <- gcd.hf(deg2rad(lat[i]), deg2rad(long[i]), deg2rad(lat[i]), deg2rad(long[k]))/100
    }
  }
  
  #get distance of each station from station 1 in x and y directions
  x <- dist_x[1, ]
  y <- dist_y[1, ]
  
  
  #dat.cuta frame
  glm.spl <- data.frame(dat.cut$l.fluoro, dat.cut$depth, as.factor(dat.cut$stn), rep(x, 1, each = length(unique(dat.cut$depth))), rep(y, 1, each = length(unique(dat.cut$depth))), dat.cut$temp, dat.cut$par, dat.cut$sal, dat.cut$oxygen, dat.cut$ice, as.factor(dat.cut$wm))
  names(glm.spl) <- c("l.obs", "z", "stn", "x", "y", "temp", "par", "sal", "oxy", "ice", "wm")
  glm.spl$z.fact <- as.factor(as.integer(glm.spl$z))
  glm.spl$x.fact <- as.factor(glm.spl$x)
  glm.spl$y.fact <- as.factor(glm.spl$y)
  glm.spl <- glm.spl[order(glm.spl$z, glm.spl$x, glm.spl$y), ] #sort by order of rcov structure
  glm.spl$l.obs[glm.spl$l.obs == -Inf] <- NA
  
  #centre and scale covariates to mean = 0 and sd = 1
  #this is required if using na.method = "include" since this sets the missing values to 0
  glm.spl$temp <- scale(glm.spl$temp)
  glm.spl$par  <- scale(glm.spl$par)
  glm.spl$sal  <- scale(glm.spl$sal)
  glm.spl$oxy  <- scale(glm.spl$oxy)
  glm.spl$ice  <- scale(glm.spl$ice)
  glm.spl$oxy  <- scale(glm.spl$oxy)
  
  
  #------------------------------- FIT ASREML MODELS -----------------------------------#
  
  #fit asreml model
  asreml.fit <- asreml(fixed = l.obs ~ z + par + temp:diag(wm) + ice + oxy, random =~ spl(z, 10) + spl(par, 10) + 
                         spl(temp, 10):diag(wm) +  spl(ice, 10) + spl(oxy, 10) + stn, 
                       dat = glm.spl, rcov=~ ar1(z.fact):agau(x.fact, y.fact),
                       na.method.X = "include", workspace = 50000000)
  asreml.fit <- update(asreml.fit)
    
  
  #----------------------- PREDICT VALUES AT OTHER STATIONS ---------------------#
  
  #predict at each point using parameter values for other stations to test extrapolation
  #using depth, par, temperature, ice and watermass, predict l.fluoro using the model
  extra.stn <- setdiff(2:120, runs[s][[1]])
  extra.dat <- dat[dat$stn %in% extra.stn, ]
  
  
  extra.dat$temp <- scale(extra.dat$temp)
  extra.dat$ice <- scale(extra.dat$ice)
  extra.dat$par <- scale(extra.dat$par)
  
  
  pval <- 0
  se <- 0
  for(i in 1:nrow(extra.dat)) {
    
    capture.output({
      pred <- predict(asreml.fit, classify = "wm:ice:par:temp:z", levels = list("z" = extra.dat$depth[i], "temp" = extra.dat$temp[i], "ice" = extra.dat$ice[i], "par" = extra.dat$par[i], "wm" = extra.dat$wm[i]))
    }, file = tempfile())
    pval[i] <- pred$predictions$pvals["predicted.value"]$predicted.value
    se[i] <- pred$predictions$pvals["standard.error"]$standard.error
    
    if(i %% 100 == 0) print(paste(Sys.time(), i))
    
  }  
  
  #write out predicted values and standard errors
  write.csv(pval, "")
  se_run <- list(se_run, se)
  
  print(paste("Finished run", s))
  
}


#plot fitted against observed by station
pval[is.na(extra.dat$l.fluoro)] <- NA
lat.plot <- xyplot(extra.dat$l.fluoro + pval ~ extra.dat$depth | extra.dat$stn, outer = FALSE, type = "l")
update(lat.plot, par.settings = simpleTheme(lwd = c(2, 1), col = c("dodgerblue", "red")))




