#fits asreml model developed in simulation to BROKE-West data
#uses the distance between stations rather than latitude and longitude
#runs cross-validation by dropping one station at a time
#author: Lisa-Marie Harrison
#date: 11/11/2014

if (Sys.info()[4] == "SCI-6246") {
  setwd(dir = "C:/Users/43439535/Documents/Lisa/phd/Mixed models")
} else {
  setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models")
}

dat <- read.csv(file = "Data/procCTD.csv", header= T)
library(asreml)
library(nlme)
library(lattice)

names(dat) <- c("survey", "stn", "lat", "long", "start.time", "end.time", "depth", "transmittance", "cond", "temp", "sal", "par", "oxygen", "fluoro", "x2", "ice", "wm")

#source required functions
function_list <- c("distFromStn1.R")

for (f in function_list) {
  source(paste("R code/R-functions-southern-ocean/", f, sep = ""))
}

#remove null values
dat$sal[dat$sal == -9] <- NA
dat$temp[dat$temp == -9] <- NA
dat$par[dat$par == -9] <- NA
dat$fluoro[dat$fluoro == -9] <- NA

#compute log transformed fluoro values
dat$l.fluoro <- log(dat$fluoro)
dat$l.fluoro[is.nan(dat$l.fluoro)] <- NA

#plot latitude and longitude of each station
plot(dat$long[seq(1, nrow(dat), 125)], dat$lat[seq(1, nrow(dat), 125)], col = "white", xlab = "Longitude", ylab = "Latitude")
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
  
  #get distance of each station from station 1 in x and y directions
  distance <- distFromStn1(lat, long)
  x <- distance$x
  y <- distance$y
  
  #dat.cuta frame
  glm.spl <- data.frame(dat.cut$l.fluoro, dat.cut$depth, as.factor(dat.cut$stn), rep(x, 1, each = length(unique(dat.cut$depth))), rep(y, 1, each = length(unique(dat.cut$depth))), dat.cut$temp, dat.cut$par, dat.cut$sal, dat.cut$oxygen, dat.cut$ice, as.factor(dat.cut$wm))
  names(glm.spl) <- c("l.obs", "z", "stn", "x", "y", "temp", "par", "sal", "oxy", "ice", "wm")
  glm.spl$z.fact <- as.factor(as.integer(glm.spl$z))
  glm.spl$x.fact <- as.factor(glm.spl$x)
  glm.spl$y.fact <- as.factor(glm.spl$y)
  glm.spl <- glm.spl[order(glm.spl$z, glm.spl$x, glm.spl$y), ] #sort by order of rcov structure
  glm.spl$l.obs[glm.spl$l.obs == -Inf] <- NA
  
  #centre and scale [temp, par, sal, oxy, ice] to mean = 0 and sd = 1
  #this is required if using na.method = "include" since this sets the missing values to 0
  glm.spl <- cbind(glm.spl[, c(1:5, 11:14)], apply(glm.spl[, 6:10], 2, scale))
  
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
  se   <- 0
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




