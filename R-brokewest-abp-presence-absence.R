#presence and absence of broke-west swarms linked to predator presence-absence
#Author: Lisa-Marie Harrison
#Date: 05/02/2015

setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/BROKE-West raw data/Echoview/Extracted data/schools detection")
broke_transects <- unique(substr(list.files(), start = 1, stop = 13))
library(chron)

for (i in 1:length(broke_transects)) {
  
  #read in predator data
  pred <- read.csv("C:/Users/Lisa/Documents/phd/southern ocean/BROKE-West raw data/Air breathing predator/broke_seabirds.csv")
  
  #create empty multi level list
  level_2 <- list("38khz" = c(), "120khz" = c())
  swarm <- list("low" = level_2, "medium" = level_2, "high" = level_2)
    
  for (j in names(swarm)) {    
    for (k in names(swarm[[j]])) {    
      
      dat <- paste(broke_transects[i], "_", j, "_aggregations_by_region_", k, ".csv", sep = "")
      
      if (file.exists(dat)) {
        swarm[[j]][[k]] <- read.csv(dat, header = T)
      } else {
        swarm[[j]] <- NULL
      }
    }  
  }
  
  #dB difference to find krill
  for (j in names(swarm)) {  
    
    db_diff <- swarm[[j]]$"120khz"$Sv_mean - swarm[[j]]$"38khz"$Sv_mean
    swarm[[j]]$"120khz"$Sv_mean[2.5 > db_diff | db_diff > 14.7] <- NA
    
  } 
  
  for (j in names(swarm)) {
    
    #exclude depths > 100m
    swarm[[j]]$"120khz" <- swarm[[j]]$"120khz"[swarm[[j]]$"120khz"$Depth_mean <= 100, ]
    swarm[[j]]$"38khz" <- swarm[[j]]$"38khz"[swarm[[j]]$"38khz"$Depth_mean <= 100, ]
    
    #convert krill times to a chron object
    swarm[[j]]$"120khz"$dt_start <- chron(chron(dates. = as.character(swarm[[j]]$"120khz"$Date_S), times. = swarm[[j]]$"120khz"$Time_S, format = list(dates. = "ymd", times. = "h:m:s"), out.format = list(dates. = "d/m/y", times. = "h:m:s")))
    swarm[[j]]$"120khz"$dt_end <- chron(chron(dates. = as.character(swarm[[j]]$"120khz"$Date_E), times. = swarm[[j]]$"120khz"$Time_E, format = list(dates. = "ymd", times. = "h:m:s"), out.format = list(dates. = "d/m/y", times. = "h:m:s")))
    
  }
    
  #plot krill swarm location at each density by depth and time  
  plot(c(swarm$"low"$"120khz"$dt_start[1], swarm$"low"$"120khz"$dt_end[1]), c(swarm$"low"$"120khz"$Depth_mean[1], swarm$"low"$"120khz"$Depth_mean[1]), type = "l", xlim = c(min(swarm$"low"$"120khz"$dt_start), max(swarm$"low"$"120khz"$dt_end)), lwd = 2, ylim = c(100, 0), xlab = "time", ylab = "swarm", col = "white")
  title(broke_transects[i])  
  
  for (j in names(swarm)) {
    
    line_col <- "black"
    if (j == "medium") line_col <- "darkorange"
    if (j == "high") line_col <- "red"
    
    for (k in 1:nrow(swarm[[j]][["120khz"]])) {
      
      points(c(swarm[[j]][["120khz"]]$dt_start[k], swarm[[j]][["120khz"]]$dt_end[k]), c(swarm[[j]][["120khz"]]$Depth_mean[k], swarm[[j]][["120khz"]]$Depth_mean[k]), type = "l", lwd = 7, col = line_col)
    }
        
  }

  #subset predator data to include only the correct date and Adelie Penguins
  pred <- pred[pred$date %in% unique(swarm[["low"]][["120khz"]]$Date_S), ]
  pred <- pred[pred$species == "Pygoscelis adeliae (Hombron and Jacquinot,1841) (Adelie Penguin)", ]
  pred$t <- chron(times. = pred$t, format = "h:m:s")
  
  if (nrow(pred) != 0) {
    for (j in 1:nrow(pred)) {
      
      points(pred$t[j], 4, pch = "x", cex = 2, col = "blue")
      
    }
  }
  
}
