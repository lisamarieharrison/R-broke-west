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
    swarm[[j]]$"120khz"$Time_S <- chron(times. = swarm[[j]]$"120khz"$Time_S, format = "h:m:s")
    swarm[[j]]$"120khz"$Time_E <- chron(times. = swarm[[j]]$"120khz"$Time_E, format = "h:m:s")
    
  }
    
  #plot krill swarm location at each density by depth and time  
  plot(c(swarm[["low"]][["120khz"]]$Time_S[1], swarm[["low"]][["120khz"]]$Time_E[1]), c(swarm[["low"]][["120khz"]]$Depth_mean[1], swarm[["low"]][["120khz"]]$Depth_mean[1]), type = "l", xlim = c(min(swarm[["low"]][["120khz"]]$Time_S), max(swarm[["low"]][["120khz"]]$Time_E)), lwd = 2, xaxt = "n", ylim = c(100, 0), yaxt = "n", xlab = "time", ylab = "swarm", col = "white")
  
  for (j in names(swarm)) {
    
    line_col <- "black"
    if (j == "medium") line_col <- "darkorange"
    if (j == "high") line_col <- "red"
    
    for (k in 1:nrow(swarm[[j]][["120khz"]])) {
      
      points(c(swarm[[j]][["120khz"]]$Time_S[k], swarm[[j]][["120khz"]]$Time_E[k]), c(swarm[[j]][["120khz"]]$Depth_mean[k], swarm[[j]][["120khz"]]$Depth_mean[k]), type = "l", lwd = 7, col = line_col)
    }
        
  }
  
  
  title(broke_transects[i])  
  axis(1, at = seq(from = min(swarm[["low"]][["120khz"]]$Time_S), to = max(swarm[["low"]][["120khz"]]$Time_E), length.out = 10), 
       labels = chron(times. = seq(from = chron(times. = min(krill["low"][[1]]$Time_S), format = "h:m:S"), to = chron(times. = max(krill["low"][[1]]$Time_E), format = "h:m:s"), length.out = 10)))
  axis(2, at = seq(100, 0, by = -10), labels = seq(100, 0, by = -10))
  
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
