#assume all krill in the survey slice are perfectly detected

calcKrillBiomass <- function(survey_area_width, survey_area_length, 
                             survey_area_depth, edsu_width, detected_width) {
  
  #survey_area_width = survey width in m
  #survey_area_length = survey length in m
  #survey_area_depth = survey depth in m
  #edsu_width = width of each edsu in m
  #detected_width = width of the survey beam
  
  n_edsu <- (survey_area_width/detected_width) * (survey_area_length/edsu_width)
  
  edsu <- matrix(0, ncol = (survey_area_width/detected_width), nrow = (survey_area_length/edsu_width))
  
  #add random krill to cells
  n_cells <- round(ncol(edsu)*nrow(edsu), 0)
  cells <- matrix(c(runif(n_cells, 1, nrow(edsu)), runif(n_cells, 1, ncol(edsu))), ncol = 2)
  edsu[cells] <- rexp(n_cells, 1/20)
  
  #choose random locations for 1000 krill swarms (20% large and 80% medium)
  krill_col <- round(runif(1000, 1, ncol(edsu)))
  krill_row <- round(runif(1000, 1, nrow(edsu)))
  
  #place krill into survey area
#   for (i in 1:200) {
#     edsu[krill_row[i], krill_col[i]] <- rnorm(1, 50)
#   }
#   
#   for (i in 210:1000) {
#     edsu[krill_row[i],  krill_col[i]] <- rnorm(1, 20)
#   }
  
  krill_biomass <- sum(edsu*detected_width*edsu_width)/1000/1000 #actual biomass in Mt
  
  calculated_biomass <- mean(edsu[, 1])*survey_area_width*survey_area_length/1000/1000 #in Mt 
  
  return(list(krill_biomass = krill_biomass, calculated_biomass = calculated_biomass))
  
}


krill_biomass <- 0
calculated_biomass <- 0
for (i in 1:200) {
  
  calc_biomass <- calcKrillBiomass(survey_area_width = 2e+05, 
                                   survey_area_length = 6e+05, survey_area_depth = 250, 
                                   detected_width = 50, edsu_width = 2000)
  
  calculated_biomass[i] <- calc_biomass$calculated_biomass
  krill_biomass[i] <- calc_biomass$krill_biomass
  
  if (i %% 100 == 0) print(i)
  
}

hist(krill_biomass - calculated_biomass, main = "Error in calc (Mt)")








