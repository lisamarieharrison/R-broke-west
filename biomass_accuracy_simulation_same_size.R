#function for simulation for krill biomass estimation
#author: Lisa-Marie Harrison


#assume all krill in the survey slice are perfectly detected

calcKrillBiomass <- function(survey_area_width, survey_area_length, 
                             survey_area_depth, edsu_width, detected_width) {
  
  #survey_area_width = survey width in m
  #survey_area_length = survey length in m
  #survey_area_depth = survey depth in m
  #edsu_width = width of each edsu in m
  #detected_width = width of the survey beam
  
  n_edsu <- (survey_area_width/detected_width) * (survey_area_length/edsu_width)
  
  krill <- matrix(rexp(n_edsu, 1/20), ncol = (survey_area_width/detected_width), nrow = survey_area_length/10)
  
  n_average <- edsu_width/10
  
  edsu <- matrix(0, ncol = 1, nrow = nrow(krill)/n_average)
  for (i in 1:(nrow(krill)/n_average)) {
    edsu[i, 1] <- mean(krill[((i-1)*n_average+1):(i*n_average), 1])
  }
    
  krill_biomass <- sum(krill*50*10)/1000/1000 #actual biomass in Mt
  
  calculated_biomass <- mean(edsu)*survey_area_width*survey_area_length/1000/1000 #in Mt 
  
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








