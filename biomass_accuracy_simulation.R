#assume krill is homogenously distributed
#assume all krill in the survey slice are perfectly detected

calcKrillBiomass <- function(krill_biomass, survey_area_width, survey_area_length, 
                             survey_area_depth, detected_width) {
  
  #krill_biomass = total krill biomass in survey area in kg
  #survey_area_width = survey width in m
  #survey_area_length = survey length in m
  #survey_area_depth = survey depth in m
  #detected_width = width of the survey beam
  
  #apparent kg/m2
  visible_kg <- krill_biomass/(survey_area_width/detected_width)
  apparent_kg_m2 <- visible_kg/(survey_area_length*survey_area_depth)
  
  #calculated krill biomass in kg
  calculated_biomass <- apparent_kg_m2*survey_area_width*survey_area_length
  
  return(list(detected_width = detected_width, calculated_biomass = calculated_biomass))
  
}

detected_width <- 0
percentage_overestimated <- 0
calculated_biomass <- 0
for (i in 1:1000) {
  
  calc_biomass <- calcKrillBiomass(krill_biomass = 5e5, survey_area_width = 1e3, 
                                   survey_area_length = 5e5, survey_area_depth = 250, 
                                   detected_width = i)
  
  detected_width[i] <- calc_biomass$detected_width
  percentage_overestimated[i] <- calc_biomass$calculated_biomass/5e5
  calculated_biomass[i] <- calc_biomass$calculated_biomass
  
}

plot(detected_width, percentage_overestimated, type = "l")
abline(h = 1, col = "red")
abline(v = 250, col = "blue")


#calculation of maximum beam width in m (detected_width) based on angle in degrees
deg <- 30
250 * tan((deg/2)*pi/180 )*2

