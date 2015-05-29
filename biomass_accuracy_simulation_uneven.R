#assume all krill in the survey slice are perfectly detected

calcKrillBiomass <- function(survey_area_width, survey_area_length, 
                             survey_area_depth, edsu_width, detected_width) {
  
  #survey_area_width = survey width in m
  #survey_area_length = survey length in m
  #survey_area_depth = survey depth in m
  #edsu_width = width of each edsu in m
  #detected_width = width of the survey beam
  
  n_edsu <- (survey_width/detected_width) * (survey_length/edsu_width)

  edsu <- matrix(0, ncol = (survey_width/detected_width), nrow = (survey_length/edsu_width))
  
  #add random krill to cells
  n_cells <- round(ncol(edsu)*nrow(edsu), 0)
  rows <- c(round(runif(n_cells, 1, nrow(edsu)), 0))
  cols <- c(round(runif(n_cells, 1, ncol(edsu)), 0))
  for (i in 1:n_cells) {
    edsu[rows[i], cols[i]] <- rexp(1, 1/20)
  }
  
  #choose random locations for 100 krill swarms
  krill_col <- round(runif(100, 1, ncol(edsu)))
  krill_row <- round(runif(100, 1, nrow(edsu)))
  
  #place krill into survey area
  for (i in 1:30) {
    cols <- c(krill_col[i] - 2, krill_col[i] - 1, krill_col[i], krill_col[i] + 1, krill_col[i] - 2)
    rows <- c(krill_row[i] - 2, krill_row[i] - 1, krill_row[i], krill_row[i] + 1, krill_row[i] + 2)
    if (any(cols < 1)) cols[cols < 1] <- 1
    if (any(cols > ncol(edsu))) cols[cols > ncol(edsu)] <-  ncol(edsu)
    if (any(rows < 1)) rows[rows < 1] <- 1
    if (any(rows > nrow(edsu))) rows[rows > nrow(edsu)] <- nrow(edsu)
    edsu[rows, cols] <- rnorm(25, 50)
  }
  
  for (i in 31:100) {
    cols <- c(krill_col[i] - 1, krill_col[i], krill_col[i] + 1)
    rows <- c(krill_row[i] - 1, krill_row[i], krill_row[i] + 1)
    if (any(cols < 1)) cols[cols < 1] <- 1
    if (any(cols > ncol(edsu))) cols[cols > ncol(edsu)] <-  ncol(edsu)
    if (any(rows < 1)) rows[rows < 1] <- 1
    if (any(rows > nrow(edsu))) rows[rows > nrow(edsu)] <- nrow(edsu)
    edsu[rows, cols] <- rnorm(9, 20)
  }
  
  krill_biomass <- sum(edsu*detected_width*edsu_width)/1000/1000 #actual biomass in Mt
  
  calculated_biomass <- mean(edsu[, 1])*survey_width*survey_length/1000/1000 #in Mt 
  
  return(list(krill_biomass = krill_biomass, calculated_biomass = calculated_biomass))
  
}


krill_biomass <- 0
calculated_biomass <- 0
improved_calculated_biomass <- 0
for (i in 1:1000) {
  
  calc_biomass <- calcKrillBiomass(survey_area_width = 5e3, 
                                   survey_area_length = 2e4, survey_area_depth = 250, 
                                   detected_width = 100, edsu_width = 100)
  
  calculated_biomass[i] <- calc_biomass$calculated_biomass
  krill_biomass[i] <- calc_biomass$krill_biomass
  improved_calculated_biomass[i] <- calc_biomass$improved_calculated_biomass
  
}

hist(krill_biomass - calculated_biomass, main = "Error in calc (Mt)")









