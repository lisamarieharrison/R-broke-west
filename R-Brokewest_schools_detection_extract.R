#For each day of the BROKE-West survey export schools detection in Echoview
#Author: Lisa-Marie Harrison
#Date: 28/01/2015
library(RDCOMClient)
library(EchoviewR)

#create an Echoview object and open Echoview
EVAppObj=COMCreate('EchoviewCom.EvApplication')

#get list of EV files
file.dates <- list.files(path = "C:/Users/43439535/Documents/Lisa/BROKE-West/EV", full.names = T, pattern = paste(".*\\.ev$", sep = ""))
transect <- substr(list.files(path = "C:/Users/43439535/Documents/Lisa/BROKE-West/EV"), 1, 13)

#loop over each date and export integration by cells separately
for (i in 1:length(file.dates)) {
  
  
  #read in an Echoview file
  EVFile= EVOpenFile(EVAppObj, fileName = file.dates[i])$EVFile
  
  #add calibration file
  EVAddCalibrationFile(EVFile, "38H-120H-200H", "C:/Users/43439535/Documents/Lisa/BROKE-West/SimradEK60DSRII2010.ecs")
  
  #read raw data names in fileset
  file.names <- EVFilesInFileset(EVFile = EVFile, filesetName = '38H-120H-200H')
  
  #clear raw data files
  EVClearRawData(EVFile = EVFile, filesetName = '38H-120H-200H')
  
  #add multiple raw data file to EV test
  raw.files <- paste("C:/Users/43439535/Documents/Lisa/BROKE-West/RAW/", file.names, sep = "")
  EVAddRawData(EVFile = EVFile, filesetName = "38H-120H-200H", dataFiles = raw.files)
  
  #change grid of Krill aggregations
  varObj = EVAcoVarNameFinder(EVFile, acoVarName = "120H hri 7x7 convolution")$EVVar
  EVChangeVariableGrid(EVFile = EVFile, acousticVar = varObj, verticalType = 4, verticalDistance = 7, horizontalType = 1, horizontalDistance = 1)
  
  varObj = EVAcoVarNameFinder(EVFile, acoVarName = "38H hri 7x7 convolution")$EVVar
  EVChangeVariableGrid(EVFile = EVFile, acousticVar = varObj, verticalType = 4, verticalDistance = 7, horizontalType = 1, horizontalDistance = 1)
  
  #create new region classes for medium and high aggregations
  EVNewRegionClass(EVFile, "medium aggregations")
  EVNewRegionClass(EVFile, "high aggregations")
  EVNewRegionClass(EVFile, "aggregations")
  
  #Run schools detection on 120 7x7 convolution for each type of aggregation (low, medium and high)
  schDet<-EVSchoolsDetect(EVFile = EVFile,
                          acoVarName='120H hri 7x7 convolution',
                          outputRegionClassName='aggregations',
                          deleteExistingRegions = TRUE,
                          distanceMode="GPS distance",
                          maximumHorizontalLink=15,#m
                          maximumVerticalLink=5,#m
                          minimumCandidateHeight=1,#m
                          minimumCandidateLength=10,#m
                          minimumSchoolHeight=2,#m
                          minimumSchoolLength=15, #m
                          dataThreshold= -80)
  
  medSchDet<-EVSchoolsDetect(EVFile = EVFile,
                             acoVarName='120H hri 7x7 convolution',
                             outputRegionClassName='medium aggregations',
                             deleteExistingRegions = TRUE,
                             distanceMode="GPS distance",
                             maximumHorizontalLink=15,#m
                             maximumVerticalLink=5,#m
                             minimumCandidateHeight=1,#m
                             minimumCandidateLength=10,#m
                             minimumSchoolHeight=2,#m
                             minimumSchoolLength=15, #m
                             dataThreshold= -65)
  
  highSchDet<-EVSchoolsDetect(EVFile = EVFile,
                              acoVarName='120H hri 7x7 convolution',
                              outputRegionClassName='high aggregations',
                              deleteExistingRegions = TRUE,
                              distanceMode="GPS distance",
                              maximumHorizontalLink=15,#m
                              maximumVerticalLink=5,#m
                              minimumCandidateHeight=1,#m
                              minimumCandidateLength=10,#m
                              minimumSchoolHeight=2,#m
                              minimumSchoolLength=15, #m
                              dataThreshold= -58)
  
  #reset the 120 7x7 convolution threshold
  varObj <- EVAcoVarNameFinder(EVFile, acoVarName = "120H hri 7x7 convolution")$EVVar
  EVminThresholdSet(varObj, -80)
  
  #export integration by regions for 120kHz
  return.wd <- "C:/Users/43439535/Documents/Lisa/BROKE-West/Extracted data/schools detection"
  EVIntegrationByRegionsExport(EVFile = EVFile, acoVarName = "120H hri 7x7 convolution", regionClassName = "aggregations", exportFn = paste(return.wd, "/", transect[i], "_low_aggregations_by_region_120khz.csv", sep = ""))
  EVIntegrationByRegionsExport(EVFile = EVFile, acoVarName = "120H hri 7x7 convolution", regionClassName = "medium aggregations", exportFn = paste(return.wd, "/", transect[i], "_medium_aggregations_by_region_120khz.csv", sep = ""))
  EVIntegrationByRegionsExport(EVFile = EVFile, acoVarName = "120H hri 7x7 convolution", regionClassName = "high aggregations", exportFn = paste(return.wd, "/", transect[i], "_high_aggregations_by_region_120khz.csv", sep = ""))
  
  #export integration by regions for 38kHz
  EVIntegrationByRegionsExport(EVFile = EVFile, acoVarName = "38H hri 7x7 convolution", regionClassName = "aggregations", exportFn = paste(return.wd, "/", transect[i], "_low_aggregations_by_region_38khz.csv", sep = ""))
  EVIntegrationByRegionsExport(EVFile = EVFile, acoVarName = "38H hri 7x7 convolution", regionClassName = "medium aggregations", exportFn = paste(return.wd, "/", transect[i], "_medium_aggregations_by_region_38khz.csv", sep = ""))
  EVIntegrationByRegionsExport(EVFile = EVFile, acoVarName = "38H hri 7x7 convolution", regionClassName = "high aggregations", exportFn = paste(return.wd, "/", transect[i], "_high_aggregations_by_region_38khz.csv", sep = ""))
  
  #save the open .EV file
  EVSaveFile(EVFile = EVFile)
  
  #close the current file
  EVCloseFile(EVFile = EVFile)
  
  message("Finished processing", transect[i])
  
}

