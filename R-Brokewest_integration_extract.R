#For each EV transect file of the BROKE-West survey export integration by cells in Echoview 
#Exports 38kHz and 120kHz on a specified grid for high resolution data
#Author: Lisa-Marie Harrison
#Date: 28/01/2015
library(RDCOMClient)
library(EchoviewR)

#create an Echoview object and open Echoview
EVAppObj=COMCreate('EchoviewCom.EvApplication')

#get list of EV files
file.dates <- list.files(path = "H:/phd/southern ocean/BROKE-West raw data/Echoview/EV", full.names = T, pattern = paste(".*\\.ev$", sep = ""))
transect <- substr(list.files(path = "H:/phd/southern ocean/BROKE-West raw data/Echoview/EV"), 1, 13)


#loop over each date and export integration by cells separately
for (i in 1:length(file.dates)) {
  
  
  #read in an Echoview file
  EVFile= EVOpenFile(EVAppObj, fileName = file.dates[i])$EVFile
  
  #add calibration file
  EVAddCalibrationFile(EVFile, "38H-120H-200H", "C:/Users/Lisa/Documents/phd/southern ocean/BROKE-West/Echoview/calibration/SimradEK60DSRII2010.ecs")
  
  #read raw data names in fileset
  file.names <- EVFilesInFileset(EVFile = EVFile, filesetName = '38H-120H-200H')
  
  #clear raw data files
  EVClearRawData(EVFile = EVFile, filesetName = '38H-120H-200H')
  
  #add multiple raw data file to EV test
  raw.files <- paste("H:/BROKE-West data/BROKE-West raw data/", file.names, sep = "")
  EVAddRawData(EVFile = EVFile, filesetName = "38H-120H-200H", dataFiles = raw.files)
  
  #change grid of 120kHz and 38kHz variables to 10m x 50pings
  varObj = EVAcoVarNameFinder(EVFile, acoVarName = "120 H hrp 0 to 250 m")$EVVar
  EVChangeVariableGrid(EVFile = EVFile, acousticVar = varObj, verticalType = 5, verticalDistance = 2000, horizontalType = 1, horizontalDistance = 250)
  EVSetAcoVarDisplayDepth(EVFile, varObj, 0, 250)
  
  varObj = EVAcoVarNameFinder(EVFile, acoVarName = "38 H hrp 0 to 250 m")$EVVar
  EVChangeVariableGrid(EVFile = EVFile, acousticVar = varObj, verticalType = 5, verticalDistance = 2000, horizontalType = 1, horizontalDistance = 250)
  EVSetAcoVarDisplayDepth(EVFile, varObj, 0, 250)
  
  #export integration by cells for 38kHz and 120kHz
  return.file.38 <- paste("C:/Users/Lisa/Documents/phd/southern ocean/BROKE-West/Echoview/Extracted data/250x2000 integration/", transect[i], "_250m_2000m_38kHz.csv", sep = "")
  return.file.120 <- paste("C:/Users/Lisa/Documents/phd/southern ocean/BROKE-West/Echoview/Extracted data/250x2000 integration/", transect[i], "_250m_2000m_120kHz.csv", sep = "")
  EVExportIntegrationByRegionByCells(EVFile, "38 H hrp 0 to 250 m", "export_region_250m", return.file.38)
  EVExportIntegrationByRegionByCells(EVFile, "120 H hrp 0 to 250 m", "export_region_250m", return.file.120)
  
  EVCloseFile(EVFile)
  
}


