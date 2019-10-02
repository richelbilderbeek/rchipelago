###################################################################
### ============================================================
### PIAC Run  
### ============================================================
###################################################################

#   This script contains the workflow for preparing the data underlying
#   the study 'A global spatially explicit database of changes in island paleo-area and
#   archipelago configuration during the late Quaternary'
#   by Norder et al. 
#   submitted to Global Ecology and Biogeography

### ============================================================
### Initialisation  
### ============================================================

rm(list=ls()) # clear workspace
gc() # garbage collector to clear RAM

# Load sea levels, choose either Cutler et al., (2003) or Lambeck et al., (2014)
slcs <- c("Lambeck", "Cutler")
slc <- slcs[2]

# Load function library
source("PIAC_functions.r")
# Libraries
library(rgdal)
library(xlsx) # make sure correct version of Java is installed on computer.
library(ncdf4)
library(raster)
library(sp)

# List of archipelagos (each archipelago has a corresponding KML file in input folder)
Arclist <- c("Aldabra","Azores","Balearic","California Channel","Canaries","Cape Verde","Comoros","Cook","Crozet","Dutch Caribbean","Galapagos","Gulf of Guinea","Hawaii","Inner Seychelles","Juan Fernandez","Kuriles", "Madeira","Marianas","Marquesas","Mascarenes","Phoenix","Pitcairn","Prince Edward","Revillagigedo","Samoa","Society","Tristan da Cunha")

# Define location of input, intermediate, and outputfolder
inputfolder<-file.path(getwd(),"input") 
intermedfolder<-file.path(getwd(),"intermediate")
outputfolder<-file.path(getwd(),"output",slc)

# Load sea level curve
time_level<-readRDS(paste(inputfolder,"/GLOBAL/","time_level_",slc,".rds",sep=""))

# Load global bathymetry DEM
GEBCO <- paste(inputfolder,"/GLOBAL/",'GEBCO','/GEBCO_2014_1D.nc',sep="")
if(!file.exists(GEBCO)){stop("The file GEBCO_2014_1D.nc could not be found in the following path: 'PIAC/input/GLOBAL/GEBCO/'. Please download the  GEBCO_2014 Grid from https://www.gebco.net/data_and_products/gridded_bathymetry_data/ and place it in the correct folder.")}
DEMworld<-raster(GEBCO)

if(!file.exists('C:\\OSGeo4W64\\OSGeo4W.bat')){warning("OSGeo4W.bat could not be found in 'C:\\OSGeo4W64\\'. This is not necessarily problematic because the function 'polygonizer' will attempt to determine its location in your file system. However, if OSGeo4W is not yet installed on your computer, it should be installed (OSGeo4W can be downloaded from https://trac.osgeo.org/osgeo4w/")}

# Define ESRI World Sinusoidal projected coordinate system (equal area)
# CRS("+init=ESRI:54008")@projargs
sinus.csy <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0"

# Google CRS, alternatively use: CRS("+init=epsg:4326")
g.csy <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"

# Optional: Select island names from KML files and store as csv.
islandnames<-IB.datastructure(IB.inputfolder=inputfolder, IB.Arclist=Arclist)
saveRDS(islandnames, file=paste(getwd(), "/islands overview",".rds",sep=""))
write.csv(islandnames, file=paste(getwd(),"/islands overview",".csv",sep=""))

### ============================================================
### PALEO AREA: Create shapefiles for various depths 
### ============================================================

# Create DEM for each archipelago: projection system is assigned.
# Use IB.WMcrop.proj.r for high resolution
# par(mar = rep(2, 4)) # Use this command if IB.plot=T returns an error.

for(Arcname in Arclist){ 
  IB.WMcrop.proj.r(IB.Arcname=Arcname, IB.projname="DEMsin", IB.proj=sinus.csy, 
                   IB.WM=DEMworld, resfact = 1, IB.inputfolder=inputfolder, 
                   IB.outputfolder=intermedfolder, IB.plot=F)
}
 
# Create separate folder for each archipelago in outputfolder. 
for(Arcname in Arclist){ 
if(!dir.exists(paste(outputfolder, "/", Arcname, sep=""))){dir.create(paste(outputfolder, "/", Arcname ,sep=""))}
}

# Create shapefiles and assign areas in area.depth
for(Arcname in Arclist){ 
  area.depth<-IB.areachange.pol(IB.Arcname=Arcname, IB.projname="DEMsin", 
              IB.proj=sinus.csy, IB.inputfolder = inputfolder, 
              IB.intermedfolder = intermedfolder, IB.outputfolder = outputfolder, 
                              time.ka=time_level$t1000y[],
              time.sl=time_level$sea_m[])
  saveRDS(area.depth, paste(outputfolder, "/", Arcname, "/", "area_", Arcname, ".rds", sep=""))
  write.csv(area.depth, file=paste(outputfolder, "/", Arcname, "/", "area_", Arcname, ".csv", sep=""))
  }

# Merge files that were just created 
  island.area <- IB.areamerge(IB.Arclist = Arclist, IB.outputfolder = outputfolder)
  island.area[,"island"]<-as.character(island.area[,"island"])
  # Save palaeo areas of all archipelagos in a single csv file
  write.csv(island.area,file = paste(outputfolder, "/","island.parea",".csv",sep = ""), row.names = FALSE)

# # Optional: project shapefiles of t=0 and save as new shapefile:
IB.projpol(IB.inputfolder=inputfolder,
             IB.outputfolder=outputfolder, IB.Arclist=Arclist, IB.projcrs=g.csy)

### ============================================================
###   Clean up dataset: rename some columns, and sort per archipelago
### ============================================================

## paleo area 
  island.parea<-read.csv(file=paste(outputfolder, "/","island.parea",".csv",sep = ""))
  colnames(island.parea)[2]<-"archipelago"

# Sort 
  island.merge.order <- island.parea[order(island.parea[,"archipelago"], 
    island.parea[,"island"]), ]
  
# Save
  write.csv(island.merge.order,file = paste(outputfolder, "/", "area_global", ".csv", sep = ""), 
            row.names = FALSE)
  write.xlsx(island.merge.order,file = paste(outputfolder, "/", "area_global", ".xlsx", sep = ""),
             row.names = FALSE)  
  