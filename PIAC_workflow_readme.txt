		# WORKFLOW FOR PRODUCING THE PIAC DATABASE

The workflow for producing the Paleo Islands and Archipelago Configuration (PIAC) database is described in the data paper 'A global spatially explicit database of changes in island paleo-area and   archipelago configuration during the late Quaternary' by Norder et al. in Global Ecology and Biogeography. Please refer to this article when using the data. 

To enable users to make reconstructions of paleo configuration for other archipelagos we also share the R scripts (PIAC_run.R and PIAC_functions.R) to generate the data. More details about the R scripts and directory structure are provided below.

		# DIRECTORY STRUCTURE AND SCRIPTS

-The input folder contains:
	* A folder for each archipelago containing a KML file with island center points. 
	* A folder 'GLOBAL' containing the sea level curves from Cutler et al. (2003) and Lambeck et al. (2014). In case users wish to run the script for another archipelago, the  GEBCO_2014 Grid should be downloaded to this folder from 'https://www.gebco.net/data_and_products/gridded_bathymetry_data/'. 
	* Archipelago_pnt.KML, an empty kml file for users who wish to calculate area change and archipelago configurations for an archipelago which is not in the PIAC database.

-The intermediate folder contains:
	* A folder for each archipelago with a bathymetry DEM (.tif format). 

-The output folder contains:
	* A folder for each archipelago containing: a folder (archipelago_shp) with polygon shapefiles for each time interval, a shapefile of the current archipelago configuration ('archipelagoAcpol_ll.shp'), the area in km2 of all islands within that archipelago for each 1 ka time interval (area_archipelago.CSV and area_archipelago.RDS).
	* area_global.CSV and area_global.xlsx, the area in km2 of all islands in the PIAC database for each 1 ka time interval.
	
-islands overview.CSV and islands overview.RDS - an overview of the names of islands in the PIAC database. 
-PIAC_run.R - R script containing the workflow for preparing the data.
-PIAC_functions.R - R script containing the functions for preparing the data. To be able to use the function 'polygonizer' (in PIAC_functions.R), OSGeo4W should already be installed on your computer (preferrably in C:/OSGeo4W64/OSGeo4W.bat). If it is not yet installed, it can be downloaded from https://trac.osgeo.org/osgeo4w/
