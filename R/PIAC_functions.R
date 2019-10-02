###################################################################
### ============================================================
### PIAC Functions
### ============================================================
###################################################################

#   This script contains the functions for preparing the data underlying
#   the study 'A global spatially explicit database of changes in island paleo-area and
#   archipelago configuration during the late Quaternary'
#   by Norder et al.
#   submitted to Global Ecology and Biogeography

### ============================================================
### IB.areachange.pol
### ============================================================

#' IB.areachange.pol creates polygon shapefiles per archipelago at each sea level stand from a sea level curve and a bathymetry DEM. It attaches a name to each polygon based on the names provided by the user in a separate KML file.
#' @param IB.Arcname, the name of the archipelago
#' @param IB.projname, name of the archipelago bathymetry DEM which was produced in IB.WMcrop.proj.r (without name of archipelago, and without file extention).
#' @param IB.proj, coordinate reference system (crs) to which the shapefiles should be projected
#' @param IB.inputfolder, file path to input folder
#' @param IB.intermedfolder, file path to intermediate folder
#' @param IB.outputfolder, file path to output folder
#' @param time.ka, vector with time periods corresponding to time.sl
#' @param time.sl, vector with sea level stands corresponding to time.ka
#' @author Sietze Norder
#' @export
IB.areachange.pol <- function(IB.Arcname, IB.projname, IB.proj, IB.inputfolder,
                             IB.intermedfolder, IB.outputfolder, time.ka, time.sl){

  require(rgdal)
  require(plyr)
  inputarcfolder <- paste(IB.inputfolder, "/", IB.Arcname, "/", sep="")
  outputarcfolder <- paste(IB.outputfolder, "/", IB.Arcname, "/", sep="")
  intermedarcfolder <- paste(IB.intermedfolder, "/", IB.Arcname, "/", sep="")

  # import KML of archipelago
  pntkml <- readOGR(dsn=paste(inputarcfolder, IB.Arcname,"_pnt",".kml",sep=""), layer="ID")
  islandpoints <- spTransform(pntkml, CRS=IB.proj)

  area.depth <- data.frame(row.names = islandpoints@data$Name)
  area.depth[,"grp"]<-IB.Arcname
  for(x in 1:length(time.ka)){
    area.depth[,paste("ka", time.ka[x], sep="")]<-NA
  }

  # Create shapefiles of paleo island extent
  for(i in 1:length(time.sl)){
    ka<-paste("ka", time.ka[i], sep="")
    cat('Calculating area for t=', ka, '\n')
    # loading raster inside loop is faster than loading before the loop and assigning.
    # Load raster
    DEMmask <- raster(x=paste(intermedarcfolder,IB.Arcname,IB.projname,".tif",sep=""))
    # create picture mask (assign NA)
    # assign NA to areas below sea-level (assigning NA reduces computation time).
    DEMmask[DEMmask <= time.sl[i]] <- NA
    DEMmask[DEMmask > time.sl[i]] <- time.ka[i]
    if(suppressWarnings(cellStats(DEMmask,'max')==time.ka[i])){
    temp.arc<- polygonizer(DEMmask, fillholes=TRUE, aggregate=FALSE)
    tmppol <- spTransform(temp.arc, CRS=IB.proj)
    # a data.frame with rows from tmppol and corresponding index of points in islandpoints is returned.
    indexpnt<-over(tmppol, as(islandpoints,"SpatialPoints"))
    # Label polygons based on names in islandpoints
    tmppol@data$iname<-as.character(islandpoints@data$Name[indexpnt])
    #plot(tmppol)
    #text(coordinates(tmppol), tmppol@data$iname,cex=0.9)
    writeOGR(tmppol, paste(outputarcfolder,IB.Arcname,"_shp",sep=""), layer=paste(IB.Arcname,time.ka[i], sep=""), driver="ESRI Shapefile", overwrite_layer = TRUE)
    # write area of each polygon to area.depth
    name.island<-tmppol@data$iname[which(!is.na(tmppol@data$iname))]
    if(!identical(name.island,character(0))){
      area.depth[name.island,ka]<-sapply(slot(tmppol[which(!is.na(tmppol@data$iname)),], "polygons"), slot, "area")/1e6
    }

    ### START of paleo polygons: In case multiple current islands were merged into a paleo island,
    ##  the paleo area should be assigned to all islands in area.depth.
    # a data.frame with rows from pnt.kml.sin and corresponding index of polygons in tmppol is returned.
    indexpol<-over(islandpoints,as(tmppol,"SpatialPolygons"))
    # Check number of occurrences. Only continue if number of paleo islands is larger than 0.
    paleocount<-count(indexpol)
    # Check for each polygon with how many points it overlaps.
    # Return ID's of paleo polygons (Polygons that consist of two or more present-day islands).
    paleopol.id<-paleocount$x[paleocount$freq>1]
    if(any(paleopol.id[!is.na(paleopol.id)])){
      paleopol.id<-paleopol.id[!is.na(paleopol.id)]
      # In theory there can be several of these paleo-islands within one archipelago.
      for (p in 1:length(paleopol.id)){
        paleopol<-paleopol.id[[p]]
        # Get area of paleo polygon
        paleopol.area <- slot(tmppol@polygons[[paleopol]],"area")/1e6
        # Get names of points that fall in same paleo polygon.
        paleopol.names<-islandpoints@data$Name[indexpol==paleopol]
        paleopol.names<-paleopol.names[!is.na(paleopol.names)]
        # Assign area of paleo polygon to area.depth
        area.depth[as.character(paleopol.names),ka]<-paleopol.area
      }
    }
    }
  }
  return(area.depth)
}

### END of function IB.areachange.pol

### ============================================================
### IB.datastructure
### ============================================================

#' IB.datastructure selects island names from KML files per archipelago and stores island names of all archipelagos in a csv file. This is just to provide an overview of which islands are contained in the dataset.
#' @param IB.inputfolder, file path to input folder
#' @param IB.Arclist, list with names of archipelagos
#' @author Sietze Norder
#' @export
IB.datastructure <- function (IB.inputfolder, IB.Arclist){
  IB.islandnames <- data.frame(cbind("island"=as.character(), "archipelago"=as.character()))
  for (IB.Arcname in IB.Arclist){
     inputarcfolder<-file.path(IB.inputfolder,IB.Arcname)
      pntkml <- readOGR(dsn=paste(inputarcfolder,"/",IB.Arcname,"_pnt",".kml",sep=""), layer="ID")
      tmp.islandnames<-cbind("island"=as.character(pntkml@data$Name), "archipelago"=IB.Arcname)
      IB.islandnames<-rbind(IB.islandnames, tmp.islandnames)
      rm(tmp.islandnames)
    }
  return(IB.islandnames)
}

## END of function IB.datastructure

### ============================================================
### IB.projpol
### ============================================================

#' IB.projpol projects shapefiles of t=0 to a user-defined projection and saves it as new shapefile:
#' @param IB.intermedfolder, file path to intermediate folder
#' @param IB.Arclist, list with names of archipelagos
#' @param IB.projcrs, crs to which the shapefile should be projected.
#' @author Sietze Norder
#' @export
IB.projpol <- function (IB.inputfolder, IB.outputfolder, IB.Arclist, IB.projcrs){

  for (IB.Arcname in IB.Arclist){
    inputarcfolder <- paste(IB.inputfolder, "/", IB.Arcname, "/", sep="")
    outputarcfolder <- paste(IB.outputfolder, "/", IB.Arcname, "/", sep="")
    Arcpol <- readOGR(dsn=paste(outputarcfolder,IB.Arcname,"_shp/",IB.Arcname,"0",".shp",sep=""), layer=paste(IB.Arcname,"0", sep=""))
    Arcpol.proj <- spTransform(Arcpol,IB.projcrs)
    writeOGR(obj= Arcpol.proj, dsn=paste(outputarcfolder,IB.Arcname,"Acpol_ll.shp",sep=""), layer=paste(IB.Arcname,"Acpol_ll", sep=""), driver="ESRI Shapefile", overwrite_layer = TRUE)
  }
}

### END of function IB.projpol

### ============================================================
### IB.WMcrop.proj.r
### ============================================================

#' IB.WMcrop.proj.r crops the worldmap to an extent which is slightly larger than the bounding box of an archipelago. It also reprojects the resulting raster to a user-defined crs.
#' @param IB.Arcname, the name of the archipelago
#' @param IB.projname, name of the archipelago bathymetry DEM to be stored (without name of archipelago, and without file extention). The name provided here should correspond to the name provided in the function IB.areachange.pol
#' @param IB.proj, coordinate reference system (crs) to which the raster should be projected
#' @param IB.WM, file path where global bathymetry DEM is stored locally.
#' @param resfact, amount of disaggregation (expressed as number of cells) to create a raster in higher resolution
#' @param IB.inputfolder, file path to input folder
#' @param IB.outputfolder, file path to output folder
#' @param IB.plot, should output be plotted on screen?
#' @param IB.mask, if TRUE: assign NA to areas below sea-level of -150m MSL (this reduces computation time)
#' @author Sietze Norder
#' @export
IB.WMcrop.proj.r <- function(IB.Arcname, IB.projname, IB.proj, IB.WM, resfact, IB.inputfolder, IB.outputfolder, IB.plot=TRUE, IB.mask=FALSE){

  require(raster)
  require(rgdal)

  inputarcfolder<-paste(IB.inputfolder, "/", IB.Arcname, "/",sep="")
  outputarcfolder<-paste(IB.outputfolder, "/", IB.Arcname, "/",sep="")
  # Create outputarcfolder if it doesn't exist yet.
  ifelse(!dir.exists(file.path(outputarcfolder)), dir.create(file.path(outputarcfolder)), FALSE)

  # import KML
  pntkml <- readOGR(dsn=paste(inputarcfolder,IB.Arcname,"_pnt",".kml",sep=""), layer="ID")
  pntkml.proj <- spTransform(pntkml, CRS=IB.proj)

  # Cut worldmap to the extent of bounding box of islands
  #Get extent
  e <- extent(pntkml)
  # Define crop extent (`xmin`, `xmax`, `ymin` , `ymax`)
  cropbox <- c(xmin(e)-1,xmax(e)+1,ymin(e)-1,ymax(e)+1)
  #crop the raster
  Distcrop <- crop(IB.WM, cropbox)
  if(IB.mask){
  # assign NA to areas below sea-level (assigning NA reduces computation time).
  # create picture mask (assign NA)
  Distcrop[Distcrop < -150] <- NA
  }
  # Time taken to project and resample Hawaii with NA assigned: 1.528078 mins. If NA is not assigned: 1.765623 mins (quicker than resample and subsequently project, that takes 3.392157 mins)
  # project
  # start.time <- Sys.time()
  Distcrop.proj<-projectRaster(Distcrop,crs=IB.proj)
  # Resample using bilinear interpolation
  Distcrop.projr <- disaggregate(Distcrop.proj, fact=resfact, method='bilinear')
  writeRaster(x=Distcrop.projr, filename=paste(outputarcfolder,IB.Arcname,IB.projname,".tif",sep=""),overwrite=TRUE)
  # end.time <- Sys.time()
  # time.taken <- end.time - start.time
  # time.taken

  # plot
  if(IB.plot){
    plot(Distcrop.proj)
    plot(pntkml.proj, bg="transparent", add=TRUE)
  }
  print(paste("Raster saved as: ",outputarcfolder,IB.Arcname,IB.projname,".tif",sep=""))
  print(IB.proj)
}

### END of function IB.WMcrop.proj.r

### ============================================================
### IB.areamerge
### ============================================================

#' IB.areamerge combines the csv files of separate archipelagos into a single file.
#' @param IB.Arclist, list with names of archipelagos
#' @param IB.outputfolder, file path to output folder
#' @author Sietze Norder
#' @export
IB.areamerge <- function(IB.Arclist, IB.outputfolder){
area.depth <- data.frame("island"=as.character())

for(IB.Arcname in IB.Arclist){
  tmp.area.depth <- read.csv(file=paste(IB.outputfolder, "/", IB.Arcname, "/", "area_", IB.Arcname,".csv",sep=""))
  colnames(tmp.area.depth)[1]<-"island"
  area.depth<-rbind.fill(area.depth, tmp.area.depth)
}
return(area.depth)
}

### END of function IB.areamerge

### ============================================================
### polygonizer
### ============================================================

#' Polygonizer creates a polygon shapefile from a raster layer
#' @param x: an R Raster layer, or the file path to a raster file recognised by GDAL
#' @param outshape: the path to the output shapefile (if NULL, a temporary file will
#'           be created)
#' @param pypath: the path to gdal_polygonize.py or OSGeo4W.bat (if NULL, the function
#'         will attempt to determine the location). OSGeo4W can be downloaded from
#'         https://trac.osgeo.org/osgeo4w/
#' @param readpoly: should the polygon shapefile be read back into R, and returned by
#'           this function? (logical)
#' @param fillholes: should holes be deleted (i.e., their area added to the containing
#'            polygon)
#' @param aggregate: should polygons be aggregated by their associated raster value?
#' @param quietish: should (some) messages be suppressed? (logical)
#' @author Sietze Norder
#' @export
polygonizer <- function(x, outshape=NULL, pypath=NULL, readpoly=TRUE,
                        fillholes=FALSE, aggregate=FALSE,
                        quietish=TRUE) {
  if (isTRUE(readpoly) || isTRUE(fillholes)) require(rgdal)
  if (is.null(pypath)) {
    cmd <- Sys.which('C:\\OSGeo4W64\\OSGeo4W.bat')
    pypath <- 'gdal_polygonize'
    if(cmd == "") {
      cmd <- "python"
      pypath <- Sys.which('gdal_polygonize.py')
      if (!file.exists(pypath))
        stop("Could not find gdal_polygonize.py or OSGeo4W on your system.")
    }
  }
  if (!is.null(outshape)) {
    outshape <- sub('\\.shp$', '', outshape)
    f.exists <- file.exists(paste(outshape, c('shp', 'shx', 'dbf'), sep='.'))
    if (any(f.exists))
      stop(sprintf('File already exists: %s',
                   toString(paste(outshape, c('shp', 'shx', 'dbf'),
                                  sep='.')[f.exists])), call.=FALSE)
  } else outshape <- tempfile()
  if (is(x, 'Raster')) {
    require(raster)
    writeRaster(x, {f <- tempfile(fileext='.tif')})
    rastpath <- normalizePath(f)
  } else if (is.character(x)) {
    rastpath <- normalizePath(x)
  } else stop('x must be a file path (character string), or a Raster object.')

  system2(cmd, args=(
    sprintf('"%s" "%s" %s -f "ESRI Shapefile" "%s.shp"',
            pypath, rastpath, ifelse(quietish, '-q ', ''), outshape)))

  if(isTRUE(aggregate)||isTRUE(readpoly)||isTRUE(fillholes)) {
    shp <- readOGR(dirname(outshape), layer=basename(outshape),
                   verbose=!quietish)
  } else return(NULL)

  if (isTRUE(fillholes)) {
    poly_noholes <- lapply(shp@polygons, function(x) {
      Filter(function(p) p@ringDir==1, x@Polygons)[[1]]
    })
    pp <- SpatialPolygons(mapply(function(x, id) {
      list(Polygons(list(x), ID=id))
    }, poly_noholes, row.names(shp)), proj4string=CRS(proj4string(shp)))
    shp <- SpatialPolygonsDataFrame(pp, shp@data)
    if(isTRUE(aggregate)) shp <- aggregate(shp, names(shp))
    writeOGR(shp, dirname(outshape), basename(outshape),
             'ESRI Shapefile', overwrite=TRUE)
  }
  if(isTRUE(aggregate) & !isTRUE(fillholes)) {
    shp <- aggregate(shp, names(shp))
    writeOGR(shp, dirname(outshape), basename(outshape),
             'ESRI Shapefile', overwrite=TRUE)
  }
  ifelse(isTRUE(readpoly), return(shp), return(NULL))
}

### END of function polygonizer
