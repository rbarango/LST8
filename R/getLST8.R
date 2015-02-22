#######################################################################################
#' @author Rodolfo de Benito Arango, \email{rbarango@@gmail.com}
#' @references \url{https://landsat.usgs.gov/Landsat8_Using_Product.php}
#' @references \url{http://landsat.usgs.gov/band_designations_landsat_satellites.php}
#' @references \url{http://earthexplorer.usgs.gov}
#' @keywords Landsat8 OLI TIRS TOA
#'
#' @title
#' Retrieves OLI & TIRS raw values and calculates TOA and BRIGHTNESS TEMPERATURE
#' 
#' @description
#' \code{getLST8} Process Level 1 GeoTIFF Data Product of Landsat 8 and converts
#' to REFLECTANCE TOA, RADIANCE TOA and BRIGHTNESS TEMPERATURE.
#' 
#' # Params
#' @param dir The folder or directory where the Landsat 8 GEOTIF files are stored.
#' @param geometry The geometry object of the geografical zone to retrieve the data from.
#' 
#' @details
#' The procedure considers as input (1) the folder where the GEOTIF data products are
#' stored taking into account that each data product must be decompress in one folder 
#' and (2) a geometry object corresponding to the geografical zone to retrieve the data from
#'  (for instance a SpatialPolygonsDataFrame)
#' The output of this procedure is a CSV file including the spatial data points of the extent with the
#' raw data of the OLI (Band B8 panchromatic is not included) and TIRS bands and their corresponding reflectance 
#' and radiance TOA values for the spectral bands and, for thermal bands, their brightness temperature. 
#' The name of the CSV file is LST8_startDate_endDate.
#' This procedure has dependences on rgdal, raster and sp R packages. 
#' The procedure do not download the data products. The data products can be downloaded from
#' http://earthexplorer.usgs.gov/ 
#' 
#' The structure of the CSV file is the following:
#' 
#' x: Coordinate x of the data point in the CRS of the data product
#' y: Coordinate y of the data point in the CRS of the data product
#' crs_proj: CRS projection string of the data product
#' date: Year + Day number of the year in the format YYYYddd
#' 30_B1: OLI band 1
#' 30_B2: OLI band 2
#' 30_B3: OLI band 3
#' 30_B4: OLI band 4
#' 30_B5: OLI band 5
#' 30_B6: OLI band 6
#' 30_B7: OLI band 7
#' 30_B9: OLI band 9
#' 30_B10: TIRS band 10
#' 30_B11: TIRS band 11
#' 30_QA: Quality indicator of the measures
#' 30_B1_RADIANCE: RADIANCE (TOA) for band 1
#' 30_B1_REFLECTANCE: REFLECTANCE (TOA) for band 1
#' 30_B2_RADIANCE: RADIANCE (TOA) for band 2
#' 30_B2_REFLECTANCE: REFLECTANCE (TOA) for band 2
#' 30_B3_RADIANCE: RADIANCE (TOA) for band 3
#' 30_B3_REFLECTANCE: REFLECTANCE (TOA) for band 3
#' 30_B4_RADIANCE: RADIANCE (TOA) for band 4
#' 30_B4_REFLECTANCE: REFLECTANCE (TOA) for band 4
#' 30_B5_RADIANCE: RADIANCE (TOA) for band 5
#' 30_B5_REFLECTANCE: REFLECTANCE (TOA) for band 5
#' 30_B6_RADIANCE: RADIANCE (TOA) for band 6
#' 30_B6_REFLECTANCE: REFLECTANCE (TOA) for band 6
#' 30_B7_RADIANCE: RADIANCE (TOA) for band 7
#' 30_B7_REFLECTANCE: REFLECTANCE (TOA) for band 7
#' 30_B9_RADIANCE: RADIANCE (TOA) for band 9
#' 30_B9_REFLECTANCE: REFLECTANCE (TOA) for band 9
#' 30_B10_TEMPERATURE: BRIGHTNESS TEMPERATURE (TOA) for band 11
#' 30_B11_TEMPERATURE: BRIGHTNESS TEMPERATURE (TOA) for band 12


getLST8  <- function(dir, geometry){
    
    # Check parameters
    if (missing(dir))
        stop("Provide the directory where the GEOTIF files are located")
    if(missing(geometry))
        stop("Provide the geometry of the area to retrieve")
    
    dirBase  <- dir
        
    layerNames  <- c("B1","B2","B3","B4","B5","B6","B7","B8","B9","B10","B11","BQA")
    
    # Initialization of variables 
    for (i in layerNames[1:9]) {
        assign(paste0("S_",i), NULL)
        assign(paste0("S_",i,".RADIANCE"), NULL)
        assign(paste0("S_",i,".REFLECTANCE"), NULL)
    }
    
    QA <- NULL
    S_B10 <- NULL
    S_B11 <- NULL
    S_B10.RADIANCE <- NULL
    S_B11.RADIANCE <- NULL
    S_B10.TEMPERATURE <- NULL
    S_B11.TEMPERATURE <- NULL
    
    coordinates <- NULL
    crs_proj  <- NULL
    dataProduct  <- "LST8"
    
    dirNames <- list.dirs(dirBase)
    
    if (length(dirNames)==1) 
        stop(paste(folder,"do not contain any directory"))
    
    dirNames  <- dirNames[2:length(dirNames)]
    
    dateList  <- NULL
    
    # Process each folder (one folder corresponds to one day)
    for (folder in dirNames){
        df <- NA
        day  <- substr(folder,nchar(folder)-11,nchar(folder)-5)
        
        # Reads metadata file
        metafile  <- list.files(path=folder, pattern="_MTL.txt")
        metadata<-read.table(file.path(folder,metafile), sep=c("="), header=F, stringsAsFactors=F, fill=T, strip.white=T) 
        
        
        # Reads projection of the spatial data points
        map_projection  <- metadata[metadata$V1=="MAP_PROJECTION",][[2]]
        datum  <- metadata[metadata$V1=="DATUM",][[2]]
        utm_zone  <- metadata[metadata$V1=="UTM_ZONE",][[2]]
        lst8_proj  <- paste0("+proj=",tolower(map_projection)," +zone=",utm_zone," +datum=",datum)
        
        # Get the constant values in metadata file for the conversions
        for (i in 1:11){
            assign(paste0("RADIANCE_MULT_BAND_",i), as.numeric(metadata[metadata$V1==paste0("RADIANCE_MULT_BAND_",i),][[2]]))
            assign(paste0("RADIANCE_ADD_BAND_",i), as.numeric(metadata[metadata$V1==paste0("RADIANCE_ADD_BAND_",i),][[2]]))
        }
        
        for (i in 1:9){
            assign(paste0("REFLECTANCE_MULT_BAND_",i), as.numeric(metadata[metadata$V1==paste0("REFLECTANCE_MULT_BAND_",i),][[2]]))
            assign(paste0("REFLECTANCE_ADD_BAND_",i), as.numeric(metadata[metadata$V1==paste0("REFLECTANCE_ADD_BAND_",i),][[2]]))
        }
        
        K1_CONSTANT_BAND_10 <- as.numeric(metadata[metadata$V1=="K1_CONSTANT_BAND_10",][[2]])
        K1_CONSTANT_BAND_11 <- as.numeric(metadata[metadata$V1=="K1_CONSTANT_BAND_11",][[2]])
        K2_CONSTANT_BAND_10 <- as.numeric(metadata[metadata$V1=="K2_CONSTANT_BAND_10",][[2]])
        K2_CONSTANT_BAND_11 <- as.numeric(metadata[metadata$V1=="K2_CONSTANT_BAND_11",][[2]])
        
        SUN_ELEVATION  <- as.numeric(metadata[metadata$V1=="SUN_ELEVATION",][[2]])
        
        
        # Process each layer of the TIFF file
        for(layer in layerNames) {
            
            # Find names of all the tiff files 
            files  <-  list.files(path = folder,  pattern=paste0("_",layer,"\\.TIF"), recursive=F)
            
            # Create a raster stack for each band and day
            dataproduct= stack(file.path(folder,files))
            
            # Convert input geometry to the projection of LST8
            geometry_lst8  <- spTransform(geometry, CRS(lst8_proj))
            
            # Crop out the study region
            dataproduct  <-  crop(dataproduct,extend(extent(geometry_lst8),50))
            
            #dataproduct  <- projectRaster(dataproduct_LST8, crs = crs(proj4string(geometry)))
            
            message(paste0("Performing ",dataProduct,"_",layer,"_",day))
            
            # Save data as an R object
            assign(layer,dataproduct)
        }
        
        # CONVERSION TO RADIANCE (TOA) 
        for (i in 1:11){
            assign(paste0("B",i,".RADIANCE"), eval(parse(text=paste0("RADIANCE_MULT_BAND_"
                                                                     ,i, " * values(B",i,") + RADIANCE_ADD_BAND_",i))))
        }
        
        # CONVERSION TO REFLECTANCE TOA WITH ANGULAR CORRECTION
        for (i in c(1:7,9)){
            assign(paste0("B",i,".REFLECTANCE"), eval(parse(text=paste0("(REFLECTANCE_MULT_BAND_"
                                                                        ,i, " * values(B",i,") + REFLECTANCE_ADD_BAND_",i,
                                                                        ") / sin(SUN_ELEVATION)"))))
        }
                
        # CONVERSION TO BRIGHTNESS TEMPERATURE
        B10.TEMPERATURE  <-  (K2_CONSTANT_BAND_10 / log((K1_CONSTANT_BAND_10/B10.RADIANCE) +1)) - 273.15
        B10.TEMPERATURE  <-  matrix(ncol=1,data=B10.TEMPERATURE)
        
        B11.TEMPERATURE  <-  (K2_CONSTANT_BAND_11 / log((K1_CONSTANT_BAND_11/B11.RADIANCE) +1)) - 273.15
        B11.TEMPERATURE  <-  matrix(ncol=1,data=B11.TEMPERATURE)
        
        # Saving the data 
        coordinates <- rbind(coordinates, coordinates(B1))
        S_B1 <- rbind(S_B1, values(B1))
        S_B2 <- rbind(S_B2, values(B2))
        S_B3 <- rbind(S_B3, values(B3))
        S_B4 <- rbind(S_B4, values(B4))
        S_B5 <- rbind(S_B5, values(B5))
        S_B6 <- rbind(S_B6, values(B6))
        S_B7 <- rbind(S_B7, values(B7))
        #S_B8 <- rbind(S_B8, values(B8))
        S_B9 <- rbind(S_B9, values(B9))
        S_B10 <- rbind(S_B10, values(B10))
        S_B11 <- rbind(S_B11, values(B11))
        QA <- rbind(QA, values(BQA))
        S_B1.RADIANCE <- rbind(S_B1.RADIANCE, B1.RADIANCE)
        S_B1.REFLECTANCE <- rbind(S_B1.REFLECTANCE, B1.REFLECTANCE)
        S_B2.RADIANCE <- rbind(S_B2.RADIANCE, B2.RADIANCE)
        S_B2.REFLECTANCE <- rbind(S_B2.REFLECTANCE, B2.REFLECTANCE)
        S_B3.RADIANCE <- rbind(S_B3.RADIANCE, B3.RADIANCE)
        S_B3.REFLECTANCE <- rbind(S_B3.REFLECTANCE, B3.REFLECTANCE)
        S_B4.RADIANCE <- rbind(S_B4.RADIANCE, B4.RADIANCE)
        S_B4.REFLECTANCE <- rbind(S_B4.REFLECTANCE, B4.REFLECTANCE)
        S_B5.RADIANCE <- rbind(S_B5.RADIANCE, B5.RADIANCE)
        S_B5.REFLECTANCE <- rbind(S_B5.REFLECTANCE, B5.REFLECTANCE)
        S_B6.RADIANCE <- rbind(S_B6.RADIANCE, B6.RADIANCE)
        S_B6.REFLECTANCE <- rbind(S_B6.REFLECTANCE, B6.REFLECTANCE)
        S_B7.RADIANCE <- rbind(S_B7.RADIANCE, B7.RADIANCE)
        S_B7.REFLECTANCE <- rbind(S_B7.REFLECTANCE, B7.REFLECTANCE)
        S_B9.RADIANCE <- rbind(S_B9.RADIANCE, B9.RADIANCE)
        S_B9.REFLECTANCE <- rbind(S_B9.REFLECTANCE, B9.REFLECTANCE)
        S_B10.RADIANCE <- rbind(S_B10.RADIANCE, B10.RADIANCE)
        S_B11.RADIANCE <- rbind(S_B11.RADIANCE, B11.RADIANCE)
        S_B10.TEMPERATURE <- rbind(S_B10.TEMPERATURE, B10.TEMPERATURE)
        S_B11.TEMPERATURE <- rbind(S_B11.TEMPERATURE, B11.TEMPERATURE)
        
        for (n in values(B1)){
            dateList  <- rbind(dateList,day)
            crs_proj  <- rbind(crs_proj, lst8_proj)
        }
    }
    
    df.sur_refl <- suppressWarnings(data.frame(coordinates=coordinates, crs_proj, date=dateList, b1=S_B1, 
                              B2=S_B2, B3=S_B3, B4=S_B4,
                              B5=S_B5, B6=S_B6, B7=S_B7,
                              #B8=S_B8,
                              B9=S_B9, B10=S_B10, B11=S_B11,
                              QA=QA,
                              B1.RADIANCE=S_B1.RADIANCE, B1.REFLECTANCE=S_B1.REFLECTANCE,
                              B2.RADIANCE=S_B2.RADIANCE, B2.REFLECTANCE=S_B2.REFLECTANCE,
                              B3.RADIANCE=S_B3.RADIANCE, B3.REFLECTANCE=S_B3.REFLECTANCE,
                              B4.RADIANCE=S_B4.RADIANCE, B4.REFLECTANCE=S_B4.REFLECTANCE,
                              B5.RADIANCE=S_B5.RADIANCE, B5.REFLECTANCE=S_B5.REFLECTANCE,
                              B6.RADIANCE=S_B6.RADIANCE, B6.REFLECTANCE=S_B6.REFLECTANCE,
                              B7.RADIANCE=S_B7.RADIANCE, B7.REFLECTANCE=S_B7.REFLECTANCE,
                              B9.RADIANCE=S_B9.RADIANCE, B9.REFLECTANCE=S_B9.REFLECTANCE,
                              B10.RADIANCE=S_B10.RADIANCE, B11.RADIANCE=S_B11.RADIANCE,
                              B10.TEMPERATURE=S_B10.TEMPERATURE, B11.TEMPERATURE=S_B11.TEMPERATURE
    ))
    
    
    
    names(df.sur_refl)  <- c("x","y", "crs_proj","date", "30_B1", "30_B2","30_B3", "30_B4",
                             "30_B5", "30_B6","30_B7", 
                             #"15_refl_b08",
                             "30_B9", "30_B10","30_B11", 
                             "30_QA",
                             "30_B1_RADIANCE","30_B1_REFLECTANCE",
                             "30_B2_RADIANCE","30_B2_REFLECTANCE",
                             "30_B3_RADIANCE","30_B3_REFLECTANCE",
                             "30_B4_RADIANCE","30_B4_REFLECTANCE",
                             "30_B5_RADIANCE","30_B5_REFLECTANCE",
                             "30_B6_RADIANCE","30_B6_REFLECTANCE",
                             "30_B7_RADIANCE","30_B7_REFLECTANCE",
                             "30_B9_RADIANCE","30_B9_REFLECTANCE",
                             "30_B10_RADIANCE","30_B11_RADIANCE",
                             "30_B10_TEMPERATURE","30_B11_TEMPERATURE"
    )
    
    startDate  <- min(dateList)
    endDate  <- max(dateList)
    message(paste0("Writing CSV file..."))
    write.csv(df.sur_refl, file=paste0(dataProduct,"_",startDate,"_",endDate,".csv"), row.names=F)
}