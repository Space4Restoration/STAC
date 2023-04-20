##  STAC_functions.R
## Helper functions for interacting with STAC 
##  Author: Jasper Van doninck
##

sign_backup <- function(item, pcKey=''){
  #   Sign STACItemCollection 
  #
  #   Description
  #     Alternative method for singing a STACItemCollection in Planetary Computer
  #     Can be used if the functions 'items_sign' and 'sign_planetary_computer' from the 'rstac' package throw an error
  #
  #   Usage
  #     sign_backup(item, collection, pcKey=NULL)
  #
  #   Arguments
  #     item        a STACItemCollection object representing the results of /stac/search
  #     pcKey       a character representing a personal Planetary Computer key. Not required
  #
  #   Value
  #     a STACItemCollection with downloadable urls

  library(httr)
  
  collection <- sapply(item$features, function(x){x$collection}) %>%
    as.factor() %>%
    levels()
  if(length(collection)>1) stop('Multiple collection IDs in item')

  if(nchar(pcKey)==0) {
    tok_ret <- httr::GET(paste0('https://planetarycomputer.microsoft.com/api/sas/v1/token/', collection)) %>%         
      httr::content() 
  } else {
    tok_ret <- httr::GET(paste0('https://planetarycomputer.microsoft.com/api/sas/v1/token/', collection,
                                '?subscription-key=', pcKey)) %>%         
      httr::content() 
  }
  
  expireTime <- as.POSIXlt(tok_ret$`msft:expiry`, format="%Y-%m-%dT%H:%M:%S", tz = "UTC") - 
    as.POSIXlt(Sys.time(), tz = "UTC")
  cat('Token expires in', round(expireTime), units(expireTime),'\n')
  token <- tok_ret$token
  
  feats <- item$features
  signFeature <- function(feat){
    assets <- feat$assets
    assets_signed <- lapply(assets, function(x){
      href <- x$href
      href_signed <-  paste0(href,'?',token)
      x$href <- href_signed
      return(x)
    })
    feat$assets <- assets_signed
    return(feat)
  }
  feats_signed <- lapply(feats, signFeature)
  item$features <- feats_signed
  return(item)
}

assets2rast <- function(feature, assets, as.list=FALSE){
  #   Read STAC assets as spatRaster
  #
  #   Description
  #    Reads selected assets from a STAC features as a spatRaster object
  #
  #   Usage
  #     assets2rast(feature, assets)
  #
  #   Arguments
  #     feature     a STAC feature: list element of STACItemCollection$features
  #     assets      character vector of assets names to extract
  #     as.list     boolean defining whether assets should be read as list of spatRasters or as multilayer spatRaster 
  #   
  #   Details
  #     Assets to read must be SpatRaster format of same extent and dimensions if read into mulitlayer SpatRaster.
  #     Set 'as.list=TRUE' to deal with assets with different extents or dimensions
  #
  #   Value
  #     a (list of) spatRaster object(s)

  library(httr)
  library(terra)
  fAssets <- feature$assets[which(names(feature$assets) %in% assets)]
  fURLs <- lapply(fAssets, function(x){paste0('/vsicurl/',x$href)})
  ras <- lapply(fURLs, rast)
  if(!isTRUE(as.list)) ras <- rast(ras)
  return(ras)
}

cropSTACFeature <- function(feature, 
                            assets,
                            bbox,
                            extend=TRUE,
                            as.list=FALSE,
                            commonRes="high"){

  #   Crop spatRraster from STAC feature
  #
  #   Description
  #    Extract spatial subset (crop) from a STAC feature assets
  #
  #   Usage
  #     cropSTACFeature(feature, assets=NULL, bbox, collection, commonRes="low")
  #
  #   Arguments
  #     feature     a STAC feature: list element of STACItemCollection$features
  #     assets      character vector of assets names to extract
  #     bbox        numeric vector representing bounding box (xmin, ymin, xmax, ymax) in geographic coordinates
  #     commonRes   numeric or character. Defines the common resolution in case speciefied assets have multiple resolutions and as.list=FALSE. 
  #                   numeric: spatial resolution in x and y, in units of feature
  #                   character: one of c("low", "high"); "low" aggregating to coarsest resolution, "high" disaggregating to finest resolution 
  #     as.list     boolean. If multiple assets specifiedm should they be returned as a list of SpatRaster instead of a multilayer SpatRaster
  #     extend      boolean. Should the cropped area be extended if the bounding box goes beyond the STAC asset's geometry    
  #
  #   Details
  #     Warning! If spatial resolutions of assets are not all multiples of each others, returning as multilayer spatRaster will result in error
  #
  #   Value
  #     a spatRaster object
  #
  #   To do
  #     handling error when no overlap between assets and bbox
  #
  
  library(terra)
  terraOptions(progress=0)

  ##  Check input parameters, define default values
    # feature
  if(missing(feature)) stop("feature must be provided")
  if(!"list" %in% class(feature)) stop("feature must be list of type Feature")
  if(feature$type!="Feature") stop("feature must be list of type Feature")
    #assets
    #If argument assets is missing, use all assets of type "image/tiff; application=geotiff; profile=cloud-optimized"
  if(missing(assets)){
    warning("assets not provided, all available assets selected")
    assets <- names(feature$assets)[lapply(feature$assets, function(x) x$type)=="image/tiff; application=geotiff; profile=cloud-optimized"]
  }
    #bbox
  if(missing(bbox)) bbox <- feature$bbox
  if(class(bbox)=="SpatExtent") bbox <- bbox[c(1,3,2,4)]  #Check if bbox is already a SpatExtent object

  ##  Read assets as spatRaster object
  featBands <- assets2rast(feature, assets, TRUE)
  
  ##  Reproject bounding box to crs of feature
  featCRS <- crs(paste0('EPSG:',feature$properties$`proj:epsg`))
  cropExt <- terra::project(ext(bbox, xy=TRUE), crs("+proj=longlat +datum=WGS84"), featCRS)

  ##  Adjust crop extent to match asset with coarsest resolution
  coarseRast <- sapply(featBands, function(x) res(x)[1]) %>% 
    which.max()
  cropExt <- align(cropExt, featBands[[coarseRast]], snap="out")

  ##  Crop assets
  cropFun <- function(b, cropExt, extend){
    #TO DO: Add extra check if extents of b and cropExt overlap
    # If no overlap, this will throw an error in the current version of terra. This should be replaced by a warning and returning an empty raster (with correct properties) 
    cr <- crop(b,cropExt, snap="near")
    if(isTRUE(extend)) cr <- extend(cr, cropExt) #Crop with extend=TRUE is bugged in older versions of terra, so using this approach for now. Should be updated
    time(cr) <- rep(as.POSIXlt(feature$properties$datetime, format="%Y-%m-%dT%H:%M:%S", tz="GMT"), nlyr(cr))
    return(cr)
  }
  croppedBands <- lapply(featBands, cropFun, cropExt, extend)
  
  if(!isTRUE(as.list)){
    # Check spatial resolutions of assets
    resolutions <- sapply(featBands, function(x) res(x)[1])
    resLevels <- as.numeric(levels(as.factor(resolutions)))
    if(length(resLevels)>1){
      #Assets with different resolution: harmonize assets to common resolution (assuming multiples)
      # Target resolution: either specified or highest/lowest
      if(is.numeric(commonRes)){
        targetRes <- commonRes
      } else {
        targetRes <- switch(commonRes,
                            low=max(resLevels),
                            high=min(resLevels))
      }
      fun <- function(x){
        if(res(x)[1]==targetRes){
          return(x)
        } else if(res(x)[1]>targetRes){
          #disaggregate
          return(disagg(x, fact=res(x)[1]/targetRes, method="near"))
        } else if(res(x)[1]<targetRes){
          #aggregate
          return(aggregate(x, fact=targetRes/res(x)[1], fun="mean"))
        }
      }
      croppedBands <- lapply(croppedBands, fun)

    }
  }
  #Extents may differ in case of different resolutions, add check to avoid errors  
  if(!isTRUE(do.call(compareGeom, c(unname(croppedBands),stopOnError=FALSE)))){
    croppedBands <- lapply(croppedBands, crop, cropExt)
  } 
  croppedBands <- rast(croppedBands)
  return(croppedBands)
}

cropSTACFeatureCollection <- function(FeatureCollection,
                                      assets,
                                      bbox,
                                      extend=TRUE,
                                      as.list=FALSE,
                                      commonRes="high",
                                      progress=TRUE){
  # Read STAC Feature Collection
  #
  # Description
  #   Read STAC FeatureCollection as list STAC features
  #
  # Arguments
  #   stacFeatures    list. STAC FeatureCollection
  #   ...             additional arguments passed to cropSTACFeature (assets, bbox, extend, as.list, commonRes)
  #
  # Value
  #   list of spatRaster objects
  #
  
  library(terra)
  library(httr)
  
  ##  Check input parameters
  if(missing(FeatureCollection)) stop("FeatureCollection must be provided")
  if(!"list" %in% class(FeatureCollection)) stop("FeatureCollection must be list of type FeatureCollection")
  if(FeatureCollection$type!="FeatureCollection") stop("feature must be list of type FeatureCollection")
  
  nFeat <- length(FeatureCollection$features)
  rastFeatures <- vector(mode='list', length=nFeat)
  if(progress) pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                                    max = nFeat,  # Maximum value of the progress bar
                                    style = 3,    # Progress bar style (also available style = 1 and style = 2)
                                    width = 50,   # Progress bar width. Defaults to getOption("width")
                                    char = "=")   # Character used to create the bar
  for(f in 1:nFeat){  # ITeratively get the features
    rastFeatures[[f]] <- tryCatch(cropSTACFeature(FeatureCollection$features[[f]],
                                                  assets, bbox, extend, as.list, commonRes), 
                                  error=function(e) NULL)
    if(progress) setTxtProgressBar(pb, f)
  }
  if(progress) close(pb)
  failed <- sapply(rastFeatures, is.null) #check for failed downloads
  if(sum(failed) > 0) {                   #remove from list
    warning(paste0(sum(failed)," failed download(s), feature(s) ignored"))
    rastFeatures <- rastFeatures[!failed]
  }
  return(rastFeatures)
}

