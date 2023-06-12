##  STAC_functions.R
##  Various functions for interacting with STAC 
##  Author: Jasper Van doninck
##
##  With focus on Landsat/Sentinel-2 


stac_endpoint <- function(name){
  #Some endpoints of STAC APIs because I'm too lazy to remember/search
  #Defaults to Planetary Computer
  if(missing(name)) {
    endpoint <- "https://planetarycomputer.microsoft.com/api/stac/v1"
  } else {
    endpoint <- switch(tolower(name),
                       "lpdaac" = "https://cmr.earthdata.nasa.gov/stac/LPCLOUD",
                       "planetarycomputer" = "https://planetarycomputer.microsoft.com/api/stac/v1",
                       stop("Endpoint not implemented")
    )
  }
  return(endpoint)
}


stac_searchWrapper <- function(endpoint,
                               collection=NULL, ids=NULL, bbox=NULL, datetime=NULL, intersects=NULL, limit=NULL, #arguments for stac_search
                               range=NULL, years=NULL, months=NULL, days=NULL, #alternative arguments to datetime to search/filter by dates
                               eocc=NULL, #argument for filtering on 'eo:cloud_cover' field (max)
                               ... #additional arguments to ?
){
  
  # STAC search wrapper
  #
  # Description
  #   Wrapper for formatting search arguments and performing STAC search. Additional filtering for selected months/days
  #
  # Arguments
  #   
  #   ... additional arguments passed to stac_format
  #
  #
  # Returns
  #   STACItemCollection
  #
  # Warning: signing of objects only implemented for Planetary Computer. 
  
  
  require(rstac)
  
  ##Format the arguments for stac_search
  searchArgs <- list()
  #q
  searchArgs$q <- stac(endpoint) #, ...)
  #collection
  searchArgs$collection <- collection
  #ids
  searchArgs$ids <- ids
  #bbox - added check if is spatExtent instead of bbox and correct if needed
  if(!is.null(bbox)) if(class(bbox)=="SpatExtent") bbox <- bbox[c(1,3,2,4)]
  searchArgs$bbox <- bbox
  #datetime - constructed from "datetime", "range" , or "years" (in that order) 
  if(!is.null(datetime)){
    searchArgs$datetime <- datetime
  } else if (!is.null(range)){
    searchArgs$datetime <- paste(format(range, "%Y-%m-%d"), collapse="/")
  } else if (!is.null(years)){
    searchArgs$datetime <- paste(c(min(years), max(years)), c("01", "12"), c("01", "31"), sep="-", collapse=("/"))
    if(!is.null(months)){
      dateMin <- as.Date(paste(min(years), min(months), "01", sep="-"))
      if(max(months)==12){
        dateMax <- as.Date(paste(max(years), max(months), "31", sep="-"), format="%Y-%m-%d") 
      } else{
        dateMax <- as.Date(paste(max(years), max(months)+1, "01", sep="-"), format="%Y-%m-%d")-1
      }
      searchArgs$datetime <- paste(c(dateMin, dateMax), collapse="/")
    }
  }
  #intersects
  searchArgs$intersects <- intersects
  
  if(endpoint=="https://planetarycomputer.microsoft.com/api/stac/v1"){
    #limit
    if(is.null(limit)) limit <- 1000 #max limit for PlanetaryComputer
    if(limit>1000){
      limit <- 1000
      warning("Limit set to 1000 (maximum for Planetary Computer)")
    }
    searchArgs$limit <- limit
    
    query <- do.call(stac_search, searchArgs)
    
    ##  Query extension (rework this to make it more flexible?)
    #Cloud cover:
    if(!is.null(eocc)) query <- ext_query(query, "eo:cloud_cover" < eocc)
    
    ## Post the request
    items <- post_request(query)
    if(length(items$features)==limit) warning(paste0("Limit of ",limit," items reached."))
  }
  
  if(endpoint=="https://cmr.earthdata.nasa.gov/stac/LPCLOUD"){
    searchArgs$limit <- NULL
    query <- do.call(stac_search, searchArgs)
    items <- post_request(query) %>%
      items_fetch()
    if(!is.null(eocc)) items <- items_filter(items, properties$`eo:cloud_cover` < eocc)
  }
  
  ##  Filter based on additional parameters (WARNING: this is after truncation of results by 'limit', might be better to call search several times and merge results)
  ## still searching for a way to do this as ext_query or similar
  if(!is.null(months) & items_length(items)>0){
    s <- items_datetime(items) %>%
      unlist() %>%
      as.Date %>%
      format(format="%m") %>%
      as.numeric() %>%
      '%in%'(months)
    items$features <- items$features[s]
  }
  if(!is.null(days) & items_length(items)>0){
    s <- items_datetime(items) %>%
      unlist() %>%
      as.Date %>%
      format(format="%d") %>%
      as.numeric() %>%
      '%in%'(days)
    items$features <- items$features[s]
  }
  
  return(items)
}



sign_backup <- function(item, pcKey='', silent=TRUE){
##DEPRECATED - REPLACED BY sign_feature_pc  
  
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
  #     silent      boolean indicating whether validity period of key is returned
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
  if(!silent) cat('Token expires in', round(expireTime), units(expireTime),'\n')
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

sign_feature_pc <- function(feature, key=''){
  #   Sign a Microsoft Planetary Computer STAC feature to enable dowmload
  
  library(httr)
  
  if(nchar(pcKey)==0) {
    tok_ret <- httr::GET(paste0('https://planetarycomputer.microsoft.com/api/sas/v1/token/', feature$collection))
  } else {
    tok_ret <- httr::GET(paste0('https://planetarycomputer.microsoft.com/api/sas/v1/token/', feature$collection,
                                '?subscription-key=', pcKey))
  }
  token <- httr::content(tok_ret)$token
  
  feature_assets <- feature$assets
  feature_assets_signed <- lapply(feature_assets, function(x){
    href <- x$href
    href_signed <-  paste0(href,'?',token)
    x$href <- href_signed
    return(x)
  })
  feature$assets <- feature_assets_signed
  return(feature)
}

assets2rast <- function(feature, assets, as_list=FALSE, mapBands = FALSE, sign_feature=NULL, ...){
  #   Read STAC assets as spatRaster
  #
  #   Description
  #    Reads selected assets from a STAC features as a spatRaster object
  #
  #   Usage
  #     assets2rast(feature, assets)
  #
  #   Arguments
  #     feature       a STAC feature: list element of STACItemCollection$features
  #     assets        character vector of assets names to extract
  #     as_list       boolean defining whether assets should be read as list of spatRasters or as multilayer spatRaster 
  #     sign_feature  function used to sign feature (e.g., sign_feature_pc), or NULL
  #     ...           Additional arguments to function in sign_feature
  #   
  #   Details
  #     Assets to read must be SpatRaster format of same extent and dimensions if read into mulitlayer SpatRaster.
  #     Set 'as_list=TRUE' to deal with assets with different extents or dimensions
  #
  #   Value
  #     a (list of) spatRaster object(s)

  library(terra)
  if(!is.null(sign_feature)) feature <- sign_feature(feature, ...)
  if(mapBands) assets <- mapBand(assets, feature$collection)
  fAssets <- feature$assets[which(names(feature$assets) %in% assets)]
  fURLs <- lapply(fAssets, function(x){paste0('/vsicurl/',x$href)})
  ras <- lapply(fURLs, rast)
  ras <- ras[assets]
  if(!isTRUE(as_list)) ras <- rast(ras)
  return(ras)
}

cropSTACFeature <- function(feature, 
                            assets,
                            bbox,
                            extend=TRUE,
                            commonRes="high",
                            as_list=FALSE,
                            ...
                            ){

  #   Crop SpatRaster from STAC feature
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
  #     extend      boolean. Should the cropped area be extended if the bounding box goes beyond the STAC asset's geometry    
  #     commonRes   numeric or character. Defines the common resolution in case speciefied assets have multiple resolutions and as.list=FALSE. 
  #                   numeric: spatial resolution in x and y, in units of feature
  #                   character: one of c("low", "high"); "low" aggregating to coarsest resolution, "high" disaggregating to finest resolution 
  #     as_list     boolean. If multiple assets specified should they be returned as a list of SpatRaster instead of a multilayer SpatRaster
  #     ...         Additional arguments passed to assets2rast (for signing feature)
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
  library(magrittr)
  terraOptions(progress=0)

  ##  Check input parameters, define default values
    # feature
  if(missing(feature)) stop("feature must be provided")
  # if(!"list" %in% class(feature)) stop("feature must be list of type Feature")
  # if(feature$type!="Feature") stop("feature must be list of type Feature")
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
  featBands <- assets2rast(feature, assets, as_list=TRUE, ...)
  
  ##  Reproject bounding box to crs of feature
  cropExt <- terra::project(ext(bbox, xy=TRUE), crs("+proj=longlat +datum=WGS84"), crs(featBands[[1]]))

  ##  Adjust crop extent to match asset with coarsest resolution
  coarseRast <- sapply(featBands, function(x) res(x)[1]) %>% 
    which.max()
  cropExt <- align(cropExt, featBands[[coarseRast]], snap="out")

  ##  Crop assets
  cropFun <- function(b, cropExt, extend){
    #TO DO: Add extra check if extents of b and cropExt overlap # If no overlap, this will throw an error in the current version of terra. This should be replaced by a warning and returning an empty raster (with correct properties) 
    cr <- crop(b,cropExt, snap="near")
    if(isTRUE(extend)) cr <- extend(cr, cropExt) #Crop with extend=TRUE is bugged in older versions of terra, so using this approach for now. Should be updated
    
    #set "time" field
    if(!is.null(feature$properties$datetime)){
      time(cr) <- rep(as.POSIXlt(feature$properties$datetime, format="%Y-%m-%dT%H:%M:%S", tz="GMT"), nlyr(cr))
    }
    return(cr)
  }
  croppedBands <- lapply(featBands, cropFun, cropExt, extend)
  
  if(!isTRUE(as_list)){
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
    
    #Extents may differ in case of different resolutions, add check to avoid errors  
    if(length(croppedBands)>1){
      if(!isTRUE(do.call(compareGeom, c(unname(croppedBands),stopOnError=FALSE)))){
        croppedBands <- lapply(croppedBands, crop, cropExt)
      } 
    }
    croppedBands <- rast(croppedBands)
  }
  return(croppedBands)
}

cropSTACFeatureCollection <- function(FeatureCollection,
                                      assets,
                                      bbox,
                                      progress=TRUE,
                                      ...){
  # Read STAC Feature Collection
  #
  # Description
  #   Read STAC FeatureCollection as list STAC features
  #
  # Arguments
  #   stacFeatures    list. STAC FeatureCollection
  #   assets      character vector of assets names to extract
  #   bbox        numeric vector representing bounding box (xmin, ymin, xmax, ymax) in geographic coordinates
  #   
  #   ...             additional arguments passed to cropSTACFeature (assets, bbox, extend, as.list, commonRes)
  #
  # Value
  #   list of SpatRaster objects
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
    rastFeatures[[f]] <- tryCatch(cropSTACFeature(FeatureCollection$features[[f]], assets, bbox, ...), 
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



#not strictly stac-related
mapBand <- function(bands, collection){
  #Re-assign band names to combine Landsat/Sentinel-2 more easily
  
  mapFun <- function(band, collection){
    
    #Landsat (for planetary computer, check for other catalogs)
    if(collection=="landsat-c2-l2"){
      assetName <- switch(band,
                          "coastal" = "coastal", "coastal aerosol" = "coastal",
                          "blue" = "blue",
                          "green" = "green",
                          "red" = "red",
                          "nir" = "nir08",
                          "swir1" = "swir16",
                          "swir2" = "swir22",
                          "qa" ="qa_pixel",
                          band)
    }
    
    #HLS Landsat (for LPDAAC, check for other catalogs)
    if(collection=="HLSL30.v2.0"){
      assetName <- switch(band,
                          "B01 "= "B01", "coastal" = "B01", "coastal aerosol" = "B01",
                          "B02" = "B02", "blue" = "B02", 
                          "B03" = "B03", "green" = "B03",
                          "B04" = "B04", "red" = "B04", 
                          "B05" = "B05", "nir" = "B05",
                          "B06" = "B06", "swir1" = "B06",
                          "B07" = "B07", "swir2" = "B07",
                          "B09" = "B09", "cirrus" = "B09",
                          "B10" = "B10", "thermal infrared 1" = "B10", "thermal1" = "B10",
                          "B11" = "B11", "thermal" = "B11", "thermal2" = "B11",
                          band)
    } 
    
    #Sentinel (PC), HLS Sentinel (LPDAAC)  
    if(collection=="sentinel-2-l2a" | collection=="HLSS30.v2.0"){
      assetName <- switch(band,
                          "B01 "= "B01", "coastal" = "B01", "coastal aerosol" = "B01",
                          "B02" = "B02", "blue" = "B02", 
                          "B03" = "B03", "green" = "B03",
                          "B04" = "B04", "red" = "B04", 
                          "B05" = "B05", "vre1" = "B05", "red-edge 1" = "B05",
                          "B06" = "B06", "vre2" = "B06", "red-edge 2" = "B06",
                          "B07" = "B07", "vre3" = "B07", "red-edge 3" = "B07",
                          "B08" = "B08", "nir" = "B08", "nir broad" = "B08",
                          "B8A" = "B8A", "nir narrow" = "B8A",
                          "B09" = "B09", "water vapor" = "B09", "water vapour" = "B09",
                          "B10" = "B10", "cirrus" = "B10",
                          "B11" = "B11", "swir1" = "B11", "swir 1" = "B11",
                          "B12" = "B12", "swir2" = "B12", "swir 2" = "B12",
                          "SCL" = "SCL", "qa" = "SCL",
                          band)
    }
    return(assetName)
  }
  
  assets <- sapply(bands, mapFun, collection)
  return(assets)
}

