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
  #     Assets to download must be raster format of with same extent and dimensions if read into mulitlayer spatRaster.
  #     Set 'as.list=TRUE' to deal with assets with different extents or dimensions
  #
  #   Value
  #     a spatRaster object

  library(httr)
  library(terra)
  fAssets <- feature$assets[which(names(feature$assets) %in% assets)]
  fURLs <- lapply(fAssets, function(x){paste0('/vsicurl/',x$href)})
  ras <- lapply(fURLs, rast)
  if(!isTRUE(as.list)) ras <- rast(ras)
  return(ras)
}

cropFromSTAC <- function(feature, 
                         bbox=feature$bbox, 
                         assets=NULL,
                         collection=feature$collection,
                         commonRes="low",
                         extend=FALSE){
  #   Crop raster from STAC feature
  #
  #   Description
  #    Extract spatial subset (crop) from a STAC feature assets
  #
  #   Usage
  #     cropFromSTAC(feature, assets=NULL, bbox, collection, commonRes="low")
  #
  #   Arguments
  #     feature     a STAC feature: list element of STACItemCollection$features
  #     assets      character vector of assets names to extract
  #     bbox        numeric vector representing bounding box (xmin, ymin, xmax, ymax) in geographic coordinates
  #     collection  character
  #     commonRes   character. If assets have multiple resolutions: "low" aggregating to coarsest resolution, "high" for disaggregating to finest resolution 
  #   
  #   Details
  #     
  #
  #   Value
  #     a spatRaster object
  
  library(terra)
  
  ##  Set default assets for Landsat/Sentinel-2, if not provided as arguments
  if(is.null(assets)){
    if(collection=='landsat-c2-l2'){
      assets <- c("blue" , "green", "red", "nir08", "swir16", "swir22")
    } else if(collection=='sentinel-2-l2a'){
      assets <- c("B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B11", "B12")
    } else {
      stop("Only default assets for Landsat/Sentinel-2 implemented")
    }
  }
  
  ##  Reproject bounding box to crs of feature
  featCRS <- crs(paste0('EPSG:',feature$properties$`proj:epsg`))
  cropExt <- terra::project(ext(bbox, xy=TRUE), crs("+proj=longlat +datum=WGS84"), featCRS)
  
  ##  Read assets as spatRaster object - crop - harmonize resolutions (agg/disagg) if required 
  featBands <- assets2rast(feature, assets, TRUE)
  resolutions <- sapply(featBands, function(x) res(x)[1])
  resLevels <- as.numeric(levels(as.factor(resolutions)))
  
  if(length(resLevels)==1){
    #Simple case in which all assets have same resolution
    featBands <- rast(featBands) %>%
      crop(cropExt, extend=extend)
  } else {
    #split bands in list by resolution and crop
    featBands <- lapply(resLevels, function(r) rast(featBands[resolutions==r]) %>% crop(cropExt, snap="in", extend=extend))
    
    #crop all resolutions to extent of lowest spatial resolution to make sure they allign
    #TODO: Should search for a better way to do this by first defining the commom extent, then cropping
    commonExt <- ext(featBands[[which.max(resLevels)]])
    featBands <- lapply(featBands, crop, commonExt, extend=extend)
    
    if(commonRes=="low"){
      suppressWarnings( #will give a warning for resolution where aggregation is not necessary
        featBands <- lapply(featBands, function(x) aggregate(x, fact=max(resLevels)/res(x)[1], fun="mean")) %>%
          rast()
      )        
    } else if (commonRes=="high") {
      suppressWarnings(
        featBands <- lapply(featBands, function(x) disagg(x, fact=res(x)[1]/min(resLevels), method="near")) %>%
          rast()
      )
    } else {
      stop('Argument "commonRes" must be either "low" or "high"')
    }
  } 
  
  ##  Apply mask
  #add if() with argument to make this optional?
  if(collection=='landsat-c2-l2'){
    #Apply mask based on Landsat 'pa_pixel' bits
    numToBinary <- function(nums, nBits){
      sapply(nums, function(x){as.integer(intToBits(x)[1:nBits])})
    }
    metaBits <- assets2rast(feature, "qa_pixel") %>%
      crop(cropExt, extend=extend) %>%
      app(numToBinary, nBits=8)
    featBands <- terra::mask(featBands, app(terra::subset(metaBits, 7, negate=TRUE), sum)==0, maskvalues=FALSE, updatevalue=NA)
  } else if(collection=='sentinel-2-l2a'){
    #Apply mask based on "scene classification map" layer
    SCL <- assets2rast(feature, "SCL") %>%
      crop(commonExt, extend=extend)
    if(res(SCL)[1] > res(featBands)[1]) SCL <- disagg(SCL, fact=res(SCL)[1]/res(featBands)[1], method="near")
    SCL <- classify(SCL, matrix(data=c(
      0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
      0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0,  0), ncol=2))
    featBands <- terra::mask(featBands, SCL, maskvalues=0, updatevalue=NA)
  }
  
  ##  Attach time to output spatRast
  time(featBands) <- rep(as.POSIXlt(feature$properties$datetime, format="%Y-%m-%dT%H:%M:%S", tz="GMT"), nlyr(featBands))
  
  ##  Re-order layers to order of input assets
  featBands <- featBands[[assets]]
  
  return(featBands)
}

  

  
