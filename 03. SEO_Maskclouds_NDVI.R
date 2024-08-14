#Masking clouds and average NDVI

library(terra)
library(RStoolbox)

#########
setwd("C:/Users/jdat9/OneDrive/Documents/GitHub/SEO_RSF_Before")
image_list <- list.files("./data_S2data_SantaCruz_openeo/", full.names = TRUE)
bandnames <- readRDS("bandnames.RDS")
dir.create("NDVI", showWarnings = FALSE)

lapply(image_list, function(i){
  # i = image_list[1]
  img <- rast(i)
  # add bandnames
  names(img) <- bandnames
  ### mask clouds
  # create cloud mask
  cld_msk <- terra::ifel(img[['SCL']] %in% c(7, 8, 9), NA, 1)
  terra::plot(cld_msk, main= i) #this is the last edit, make sure it works when doing the masking and plotting of everything
  
  # mask clouds in bands B01 to B12
  # define the bands you want to mask
  bands_to_mask <- bandnames[1:12]
  bands_to_mask
  
  # apply the mask 
  img_msked <- terra::mask(img[[bands_to_mask]], cld_msk)
  
  ### NDVI
  r_norm <- normImage(img_msked) # center and scaling to sd
  r_resc <- rescaleImage(r_norm, ymin=0, ymax=1) # reflectance range 0 to 1
  # names(r_resc)
  # calculate NDVI for Sentinel-2:
  # returns a layer stack if multiple indices are calculated
  r_indices <- spectralIndices(
    r_resc,
    red = 4,
    nir = 8,
    indices = c("NDVI")
  )
  # terra::plot(r_indices)
  img_name <- gsub(basename(i), pattern = ".tif", replacement = "_NDVI.tif")
  img_path <- paste0("./NDVI/", img_name)
  writeRaster(r_indices, img_path,overwrite=T)
})

#Average NDVI

image_list <- list.files("./NDVI/", full.names = TRUE)

vrt_file <- "./NDVI/NDVI.vrt"

NDVI_vrt <- terra::vrt(image_list, vrt_file, "-separate", overwrite=TRUE)

mean_NDVI <- mean(NDVI_vrt, na.rm = TRUE)

ndvi_scale <- ndvi_colors <- colorRampPalette(c("#a60027", "#fffebe", "#036e3a"))(100)

terra::plot(mean_NDVI, col=ndvi_scale)

writeRaster(mean_NDVI, "./Rasters/mean_2023Jan_2024March_NDVI_SantaCruz.tif")

mean_NDVI <- terra::rast("./Rasters/mean_2023Jan_2024March_NDVI_SantaCruz.tif")

terra::plot(mean_NDVI, col=ndvi_scale)
terra::plot(owls_tele, 
     max.plot=1, pch=16, cex=0.5, axes=T, err=F,add=T)
