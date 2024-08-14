# load/install necessary libraries
library(openeo)
library(sf)
library(terra)
library(RStoolbox)
library(ggplot2)
library(dplyr)

setwd("C:/Users/jdat9/OneDrive/Documents/GitHub/SEO_RSF_Before")

# temporal extent
t <- c("2023-12-09", "2024-03-21") # 14 months to extract info during time we have information for the owls

# connect to the back-end
connection <- openeo::connect(host = "https://openeo.dataspace.copernicus.eu")
login()

# assign the graph-building helper object to p for easy access to all openEO processes
# see > ?processes()
p <- processes()

# Define the bounding box for the Galapagos Islands
galapagos_bbox_coords <- matrix(c(
  -92.2, -1.7,
  -89.1, -1.7,
  -89.1, 0.7,
  -92.2, 0.7,
  -92.2, -1.7
), ncol = 2, byrow = TRUE)

# Get OSM data for the Galapagos Floreana
# load S2 footprint

foot_wgs84 <- st_read('Limit SantaCruz.gpkg')
bbox_foot <- st_bbox(foot_wgs84)


ext <- list(west = bbox_foot[[1]], south = bbox_foot[[2]], 
            east = bbox_foot[[3]], north = bbox_foot[[4]])

# load first datacube for defined spatial and temporal extent
# Sentinel-2 Level 2A
cube_s2 <- p$load_collection(
  id = "SENTINEL2_L2A",
  spatial_extent = ext,
  temporal_extent = t #, bands = c('Bandname1', 'Bandname2')
)

s2 <- describe_collection("SENTINEL2_L2A") # or use the collection entry from the list, e.g. collections$`COPERNICUS/S2`
print(s2)

# extract bandnames from s2
bandnames <- unlist(s2$`cube:dimensions`$bands$values)
bandnames

#saveRDS(bandnames,"bandnames.RDS")

# list file formats for saving the result
formats <- list_file_formats()

# save using save_result, give format via list
res <- p$save_result(data = cube_s2, format = formats$output$GTiff)

# send job to back-end
job <- create_job(graph = res, title = "S2_GLP")

##################################
# start job 
##################################
start_job(job = job)
# takes some time

##################################
# list job status
###################################
jobs_status <- list_jobs()

jobs_ls <- lapply(1:length(jobs_status), function(i){
  df <- do.call(cbind.data.frame, jobs_status[[i]]) %>%  
    dplyr::select(created, id, status, updated)
  return(df)
})
jobs_df <- do.call(rbind, jobs_ls) %>% 
  arrange(created)

jobs_df

##################################
# list the processed resultest
##################################
list_results(job = job)

#########################################################
# download all the files into a folder on your computer
#########################################################

# dir.create("../data/S2data_Radolfzell_openeo", showWarnings = FALSE)
download_results(job = job, folder = "data_S2data_SantaCruz2_openeo")

#Best ones so far: 03-17, 03-02, 02-15, 02-10

S2_openeo <- rast('data_S2data_SantaCruz_openeo/openEO_2023-07-15Z.tif')
S2_openeo
bandnames

# add bandnames
# !!! careful, this works only if you actually downloaded all bands in the collection
names(S2_openeo)
names(S2_openeo) <- bandnames
names(S2_openeo)

# plot RGB image
ggRGB(S2_openeo,r = 'B04', g = 'B03', b = 'B02', geom_raster = TRUE,
      stretch = 'lin')

#Code to mask clouds individually but there is code for all of them at once. 

# ####################################
# # mask clouds using SCL layer
# ####################################
# 
# # plot scene classification layer
# plot(S2_openeo[['SCL']])
# 
# # scl band (produced by Sen2Cor): scene classification layer can give you information on clouds, 
# # cloud shadow, snow etc.
# # see website for further information: 
# # https://www.sentinel-hub.com/faq/how-get-s2a-scene-classification-sentinel-2/
# # 7,8,9 = low, medium and high probability clouds 
# 
# # create cloud mask
# cld_msk <- terra::ifel(S2_openeo[['SCL']] %in% c(7, 8, 9), NA, 1)
# plot(cld_msk)
# 
# # mask clouds in bands B01 to B12
# # define the bands you want to mask
# bands_to_mask <- bandnames[1:12]
# bands_to_mask
# 
# # apply the mask 
# s2_msked <- terra::mask(S2_openeo[[bands_to_mask]], cld_msk)
# 
# # check if worked
# plot(S2_openeo[['B04']])
# plot(s2_msked[['B04']])
# 


