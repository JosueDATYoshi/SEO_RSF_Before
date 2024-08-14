#### This code was created only to get the rasters necessary for the analysis ####

### Still missing to get is the NDVI, human scape (land cover) and 

# Install and load necessary packages
#install.packages(c("elevatr", "sf", "raster", "ggplot2", "osmdata", "terra"))
library(elevatr)
library(sf)
library(raster)
library(ggplot2)
library(osmdata)
library(terra)

setwd("C:/Users/jdat9/OneDrive/Documents/GitHub/SEO_RSF_Before/Rasters/")

# Define the bounding box for the Galapagos Islands
galapagos_bbox_coords <- matrix(c(
  -92.2, -1.7,
  -89.1, -1.7,
  -89.1, 0.7,
  -92.2, 0.7,
  -92.2, -1.7
), ncol = 2, byrow = TRUE)

# Get OSM data for the Galapagos Islands
bbox <- c(-92.2, -1.7, -89.1, 0.7)

# Fetch buildings data
buildings <- opq(bbox = bbox) %>%
  add_osm_feature(key = 'building') %>%
  osmdata_sf() %>%
  .$osm_polygons

# # Fetch roads data
# roads <- opq(bbox = bbox) %>%
#   add_osm_feature(key = 'highway') %>%
#   osmdata_sf() %>%
#   .$osm_lines
# 
# # Fetch land use data
# landuse <- opq(bbox = bbox) %>%
#   add_osm_feature(key = 'landuse') %>%
#   osmdata_sf() %>%
#   .$osm_polygons

# # Ensure the combined data has the same CRS as the elevation data
# buildings <- st_transform(buildings, crs = st_crs(galapagos_elevation))
# roads <- st_transform(roads, crs = st_crs(galapagos_elevation))
# landuse <- st_transform(landuse, crs = st_crs(galapagos_elevation))
# 
# # Convert the sf objects to SpatVector objects
# buildings_vect <- terra::vect(buildings)
# roads_vect <- terra::vect(roads)
# landuse_vect <- terra::vect(landuse)

# # Rasterize buildings, roads, and land use data
# buildings_raster <- terra::rasterize(buildings_vect, galapagos_elevation, field = 1, background = NA)
# roads_raster <- terra::rasterize(roads_vect, galapagos_elevation, field = 1, background = NA)
# landuse_raster <- terra::rasterize(landuse_vect, galapagos_elevation, field = 1, background = NA)
# 
# # Combine all human activity rasters
# human_activity_raster <- merge(buildings_raster, roads_raster, landuse_raster, fun = "sum", na.rm = TRUE)

# Fetch coastline data
coastlines <- opq(bbox = bbox) %>%
  add_osm_feature(key = 'natural', value = 'coastline') %>%
  osmdata_sf() %>%
  .$osm_lines

# Create an island outline by converting the coastline lines to a polygon
coastlines <- st_transform(coastlines, crs = st_crs(buildings))

# Create a polygon from the coastlines
islands_outline <- st_union(st_cast(coastlines, "MULTILINESTRING")) %>%
  st_cast("POLYGON")

# Plot to verify the islands outline
plot(st_geometry(islands_outline))

# Create a polygon from the bounding box coordinates
galapagos_bbox_polygon <- st_polygon(list(galapagos_bbox_coords))

# Convert the polygon to an sf object with WGS 84 projection
galapagos_bbox_sf <- st_sf(geometry = st_sfc(galapagos_bbox_polygon, crs = 4326))

# Fetch elevation data
galapagos_elevation <- get_elev_raster(
  locations = galapagos_bbox_sf,
  z = 10, # Adjust zoom level for desired resolution
  prj = "EPSG:4326" # Set projection to center of individuals
)

# Ensure the elevation data is a SpatRaster (terra object)
galapagos_elevation <- terra::rast(galapagos_elevation)

# Plot the elevation data to verify
plot(galapagos_elevation, main = "Galapagos Islands Elevation (10m resolution)")

# Define the target projection string
target_proj <- "+proj=tpeqd +lat_1=-1.30589644452226 +lon_1=-90.4477761048985 +lat_2=-1.26375754076642 +lon_2=-90.4270820923035 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# Transform the elevation data to the specified projection
galapagos_elevation_utm <- terra::project(galapagos_elevation, projection(owls_tele))

# Plot the transformed elevation data to verify
raster::plot(galapagos_elevation_utm, main = "Galapagos Islands Elevation in UTM (10m resolution)")

# Ensure the islands outline has the same CRS as the raster
islands_outline <- st_transform(islands_outline, crs = st_crs(galapagos_elevation))
islands_outline_utm <- st_transform(islands_outline, crs = st_crs(galapagos_elevation_utm))



# Convert the islands outline to a SpatVector object (needed for masking)
islands_vect <- terra::vect(islands_outline)
islands_vect_utm <- terra::vect(islands_outline_utm)

# Mask the raster using the islands outline to keep only the islands
islands_elevation <- terra::mask(galapagos_elevation, islands_vect)
islands_elevation_utm <- terra::mask(galapagos_elevation_utm, islands_vect_utm)

#errase any negative values
islands_elevation[islands_elevation < 0] <- NA

islands_elevation_utm[islands_elevation_utm < 0] <- NA

# Create a binary mask (1 for islands, NA for ocean)
islands_mask <- islands_elevation
islands_mask[!is.na(islands_mask)] <- 1

islands_mask_utm <- islands_elevation_utm
islands_mask_utm[!is.na(islands_mask_utm)] <- 1

plot(islands_elevation, main = "Galapagos Islands Elevation (10m resolution)")
plot(islands_elevation_utm, main = "Galapagos Islands Elevation in UTM (10m resolution)")

# mask for Galapagos Floreana

foot_wgs84 <- st_read('Limit Floreana.gpkg')
bbox_foot <- st_bbox(foot_wgs84)


ext <- list(west = bbox_foot[[1]], south = bbox_foot[[2]], 
            east = bbox_foot[[3]], north = bbox_foot[[4]])

floreana_elev <- terra::crop(islands_elevation, foot_wgs84)

floreana_elev <- terra::extend(floreana_elevation, ext)



terra::plot(floreana_elev)

summary(floreana_elev)

# Save the all files

# Save and reload the raster data to ensure it's in the correct format
terra::writeRaster(galapagos_elevation, "galapagos_elevation.tif", overwrite = TRUE)
terra::writeRaster(galapagos_elevation_utm, "galapagos_elevation_utm.tif", overwrite = TRUE)

terra::writeRaster(islands_elevation_utm, "galapagos_islands_elevation_10m_utm.tif", overwrite = TRUE)
terra::writeRaster(islands_elevation, "galapagos_islands_elevation_10m.tif", overwrite = TRUE)

terra::writeRaster(islands_mask_utm, "galapagos_islands_mask_utm.tif", overwrite = TRUE)
terra::writeRaster(islands_mask, "galapagos_islands_mask.tif", overwrite = TRUE)

st_write(islands_outline,"Islands_outline.shp")
st_write(islands_outline_utm,"Islands_outline_utm.shp")


# Save the human activity raster
# terra::writeRaster(human_activity_raster, "galapagos_human_activity.tif", overwrite = TRUE)

# # Plot the human activity raster to verify
# plot(human_activity_raster, main = "Galapagos Islands Human Activity")

### 2.Extract environmental data rasters and annotations ####

# temporal extent
t <- c("2023-01-20", "2023-09-20") # 14 months to extract info during time

# connect to the back-end
connection <- openeo::connect(host = "https://openeo.dataspace.copernicus.eu")
login()

# assign the graph-building helper object to p for easy access to all openEO processes
# see > ?processes()
p <- processes()

bbox_foot <- st_bbox(galapagos_bbox)

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
# list file formats for saving the result
formats <- list_file_formats()

# save using save_result, give format via list
res <- p$save_result(data = cube_s2, format = formats$output$GTiff)

# send job to back-end
job <- create_job(graph = res, title = "S2_Floreana")

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
download_results(job = job, folder = "data_S2data_Galapagos_openeo")

# dir.create("../data/S2data_SouthAmerica_openeo", showWarnings = FALSE)
# download_results(job = job, folder = "../data/S2data_SouthAmerica_openeo")

####################################
# load and and inspect result
####################################

S2_openeo <- rast('../data/S2data_Radolfzell_openeo/openEO_2023-06-13Z.tif')
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

