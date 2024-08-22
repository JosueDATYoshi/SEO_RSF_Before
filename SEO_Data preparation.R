#_________________________________
#Short-eared Owl Data Cleaning and preparation ####


#loading packages

library(animove)
library(move)
library(move2)
library(dplyr)
library(sf)
library(vroom)
library(ggplot2)
library(rnaturalearth) #needs rnaturalearthhires also installed
library(ggmap)
library(ggspatial) ## needs prettymapr also installed
library(gganimate)
library(plotly)
library(mapview)
library(gganimate)
library(units)
library(lubridate)
library(moveVis)
library(openeo)
library(RStoolbox)
library(ctmm)
library(sf)
library(mvtnorm)
library(terra)
library(move2)
library(ctmm)
library(sf)

#working directory

setwd("C:/Users/jdat9/OneDrive/Documents/GitHub/SEO_RSF_Before")


### 1. Data MoveBank ####

#Login to MoveBank

#movebank_store_credentials("Josue_YOSHI_Arteaga", key_name = "Personal Account")
#keyring::key_list()
        #study_id=2568569325        "Rewilding Galapagos Short-eared Owls"

movebank_store_credentials("Josue_YOSHI_Arteaga", "AmW6*ptdfmubaH4")
studyID <- movebank_download_study_info(study_id=2568569325)$id

movebank_download_study_info(study_id = studyID) %>%
  print(width = Inf) #all cols

##### a. Extracting information ####

deploy_db <- movebank_download_deployment(study_id = studyID)
individuals <- levels(deploy_db$individual_local_identifier) #to see all the individuals
transmitter <- levels(deploy_db$model) #To have direct access information for the tag used
islands <- unique(deploy_db$study_site) #To have an easy reference to the island to select
                                        #[1]=NA which are not deployed tags
                                        #[2]=Santa Cruz
                                        #[3]=Floreana
movebank_download_study_info(study_id = studyID)$sensor_type_ids

owls <- movebank_download_study(study_id = studyID,
                                sensor_type_ids="GPS", 
                                )

#owls <- mt_read("SEO_MoveBank_Download_20240624.csv") #Includes all the records from ACC and GPS
#owls_gps <- mt_read("SEO_MoveBank_Download_20240624_GPS.csv") #Pre-filtered GPS

names(owls)

mt_is_track_id_cleaved(owls) # see if it grouped tracks
mt_is_time_ordered(owls) # see if it organized tracks

#Only necessary when using full database
#owls_gps <- dplyr::filter(owls, !sf::st_is_empty(owls)) #Removing ACC data lines using empty Geographic data
            
  
mt_has_unique_location_time_records(owls_gps) #check for duplicated points

#remove duplicated points
owls_gps <- owls_gps %>%
  mutate(n_na = rowSums(is.na(pick(everything())))) %>%
  dplyr::arrange(n_na) %>%
  mt_filter_unique(criterion='first') %>% # this always needs to be "first" because the duplicates get ordered according to the number of columns with NA. 
  dplyr::arrange(mt_track_id()) %>% ## data needs be reordered again
  dplyr::arrange(mt_track_id(),mt_time())

mt_has_unique_location_time_records(owls_gps) #check for duplicated points

mt_time_column(owls_gps)
mt_track_id_column(owls_gps)

sf::st_coordinates(owls_gps) %>% head()
sf::st_crs(owls_gps)

### Save the move2 object as RData file
# save(owls_gps, file="SEO_gps_clean.Rdata")
# load("SEO_gps_clean.Rdata") #cannot be assigned a new name, but can contain multiple objects

#sf::st_transform(owls_gps, crs = "EPSG:4326") %>% sf::st_crs()

#movebank_get_vocabulary(owls_gps) #Information of each of the columns for the object

## this script contains translators (mostly both directions) between the move2 object and the:
# - telemetry class of ctmm package: move2_TO_telemetry() (mt_as_move2() reads in telemetry & telemetry lists)
# - track_xyt from amt: move2_TO_track.xyt() (mt_as_move2() reads in track.xyt objects)
# - binClstPath and binClstStck from EMbC: move2_TO_embctrk() & embctrk_TO_move2()
# - ltraj from adehabitatLT: move2_TO_ltraj() & ltraj_TO_move2()
# - data.frame for recurse functions: move2_TO_recurseDF()
# - data.frame containing all event and track attributes: move2_TO_df() (mt_as_move2() reads in data.frame)

# Refers to Anne R Script for translation functions

# Function necessary to transform move2 object to telemetry---------------------------------------------
move2_TO_telemetry <- function(mv2) {
  # needed columns: individual.local.identifier (or tag.local.identifier), timestamp, location.long and location.lat
  mv2 <- mt_as_event_attribute(mv2, names(mt_track_data(mv2)))
  mv2 <- dplyr::mutate(mv2, location.long=sf::st_coordinates(mv2)[,1],
                       location.lat=sf::st_coordinates(mv2)[,2])
  
  mv2df <- data.frame(mv2)
  ## as.telemetry expects the track id to be called "individual.local.identifier" this is a quick fix, it might need some more thought to it to make it nicer. HOPE THIS IS FIXED ONCE ctmm INTEGRATES READING IN move2
  # fix: idtrack colum gets the prefix "track_id:", individual.local.identifier gets the sufix "_original" to maintain this original information
  colnames(mv2df)[colnames(mv2df)%in%make.names(mt_track_id_column(mv2))] <- paste0("track_id:",make.names(mt_track_id_column(mv2)))
  colnames(mv2df)[colnames(mv2df)%in%c("individual.local.identifier","individual_local_identifier","individual-local-identifier")] <- paste0(colnames(mv2df)[colnames(mv2df)%in%c("individual.local.identifier","individual_local_identifier","individual-local-identifier")],"_original")
  mv2df$individual_local_identifier <-mt_track_id(mv2)
  mv2df$timestamp <- mt_time(mv2) # ensuring used timestamps are in the column "timestamp" as expected by as.telemetry()
  telem <- as.telemetry(mv2df,
                        timezone=tz(mt_time(mv2)),
                        projection= if(st_is_longlat(mv2)){NULL}else{projection(mv2)},
                        na.rm= "col",
                        keep=T)
  return(telem)
}

# telemetry_TO_move2 => mt_as_move2 reads in a telemetry /list telemetry object

owls_tele<- move2_TO_telemetry(owls_gps)


####b. Creating a Square reference for Area of Interest AOI from Google maps ####

#### world map to have a reference layer for the background
worldMap <- rnaturalearth::ne_countries(returnclass = "sf", scale = "large")

#Manually creating a box for all the islands using bonuding box 
galapagos_bbox_coords <- matrix(c(
  -92.2, -1.7,
  -89.1, -1.7,
  -89.1, 0.7,
  -92.2, 0.7,
  -92.2, -1.7
), ncol = 2, byrow = TRUE)

bbox <- c(-92.2, -1.7, -89.1, 0.7)

# Create a polygon from the bounding box coordinates
galapagos_bbox_polygon <- st_polygon(list(galapagos_bbox_coords))

# Convert the polygon to an sf object with WGS 84 projection
galapagos_bbox <- st_sf(geometry = st_sfc(galapagos_bbox_polygon, crs = 4326))

  ######A.Loading rasters and Shapes for the islands ####

#Loading raster for elevation
# galapagos_elevation <- terra::rast('Rasters/galapagos_elevation.tif')
islands_elevation <- terra::rast('Rasters/galapagos_islands_elevation_10m.tif')
islands_elevation_utm <- terra::rast('Rasters/galapagos_islands_elevation_10m_utm.tif')

#Loading raster for the outline of the islands
islands_mask <- terra::rast('Rasters/galapagos_islands_mask.tif')
islands_mask_utm <- terra::rast('Rasters/galapagos_islands_mask_utm.tif')


options(viewer = NULL)

raster::plot(islands_elevation, ext=owls_gps)
#### Basic plot by islands MAYBE EVENT ID is classifying by island
raster::plot(owls_gps, 
     max.plot=1, pch=16, cex=0.5, axes=T,add=T,err=F)


raster::plot(islands_mask, col="Dark Green", ext=owls_gps)
#### Basic plot by individuals
raster::plot(owls_gps["individual-tag-local-identifier"], 
     max.plot=1, pch=16, cex=0.5, axes=T,add=T,err=F)
raster::plot(owls_gps, ### WHEN I DON'T SPECIFY, SEEMS LIKE IT IS CLASSIFYING BY ISLAND BUT NOT SURE WHY
             max.plot=1, pch=16, cex=0.5, axes=T,add=T,err=F)


hist(islands_elevation)

terra::plot(islands_elevation_utm, ext=length(PS16)*2)
#### Basic plot by islands 
plot(PS16, 
     max.plot=1, pch=16, cex=0.5, axes=T, add=T,err=F)


hist(islands_elevation_utm)


#Loading raster mask of island
islands_mask <- rast('Rasters/galapagos_islands_mask.tif')
islands_mask_utm <- rast('Rasters/galapagos_islands_mask_utm.tif')


terra::plot(islands_mask)
  
class(owls_gps)

### Dynamic Map: Mapview ####

# by using mapview, after changing the class so that it is recognised as an SF object
owlsSF <- owls_gps
class(owlsSF) <- class(owlsSF) %>% setdiff("move2") # remove class "move2" from object
options(viewer = NULL)
mapview::mapView(owlsSF, zcol="individual-tag-local-identifier", legend=F) #as points
# mapview::mapView(mt_track_lines(owls_gps), zcol="individual-tag-local-identifier", legend=F) #as lines
######????? mt_track_lines not working

