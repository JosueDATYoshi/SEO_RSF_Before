#' # Resource selection functions
#' Modified from Code by Björn Reineking, 2024-06-20 ANIMOVE 2024
#'
#' Model the relative density of animals (also called range distribution or utilisation distribution) as a function of environmental predictors.
#'

#' 
#' Loading packages
library(animove)
library(ctmm)
library(dplyr)
library(future.apply) # for parallel processing
library(sf)
library(viridis)

plan("multisession", workers = min(6, parallelly::availableCores()))

#' ## Load buffalo data
# data("buffalo") #example
# Loading Owl Data
load("SEO_gps_clean.Rdata") #Move2 object in GPS form

#' ## Load owls in telemetry form

#Important function to transform data from Move2 Object into Telemetry
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

owls_tele<- move2_TO_telemetry(owls_gps) #transforming information 
owls_tele <- owls_tele[-44] #to remove empty telemetry

#' ## Environmental data: topography and NDVI

# data(buffalo_env) #example
NDVI_flor <- terra::rast("./Rasters/mean_2023Jan_2024March_NDVI_Floreana.tif")
NDVI_SC <- terra::rast("./Rasters/mean_2023JanDec_NDVI_SantaCruz.tif")
NDVI <- terra::merge(NDVI_flor, NDVI_SC)
elev_galap <- raster::readAll(raster("./Rasters/galapagos_elevation_30m.tif"))
elev_AOI <- raster::readAll(raster("./Rasters/AOI_elevation_30m.tif"))

ndvi_scale <- ndvi_colors <- colorRampPalette(c("#a60027", "#fffebe", "#036e3a"))(100) #Scale for plotting

#' Project buffalo data to projection of rasters 
# needed later when doing akde with a provided grid

ctmm::projection(owls_tele) <- raster::crs(NDVI) # working for NDVI raster but not for elevation


##### Plotting code to validate information #####

terra::plot(elev_AOI)
terra::plot(owls_gps, 
            max.plot=1, pch=16, cex=0.75, axes=T,err=F, add=T)

terra::plot(NDVI, col=ndvi_scale) # alternative color blind 
# terra::plot(NDVI, col = viridis(256))
terra::plot(owls_tele, 
            max.plot=1, pch=16, cex=0.75, axes=T,err=F, add=T)

# terra::plot(NDVI_SC, col=ndvi_scale)
# terra::plot(owls_tele, 
#             max.plot=1, pch=16, cex=0.75, axes=T,err=F, add=T)
# 
# terra::plot(NDVI_flor, col=ndvi_scale)
# terra::plot(owls_tele, 
#             max.plot=1, pch=16, cex=0.75, axes=T,err=F, add=T)


#' Create named list of rasters
# buffalo_env_list <- list("elev" = raster::raster(buffalo_env, "elev"),
#            "slope" = raster::raster(buffalo_env, "slope"),
#            "var_NDVI" = raster::raster(buffalo_env, "var_NDVI"))


#' Create a list of ctmm models
owls_guess <- lapply(owls_tele, ctmm.guess, 
                         CTMM = ctmm(isotropic = TRUE), 
                         interactive = FALSE)
owls_ctmm <- future_mapply(ctmm.select, owls_tele, owls_guess, 
                               SIMPLIFY = FALSE, future.seed = TRUE)

#' Calling akde on the list of buffalos and their ctmms ensures 
#' that a common grid is used; here we explicitly provide a common grid
owls_akde <- akde(owls_tele, owls_ctmm) # , grid = NDVI)
### owls_akde <- akde(owls_tele, owls_ctmm, grid = elev_AOI) #Maybe a secondary option for all the area of study

#' See ctmm_3_meta.R from Inês and Chris for more information

#' Calculate rsf.fit for each individual
owls_rsf <- future_mapply(rsf.fit, owls_tele, owls_akde, 
                       MoreArgs = list(R = NDVI, #would need to be change for the analysis with elevation
                                       integrator = "Riemann"), 
                       SIMPLIFY = FALSE, future.seed = TRUE)

#' Have a look at the individual rsf selection estimates
selection_coef <- sapply(owls_rsf, "[[", 3)
colnames(selection_coef) <- names(owls_tele)
selection_coef
#' Except for Queen, buffalos have negative estimates for elevation 
#' (not necessarily statistically significant)

#' Calculation of population-level rsf
mean_rsf <- mean(owls_rsf)
summary(mean_rsf)

#' Map the suitability for the "average" animal
suitability_mean <- ctmm::suitability(mean_rsf, 
                                      R = elev_Galap, #HERE I Should have all the rasters in stack
                                      grid = elev_Galap)
raster::plot(suitability_mean)



