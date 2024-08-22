#_________________________________
#Short-eared Owl RSF by individual using loop ####

#load data ####

#libraries
library(animove)
library(ctmm)
library(dplyr)
library(future.apply) # for parallel processing
library(sf)
library(mvtnorm)
library(terra)
library(viridis)

#code to use more power from the computer
plan("multisession", workers = min(6, parallelly::availableCores()))

#working directory 
setwd("C:/Users/jdat9/OneDrive/Documents/GitHub/SEO_RSF_Before")

#loading owl data
load("SEO_gps_clean.Rdata") #Move2 object in GPS form

#loading raster environmental data
NDVI_flor <- raster::readAll(raster("./Rasters/mean_2023Jan_2024March_NDVI_Floreana.tif"))
NDVI_SC <- raster::readAll(raster("./Rasters/mean_2023JanDec_NDVI_SantaCruz.tif"))
NDVI <- terra::merge(NDVI_flor, NDVI_SC)
elev_galap <- raster::readAll(raster("./Rasters/galapagos_elevation_30m.tif"))
elev_AOI <- raster::readAll(raster("./Rasters/AOI_elevation_30m.tif"))
ndvi_scale <- ndvi_colors <- colorRampPalette(c("#a60027", "#fffebe", "#036e3a"))(100)

#_________________________________
#cleaning data and transformation ####

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
#list of problematic island-hopppers individuals PS0002, PS0006, PS0011, PS0021, PS0026
e <- c(2,6,12,28,44) #empty telemetry PS0100_OWL09
owls_tele <- owls_tele[-e] #to remove empty telemetry

#Reproject information
ctmm::projection(owls_tele) <- raster::crs(NDVI) #Projecting telemetry into NDVI CRS
owls_mv <- move2::mt_as_move2(owls_tele) #Correct telemetry info back into move2 for plots
crs(elev_galap) <- crs(NDVI) #Defining coordinate system for elevation
crs(elev_AOI) <- crs(NDVI) #Defining coordinate system for elevation

#_________________________________
#loop to run AKDE RSF and AKDE+RSF for each owl independently ####

#creating list of owls to be processed
owls_names <- names(owls_tele)

for (owl_name in owls_names) {

  #owl_name will be each individual string from owls_names
  print(paste("Processing owl:", owl_name))  #print each owl name to check progress/errors

# Defining individual to work with
id_owl_mv<- subset(owls_mv, track == owl_name) #Select move2 object for each individual
id_owl<-owls_tele[[owl_name]] #Select telemetry for each individual
name <- owl_name #creating a vector name for saving the records of each individual

#cropping raster to extent of owls observations
r <- as(extent(id_owl_mv)*2,'SpatialPolygons') #to crop raster to individuals extent
crs(r) <- crs(NDVI) #making sure crs of the crop length is the same as the whole raster
rast_owl <- crop(NDVI, r) #croping raster

#plot to check owls distribution 
plot(rast_owl) 
plot(st_geometry(id_owl_mv), add = TRUE, type = "p", 
     col = "red", pch = 19, cex = 0.5)

#' Create named list of rasters
#' rsf.fit can deal with raster layers that differ in resolution or even projection
owl_env_list <- list("Raster_Individual" = rast_owl) #we only have one working raster

#' Definition of reference grid the akde predictions should get aligned to
reference_grid <- rast_owl

#_________________________________
#fitting ctmm model ####
#isotropic ctmm model, because this is a requirement for rsf.fit
id_owl_guess <- ctmm.guess(id_owl, CTMM = ctmm(isotropic = TRUE), 
                          interactive = FALSE)

#creating a new path to save ctmm model
model_name <- paste0("./outputs/models/fit_ctmm_", name,".RDS")

#Conditional to run a fit model if not already done
if (file.exists(model_name)) {
  owl_fit <- readRDS(model_name)
} else {
  owl_fit <- ctmm.select(id_owl, id_owl_guess, cores = -1)
  saveRDS(owl_fit, model_name)
}

summary(owl_fit) # show results from ctmm model

#_________________________________
#akde model #### 
#for the initial range distribution without considering the environment
owl_akde <- akde(id_owl, owl_fit, grid=reference_grid) 

#creating a new path to save akde model
AKDE_name <- paste0("./AKDEs/AKDE_", name,".RData")

#saving akde model 
save(owl_akde, file = here::here("outputs", AKDE_name))

#showing the initial akde
plot(owl_akde) #just plot with 95CI
pal <- color(id_owl, by = "time") #independent colors for the owl points
plot(id_owl, col = pal, main = "id_owl",err=F,add=T) #plotting the points within the akde

plot(rast_owl) #plotting NDVI raster 
plot(st_geometry(id_owl_mv), add = TRUE, type = "p", 
     col = "red", pch = 19, cex = 0.5) #plotting individual locations
plot(ctmm::raster(owl_akde, DF = "PMF")) #plotting akde in raster mode

#_________________________________
#rsf model ####
owl_rsf_riemann <- rsf.fit(id_owl, owl_akde, 
                             R = owl_env_list, integrator = "Riemann") #not sure why
#path for saving rsf model
RSF_name <- paste0("./RSFs/RSF_", name,".RData")

#saving rsf model
save(owl_rsf_riemann, file = here::here("outputs", RSF_name))


#' A suitability map - with 95% confidence interval using RSF model
suitability_riemann <- ctmm::suitability(owl_rsf_riemann, R = owl_env_list, 
                                         grid = reference_grid)
terra::plot(suitability_riemann)

#' Range distribution (includes the ranging behaviour)
agde_owl <- agde(CTMM = owl_rsf_riemann, R = owl_env_list, 
                   grid = reference_grid)
plot(agde_owl)
# Plot raster of range distribution
agde_raster <- ctmm::raster(agde_owl, DF = "PMF")
plot(agde_raster)

#_________________________________
# Selection-informed akde (AKDE+RSF)####
akde_rsf <- akde(id_owl, CTMM = owl_rsf_riemann, R = owl_env_list, 
                 grid = reference_grid)
plot(ctmm::raster(akde_rsf, DF = "PMF"))

#path to save the akde+rsf
output_name <- paste0("./AKDE+RSF/AKDE+RSF_", name,".RData")

#Saving akde+rsf
save(akde_rsf, file = here::here("outputs", output_name))
#load(here::here("outputs", output_name))

#plotting to compare all three rasters resulting from the models
pmf <- raster::stack(ctmm::raster(owl_akde, DF = "PMF"),
                     ctmm::raster(agde_owl, DF = "PMF"),
                     ctmm::raster(akde_rsf, DF = "PMF"))
names(pmf) <- c("akde", "rsf.fit", "akde_rsf.fit")
plot(pmf, zlim=c(0, max(raster::getValues(pmf)))) # use zlim to force same color scale

}
