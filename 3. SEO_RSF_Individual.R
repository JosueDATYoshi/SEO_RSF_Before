#_________________________________
#Short-eared Owl RSF by individual ####

#Load libraries
library(animove)
library(ctmm)
library(sf)
library(mvtnorm)
library(terra)

#Loading Raster data ####
NDVI_flor <- terra::rast("./Rasters/mean_2023Jan_2024March_NDVI_Floreana.tif")
NDVI_SC <- terra::rast("./Rasters/mean_2023JanDec_NDVI_SantaCruz.tif")
NDVI <- terra::merge(NDVI_flor, NDVI_SC)
elev_galap <- raster::readAll(raster("./Rasters/galapagos_elevation_30m.tif"))
elev_AOI <- raster::readAll(raster("./Rasters/AOI_elevation_30m.tif"))

#reproject telemetry data 

ctmm::projection(owls_tele) <- raster::crs(NDVI) # working for NDVI raster but not for elevation

ndvi_scale <- ndvi_colors <- colorRampPalette(c("#a60027", "#fffebe", "#036e3a"))(100)


# Defining individual to work with
PS02<-owls_tele$PS0002_OWL11 #Select manually which owl to work with

id_owl <- PS02 #identifying which owl you want to work with
name <- "PS0002"



#' Fit ctmm model
#' We start with a guesstimate of parameter values; we demand a
#' isotropic ctmm model, because this is a requirement for rsf.fit
id_owl_guess <- ctmm.guess(id_owl, CTMM = ctmm(isotropic = TRUE), 
                          interactive = FALSE)

# Creating a new path for the models since it takes time to run the models
model_name <- paste0("./models/fit_pHREML_", name,".RDS")


#Conditional to run a fit model if not already done
if (file.exists(model_name)) {
  owl_fit <- readRDS(model_name)
} else {
  owl_fit <- ctmm.select(id_owl, id_owl_guess, cores = -1)
  saveRDS(owl_fit, model_name)
}

summary(owl_fit)

#' Create named list of rasters
#' rsf.fit can deal with raster layers that differ in resolution or even projection

owl_env_list <- list("NDVI" = NDVI)


#' Definition of reference grid the akde predictions should get aligned to
reference_grid <- crop(NDVI)

owl_akde <- akde(id_owl, owl_fit, grid=NDVI)
 
# Error in seq.default(EXT[1, i] - dr[i], EXT[2, i] + dr[i], length.out = 1 +  : 
#                        'length.out' must be a non-negative number

plot(owl_akde)
pal <- color(id_owl, by = "time")
plot(id_owl, col = pal, main = "id_owl",err=F,add=T)

##### Not working #####

plot(ctmm::raster(owl_akde, DF = "PMF")) 

owl_rsf_riemann <- rsf.fit(id_owl, owl_akde, 
                             R = elev_raster_utm, integrator = "Riemann")
class()
