#_________________________________
#Short-eared Owl aKDE by individual ####

#' This code needs ctmm version 1.1.1 (or newer)
#' remotes::install_github("ctmm-initiative/ctmm")

#loading packages

library(animove)
library(ctmm)
library(sf)
library(mvtnorm)
library(terra)
library(move2)
library(ctmm)
library(sf)

# Convert from raster to terra object. Terra is the recommended package.
#islands_elevation <- terra::rast('Raster Files/galapagos_islands_elevation_10m.tif')

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

ctmm::projection(owls_tele) <- crs(islands_elevation)

plot(islands_elevation)
plot(owls_gps["individual-tag-local-identifier"], 
         +      max.plot=1, pch=16, cex=0.5, axes=T,add=T)

class(owls_tele)
crs(owls_tele)
projection(owls_tele)
crs(islands_elevation)
summary(owls_tele)


# Owl with longest tracking info in Floreana PS0016_ASI08
# Owl with average tracking info in Floreana PS0008_OWL06
# Owl with shortest tracking info in Floreana NR87_OWL01
# Owl with longest tracking info in Santa Cruz PS0100_OWL12
# Owl with average tracking info in Santa Cruz PS0098_ASI11
# Owl with shortest tracking info in Santa Cruz PS0095_ASI12

PS16<-owls_tele$PS0016_ASI08
PS08<-owls_tele$PS0008_OWL06
NR87<-owls_tele$NR87_OWL01

PS100<-owls_tele$PS0100_OWL12
PS98<-owls_tele$PS0098_ASI11
PS95<-owls_tele$PS0095_ASI12

id_owl <- PS16 #Change this line to select individuals

head(id_owl)

# Plotting locations:
mapview::mapView(mt_track_lines(owls_gps$PS0016_ASI08), zcol="", legend=F) #as points

pal <- color(id_owl, by = "time")
plot(id_owl, col = pal, main = "id_owl",err=F,add=T)

# time sampling plot to check on gaps in the information:
dt.plot(id_owl)

# Range residency assumption: ---------------------------------------------

level <- 0.95 # we want to display 95% confidence intervals
xlim <- c(0,1 %#% "day") # to create a window of one day  #### Chris?

# Checking for the range residency assumption:
SVF <- variogram(id_owl)
par(mfrow = c(1,2))
plot(SVF, fraction = 0.5, level = level)
abline(v = 1, col = "red", lty = 2) # adding a line at 1 month
plot(SVF, xlim = xlim, level = level)
par(mfrow = c(1,1))

# Model selection ---------------------------------------------------------
# Selecting the best-fit movement model through model selection:

# Calculate an automated model guesstimate:
GUESS1 <- ctmm.guess(id_owl, interactive = FALSE)
summary(GUESS1)

# Automated model selection, starting from GUESS1:
start_time <- Sys.time()
FIT1_ML <- ctmm.select(id_owl, GUESS1, 
                       method = "ML", verbose = TRUE)
Sys.time() - start_time # Time difference of 3.745824 mins
summary(FIT1_ML)

plot(SVF, CTMM = FIT1_ML[[1]],
     units = TRUE, fraction = 0.5, level = c(0.95, 0.50), 
     col = "black", col.CTMM = "red")

start_time <- Sys.time()
FIT1_pHREML <- ctmm.select(id_owl, GUESS1,
                           method = "pHREML", verbose = TRUE)
## reminder: it will default to pHREML if no method is specified.
Sys.time() - start_time # Time difference of 2.435679 mins
summary(FIT1_pHREML)

plot(SVF, CTMM = FIT1_pHREML[[1]],
     units = TRUE, fraction = 0.5, level = c(0.95, 0.50), 
     col = "black", col.CTMM = "red")

summary(FIT1_ML[[1]]) # best-fit model only
summary(FIT1_pHREML[[1]])

# fit_phREML_PS16 <- readRDS( "models/fit_pHREML_PS16.RData")

saveRDS(FIT1_pHREML, "models/fit_pHREML_PS16.RDS")




# Home range estimator ----------------------------------------------------
# Feeding a movement model into the home range estimator

# Run an area-corrected AKDE (default):
AKDE1_ML <- akde(id_owl, FIT1_ML[[1]], debias = TRUE)
AKDE1_pHREML <- akde(id_owl, FIT1_pHREML[[1]], debias = TRUE)

summary(AKDE1_ML, level.UD = 0.95)$CI # 95% home range area
summary(AKDE1_pHREML, level.UD = 0.95)$CI

( 1 - summary(AKDE1_ML)$CI[1,2] / summary(AKDE1_pHREML)$CI[1,2] ) * 100
# ML overestimates by ~6%

m.iid <- ctmm.fit(id_owl) # IID
KDE1 <- akde(id_owl, m.iid)

# Creating an extent that includes both UDs at the 95% CI level:
newEXT <- extent(list(AKDE1_pHREML, KDE1))

# Plotting KDE and AKDE side-by-side:
par(mfrow = c(1,2))
plot(id_owl, UD = KDE1, ext = newEXT)
title(expression("KDEc"))
plot(id_owl, UD = AKDE1_pHREML, ext = newEXT)
title(expression("AKDEc"))
par(mfrow = c(1,1))

# Mitigation measures -----------------------------------------------------
# Evaluating additional biases, applying mitigation measures

## Irregular representation in time: ONLY USE WHEN BROKEN SAMPLING TIMES AVAILABLE --------------------------------------

plot(id_owl, lwd = 3)

## Sample sizes:
summary(AKDE1_pHREML)$DOF["area"] # effective sample size of animal1
nrow(id_owl) # absolute sample size

# plot all sampling intervals
dt.plot(id_owl) # id_owl (buffalo[[4]])
abline(h = 1 %#% "hours", col = "red")

pal <- "hr" %#% diff(id_owl$t)
# minimum adjacent sampling interval
pal <- pmin(c(Inf, pal), c(pal, Inf))
# sampling intervals under 1.5 hours
pal <- (pal < 1.5)
# red (low-frequency) or yellow (high-frequency)
pal <- grDevices::rgb(1, pal, 0)
plot(id_owl, col = pal, lwd = 2)

# minimum sampling interval
"minutes" %#% min(diff(id_owl$t))

# Calculate wAKDE:
start_time <- Sys.time()
wAKDE1_pHREML <- akde(id_owl,
                      CTMM = FIT1_pHREML[[1]],
                      weights = TRUE)
# you only need this with irregular sampling (can be slow!)
# unweighted AKDE places too much density on oversampled times
Sys.time() - start_time
summary(wAKDE1_pHREML)$CI # 95% home range area (weighted)
## Time difference of 1.66443 mins

plot(id_owl, UD = wAKDE1_pHREML)

EXT <- extent(list(AKDE1_ML, AKDE1_pHREML, wAKDE1_pHREML), level = 0.95)

# Plotting pHREML (with and without weights) side-by-side:
par(mfrow = c(1,2))
plot(id_owl, UD = AKDE1_pHREML, ext = EXT)
title(expression("pHREML AKDE"["C"]))
plot(id_owl, UD = wAKDE1_pHREML, ext = EXT)
title(expression("pHREML wAKDE"["C"]))
par(mfrow = c(1,1))

( 1 - summary(AKDE1_pHREML)$CI[1,2] / summary(wAKDE1_pHREML)$CI[1,2] ) * 100
# Unweighted AKDE underestimates by 3%



## Low effective sample sizes: ONLY USE WHEN INDIVIDUALS HAVE LOW EFFECTIVE SAMPLING--------------------------------------------

# For a target bias of O(5%)
# ctmm.fit method="ML" requires DOF[area]>=20       (CONVENTIONAL)
# ctmm.fit method="pHREML" requires DOF[area]>=4-5  (DEFAULT)
# ctmm.boot method="pHREML" requires DOF[area]>=2-3 (SLOW)
# but in all cases DOF[area] is an estimate

data("gazelle")
gazelle <- gazelle[[11]]
head(gazelle)
plot(gazelle, col = "blue", lwd = 3)

GUESS2 <- ctmm.guess(gazelle, interactive = FALSE)
FIT2_pHREML <- ctmm.select(gazelle, GUESS2, method = "pHREML")
summary(FIT2_pHREML)
# note the effective sample sizes

# SI units converter
help("%#%")
1 %#% 'hr'

# approximate effective sample size
summary(gazelle)
(8.100378 %#% "month") / (5.881915 %#% "month")
summary(FIT2_pHREML)$DOF[["area"]]

UD2_pHREML <- akde(gazelle, FIT2_pHREML)
summary(UD2_pHREML)

summary(UD2_pHREML)$DOF["area"] # effective sample size
nrow(gazelle) # absolute sample size

# Expected order of pHREML bias:
1/summary(FIT2_pHREML)$DOF["area"]^2

help("ctmm.boot")
start_time <- Sys.time() # start recording running time
BOOT <- ctmm.boot(gazelle, FIT2_pHREML, 
                  error = 0.01, trace = 2, cores = -1)
# save(BOOT, file = here::here("data", "outputs", "bootstrap.RData"))
## note: this function incurs substantial computational cost, may take hours.
( total_time <- Sys.time() - start_time ) # output running time
# Time difference of 43.93944 mins

load(here::here("outputs", "bootstrap.RData"))
summary(BOOT)

# Expected order of bootstrap bias:
1/summary(BOOT)$DOF["area"]^3

UD2_bpHREML <- akde(gazelle, BOOT, weights = TRUE)

summary(UD2_pHREML)$CI
summary(UD2_bpHREML)$CI

( 1 - summary(UD2_pHREML)$CI[1,2] / summary(UD2_bpHREML)$CI[1,2] ) * 100

EXT <- extent(list(UD2_pHREML, UD2_bpHREML), level = 0.95)

# Plotting pHREML and bootstrapped-pHREML side-by-side:
par(mfrow = c(1,2))
plot(gazelle, UD = UD2_pHREML, ext = EXT)
title(expression("pHREML AKDE"["C"]))
plot(gazelle, UD = UD2_bpHREML, ext = EXT)
title(expression("Bootstrapped pHREML wAKDE"["C"]))
par(mfrow = c(1,1))

plot(gazelle, UD = UD2_pHREML, ext = EXT)
plot(gazelle, UD = UD2_bpHREML, ext = EXT)

# save(FIT1_ML,
#      FIT1_pHREML,
#      wAKDE1_pHREML,
#      FIT2_ML,
#      FIT2_pHREML,
#      BOOT,
#      file = here::here("outputs", "akde.RData"))






# RSF Analysis ------------------------------------------------------------

##### Toca crear un raster compilado asi que todavia no podemos usar un stack

env <- islands_elevation

#' ## Animals and elevation
plot(env)
owl_mv <- move2::mt_as_move2(id_owl) # Convert to move2 object for plotting
plot(owl_mv, add = TRUE)

#' # Minimal example of rsf.fit
#' Select one individual, Cilla
cilla <- buffalo$Cilla

#' Fit ctmm model
#' We start with a guesstimate of parameter values; we demand a
#' isotropic ctmm model, because this is a requirement for rsf.fit
cilla_guess <- ctmm.guess(cilla, CTMM = ctmm(isotropic = TRUE), 
                          interactive = FALSE)

# Fit ctmm model and perform model selection, using all but one core
if (file.exists("cilla_fit_rsf.rds")) {
  cilla_fit <- readRDS("cilla_fit_rsf.rds")
} else {
  cilla_fit <- ctmm.select(cilla, cilla_guess, cores = -1)
  saveRDS(cilla_fit, "cilla_fit_rsf.rds")
}

#' Create named list of rasters
#' rsf.fit can deal with raster layers that differ in resolution or even projection
env_list <- list("elev" = raster::raster(env[["elev"]]),
                         "slope" = raster::raster(env[["slope"]]),
                         "var_NDVI" = raster::raster(env[["var_NDVI"]]))




