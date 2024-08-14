#### Telemetry object to SF

# Load the required packages
library(ctmm)
library(sf)

# Load or create a telemetry object (example)
telemetry_object <- owls_tele

# Extract coordinates and time
coords <- data.frame(
  x = telemetry_object$longitude,
  y = telemetry_object$latitude,
  time = telemetry_object$timestamp
)

# Convert to sf object
owls_sf <- st_as_sf(coords, coords = c("x", "y"), crs = crs(mean_NDVI_flor), agr = "constant") #DEFINE Projection in CRS

# Print the sf object
print(id_owl_sf)

# Optionally, you can save it to a file
st_write(sf_object, "telemetry_data.shp")