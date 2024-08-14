#### Telemetry object to SF

# Load the required packages
library(ctmm)
library(sf)

# Assume owls_tele is your list of telemetry objects
# Create a directory to save the shapefiles if it doesn't exist
if(!dir.exists("./owls_shape")) {
  dir.create("./owls_shape")
}

# Initialize a list to store the individual SF objects
sf_list <- lapply(seq_along(owls_tele), function(index) {
  i <- owls_tele[[index]]
  
  # Extract coordinates and time
  coords <- data.frame(
    x = i$longitude,
    y = i$latitude,
    time = as.POSIXct(i$timestamp) # Convert timestamp to POSIXct for correct DateTime format
  )
  
  # Convert to sf object
  owls_sf <- st_as_sf(coords, coords = c("x", "y"), crs = crs(mean_NDVI_flor), agr = "constant") # Define Projection in CRS
  
  # Generate a unique filename for each owl using the index
  Opath <- paste0("./owls_shape/owl_", index, ".shp")
  
  # Check if file exists and delete if necessary to avoid append/overwrite issues
  if (file.exists(Opath)) {
    file.remove(Opath)
  }
  
  # Save each SF object separately
  st_write(owls_sf, Opath, append = FALSE)
  
  # Return the SF object
  return(owls_sf)
})

# Combine all individual SF objects into a single SF object
combined_sf <- do.call(rbind, sf_list)

# Save the combined SF object
combined_path <- "./owls_shape/combined_owls.shp"
if (file.exists(combined_path)) {
  file.remove(combined_path)
}
st_write(combined_sf, combined_path, append = FALSE)

# Print the combined SF object
print(combined_sf)
