### Code to transform shapefile downloaded from Google Earth (kml) into gpkg file so we can use as a Area of interest (AOI)

# Load the sf package
library(sf)

# Specify the file path for the KML file and the output GeoPackage
kml_file <- "Raster limit SCruzArea.kml"
gpkg_file <- "Limit SantaCruz.gpkg"

# Read the KML file
sf_data <- st_read(kml_file)

# Write the data to a GeoPackage
st_write(sf_data, gpkg_file, driver = "GPKG")

# Print a message indicating success
cat("KML file has been successfully transformed into a GeoPackage.\n")
