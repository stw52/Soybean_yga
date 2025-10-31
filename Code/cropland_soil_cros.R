#### Looking into cropland CROS for identifying soybean areas and then most dominant soil type for soybeans

library(sf)
library(raster)
library(terra)
library(CropScapeR)
library(soilDB)
library(tidyverse)
library(janitor)


years <- 2024
out_dir <- "cdl_ny_soy" 
dir.create(out_dir, showWarnings = FALSE)

for (y in years) {
  message("Requesting year ", y)
  # Get full-state CDL GeoTIFF (this calls CropScape web service)
  # GetCDLData returns a URL or downloads depending on function version; adjust if needed
  cdl_file <- GetCDLData(year = y, aoi = 36, type = "f",save_path = file.path(out_dir, paste0("CDL_NY_", y, ".tif")))
  
  # Read CDL
  r <- rast(cdl_file)
  
  # Soybeans are usually code 5 — check the CDL legend to confirm
  soy_mask <- r == 5
  
  # Apply mask (soybean pixels = 1, others = NA)
  soy_only <- mask(r, soy_mask, maskvalues = 0, updatevalue = NA)
  
  # Save soybean-only raster
  out_file <- file.path(out_dir, paste0("CDL_NY_soybeans_", y, ".tif"))
  writeRaster(soy_only, out_file, overwrite = TRUE)
}

### Get data for Wisconsin
out_dir_wi <- "cdl_wi_soy"
dir.create(out_dir_wi, showWarnings = F)

cdl_file_wi <- GetCDLData(year = 2024, aoi = 55, type = "f",save_path = file.path(out_dir_wi, paste0("CDL_WI_", 2024, ".tif")))


r_wi <- rast(cdl_file_wi)

# Soybeans are usually code 5 — check the CDL legend to confirm
soy_mask_wi <- r_wi == 5

# Apply mask (soybean pixels = 1, others = NA)
soy_only_wi <- mask(r_wi, soy_mask_wi, maskvalues = 0, updatevalue = NA)

# Save soybean-only raster
out_file_wi <- file.path(out_dir_wi, paste0("CDL_wi_soybeans_", 2024, ".tif"))
writeRaster(soy_only_wi, out_file_wi, overwrite = TRUE)


# Bring in QGIS processed joins

ny_joined <- clean_names(read_csv("QGIS_Project_APSIM_Soils/ny_soy_ssurgo_filtered.csv"))

ny_smaller <- ny_joined %>%
  dplyr::select(soy_pixelssum, comppct_r, compname) %>%
  group_by(compname) %>%
  summarise(soy_pixel_count = sum(comppct_r/100*soy_pixelssum)) %>%
  mutate(hectares = soy_pixel_count*900/10000) %>%
  slice_max(n = 20, order_by = hectares) 

wi_joined <- clean_names(read_csv("QGIS_Project_APSIM_Soils/wi_joined_ssurgo_right.csv"))


wi_smaller <- wi_joined %>%
  dplyr::select(soy_pixelssum, comppct_r, compname, drainagecl) %>%
  group_by(compname) %>%
  summarise(soy_pixel_count = sum(comppct_r/100*soy_pixelssum)) %>%
  mutate(hectares = soy_pixel_count*900/10000) %>%
  slice_max(n = 20, order_by = hectares) 

  


