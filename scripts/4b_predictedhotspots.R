require(sf)
require(tidyverse)
require(terra)
require(raster)

# Requires datasets of predicted hosts (from scripts 3a-f)
base_path <- "MAMMALS_10km_sum/" # path to mammals AOH maps (Lumbierres et al., 2022; doi.org/10.1038/s41597-022-01838-w)

# Create a blank global raster in the Equal Earth projection
equalearth <- CRS("+proj=eqearth")
blank_raster_ee <- raster(xmn = -16920565,
                          xmx = 16920565,
                          ymx = 8315130,
                          ymn = -8392928,
                          res = 10000,
                          crs = equalearth)

# Create a global bounding box in the Equal Earth projection
bb <- sf::st_union(sf::st_make_grid(
  st_bbox(c(xmin = -180,
            xmax = 180,
            ymax = 90,
            ymin = -90), crs = st_crs(4326)),
  n = 100))

bb_ee <- st_transform(bb, as.character(equalearth))

# Loop through viral groups
for (vir in c("ohv", "flv", "mmv", "ebv", "hpv", "phv", "bcov")){ 
  
  # Load the corresponding predictions file
  final_predictions <- read_csv(paste0(".../data/", vir, "_predictions_main.csv"))
  final_predictions$order <- str_to_lower(final_predictions$order)
  target_orders <- unique(final_predictions %>% pull(order))
  
  final_predictions$iucn2020_binomial <- str_replace(final_predictions$iucn2020_binomial, " ",  "_")
  
  # Filter predicted hosts
  reservoirs <- final_predictions %>% filter((real_status == "high-evidence host" | hard_class == 1) & predictable_fam == 1) %>% pull(iucn2020_binomial)
  
  raster_list <- NULL
  
  for (sp in reservoirs){
    
    print(paste("Current species:", sp, "(", which(reservoirs == sp), "/", length(reservoirs), ")"))
    
    final_predictions %>% 
      distinct(iucn2020_binomial, order, family) %>% 
      filter(iucn2020_binomial == sp) %>% 
      pull(order) -> o
    
    if (file.exists(paste0(base_path, "Mammals_", str_to_lower(o), "/", sp, ".tif"))){
      ras_r <- rast(paste0(base_path, "Mammals_", str_to_lower(o), "/", sp, ".tif")) # Get AOH raster
  
      # Reclassify the raster to binary 
      ras_r[ras_r[] < 2500] <- NA
      ras_r[ras_r[] >= 2500] <- 1
      
      # Reproject the species raster to the 10km Equal Earth grid
      ras_r_ee <- project(ras_r, as.character(equalearth), res = 10000, method = "near", origin = origin(blank_raster_ee))
      
      # Extend the raster's extent to perfectly match the others
      ras_r_ee_extended <- raster::extend(raster(ras_r_ee), as(bb_ee, "Spatial"), snap = "near")

      # Mask out any artifacts
      ras_r_ee_extended <- raster::mask(ras_r_ee_extended, as(bb_ee, "Spatial"))
      
      # Append the aligned rasters
      raster_list <- append(raster_list, ras_r_ee_extended)
    
    }
    
  }
  
  # Create multi-layer stack (and save it for Figure 3)
  stacked <- stack(raster_list)
  terra::writeRaster(stacked, paste0(".../data/predictedstack_", vir,"_2025.grd"), overwrite = TRUE)

}


viral_groups <- c("phv", "flv", "ohv", "mmv", "hpv", "ebv")

for (vir in viral_groups) {
  
  stack_path <- paste0("predictedstack_", vir, "_2025.grd")
  stacked <- brick(stack_path)
  
  # Sum the values across all layers for every pixel
  # This creates a map of known host richness
  predicted_hotspots <- calc(stacked, sum, na.rm = TRUE)
  
  predicted_hotspots[predicted_hotspots[] == 0] <- NA
  
  output_path <- paste0(".../data/", vir, "_pred_hotspots_2025.tif")
  writeRaster(predicted_hotspots, filename = output_path, format = "GTiff", overwrite = TRUE)
}

