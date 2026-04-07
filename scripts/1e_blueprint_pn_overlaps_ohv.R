setwd("...")

require(tidyverse)
require(terra)

raster::rasterOptions(tmpdir = "tmp_over1")

base_path <- "MAMMALS_10km_sum/"

# There are no known high-evidence hosts (please refer to Supplementary Dataset 1):
ohv_reservoirs <- NULL

# Checking if all of them are in COMBINE
combine <- read.csv(".../data/combine_imputed.csv") %>% filter(iucn2020_binomial != "Not recognised")
virion_m_clean <- read.csv(".../data/virion_m_clean.csv")

# Aquatic species
combine %>% filter(freshwater == "1" | marine == "1") %>% pull(phylacine_binomial) -> aquatic

# calling 'not in' function
'%!in%' <- function(x, y)!('%in%'(x, y)) # %!in%

combine <- combine %>% filter(phylacine_binomial %!in% aquatic)

# These are the species that tested positive to the target viral genus but are not high-evidence hosts (according to WOS search)
# chances are that they have not been sufficiently studied, therefore they might be unrecognised reservoirs 
# 1) they are never pseudoabsences, 2) they enter into the analysis as 1s, but with lower instance weigths

pseudo_reservoir_pool <- c(unique(virion_m_clean %>% filter(VirusGenus %in% c("orthonairovirus") & Host %!in% ohv_reservoirs) %>%
                                    pull(Host)))

positive_species <- c(pseudo_reservoir_pool, ohv_reservoirs)

# Computing PseudoNegatives:
# 1) selecting target families (families in which there is at least one positive species)
target_families <- unique(combine %>% 
                            filter(iucn2020_binomial %in% positive_species) %>% 
                            pull(family))

setdiff(positive_species, combine %>% 
          filter(iucn2020_binomial %in% positive_species) %>% pull(iucn2020_binomial))

require(raster)
require(sf)
require(stringr)

equalearth <- CRS("+proj=eqearth")

blank_wgs <- raster(xmn = -180,
                    xmx = 180,
                    ymx = 90,
                    ymn = -90,
                    res = 0.089,
                    crs = "+proj=longlat +datum=WGS84 +no_defs")

blank_wgs <- rast(blank_wgs)

blank_raster_ee <- raster(xmn = -16920565,
                          xmx = 16920565,
                          ymx = 8315130,
                          ymn = -8392928,
                          res = 10000,
                          crs = "+proj=eqearth")

e <- terra::ext(blank_raster_ee)

bb <- sf::st_union(sf::st_make_grid(
  st_bbox(c(xmin = -180,
            xmax = 180,
            ymax = 90,
            ymin = -90), crs = st_crs(4326)),
  n = 100))

bb_ee <- st_transform(bb, as.character(equalearth))

# bb <- vect(bb)
# bb_ee <- terra::project(bb, equalearth)

combine %>% 
  filter(iucn2020_binomial %in% positive_species) %>% 
  distinct(order, family) -> ord_fam_tab

# tab_overlaps <- NULL
tab_overlaps <- read.csv(".../data/tab_overlaps_ohv.csv")

target_orders <- unique(combine %>% 
                          filter(iucn2020_binomial %in% positive_species) %>% 
                          pull(order))

# Computing pseudoabsences
# (if reservoir status is unknown put 0 when the species' range overlaps (50% threshold) with that of one or more positive species belonging to the same family)
combine$iucn2020_binomial <- str_replace(combine$iucn2020_binomial, " ", "_")
combine$order <- str_to_lower(combine$order)

# loading biogeographical realms (Olson et al. 2001) to mask areas of endemicity where to take PseudoNegatives
ecoregions <- read_sf(".../data/Holt2013/Regions.shp")

ecoregions_mod <- ecoregions
ecoregions_mod <- st_set_crs(ecoregions_mod, value = "PROJCRS[\"World_Plate_Carree12\",\n    BASEGEOGCRS[\"WGS 84\",\n        DATUM[\"World Geodetic System 1984\",\n            ELLIPSOID[\"WGS 84\",6378137,298.257223563,\n                LENGTHUNIT[\"metre\",1]],\n            ID[\"EPSG\",6326]],\n        PRIMEM[\"Greenwich\",12,\n            ANGLEUNIT[\"Degree\",0.0174532925199433]]],\n    CONVERSION[\"unnamed\",\n        METHOD[\"Equidistant Cylindrical\",\n            ID[\"EPSG\",1028]],\n        PARAMETER[\"Latitude of 1st standard parallel\",0,\n            ANGLEUNIT[\"Degree\",0.0174532925199433],\n            ID[\"EPSG\",8823]],\n        PARAMETER[\"Longitude of natural origin\",0,\n            ANGLEUNIT[\"Degree\",0.0174532925199433],\n            ID[\"EPSG\",8802]],\n        PARAMETER[\"False easting\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8806]],\n        PARAMETER[\"False northing\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8807]]],\n    CS[Cartesian,2],\n        AXIS[\"(E)\",east,\n            ORDER[1],\n            LENGTHUNIT[\"metre\",1,\n                ID[\"EPSG\",9001]]],\n        AXIS[\"(N)\",north,\n            ORDER[2],\n            LENGTHUNIT[\"metre\",1,\n                ID[\"EPSG\",9001]]]]")

# convert the sf object to Equal Earth  
ecoregions_transformed_reproj <- st_transform(ecoregions_mod, crs = 4326)
endem_map <- ecoregions_transformed_reproj[ecoregions_transformed_reproj$Regions %in% c("African", "Guineo-Congolian", "Indo-Malayan", "Eurasian", "Saharo-Arabian"), ]
endem_map_vect <- terra::vect(endem_map)
endem_map_vect_ee <- terra::project(endem_map_vect, blank_raster_ee)

# Loop over orders with at least one known host
for (ord in target_orders){
  
  print(paste("Current order:", ord, "(", which(target_orders == ord), "/", length(target_orders), ")"))
  
  target_fam <- ord_fam_tab %>% filter(order == ord) %>% pull(family)
  
  # Loop over families with known hosts
  for (fam in target_fam){
    
    print(paste("Current family:", fam, "(", which(target_fam == fam), "/", length(target_fam), ")"))
    
    res <- unique(combine %>% filter(family == fam  & phylacine_binomial %in% positive_species) %>% pull(iucn2020_binomial)) # getting names of reservoirs in fam
    non_res <- unique(combine %>% filter(family == fam & phylacine_binomial %!in% positive_species) %>% pull(iucn2020_binomial)) # getting names of non-reservoir/susceptible
    
    tab_comb <- crossing(non_host = str_replace(non_res, " ", "_"), host = str_replace(res, " ", "_"))
    tab_comb$combination <- seq(1:length(tab_comb$non_host))
    
    for (n in str_replace(non_res, " ", "_")){ # for each species with unknown status
      
      # Check if the Area of Habitat (AOH) raster file exists for this species
      if (file.exists(paste0(base_path, "Mammals_", str_to_lower(ord), "/", n, ".tif"))){
        
        # Load the non-host raster
        ras_nr <- rast(paste0(base_path, "Mammals_", str_to_lower(ord), "/", n, ".tif")[1]) 
        
        # Reclassify the raster: values < 2500 become NA (background), values >= 2500 become 1 (presence)
        # 2500 is likely a specific suitability threshold from the source data
        ras_nr[ras_nr[] < 2500] <- NA
        ras_nr[ras_nr[] >= 2500] <- 1
        
        for(r in str_replace(res, " ", "_")){ # For each known host
          
          # Check if the AOH raster file exists for the host species
          if (file.exists(paste0(base_path, "Mammals_", str_to_lower(ord), "/", r, ".tif"))){
            
            # Load the host raster
            ras_r <- rast(paste0(base_path, "Mammals_", str_to_lower(ord), "/", r, ".tif")) 
            
            # Reclassify host raster using the same 25% threshold
            ras_r[ras_r[] < 2500] <- NA
            ras_r[ras_r[] >= 2500] <- 1
            
            # Project rasters to Equal Earth projection at 10km (10000m) resolution
            # Nearest neighbor for categorical data
            ras_nr_ee <- project(ras_nr, as.character(equalearth), res = 10000, method = "near", origin = origin(blank_raster_ee))
            ras_r_ee <- project(ras_r, as.character(equalearth), res = 10000, method = "near", origin = origin(blank_raster_ee))
            
            # Mask the bounding box to remove outside artifacts
            ras_nr_ee_extended <- raster::extend(raster(ras_nr_ee), as(bb_ee, "Spatial"), snap = "near")
            ras_r_ee_extended <- raster::extend(raster(ras_r_ee), as(bb_ee, "Spatial"), snap = "near")
            ras_nr_ee_extended <- raster::mask(ras_nr_ee_extended, as(bb_ee, "Spatial"))
            ras_r_ee_extended <- raster::mask(ras_r_ee_extended, as(bb_ee, "Spatial"))
            
            print(paste("Overlapping unknown", n, "with positive", r, "(", tab_comb %>% filter(non_host == n, host == r) %>% pull(combination), "/", length(str_replace(res, " ", "_"))*length(str_replace(non_res, " ", "_")), ")"))
            
            # Calculate the intersection (overlap) between the non-host and host rasters
            overlap <- terra::mask(rast(ras_nr_ee_extended), rast(ras_r_ee_extended)) 
            
            # Restrict the overlap the "endemic ecoregions"
            overlap_endem <- terra::mask(overlap, endem_map_vect_ee)
            
            # Count the number of overlapping cells
            ncell <- length(overlap_endem[!is.na(overlap_endem)])
            
            
            rast_nonres <- mask(rast(ras_nr_ee_extended), endem_map_vect_ee)
            rast_res <- mask(rast(ras_r_ee_extended), endem_map_vect_ee)
            
            # Count the total number of cells for each species
            ncell_nonres <- length(rast_nonres[!is.na(rast_nonres)])
            ncell_res <- length(rast_res[!is.na(rast_res)])
            
            # A species is a PseudoNegative if 50% of its range overlaps with that of a host
            if (ncell >= 0.5*ncell_nonres & ncell_nonres != 0){
              
              print(paste("It's a match:", "positive", r, "and unknown", n))
              
              # Append the results
              tab_overlaps <- rbind(tab_overlaps, tibble(sp = n, positive_sp = r, overlap_on_target = ncell/ncell_nonres, overlap_on_positive = ncell/ncell_res))
              
              # Write the updated results of each iteration
              write.csv(tab_overlaps, ".../data/tab_overlaps_ebv.csv", row.names = F)
              
            }
          }
        }
      }
      
      # housekeeping
      rm(ras_r)
      rm(ras_nr)
      rm(ras_r_ee)
      rm(ras_nr_ee)
      rm(ras_r_ee_extended)
      rm(ras_nr_ee_extended)
      rm(overlap)
      gc()
      
    }
  }
}  