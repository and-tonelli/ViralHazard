# Script to obtain cumulative hazard hotspots and reproduce Figure 3

require(Matrix)
require(raster)
require(terra)
require(tidyverse)
require(sf)
require(rnaturalearth)

ecoregions <- read_sf(".../data/Holt2013/Regions.shp") # Zoogeographic regions from Holt et al., 2013 (10.1126/science.12282)
equalearth <- CRS("+proj=eqearth")

# countries shapes (background)
countries <- ne_countries(scale = 50, returnclass = c("sf"))
countries_ee <- st_transform(countries, equalearth)

countries_ee_grouped <- countries_ee %>%
  group_by(`region_wb`) %>%
  summarise(n = n())

to_sparse <- function(r) {
  v <- terra::values(r, mat = TRUE)
  v[is.na(v)] <- 0
  Matrix(v, sparse = TRUE)
}

countries_ee_rast_list <- lapply(countries_ee_rast_list, rast)
countries_ee_mat_list <- lapply(countries_ee_rast_list, to_sparse)

# Compute global and regional hotspots (top 10%) for the different viral groups
for (vir in c("phv", "ohv", "flv", "ebv", "mmv", "bcov", "hpv")){ 
  
  hotspots_predicted <- raster(paste0(".../data/", vir, "_pred_hotspots_2025.tif"))
  
  hotspots_predicted[is.na(hotspots_predicted[])] <- 0
  
  # Use the bounding box to trim (mask) the raster layer (otherwise it will have weird duplicated areas in the corners)
  hotspots_predicted <- raster::mask(hotspots_predicted, as(bb_ee, "Spatial"))
  
  v <- terra::values(rast(hotspots_predicted), mat = TRUE)
  v[is.na(v)] <- 0
  Matrix(v, sparse = TRUE)
  
  # Hotspot - global
  hot <- quantile(v[v[, 1] > 0,] , 0.90)[[1]]
  
  v_hot <- ifelse(v[, 1] >= hot, TRUE, FALSE)
  
  for (i in c(2:8)) {
  
  # Mask the raster for the specific region
  maskd_hot <- (countries_ee_mat_list[[i]] != 0) & (v != 0)

  # Hotspot - regional
  hot_reg <- quantile(v[as.matrix(maskd_hot)] , 0.90)[[1]]
  v_reg <- ifelse(v[, 1] >= hot_reg & maskd_hot, TRUE, FALSE)
  
  if (i == 2){
    
    v_reg_final <- ifelse(v_reg | v_hot, TRUE, FALSE)
    
  }
  
  if (i > 2) {
    
    v_reg_final <- ifelse(v_reg_final | v_reg, TRUE, FALSE)
    
    }
  }
  
  assign(paste0("mat_hs_", vir, "_regions"), v_reg_final)
  
}

viral_groups <- c("phv", "flv", "ohv", "mmv", "hpv", "ebv")

for (vir in viral_groups) {
  
  mat_name <- paste0("mat_hs_", vir, "_regions")
  mat <- get(mat_name)
  mat <- ifelse(mat, 1, 0)

  r <- rast(nrows = 1679, ncols = 3448,
            xmin = -17240565, xmax = 17239435,
            ymin = -8394870, ymax = 8395130,
            crs = "+proj=eqearth")
  
  values(r) <- as.vector(t(mat)) # transpose to match raster row-major order

  append(paste0("r_g_", vir), r)
  
}

# Getting cumulative raster
raster_cum_global <- sum(stack(c(r_g_hpv, r_g_phv, r_g_bcov, r_g_flv, r_g_ohv, r_g_mmv, r_g_ebv)))

df_predicted_reg <- raster_cum_global %>%
  rasterToPoints %>%
  as.data.frame() %>%
  `colnames<-`(c("x", "y", "cumulative_hotspots")) %>% 
  filter(cumulative_hotspots != 0)

df_predicted_reg$cumulative_hotspots <- as.factor(df_predicted_reg$cumulative_hotspots)

# Cumulative hotspots map
ggplot() +
  # Background
  geom_sf(data = bb_ee,   
          colour = "#9B9EA0",
          linetype = 'solid',
          fill = "grey90", # lighter grey for better contrast
          size = 0.3,
          colour = NA) +
  geom_sf(data = countries_ee, # CONTINENTS
          colour = NA,
          linetype = 'solid',
          fill = "grey22", # dark grey #3B444D
  ) +
  geom_raster(
    data = df_predicted_reg,   # HOTSPOTS
    aes(
      x = x,
      y = y,
      fill = cumulative_hotspots))+
  scale_fill_manual(values = c("#3A0CA3", "#7209B7", "#B5179E", "#F72585", "#FF6D5A", "#FF9E51", "#FFC94C"))+ 
  theme_void() +
  theme(legend.position = "right",
        legend.box = "horizontal",
        legend.key.width = unit(0.5, 'cm'),
        legend.key.height = unit(1.5, 'cm'),
        legend.title = element_text(hjust = 0.5, size = 14, face = "bold"))+
  guides(fill = guide_legend(reverse = T))+
  labs(fill = "Number of\nhotspots") -> cum_hotspots_map

# Now obtain list of species within the different hotspots
## Loading the hotspots
hotspot_eap <- mask(raster_cum_global, countries_ee_grouped[countries_ee_grouped$region_wb == "East Asia & Pacific", ]) 
hotspot_eca <- mask(raster_cum_global, countries_ee_grouped[countries_ee_grouped$region_wb == "Europe & Central Asia", ]) 
hotspot_lac <- mask(raster_cum_global, countries_ee_grouped[countries_ee_grouped$region_wb == "Latin America & Caribbean", ]) 
hotspot_na <- mask(raster_cum_global, countries_ee_grouped[countries_ee_grouped$region_wb == "North America", ]) 
hotspot_sa <- mask(raster_cum_global, countries_ee_grouped[countries_ee_grouped$region_wb == "South Asia", ]) 
hotspot_ssa <- mask(raster_cum_global, countries_ee_grouped[countries_ee_grouped$region_wb == "Sub-Saharan Africa", ]) 
hotspot_mena <- mask(raster_cum_global, countries_ee_grouped[countries_ee_grouped$region_wb == "Middle East & North Africa", ]) 

eap_poly <- raster(".../data/eap_raster.tif")
eca_poly <- raster(".../data/eca_raster.tif")
lac_poly <- raster(".../data/lac_raster.tif")
mena_poly <- raster(".../data/mena_raster.tif")
na_poly <- raster(".../data/na_raster.tif")
sa_poly <- raster(".../data/sa_raster.tif")
ssa_poly <- raster(".../data/ssa_raster.tif")

lista_hot <- list(eap_poly, eca_poly, lac_poly, mena_poly, na_poly, sa_poly, ssa_poly)
names(lista_hot) <- c("eap", "eca", "lac", "mena", "na", "sa", "ssa")

# loading raster stack that store info of species
bcv_rast <- raster::stack(".../data/predictedstack_bcov_2025.grd")
phv_rast <- raster::stack(".../data/predictedstack_phv_2025.grd")
flv_rast <- raster::stack(".../data/predictedstack_flv_2025.grd")
ohv_rast <- raster::stack(".../data/predictedstack_ohv_2025.grd")
mmv_rast <- raster::stack(".../data/predictedstack_mmv_2025.grd")
hpv_rast <- raster::stack(".../data/predictedstack_hpv_2025.grd")
ebv_rast <- raster::stack(".../data/predictedstack_ebv_2025.grd")

to_sparse <- function(r) {
  v <- terra::values(r, mat = TRUE)
  v[is.na(v)] <- 0
  Matrix(v, sparse = TRUE)
}

results_list <- list()

for (vir in c("phv", "ohv", "flv", "ebv", "mmv", "bcov", "hpv")){
  
  vir_stack <- raster::stack(paste0(".../data/predictedstack_", vir, "_2025.grd"))
  
  for (n in 1:length(lista_hot)){
    
    valid_cells <- which(!is.na(getValues(lista_hot[[n]])))
    
    # Get coordinates of those cells
    coords <- xyFromCell(lista_hot[[n]], valid_cells)
    
    # Extract values from rast at those coordinates
    vals <- extract(vir_stack, coords)
    
    # Check which layers have at least one non-NA value
    has_overlap <- apply(vals, 2, function(x) any(!is.na(x)))
    
    # Get names of overlapping layers
    overlap_names <- names(vir_stack)[has_overlap]
    
    if (length(overlap_names) > 0) {
      for (sp in overlap_names) {
        results_list[[length(results_list) + 1]] <- data.frame(
          virus = vir,
          mask_id = names(lista_hot)[n],
          Species = sp,
          stringsAsFactors = FALSE
        )
      }
    }
    
  }
}

df_vir <- do.call(rbind, results_list)
df_vir <- df_vir %>% distinct(virus, mask_id, Species)

write.csv(df_vir, ".../data/species_within_hotspotregions.csv", row.names = F)

eap_hotspots_o <- df_vir_fin %>% filter(mask_id == "eap") %>% mutate(ViralGenus = fct_relevel(ViralGenus, c("Ebolavirus & Marburgvirus",  "Henipavirus", "Betacoronavirus", "Phlebovirus", "Orthonairovirus", "Flavivirus", "Mammarenavirus"))) %>% mutate(continent = "East Asia & Pacific")
eap_hotspots_o$order <- factor(eap_hotspots_o$order, levels=sort(unique(eap_hotspots_o$order)))

eca_hotspots_o <- df_vir_fin %>% filter(mask_id == "eca") %>% mutate(ViralGenus = fct_relevel(ViralGenus, c("Ebolavirus & Marburgvirus",  "Henipavirus", "Betacoronavirus", "Phlebovirus", "Orthonairovirus", "Flavivirus", "Mammarenavirus"))) %>% mutate(continent = "Europe & Central Asia")
eca_hotspots_o$order <- factor(eca_hotspots_o$order, levels=sort(unique(eca_hotspots_o$order)))

mena_hotspots_o <- df_vir_fin %>% filter(mask_id == "mena") %>% mutate(ViralGenus = fct_relevel(ViralGenus, c("Ebolavirus & Marburgvirus",  "Henipavirus", "Betacoronavirus", "Phlebovirus", "Orthonairovirus", "Flavivirus", "Mammarenavirus"))) %>% mutate(continent = "Middle East & North Africa")
mena_hotspots_o$order <- factor(mena_hotspots_o$order, levels=sort(unique(mena_hotspots_o$order)))

lac_hotspots_o <- df_vir_fin %>% filter(mask_id == "lac") %>% mutate(ViralGenus = fct_relevel(ViralGenus, c("Ebolavirus & Marburgvirus",  "Henipavirus", "Betacoronavirus", "Phlebovirus", "Orthonairovirus", "Flavivirus", "Mammarenavirus"))) %>% mutate(continent = "Latin America & Caribbean")
lac_hotspots_o$order <- factor(lac_hotspots_o$order, levels=sort(unique(lac_hotspots_o$order)))

na_hotspots_o <- df_vir_fin %>% filter(mask_id == "na") %>% mutate(ViralGenus = fct_relevel(ViralGenus, c("Ebolavirus & Marburgvirus",  "Henipavirus", "Betacoronavirus", "Phlebovirus", "Orthonairovirus", "Flavivirus", "Mammarenavirus"))) %>% mutate(continent = "North America")
na_hotspots_o$order <- factor(na_hotspots_o$order, levels=sort(unique(na_hotspots_o$order)))

sa_hotspots_o <- df_vir_fin %>% filter(mask_id == "sa") %>% mutate(ViralGenus = fct_relevel(ViralGenus, c("Ebolavirus & Marburgvirus",  "Henipavirus", "Betacoronavirus", "Phlebovirus", "Orthonairovirus", "Flavivirus", "Mammarenavirus"))) %>% mutate(continent = "South Asia")
sa_hotspots_o$order <- factor(sa_hotspots_o$order, levels=sort(unique(sa_hotspots_o$order)))

ssa_hotspots_o <- df_vir_fin %>% filter(mask_id == "ssa") %>% mutate(ViralGenus = fct_relevel(ViralGenus, c("Ebolavirus & Marburgvirus",  "Henipavirus", "Betacoronavirus", "Phlebovirus", "Orthonairovirus", "Flavivirus", "Mammarenavirus"))) %>% mutate(continent = "Sub-Saharan Africa")
ssa_hotspots_o$order <- factor(ssa_hotspots_o$order, levels=sort(unique(ssa_hotspots_o$order)))

palette_vir <- c("#71014B", # burgundy
                 "#FF5C77", # pink
                 "#FFB300", # orange
                 "#5DC122", # green
                 "#046C7A", # aqua
                 "#002A86", # blue
                 "#470CED") # violet

names(palette_vir) <- unique(levels(ssa_hotspots_o$ViralGenus))

# individual subplots:
eap_hotspots_o %>% 
  # group_by(continent, order) %>% 
  # mutate(n = n()
  #        # order = reorder_within(order, n, continent)
  # ) %>% 
  ggplot(aes(y = as.character(order), fill = `ViralGenus`))+
  geom_bar()+
  scale_fill_manual(values = palette_vir)+
  facet_wrap(~ continent, nrow = 1, ncol = 7, scales = "free_y")+
  # scale_y_reordered()+
  theme_BP()+
  scale_y_discrete(limits = rev(levels(eap_hotspots_o$order)))+
  theme(panel.spacing = unit(1.5, "lines"),
        legend.position = "none",
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 9, margin = margin(t = 0, r = 0, b = 0, l = 0)))+
  labs(y = "", x = "Number of unique associations", fill = "") -> east_asia_plot

eca_hotspots_o %>% 
  # group_by(continent, order) %>% 
  # mutate(n = n()
  #        # order = reorder_within(order, n, continent)
  # ) %>% 
  ggplot(aes(y = as.character(order), fill = `ViralGenus`))+
  geom_bar()+
  scale_fill_manual(values = palette_vir)+
  facet_wrap(~ continent, nrow = 1, ncol = 7, scales = "free_y")+
  # scale_y_reordered()+
  theme_BP()+
  scale_y_discrete(limits = rev(levels(eca_hotspots_o$order)))+
  theme(panel.spacing = unit(1.5, "lines"),
        legend.position = "none",
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 9, margin = margin(t = 0, r = 0, b = 0, l = 0)))+
  labs(y = "", x = "Number of unique associations", fill = "") -> europe_ca_plot

mena_hotspots_o %>% 
  group_by(continent, order) %>% 
  mutate(n = n()
         # order = reorder_within(order, n, continent)
  ) %>% 
  ggplot(aes(y = reorder(order, n), fill = `ViralGenus`))+
  geom_bar()+
  scale_fill_manual(values = palette_vir)+
  facet_wrap(~ continent, nrow = 1, ncol = 7, scales = "free_y")+
  # scale_y_reordered()+
  theme_BP()+
  scale_y_discrete(limits = rev(levels(mena_hotspots_o$order)))+
  theme(panel.spacing = unit(1.5, "lines"),
        legend.position = "none",
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 9, margin = margin(t = 0, r = 0, b = 0, l = 0)))+
  labs(y = "", x = "Number of unique associations", fill = "") -> middle_easte_na_plot

lac_hotspots_o %>% group_by(continent, order) %>% 
  mutate(n = n()
         # order = reorder_within(order, n, continent)
  ) %>% 
  ggplot(aes(y = reorder(order, n), fill = `ViralGenus`))+
  geom_bar()+
  scale_fill_manual(values = palette_vir)+
  facet_wrap(~ continent, nrow = 1, ncol = 7, scales = "free_y")+
  # scale_y_reordered()+
  theme_BP()+
  scale_y_discrete(limits = rev(levels(lac_hotspots_o$order)))+
  theme(panel.spacing = unit(1.5, "lines"),
        legend.position = "none",
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 9, margin = margin(t = 0, r = 0, b = 0, l = 0)))+
  labs(y = "", x = "Number of unique associations", fill = "") -> latam_carib_plot

ssa_hotspots_o %>% group_by(continent, order) %>% 
  mutate(n = n()
         # order = reorder_within(order, n, continent)
  ) %>% 
  ggplot(aes(y = reorder(order, n), fill = `ViralGenus`))+
  geom_bar()+
  scale_fill_manual(values = palette_vir)+
  facet_wrap(~ continent, nrow = 1, ncol = 7, scales = "free_y")+
  # scale_y_reordered()+
  theme_BP()+
  scale_y_discrete(limits = rev(levels(ssa_hotspots_o$order)))+
  theme(panel.spacing = unit(1.5, "lines"),
        legend.position = "none",
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 9, margin = margin(t = 0, r = 0, b = 0, l = 0)))+
  labs(y = "", x = "Number of unique associations", fill = "") -> africa_plot

sa_hotspots_o %>% group_by(continent, order) %>% 
  mutate(n = n()
         # order = reorder_within(order, n, continent)
  ) %>% 
  ggplot(aes(y = reorder(order, n), fill = `ViralGenus`))+
  geom_bar()+
  scale_fill_manual(values = palette_vir)+
  facet_wrap(~ continent, nrow = 1, ncol = 7, scales = "free_y")+
  # scale_y_reordered()+
  theme_BP()+
  scale_y_discrete(limits = rev(levels(sa_hotspots_o$order)))+
  theme(panel.spacing = unit(1.5, "lines"),
        legend.position = "none",
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 9, margin = margin(t = 0, r = 0, b = 0, l = 0)))+
  labs(y = "", x = "Number of unique associations", fill = "") -> southasia_plot

na_hotspots_o %>% group_by(continent, order) %>% 
  mutate(n = n()
         # order = reorder_within(order, n, continent)
  ) %>% 
  ggplot(aes(y = reorder(order, n), fill = `ViralGenus`))+
  geom_bar()+
  scale_fill_manual(values = palette_vir)+
  facet_wrap(~ continent, nrow = 1, ncol = 7, scales = "free_y")+
  # scale_y_reordered()+
  theme_BP()+
  scale_y_discrete(limits = rev(levels(na_hotspots_o$order)))+
  theme(panel.spacing = unit(1.5, "lines"),
        legend.position = "null",
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 9, margin = margin(t = 0, r = 0, b = 0, l = 0)))+
  labs(y = "", x = "Number of unique associations", fill = "") -> northa_plot

require(cowplot)
virus_legend <- get_legend(northa_plot)

# Assemble final figure
# northam - europe - east asia
# latam -  MAP  - southasia
# africa - middle east - legends

require(tidyverse)
require(cowplot)
require(grid)

cum_hotspots_map+
  theme(legend.position = "right",
        legend.box = "horizontal",
        legend.key.width = unit(0.3, 'cm'),
        legend.key.height = unit(0.6, 'cm'),
        legend.title = element_text(hjust = 0.5, size = 9, face = "plain")) -> map_withleg


map_legend <- get_legend(map_withleg)

ggdraw() +
  draw_plot(cum_hotspots_map+
              theme(legend.position = "none",
                    plot.margin = unit(c(15, 15, 15, 15), "pt")
                    ),
                    x = 0.22, y = 0.22, width = 0.55, height = 0.55) +
  draw_plot(northa_plot, 0.08, 0.7, 0.25, 0.3) +
  draw_grob(rectGrob(gp = gpar(col = "gray55", fill = NA, lwd = 1)), 0.08, 0.7, 0.25, 0.298) +
  geom_segment(aes(x = 0.355, y = 0.65, xend = 0.08+0.25/2, yend = 0.7), color = "gray20") +
  draw_plot(europe_ca_plot, 0.37, 0.70, 0.25, 0.3) +
  draw_grob(rectGrob(gp = gpar(col = "gray55", fill = NA, lwd = 1)), 0.37, 0.70, 0.25, 0.298) +
  geom_segment(aes(x = 0.59, y = 0.65, xend = 0.37+0.25/2, yend = 0.7), color = "gray20") +
  draw_plot(southasia_plot, 0.66, 0.70, 0.25, 0.3) +
  draw_grob(rectGrob(gp = gpar(col = "gray55", fill = NA, lwd = 1)), 0.66, 0.70, 0.25, 0.298) +
  geom_segment(aes(x = 0.6, y = 0.53, xend = 0.66+0.25/2, yend = 0.7), color = "gray20") +
  draw_plot(middle_easte_na_plot, 0.37, 0.0, 0.25, 0.3) +
  draw_grob(rectGrob(gp = gpar(col = "gray55", fill = NA, lwd = 1)), 0.37, 0.0, 0.25, 0.298) +
  geom_segment(aes(x = 0.56, y = 0.54, xend = 0.37+0.25/1.2, yend = 0.30), color = "gray20") +
  draw_plot(latam_carib_plot, 0.01, 0.35, 0.25, 0.3) +
  draw_grob(rectGrob(gp = gpar(col = "gray55", fill = NA, lwd = 1)), 0.01, 0.35, 0.25, 0.3) +
  geom_segment(aes(x = 0.39, y = 0.5, xend = 0.26, yend = 0.5), color = "gray20") +
  draw_plot(east_asia_plot, 0.73, 0.35, 0.25, 0.3) +
  draw_grob(rectGrob(gp = gpar(col = "gray55", fill = NA, lwd = 1)), 0.73, 0.35, 0.25, 0.3) +
  geom_segment(aes(x = 0.66, y = 0.5, xend = 0.73, yend = 0.5), color = "gray20") +
  draw_plot(africa_plot, 0.08, 0.0, 0.25, 0.3) +
  draw_grob(rectGrob(gp = gpar(col = "gray55", fill = NA, lwd = 1)), 0.08, 0.0, 0.25, 0.3) +
  geom_segment(aes(x = 0.50, y = 0.49, xend = 0.30, yend = 0.3), color = "gray20") +
  draw_plot(map_legend, 0.68, 0.06, 0.1, 0.20)+
  draw_plot(virus_legend, 0.75, 0.06, 0.25, 0.22) -> plot_map_and_subplots


plot_map_and_subplots
ggsave("Figure_3.jpeg", width = 13, height = 9, dpi = 600)
