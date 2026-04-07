# Script to obtain trends of population growth and deforestation rates within hotspot areas and reproduce Figure 4
require(tidyverse)
require(sf)
require(raster)
require(terra)
require(exactextractr)
require(rnaturalearth)

# Load hotspot rasters of the different regions
hotspot_eap <- rast(".../data/eap_raster.tif") %>% as.polygons(., na.rm = TRUE, dissolve = TRUE) %>% st_as_sf() # %>% st_transform("+proj=longlat +datum=WGS84 +no_defs +type=crs")
hotspot_na <- rast(".../data/na_raster.tif") %>% as.polygons(., na.rm = TRUE, dissolve = TRUE) %>% st_as_sf() 
hotspot_sa <- rast(".../data/sa_raster.tif")%>% as.polygons(., na.rm = TRUE, dissolve = TRUE) %>% st_as_sf() 
hotspot_ssa <- rast(".../data/ssa_raster.tif")%>% as.polygons(., na.rm = TRUE, dissolve = TRUE) %>% st_as_sf()
hotspot_eca <- rast(".../data/eca_raster.tif")%>% as.polygons(., na.rm = TRUE, dissolve = TRUE) %>% st_as_sf()  
hotspot_lac <- rast(".../data/lac_raster.tif") %>% as.polygons(., na.rm = TRUE, dissolve = TRUE) %>% st_as_sf() 
hotspot_mena <- rast(".../data/mena_raster.tif")%>% as.polygons(., na.rm = TRUE, dissolve = TRUE) %>% st_as_sf()

# Background regions
countries <- ne_countries(scale = 50, returnclass = c("sf"))
countries_ee <- st_transform(countries, "+proj=eqearth")
countries_ee_grouped <- countries_ee %>%
  group_by(`region_wb`) %>%
  summarise(n = n())

blank_raster_ee <- rast(xmin = -16920565,
                        xmax = 16920565,
                        ymax = 8315130,
                        ymin = -8392928,
                        resolution = 1000,
                        crs = "+proj=eqearth")

# Population data from WorldPop
raster_files <- list.files(pathPop, all.files = T)
raster_files_pop <- raster_files[-c(1, 2)]
pop_rasters <- lapply(raster_files_pop, rast)

# Stack the rasters
raster_stack <- rast(c(pop_rasters[c(2:21)]))
# Reproj to ee
stack_pop <- terra::project(raster_stack, blank_raster_ee, method = "near")

# Forest data from Hansen et al., 2013
raster_files_forest <- list.files(pathForest, all.files = T, pattern = "treecover900")
forest_rasters <- lapply(raster_files_forest, rast)
# Reproj to ee
stack_forest <- terra::project(raster_stack_forest, blank_raster_ee, method = "near")

#### POLYGON TOTAL ESTIMATES ####
# extract Population data from hotspots
hotspot_eap_pop <- exact_extract(stack_pop, hotspot_eap, fun = "mean")
hotspot_eca_pop <- exact_extract(stack_pop, hotspot_eca, fun = "mean")
hotspot_mena_pop <- exact_extract(stack_pop, hotspot_mena, fun = "mean")
hotspot_sa_pop <- exact_extract(stack_pop, hotspot_sa, fun = "mean")
hotspot_ssa_pop <- exact_extract(stack_pop, hotspot_ssa, fun = "mean")
hotspot_lac_pop <- exact_extract(stack_pop, hotspot_lac, fun = "mean")
hotspot_na_pop <- exact_extract(stack_pop, hotspot_na, fun = "mean")


rbind(tibble(Population = hotspot_eap_pop) %>% mutate("region" = "East Asia & Pacific"), 
      tibble(Population = hotspot_eca_pop) %>% mutate("region" = "Europe & Central Asia"),
      tibble(Population = hotspot_mena_pop) %>% mutate("region" = "Middle East & North Africa"),
      tibble(Population = hotspot_sa_pop) %>% mutate("region" = "South Asia"),
      tibble(Population = hotspot_ssa_pop) %>% mutate("region" = "Sub-Saharan Africa"),
      tibble(Population = hotspot_lac_pop) %>% mutate("region" = "Latin America & Caribbean"),
      tibble(Population = hotspot_na_pop) %>% mutate("region" = "North America")) -> hotspots_pop_static

hotspots_pop_total$Year <- str_replace(hotspots_pop_total$Year, "sum.ppp_", "")
hotspots_pop_total$Year <- str_replace(hotspots_pop_total$Year, "_1km_Aggregated", "")

hotspots_pop_total$`Growth Rate` <- rep(NA, nrow(hotspots_pop_total))

hotspots_pop_total_final <- NULL
for (r in unique(hotspots_pop_total$region)){
  
  hotspots_pop_total_temp <- hotspots_pop_total %>% filter(region == r) %>% arrange(Year)
  
  for(i in seq(1:20)){
    
    t0 <- hotspots_pop_total_temp[i, 3]
    t1 <- hotspots_pop_total_temp[i+1, 3]
    
    hotspots_pop_total_temp[i+1, 4] <- (t1-t0)/t0
    
  }
  
  hotspots_pop_total_final <- rbind(hotspots_pop_total_final, hotspots_pop_total_temp)
  
}

hotspots_pop_total_final$`Growth Rate` <- hotspots_pop_total_final$`Growth Rate`*100


# extract Population data from background (whole region)
bckg_ssa_pop <- exact_extract(stack_pop, countries_ee_grouped[8, ], fun = "mean")
bckg_eap_pop <- exact_extract(stack_pop, countries_ee_grouped[2, ], fun = "mean")
bckg_eca_pop <- exact_extract(stack_pop, countries_ee_grouped[3, ], fun = "mean")
bckg_lac_pop <- exact_extract(stack_pop, countries_ee_grouped[4, ], fun = "mean")
bckg_mena_pop <- exact_extract(stack_pop, countries_ee_grouped[5, ], fun = "mean")
bckg_na_pop <- exact_extract(stack_pop, countries_ee_grouped[6, ], fun = "mean")
bckg_sa_pop <- exact_extract(stack_pop, countries_ee_grouped[7, ], fun = "mean")


rbind(bckg_eap_pop %>% mutate("region" = "East Asia & Pacific"), 
      bckg_eca_pop %>% mutate("region" = "Europe & Central Asia"),
      bckg_mena_pop %>% mutate("region" = "Middle East & North Africa"),
      bckg_sa_pop %>% mutate("region" = "South Asia"),
      bckg_ssa_pop %>% mutate("region" = "Sub-Saharan Africa"),
      bckg_lac_pop %>% mutate("region" = "Latin America & Caribbean"),
      bckg_na_pop %>% mutate("region" = "North America")) %>% 
  pivot_longer(cols = 1:21, names_to = "Year", values_to = "Population") -> hotspots_pop_bck

hotspots_pop_bck$Year <- str_replace(hotspots_pop_bck$Year, "sum.ppp_", "")
hotspots_pop_bck$Year <- str_replace(hotspots_pop_bck$Year, "_1km_Aggregated", "")

hotspots_pop_bck_final <- NULL
hotspots_pop_bck$`Growth Rate` <- rep(NA, nrow(hotspots_pop_bck))

for (r in unique(hotspots_pop_bck$region)){
  
  hotspots_pop_bck_temp <- hotspots_pop_bck %>% filter(region == r) %>% arrange(Year)
  
  for(i in seq(1:20)){
    
    t0 <- hotspots_pop_bck_temp[i, 3]
    t1 <- hotspots_pop_bck_temp[i+1, 3]
    
    hotspots_pop_bck_temp[i+1, 4] <- (t1-t0)/t0
    
  }
  
  hotspots_pop_bck_final <- rbind(hotspots_pop_bck_final, hotspots_pop_bck_temp)
  
}

hotspots_pop_bck_final$`Growth Rate` <- hotspots_pop_bck_final$`Growth Rate`*100


# extract forest data from hotspots
hotspot_eap_forest <- exact_extract(stack_forest, hotspot_eap, fun = "sum")
hotspot_eca_forest <- exact_extract(stack_forest, hotspot_eca, fun = "sum")
hotspot_mena_forest <- exact_extract(stack_forest, hotspot_mena, fun = "sum")
hotspot_sa_forest <- exact_extract(stack_forest, hotspot_sa, fun = "sum")
hotspot_ssa_forest <- exact_extract(stack_forest, hotspot_ssa, fun = "sum")
hotspot_lac_forest <- exact_extract(stack_forest, hotspot_lac, fun = "sum")
hotspot_na_forest <- exact_extract(stack_forest, hotspot_na, fun = "sum")

rbind(hotspot_eap_forest %>% mutate("region" = "East Asia & Pacific"), 
      hotspot_eca_forest %>% mutate("region" = "Europe & Central Asia"),
      hotspot_mena_forest %>% mutate("region" = "Middle East & North Africa"),
      hotspot_sa_forest %>% mutate("region" = "South Asia"),
      hotspot_ssa_forest %>% mutate("region" = "Sub-Saharan Africa"),
      hotspot_lac_forest %>% mutate("region" = "Latin America & Caribbean"),
      hotspot_na_forest %>% mutate("region" = "North America")) %>% 
  pivot_longer(cols = 1:24, names_to = "Year", values_to = "Forest Cover") -> hotspots_forest_trend

hotspots_forest_trend$Year <- str_replace(hotspots_forest_trend$Year, "sum.treecover900_", "")

hotspots_forest_trend$`Forest Cover` <- hotspots_forest_trend$`Forest Cover`/100 # it was % coverage (now sqrd km)
hotspots_forest_trend$`Forest Loss` <- rep(NA, nrow(hotspots_forest_trend))
hotspots_forest_trend$`Percentage Loss` <- rep(NA, nrow(hotspots_forest_trend))

hotspots_forest_trend_final <- NULL
for (r in unique(hotspots_forest_trend$region)){
  
  hotspots_forest_trend_temp <- hotspots_forest_trend %>% filter(region == r)%>% arrange(Year)
  
  for(i in seq(1:23)){
    
    t0 <- hotspots_forest_trend_temp[i, 3]
    t1 <- hotspots_forest_trend_temp[i+1, 3]
    
    hotspots_forest_trend_temp[i+1, 4] <- t0-t1
    hotspots_forest_trend_temp[i+1, 5] <- (t0-t1)/t0
    
  }
  
  hotspots_forest_trend_final <- rbind(hotspots_forest_trend_final, hotspots_forest_trend_temp)
  
}

hotspots_forest_trend_final$`Percentage Loss` <- hotspots_forest_trend_final$`Percentage Loss`*100

#extract forest data from background
bckg_ssa_forest <- exact_extract(stack_forest, countries_ee_grouped[8, ], fun = "sum")
bckg_eap_forest <- exact_extract(stack_forest, countries_ee_grouped[2, ], fun = "sum")
bckg_eca_forest <- exact_extract(stack_forest, countries_ee_grouped[3, ], fun = "sum")
bckg_lac_forest <- exact_extract(stack_forest, countries_ee_grouped[4, ], fun = "sum")
bckg_mena_forest <- exact_extract(stack_forest, countries_ee_grouped[5, ], fun = "sum")
bckg_na_forest <- exact_extract(stack_forest, countries_ee_grouped[6, ], fun = "sum")
bckg_sa_forest <- exact_extract(stack_forest, countries_ee_grouped[7, ], fun = "sum")

rbind(bckg_eap_forest %>% mutate("region" = "East Asia & Pacific"), 
      bckg_eca_forest %>% mutate("region" = "Europe & Central Asia"),
      bckg_mena_forest %>% mutate("region" = "Middle East & North Africa"),
      bckg_sa_forest %>% mutate("region" = "South Asia"),
      bckg_ssa_forest %>% mutate("region" = "Sub-Saharan Africa"),
      bckg_lac_forest %>% mutate("region" = "Latin America & Caribbean"),
      bckg_na_forest %>% mutate("region" = "North America")) %>% 
  pivot_longer(cols = 1:24, names_to = "Year", values_to = "Forest Cover") -> hotspots_forest_bck

hotspots_forest_bck$Year <- str_replace(hotspots_forest_bck$Year, "sum.treecover900_", "")

hotspots_forest_bck$`Forest Cover` <- hotspots_forest_bck$`Forest Cover`/100 # it was % coverage (now sqrd km)
hotspots_forest_bck$`Forest Loss` <- rep(NA, nrow(hotspots_forest_bck))
hotspots_forest_bck$`Percentage Loss` <- rep(NA, nrow(hotspots_forest_bck))

# Static deforestation
hotspots_forest_static_final <- hotspots_forest_trend_final %>% 
  filter(Year %in% c(2000, 2020)) %>% 
  dplyr::select(c(region, `Forest Cover`, Year)) %>% 
  pivot_wider(values_from = `Forest Cover`, names_from = Year) %>% 
  group_by(region) %>% 
  mutate(`Percentage Loss` = 100*(`2000`-`2020`)/ `2000`)

bck_forest_static_final <- hotspots_forest_trend_bck_final %>% 
  filter(Year %in% c(2000, 2020)) %>% 
  dplyr::select(c(region, `Forest Cover`, Year)) %>% 
  pivot_wider(values_from = `Forest Cover`, names_from = Year) %>% 
  group_by(region) %>% 
  mutate(`Percentage Loss` = 100*(`2000`-`2020`)/ `2000`)

rbind(bck_forest_static_final %>% dplyr::select(region, `Percentage Loss`) %>% mutate(backg = "region"),
      hotspots_forest_static_final %>% dplyr::select(region, `Percentage Loss`) %>% mutate(backg = "main")) %>% 
  mutate(driver = "Deforestation") -> def_comp

#
rbind(bcks_pop_static %>% mutate(backg = "region") %>% mutate(driver = "Population"),
      hotspots_pop_static %>% mutate(backg = "main") %>% mutate(driver = "Population")) %>% 
  rename(mean = Population) -> pop_comp

comp_all <- rbind(pop_comp,
                  def_comp %>% rename(mean = `Percentage Loss`) %>%  mutate(mean = mean*19.65176)) %>% 
  mutate(driver = fct_relevel(driver, "Population"),
         backg = fct_relevel(backg, "main"))


ggplot() +
  geom_smooth(data = hotspots_forest_trend_bck_final %>% filter(as.numeric(as.character(Year)) < 2021), aes(y = `Percentage Loss`*3, x = as.numeric(as.character(Year)), color = "Regional deforestation rate"), size = 1, se = F) +
  geom_smooth(data = hotspots_pop_bck_final, aes(y = `Growth Rate`, x = as.numeric(as.character(Year)), color = "Regional population growth rate"), size = 1, se = F) +
  geom_smooth(data = hotspots_pop_total_final, aes(y = `Growth Rate`, x = as.numeric(as.character(Year)), color = "Population growth rate"), size = 1, se = F) +
  geom_smooth(data = hotspots_forest_trend_final %>% filter(as.numeric(as.character(Year)) < 2021), aes(y = `Percentage Loss`*3, x = as.numeric(as.character(Year)), color = "Deforestation rate"), size = 1, se = F) +
  facet_wrap(~region, nrow = 7) +
  scale_color_manual(
    name = "Trend",
    values = c("Deforestation rate" = "#012169", "Population growth rate" = "#9e1b32",
               "Regional population growth rate" = "#D8A4AD", "Regional deforestation rate" = "#99A6C3")
  ) +
  theme_BP()+
  theme(legend.position = "bottom")+
  # theme_minimal(base_size = 13) +
  scale_y_continuous(name = "Population growth rate (%)",
                     sec.axis = sec_axis(~ ./3, name = "Deforestation rate (%)"))+
  scale_x_continuous(
    breaks = seq(2000, 2020, 5),
    limits = c(2000, 2020)
  )+
  # theme_minimal(base_size = 13) +
  theme(
    axis.line.y.left = element_line(color = "black", size = 0.5),
    axis.line.y.right = element_line(color = "black", size = 0.5),
    axis.ticks.y.left = element_line(color = "black"),
    axis.ticks.y.right = element_line(color = "black"),
    axis.title.y = element_text(color = "black"),
    axis.title.y.right = element_text(color = "black"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )+  guides(fill = "none", alpha = "none", color = "none")+
  labs(x = "Year") -> left_plot

ggplot() +
  geom_bar(
    data = comp_all,
    aes(y = mean, fill = driver, x = driver, alpha = backg),
    stat = "identity",
    width = 0.9,
    position = "dodge2"
  ) +
  scale_fill_manual(values = c("Deforestation" = "#012169", "Population" = "#9e1b32")) +
  scale_alpha_manual(values = c("region" = 0.4, "main" = 1), labels = c("Hotspot", "Whole region")) +
  facet_wrap(~region, nrow = 7) +
  scale_y_continuous(sec.axis = sec_axis(~ . /19.65176, name = "Forest cover loss % (2000–2020)")) +
  theme_BP()+
  labs(
    y = expression("Population density (people/km" ~ plain("²") * ")"),
    x = "",
    alpha = "",
    fill = ""
  ) +
  theme(
    axis.line.y.left = element_line(color = "black", size = 0.5),
    axis.line.y.right = element_line(color = "black", size = 0.5),
    axis.ticks.y.left = element_line(color = "black"),
    axis.ticks.y.right = element_line(color = "black"),
    axis.title.y = element_text(color = "black"),
    axis.title.y.right = element_text(color = "black"),
    strip.text = element_text(face = "bold"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "bottom"
  ) -> right_plot

require(patchwork)
space <- ggplot() + theme_void()

label_df <- data.frame(
  region = unique(hotspots_forest_trend_bck_final$region),
  label = c("A", "B", "D", "F", "G", "C", "E"),         
  x = 2000.5,                      
  y = 3.2                       
)

left_plot +
  geom_text(
    data = label_df,
    aes(x = x, y = y, label = label),
    fontface = "bold",
    size = 5                    
  ) -> left_plot_labelled

combined_plot <- left_plot_labelled/space/right_plot
  plot_layout(heights = c(1, 1, 1),
              widths = c(1, 0.1, 1),
              ncol = 3, nrow = 1,
              guides = "collect")  &
  theme(legend.position = "bottom")


ggsave(combined_plot, "Figure_4.jpeg", width = 6, height = 11, dpi = 600)
