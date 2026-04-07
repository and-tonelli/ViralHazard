# Script to reproduce Figures 1-2 and Supplementary Figure 2

'%!in%' <- function(x, y)!('%in%'(x, y)) # %!in%

require(tidyverse)
require(rphylopic)
require(patchwork)
require(cowplot)

# Figures' theme
theme_BP <- function(){ 
  font <- "Helvetica"   #assign font family up front
  
  theme_bw() %+replace%    #replace elements we want to change
    
    theme(
      
      #grid elements
      panel.grid.major = element_line(colour = "grey95"),    #strip major gridlines
      panel.grid.minor = element_line(size = rel(0.2)),    #strip minor gridlines
      axis.ticks = element_line(size = rel(1)),          #strip axis ticks
      
      #since theme_minimal() already strips axis lines, 
      #we don't need to do that again
      
      #text elements
      plot.title = element_text(             #title
        family = font,            #set font family
        size = 20,                #set font size
        face = 'bold',            #bold typeface
        hjust = 0,                #left align
        vjust = 2),               #raise slightly
      
      plot.subtitle = element_text(          #subtitle
        family = font,            #font family
        size = 14),               #font size
      
      plot.caption = element_text(           #caption
        family = font,            #font family
        size = 9,                 #font size
        hjust = 1),               #right align
      
      axis.title = element_text(             #axis titles
        family = font,            #font family
        size = 10),               #font size
      
      axis.text = element_text(              #axis text
        family = font,            #axis famuly
        size = 9),                #font size
      
      axis.text.x = element_text(            #margin for axis text
        margin=margin(2, b = 10)),
      
      #since the legend often requires manual tweaking 
      #based on plot content, don't define it here
      
      panel.border = element_rect(fill = NA,
                                  colour = "grey20",
                                  margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")),
      
      strip.background = element_rect(fill = NA,
                                      colour = NA),
      
      strip.text = element_text(family = font,
                                size = 15,
                                margin=margin(2, b = 10)),
      
      axis.line.x.bottom = element_line(
        colour = "grey20",
        linewidth = 0.5
      ),
      
      axis.line.y.left = element_line(
        colour = "grey20",
        linewidth = 0.5
      )
      
    )
}
theme_set(theme_BP())

# Loading datasets
ebv <- read.csv(".../data/ebv_dataset_main.csv") %>% dplyr::select(c(1:4, 33)) %>% mutate(`Virus genus` = "Ebolavirus & Marburgvirus")
hpv <- read.csv(".../data/hpv_dataset_main.csv") %>% dplyr::select(c(1:4, 33)) %>% mutate(`Virus genus` = "Henipavirus")
bcv <- read.csv(".../data/bcov_dataset_main.csv") %>% dplyr::select(c(1:4, 33)) %>% mutate(`Virus genus` = "Betacoronavirus")
phv <- read.csv(".../data/phv_dataset_main.csv") %>% dplyr::select(c(1:4, 33)) %>% mutate(`Virus genus` = "Phlebovirus")
ohv <- read.csv(".../data/ohv_dataset_main.csv") %>% dplyr::select(c(1:4, 33)) %>% mutate(`Virus genus` = "Orthonairovirus")
flv <- read.csv(".../data/flv_dataset_main.csv") %>% dplyr::select(c(1:4, 33)) %>% mutate(`Virus genus` = "Flavivirus")
mmv <- read.csv(".../data/mmv_dataset_main.csv") %>% dplyr::select(c(1:4, 33)) %>% mutate(`Virus genus` = "Mammarenavirus")

blueprint_dataset <- rbind(ebv, hpv, bcv, phv, ohv, flv, mmv)
blueprint_dataset$`Virus genus` <- factor(blueprint_dataset$`Virus genus`, levels = unique(blueprint_dataset$`Virus genus`))

# Loading predictions
ebv_final_predictions <- read_csv(".../data/ebv_predictions_main.csv") %>%   dplyr::select(c(1:5, 7, 8, 9)) %>% mutate(`Virus genus` = "Ebolavirus & Marburgvirus")
hpv_final_predictions <- read_csv(".../data/hpv_predictions_main.csv") %>%   dplyr::select(c(1:5, 7, 8, 9)) %>% mutate(`Virus genus` = "Henipavirus")
bcov_final_predictions <- read_csv(".../data/bcov_predictions_main.csv") %>%   dplyr::select(c(1:5, 6, 7, 8)) %>% mutate(`Virus genus` = "Betacoronavirus")
phv_final_predictions <- read_csv(".../data/phv_predictions_main.csv") %>%   dplyr::select(c(1:5, 7, 8, 9)) %>% mutate(`Virus genus` = "Phlebovirus")
ohv_final_predictions <- read_csv(".../data/ohv_predictions_main.csv") %>%   dplyr::select(c(1:5, 7, 8, 9)) %>% mutate(`Virus genus` = "Orthonairovirus")
flv_final_predictions <- read_csv(".../data/flv_predictions_main.csv") %>%   dplyr::select(c(1:5, 7, 8, 9)) %>% mutate(`Virus genus` = "Flavivirus")
mmv_final_predictions <- read_csv(".../data/mmv_predictions_main.csv") %>%  dplyr::select(c(1:5, 7, 8, 9)) %>% mutate(`Virus genus` = "Mammarenavirus")

predictions_bp <- rbind(ebv_final_predictions, hpv_final_predictions, bcov_final_predictions, phv_final_predictions, ohv_final_predictions, flv_final_predictions, mmv_final_predictions)
predictions_bp$`Virus genus` <- factor(predictions_bp$`Virus genus`, levels = names(palette_vir))

palette_vir <- c("#71014B", # burgundy
                 "#FF5C77", # pink
                 "#FFB300", # orange
                 "#5DC122", # green
                 "#046C7A", # aqua
                 "#002A86", # blue
                 "#470CED") # violet

names(palette_vir) <- unique(blueprint_dataset$`Virus genus`)

#### Figure 1 ####
predictions_bp %>% 
  filter(predictable_fam == 1, (hard_class == 1 | real_status %in% c("high-evidence host", "low-evidence host"))) %>% 
  mutate(status = ifelse(real_status %in% c("unknown", "pseudo-negative"), "Unknown", "Observed")) %>% 
  group_by(order, `Virus genus`, status) %>%  # Group by order, Virus genus, and status
  summarize(num = n()) %>%
  mutate(num = ifelse(status == "Unknown", sum(num), num)) %>% 
  mutate(num = ifelse(num == 0, 0.01, num)) %>% 
  group_by(order, `Virus genus`) -> predictions_bp2


predictions_bp2 <- predictions_bp2 %>% 
  ungroup() %>% 
  complete(order, `Virus genus`, status, fill = list(num = as.integer(0))) #1 for plot

# phylopic IDs were previously saved

predictions_bp2 %>% 
  ggplot(aes(x = order, y = num)) +
  geom_hline(yintercept = 2, linetype = "solid", color = "gray88")+
  geom_hline(yintercept = 4, linetype = "solid", color = "gray95")+
  geom_hline(yintercept = 8, linetype = "solid", color = "gray88")+
  geom_hline(yintercept = 16, linetype = "solid", color = "gray95")+
  geom_hline(yintercept = 32, linetype = "solid", color = "gray88")+
  geom_hline(yintercept = 64, linetype = "solid", color = "gray95")+
  geom_hline(yintercept = 128, linetype = "solid", color = "gray88")+
  geom_hline(yintercept = 256, linetype = "solid", color = "gray95")+
  geom_bar(position = position_dodge(width = 0.8), aes(fill = `Virus genus`), alpha = ifelse(predictions_bp2$status == "Unknown" , 0.5, 1),
           stat = "identity", width = 0.8) +
  geom_segment(data = data.frame(x = seq(1.5, 13.5, by = 1)),
               aes(x = x, xend = x, y = 1, yend = max(predictions_bp2$num)),
               linetype = "dotted",
               color = "gray78")+
  add_phylopic(uuid = car_uuid, x = 1, y = 315, fill = "grey30", color = "grey30", ysize = 14/25)+
  add_phylopic(uuid = art_uuid, x = 2, y = 315, fill = "grey30", color = "grey30", ysize = 20/25)+
  add_phylopic(uuid = chi_uuid, x = 3, y = 335, fill = "grey30", color = "grey30", ysize = 13/30)+
  add_phylopic(uuid = did_uuid, x = 4, y = 315, fill = "grey30", color = "grey30", ysize = 10/25)+
  add_phylopic(uuid = dip_uuid, x = 5, y = 315, fill = "grey30", color = "grey30", ysize = 28/25)+
  add_phylopic(uuid = eul_uuid, x = 6, y = 315, fill = "grey30", color = "grey30", ysize = 10/25)+
  add_phylopic(uuid = lag_uuid, x = 7, y = 315, fill = "grey30", color = "grey30", ysize = 20/25)+
  add_phylopic(uuid = per_uuid, x = 8, y = 315, fill = "grey30", color = "grey30", ysize = 16/25)+
  add_phylopic(uuid = pho_uuid, x = 9, y = 315, fill = "grey30", color = "grey30", ysize = 10/25)+
  add_phylopic(uuid = pil_uuid, x = 10, y = 315, fill = "grey30", color = "grey30", ysize = 16/25)+
  add_phylopic(uuid = pri_uuid, x = 11, y = 325, fill = "grey30", color = "grey30", ysize = 24/25)+
  add_phylopic(uuid = pro_uuid, x = 12, y = 315, fill = "grey30", color = "grey30", ysize = 24/25)+
  add_phylopic(uuid = rod_uuid, x = 13, y = 335, fill = "grey30", color = "grey30", ysize = 10/25)+
  add_phylopic(uuid = sca_uuid, x = 14, y = 315, fill = "grey30", color = "grey30", ysize = 20/25)+
  scale_fill_manual(values = palette_vir) +
  theme_minimal() +
  xlab("") +
  ylab("") +
  scale_y_continuous(
    trans = 'log2', 
    breaks = scales::trans_breaks("log2", function(x) 2^x), 
    # labels = scales::trans_format("log2", scales::math_format(2^.x)), 
    labels = scales::trans_format("log2", function(x) as.character(2^x)),
    limits = c(0.1, 540)
  ) +
  coord_polar(start = 0.1, direction = 1, clip = "off") +
  theme(
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),   
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank()
  ) +
  # Add a central y-axis
  annotate("segment", x = 0.5, xend = 14.5, y = 0.05, yend = 0.05, color = "black", size = 0.8) + # Bottom axis
  annotate("text", x = 0, y = c(1, 2, 4, 8, 16, 32, 64, 128, 256), 
           label = c(1, 2, 4, 8, 16, 32, 64, 128, 256), angle = 0, 
           color = "black", size = 3)+
  theme(legend.position = "none") +
  annotate("text", x = 0, y = 540, 
           label = "Number\nof species", angle = 0, 
           color = "black", size = 4) -> main_plot


cowplot::plot_grid(main_plot, ggpubr::as_ggplot(legenda), 
                   ncol = 2,
                   rel_widths = c(25,1)) 

final_circular <- main_plot + ggpubr::as_ggplot(legenda) + plot_layout(widths = c(5, 1))

ggsave(plot = final_circular, filename = "Figure_1.jpeg", width = 10, height = 10, dpi = 600, units = "in")


#### Figure 2 ####
require(raster)
require(sf)
library(patchwork)

equalearth <- CRS("+proj=eqearth")

# countries shapes (background)
require(rnaturalearth)
countries <- ne_countries(scale = 50, returnclass = c("sf"))
countries_ee <- st_transform(countries, equalearth)

# create a bounding box for EqualEarth projection
bb <- sf::st_union(sf::st_make_grid(
  st_bbox(c(xmin = -180,
            xmax = 180,
            ymax = 90,
            ymin = -90), crs = st_crs(4326)),
  n = 100))

bb_ee <- st_transform(bb, as.character(equalearth))

# List of virus-specific color palettes
list(ebv = list(letter = "A",
                fullnames = "Ebolavirus & Marburgvirus",
                color = c("#D8D8D8FF", "#CEC2C9FF", "#C5ABB5FF", "#BB94A0FF", "#B07D8BFF", "#A56676FF", "#994F61FF", "#8B394BFF", "#7D2236FF")),#", "#71014BFF
     hpv = list(letter = "B",
                fullnames = "Henipavirus",
                color = c("#D8D8D8FF", "#EACDD2FF", "#E7B6C3FF", "#E49EB5FF", "#E187A7FF", "#DE7099FF", "#DB598BFF", "#D8417DFF", "#D52A6FFF")),
     bcov = list(letter = "C",
                 fullnames = "Betacoronavirus",
                 color = c("#D8D8D8FF", "#EADDC4FF", "#E4CFACFF", "#DEBF94FF", "#D7AF7CFF", "#D19F64FF", "#CB8F4CFF", "#C47F34FF", "#BE6F1CFF")),
     phv = list(letter = "D",
                fullnames = "Phlebovirus",
                color = c("#D8D8D8FF", "#CBD8C7FF", "#B7D1B5FF", "#A3C9A4FF", "#8FC292FF", "#7BBB81FF", "#67B36FFF", "#52AC5EFF", "#3EA44CFF")),
     ohv = list(letter = "E",
                fullnames = "Orthonairovirus",
                color = c("#D8D8D8FF", "#C2D3D4FF", "#A5C6C9FF", "#89BABEFF", "#6CAEB2FF", "#4FA2A6FF", "#33869AFF", "#176B8EFF", "#045E82FF")),# "#046C7AFF"
     flv = list(letter = "F",
                fullnames = "Flavivirus",
                color = c("#D8D8D8FF", "#C7CDE0FF", "#AEBBD3FF", "#94A8C5FF", "#7B96B8FF", "#6284AAFF", "#49629DFF", "#30408FFF",  "#002A86FF")), #"#172082FF", ##CBCBD6FF", "#B3B2C9FF"
     mmv = list(letter = "G",
                fullnames = "Mammarenavirus",
                color =  c("#D8D8D8FF", "#CDC5E7FF", "#B8ABDFFF", "#A292D8FF", "#8C78D0FF", "#775FC8FF", "#6146BFFF", "#4B2DB7FF", "#3415AFFF"))) -> listvir # "#470CEDFF"

for (vir in c("hpv", "ebv", "ohv", "flv", "mmv", "bcov", "phv")){
  
  hotspots_observed <- raster(paste0(".../data/", vir, "_observed_2025.tif"))
  hotspots_observed[hotspots_observed[] == "0"] <- NA
  hotspots_observed[hotspots_observed[] == "0"] <- NA
  
  hotspots_predicted <- raster(paste0(".../data/", vir, "_pred_hotspots_2025.tif"))
  hotspots_predicted[hotspots_predicted[] == "0"] <- NA
  
  # Use the bounding box to trim (mask) the raster layer (otherwise it will have weird duplicated areas in the corners)
  mskd_hotspots_observed <- raster::mask(hotspots_observed, as(bb_ee, "Spatial"))
  # mskd_hotspots_observed_agg <- raster::mask(hotspots_observed_agg, as(bb_ee, "Spatial"))
  mskd_hotspots_predicted <- raster::mask(hotspots_predicted, as(bb_ee, "Spatial"))
  
  # RasterLayer to DF (in order to plot it with ggplot2)
  df_observed <- mskd_hotspots_observed %>%
    rasterToPoints %>%
    as.data.frame() %>%
    `colnames<-`(c("x", "y", "richness"))
  
  df_predicted <-  mskd_hotspots_predicted %>%
    rasterToPoints %>%
    as.data.frame() %>%
    `colnames<-`(c("x", "y", "richness")) 
  
  df_obs_pred <- rbind(df_observed %>% mutate(richnesstype = "Observed"), df_predicted %>% mutate(richnesstype = "Observed & Predicted"))
  
  ggplot()+
    geom_sf(data = bb_ee,   # BACKGROUND
            colour = "gray45",
            linetype = 'solid',
            fill = "white", # light grey
            size = 0.1) +
    geom_sf(data = countries_ee, # CONTINENTS
            colour = NA,
            linetype = 'solid',
            fill = "gray45", #dark grey #3B444D
    ) +
    geom_raster(
      data = df_obs_pred,   # HOTSPOTS
      aes(
        x = x,
        y = y,
        fill = richness))+
    scale_fill_gradientn(colours = listvir[[vir]][["color"]],
                         name = "Host\nrichness",
                         breaks = c(1, max(df_predicted$richness)))+
    facet_wrap(~ richnesstype)+
    theme_void() +
    theme(legend.position = "right",
          legend.box = "vertical",
          legend.key.width = unit(0.5, 'cm'),
          legend.key.height = unit(0.7, 'cm'),
          plot.title = element_text(hjust = 0.5, vjust = 8, size = 10, face = "bold")) +
    guides(fill = guide_colorbar(title.position = "right", title.hjust = 0.5),
           size = guide_legend(title.position = "right", title.hjust = 0))+
    labs(tag = listvir[[vir]][["letter"]])+
    ggtitle(listvir[[vir]][["fullnames"]])+
    theme(plot.tag = element_text(hjust = 0.8, size = 15, face = "bold"),
          strip.background = element_blank(),
          plot.title = element_text(margin=margin(0,0,-15,0), face = "italic"),
          strip.text.x = if (vir != "ebv") element_blank() else element_text(angle = 0)# element_text(ifelse(vir == "hpv", NULL, size = 15))
          # ,
          # plot.title.position = "plot"
          ) -> virplot

  assign(paste0("plot", vir), virplot)

}

ggarrange(plotebv, plothpv, plotbcov, plotphv, plotohv, plotflv, plotmmv, ncol = 1, heights = c(1.12, 1, 1, 1, 1, 1, 1))
ggsave("Figure_2.jpeg", width = 9, height = 14, dpi = 600)


#### Figure S2 ####
Genbank_final <- read.csv(".../data/DatasetS2.csv")

unique(Genbank_final %>% 
         filter(VirusGenus %in% c("Orthoebolavirus", "Orthomarburgvirus"), 
                species %in% ebv_final_predictions$iucn2020_binomial[ebv_final_predictions$real_status %in% c("unknown", "pseudo-negative")]) %>% 
         pull(species)) -> new_ebv_hosts

unique(Genbank_final %>% 
         filter(VirusGenus %in% "Orthoflavivirus", Virus %!in% c("St. Louis encephalitis virus", "Japanese encephalitis virus"), 
                species %in% flv_final_predictions$iucn2020_binomial[flv_final_predictions$real_status %in% c("unknown", "pseudo-negative")]) %>% 
         pull(species)) -> new_flv_hosts

unique(Genbank_final %>% 
         filter(VirusGenus %in% "Betacoronavirus", 
                species %in% bcov_final_predictions$iucn2020_binomial[bcov_final_predictions$real_status %in% c("unknown", "pseudo-negative")]) %>% 
         pull(species)) -> new_bcv_hosts

unique(Genbank_final %>% 
         filter(VirusGenus %in% "Orthonairovirus", 
                species %in% ohv_final_predictions$iucn2020_binomial[ohv_final_predictions$real_status %in% c("unknown", "pseudo-negative")]) %>% 
         pull(species)) -> new_ohv_hosts

unique(Genbank_final %>% 
         filter(VirusGenus %in% "Phlebovirus", 
                species %in% phv_final_predictions$iucn2020_binomial[phv_final_predictions$real_status %in% c("unknown", "pseudo-negative")]) %>% 
         pull(species)) -> new_phv_hosts

unique(Genbank_final %>% 
         filter(VirusGenus %in% "Henipavirus", 
                species %in% hpv_final_predictions$iucn2020_binomial[hpv_final_predictions$real_status %in% c("unknown", "pseudo-negative")]) %>% 
         pull(species)) -> new_hpv_hosts

unique(Genbank_final %>% 
         filter(VirusGenus %in% "Mammarenavirus", 
                species %in% mmv_final_predictions$iucn2020_binomial[mmv_final_predictions$real_status %in% c("unknown", "pseudo-negative")]) %>% 
         pull(species)) -> new_mmv_hosts

predictions_bp %>% 
  filter(predictable_fam == 1) %>%
  mutate(new_host = case_when(`Virus genus`  == "Ebolavirus & Marburgvirus" & iucn2020_binomial %in% new_ebv_hosts ~ 2,
                              `Virus genus`  == "Mammarenavirus" & iucn2020_binomial %in% new_mmv_hosts ~ 2,
                              `Virus genus`  == "Henipavirus" & iucn2020_binomial %in% new_hpv_hosts ~ 2,
                              `Virus genus`  == "Phlebovirus" & iucn2020_binomial %in% new_phv_hosts ~ 2,
                              `Virus genus`  == "Betacoronavirus" & iucn2020_binomial %in% new_bcv_hosts ~ 2,
                              `Virus genus`  == "Orthonairovirus" & iucn2020_binomial %in% new_ohv_hosts ~ 2,
                              `Virus genus`  == "Flavivirus" & iucn2020_binomial %in% new_flv_hosts ~ 2,
                              real_status %in% c("high-evidence host", "low-evidence host") ~ 1,
                              TRUE ~ 0)) %>% 
  group_by(`Virus genus`) %>% 
  mutate(facet_label = paste0(`Virus genus`, " (n = ", sum(new_host == 2), ")")) %>% 
  ungroup() %>% 
  filter(`Virus genus`  != "Ebolavirus & Marburgvirus") %>% 
  ggplot()+
  # geom_point(aes(y = mean_pred, x = factor(new_host), fill = factor(new_host)), alpha = 0.10,  position = position_jitterdodge(jitter.width = 0.80, jitter.height = 0, dodge.width = 0.7, seed = 42), color = "steelblue4", size = 3)+
  gghalves::geom_half_violin(aes(x = factor(new_host), y = mean_pred, fill = factor(new_host)), color = NA)+
  scale_fill_manual(values = c("#364C54FF", "#94475EFF", "#E5A11FFF"), label = c("Unknown", "Positive", "New positive"))+
  stat_summary(fun.y = "mean", mapping = aes(x = factor(new_host), y = mean_pred, group = factor(new_host)), color = "black", geom = "point", size = 2)+ #size = 6
  # ggrepel::geom_text_repel(aes(label = label, y = mean_pred, x = factor(new_host)), size = 3,
  #                          box.padding = unit(3, "lines"), max.overlaps = 1000, fontface = "italic",
  #                          position = position_jitter(width = 0.55, height = 0, seed = 42))+
  ylim(c(0, 1))+
  theme_BP()+
  facet_wrap(~ facet_label)+
  labs(x = "", y = "Host status probability", fill = "")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        strip.text = element_text(size = 10),
        legend.position = "bottom" #c(0.65, 0.16)
        ) -> plot_singles


predictions_bp %>% 
  filter(predictable_fam == 1) %>%
  mutate(new_host = case_when(`Virus genus`  == "Ebolavirus & Marburgvirus" & iucn2020_binomial %in% new_ebv_hosts ~ 2,
                              `Virus genus`  == "Mammarenavirus" & iucn2020_binomial %in% new_mmv_hosts ~ 2,
                              `Virus genus`  == "Henipavirus" & iucn2020_binomial %in% new_hpv_hosts ~ 2,
                              `Virus genus`  == "Phlebovirus" & iucn2020_binomial %in% new_phv_hosts ~ 2,
                              `Virus genus`  == "Betacoronavirus" & iucn2020_binomial %in% new_bcv_hosts ~ 2,
                              `Virus genus`  == "Orthonairovirus" & iucn2020_binomial %in% new_ohv_hosts ~ 2,
                              `Virus genus`  == "Flavivirus" & iucn2020_binomial %in% new_flv_hosts ~ 2,
                              real_status %in% c("high-evidence host", "low-evidence host") ~ 1,
                              TRUE ~ 0)) %>%
  ggplot()+
  # geom_point(aes(y = mean_pred, x = factor(new_host), fill = factor(new_host)), alpha = 0.10,  position = position_jitterdodge(jitter.width = 0.80, jitter.height = 0, dodge.width = 0.7, seed = 42), color = "steelblue4", size = 3)+
  gghalves::geom_half_violin(aes(x = factor(new_host), y = mean_pred, fill = factor(new_host)), color = NA)+
  scale_fill_manual(values = c("#364C54FF", "#94475EFF", "#E5A11FFF"), label = c("Unknown", "Positive", "New positive"))+
  stat_summary(fun.y = "mean", mapping = aes(x = factor(new_host), y = mean_pred, group = factor(new_host)), color = "black", geom = "point", size = 4.5)+ #size = 6
  # ggrepel::geom_text_repel(aes(label = label, y = mean_pred, x = factor(new_host)), size = 3,
  #                          box.padding = unit(3, "lines"), max.overlaps = 1000, fontface = "italic",
  #                          position = position_jitter(width = 0.55, height = 0, seed = 42))+
  ylim(c(0, 1))+
  theme_BP()+
  labs(x = "", y = "Host status probability", fill = "")+
  theme(axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        legend.position = "null")+
  scale_x_discrete(labels=c("0" = "Unknown", "1" = "Positive",
                            "2" = "New positive")) -> plot_tot

require(patchwork)
require(ggpubr)

ggarrange(plot_tot, plot_singles, nrow = 1, labels = "AUTO")
ggsave("Fig_S2.jpeg", width = 11, height = 4, dpi = 300)

predictions_bp %>% 
  filter(predictable_fam == 1) %>%
  mutate(new_host = case_when(
    `Virus genus` == "Ebolavirus & Marburgvirus" & iucn2020_binomial %in% new_ebv_hosts ~ 2,
    `Virus genus` == "Mammarenavirus" & iucn2020_binomial %in% new_mmv_hosts ~ 2,
    `Virus genus` == "Henipavirus" & iucn2020_binomial %in% new_hpv_hosts ~ 2,
    `Virus genus` == "Phlebovirus" & iucn2020_binomial %in% new_phv_hosts ~ 2,
    `Virus genus` == "Betacoronavirus" & iucn2020_binomial %in% new_bcv_hosts ~ 2,
    `Virus genus` == "Orthonairovirus" & iucn2020_binomial %in% new_ohv_hosts ~ 2,
    `Virus genus` == "Flavivirus" & iucn2020_binomial %in% new_flv_hosts ~ 2,
    real_status %in% c("high-evidence host", "low-evidence host") ~ 1,
    TRUE ~ 0
  )) %>% 
  group_by(order) %>% 
  # 1. Create the count and the label
  mutate(n_new = sum(new_host == 2),
         facet_label = paste0(order, " (n = ", n_new, ")")) %>% 
  ungroup() %>% 
  # 2. FILTER: Only keep orders where the count of new hosts is > 0
  filter(n_new > 0) %>% 
  # 3. Rest of your plot code
  filter(`Virus genus` != "Ebolavirus & Marburgvirus") %>% 
  ggplot() +
  gghalves::geom_half_violin(aes(x = factor(new_host), y = mean_pred, fill = factor(new_host)), color = NA) +
  scale_fill_manual(values = c("#364C54FF", "#94475EFF", "#E5A11FFF"), 
                    labels = c("Unknown", "Positive", "New positive")) +
  stat_summary(fun = "mean", mapping = aes(x = factor(new_host), y = mean_pred, group = factor(new_host)), 
               color = "black", geom = "point", size = 2) + 
  ylim(c(0, 1)) +
  theme_BP() +
  facet_wrap(~ facet_label) +
  labs(x = "", y = "Host status probability", fill = "") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        strip.text = element_text(size = 10),
        legend.position = c(0.70, 0.15)
  ) -> plot_single_orders


ggsave(plot = plot_single_orders, filename = "FigureS3.tif", width = 6, height = 6, dpi = 300, units = "in")
