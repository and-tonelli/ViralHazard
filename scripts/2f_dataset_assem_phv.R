require(tidyverse)

# These are known high-evidence hosts (please refer to Supplementary Dataset 1):
phv_reservoirs <- c("Arvicanthis niloticus", "Mastomys erythroleucus", "Aethomys namaquensis", "Micaelamys namaquensis", "Erinaceus amurensis", "Atelerix albiventris")

# Checking if all of them are in COMBINE
combine <- read.csv(".../data/combine_imputed.csv")

# Aquatic species
combine %>% filter(freshwater == "1" | marine == "1") %>% pull(phylacine_binomial) -> aquatic

# calling 'not in' function
'%!in%' <- function(x, y)!('%in%'(x, y)) # %!in%

virion_m_clean <- read_csv(".../data/virion_m_clean.csv") %>% 
  filter(Host %!in% aquatic,
         #DetectionMethod != "Not specified"
  ) %>% 
  mutate(Host = str_to_sentence(Host))

# How many species are associated to phv just via Serology?
virion_m_clean %>% 
  filter(VirusGenus %in% c("phlebovirus") & Host %!in% phv_reservoirs, DetectionMethod == "Antibodies") %>% 
  pull(Host) -> antibodies

virion_m_clean %>% 
  filter(VirusGenus %in% c("phlebovirus") & Host %!in% phv_reservoirs, DetectionMethod != "Antibodies") %>% 
  pull(Host) -> species

associated_just_via_antibodies <- setdiff(antibodies, species) # 16 species

# Getting low-evidence hosts (i.e., potentially unrecognized reservoirs)
# These are the species that tested positive but are not reservoirs (according to WOS search)
# chances are that they have not been sufficiently studied, therefore they might be unrecognised reservoirs 
# What to do with them? 1) they are never pseudoabsences, 2) they enter into the analysis as 1s, but with lower weights

virion_m_clean$Host[virion_m_clean$Host == "Micaelamys namaquensis"] <- "Aethomys namaquensis"


pseudo_reservoir_pool <- c(unique(virion_m_clean %>% 
                                    filter(VirusGenus %in% c("phlebovirus") & Virus %!in% c("hughes phlebovirus", "sakhalin phlebovirus") & Host %!in% phv_reservoirs) %>%
                                    pull(Host)))
# These are the p-n species
pseudonegative_spp <- unique(str_replace(tab_overlaps %>% pull(sp), "_", " "))

# Loading predictors (Life-history, phylogeny, bioclimatic variables)
phyclimbine <- read.csv(".../data/all_predictors.csv")
dist_known_hosts <- read_csv(".../data/phvmean_phylo_dist_to_known_hosts.csv")

# Assembling the pseudonegative set
pseudonegative_set <- phyclimbine %>% filter(iucn2020_binomial %in% pseudonegative_spp) %>% 
  merge.data.frame(., dist_known_hosts, by.x = 1, by.y = 1) %>% 
  mutate(status = 0) %>%
  dplyr::select(c(phylacine_binomial, iucn2020_binomial, order, family, 
                  adult_mass_g, max_longevity_d, gestation_length_d, litter_size_n, litters_per_year_n, 
                  interbirth_interval_d, weaning_age_d, trophic_level, area.1:mean_dist)) %>% 
  filter(complete.cases(.))

length(pseudonegative_set$iucn2020_binomial)

# Assembling the low-evidence set (getting the relevant variables)
lowevidence_set <- phyclimbine %>% filter(iucn2020_binomial %in% pseudo_reservoir_pool) %>% 
  merge.data.frame(., dist_known_hosts, by.x = 1, by.y = 1) %>% 
  mutate(status = 1) %>%
  dplyr::select(c(phylacine_binomial, iucn2020_binomial, order, family, 
                  adult_mass_g, max_longevity_d, gestation_length_d, litter_size_n, litters_per_year_n, 
                  interbirth_interval_d, weaning_age_d, trophic_level, area.1:mean_dist)) %>% 
  filter(complete.cases(.))

# Assembling the high-evidence set
highevidence_set <- phyclimbine %>% filter(phylacine_binomial %in% phv_reservoirs) %>% 
  merge.data.frame(., dist_known_hosts, by.x = 1, by.y = 1) %>% 
  mutate(status = 1) %>% 
  dplyr::select(c(phylacine_binomial, iucn2020_binomial, order, family, 
                  adult_mass_g, max_longevity_d, gestation_length_d, litter_size_n, litters_per_year_n, 
                  interbirth_interval_d, weaning_age_d, trophic_level, area.1:mean_dist)) %>% 
  filter(complete.cases(.))

# Getting species-level research effort (based on reports recorded in VIRION)
species_level_re <- virion_m_clean %>% 
  filter(VirusGenus %in% c("phlebovirus")) %>% 
  group_by(Host) %>%
  summarise(records = n())

# Getting family-level bcov research effort (based on reports recorded in VIRION)
family_level_re <- virion_m_clean %>% 
  filter(VirusGenus %in% c("phlebovirus")) %>% 
  group_by(HostFamily) %>%
  summarise(records = n())

family_level_re$HostFamily <- str_to_sentence(family_level_re$HostFamily)

# Compute weights for high-evidence hosts (see Tonelli et al., 2025; doi.org/10.1111/2041-210X.14500)
highevidence_set_w <- highevidence_set %>%
  mutate(w = 1/length(highevidence_set$iucn2020_binomial),
         status = 1,
         real_status = "high-evidence host")

# Compute weights for susceptible species (see Tonelli et al., 2025; doi.org/10.1111/2041-210X.14500)
phv_lowevidence_set <- lowevidence_set %>%
  merge(., species_level_re, by.x = 2, by.y = 1) %>% 
  mutate(status = 1,
         real_status = "low-evidence host",
         w_pre = 1/log10(1+records),
         w = w_pre/sum(w_pre)) %>% 
  dplyr::select(-c(w_pre, records))

# Compute weights for pseudo-negatives (see Tonelli et al., 2025; doi.org/10.1111/2041-210X.14500)
pseudonegative_set <- pseudonegative_set %>% 
  mutate(w = 2/length(pseudonegative_set),
         real_status = "pseudo-negative")

# Putting everything together in the final dataset
phv_dataset <- rbind(highevidence_set_w[-35] %>% relocate(w, .after = last_col()),
                     phv_lowevidence_set,
                     pseudonegative_set %>% relocate(phylacine_binomial, iucn2020_binomial, order) %>% relocate(w, .after = last_col()))

# checking that sum(w of positive species) = sum(w of pseudo-negatives)
sum(phv_dataset$w[phv_dataset$status == "0"])
sum(phv_dataset$w[phv_dataset$status == "1"])

write.csv(phv_dataset, ".../data/phv_dataset_main.csv", row.names = F)