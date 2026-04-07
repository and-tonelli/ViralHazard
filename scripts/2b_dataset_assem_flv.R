require(tidyverse)

# These are known high-evidence hosts (See Dataset S1):
flv_reservoirs <- c("Apodemus flavicollis", "Clethrionomys rutilus", "Apodemus agrarius", 
                    "Bison bonasus", "Cercopithecus ascanius", "Myodes glareolus", "Myodes rutilus", "Sorex araneus", "Sorex betulina")

# Checking if all of them are in COMBINE
combine <- read.csv(".../data/combine_imputed.csv")

unique(combine %>% 
         filter(phylacine_binomial %in% flv_reservoirs) %>% 
         pull(phylacine_binomial)) # all of them are in COMBINE

# Aquatic species
combine %>% filter(freshwater == "1" | marine == "1") %>% pull(phylacine_binomial) -> aquatic

# calling 'not in' function
'%!in%' <- function(x, y)!('%in%'(x, y)) # %!in%

virion_m_clean <- read_csv(".../data/virion_m_clean.csv") %>% 
  filter(Host %!in% aquatic,
         #DetectionMethod != "Not specified"
  ) %>% 
  mutate(Host = str_to_sentence(Host))

# How many species are associated to flv just via Serology?
virion_m_clean %>% 
  filter(VirusGenus %in% c("flavivirus") & Host %!in% flv_reservoirs, DetectionMethod == "Antibodies") %>% 
  pull(Host) -> antibodies

virion_m_clean %>% 
  filter(VirusGenus %in% c("flavivirus") & Host %!in% flv_reservoirs, DetectionMethod != "Antibodies") %>% 
  pull(Host) -> species

associated_just_via_antibodies <- setdiff(antibodies, species) # 13 species

# Getting low-evidence hosts (i.e., potentially unrecognized reservoirs)
# These are the species that tested positive to the target viral genus but are not high-evidence hosts (according to WOS search)
# chances are that they have not been sufficiently studied, therefore they might be unrecognised reservoirs 
# 1) they are never pseudoabsences, 2) they enter into the analysis as 1s, but with lower instance weigths

# Also removing viruses whose reservoirs are birds (See Suppl. Table 1)
pseudo_reservoir_pool <- c(unique(virion_m_clean %>% 
                                    filter(VirusGenus %in% c("flavivirus") & Virus %!in% c("ilheus virus", "israel turkey meningoencephalomyelitis virus", "japanese encephalitis virus", "saint louis encephalitis virus", "tembusu virus", "usutu virus", "west nile virus") & Host %!in% flv_reservoirs) %>%
                                    pull(Host)))

# These are the p-n species
pseudonegative_spp <- unique(str_replace(tab_overlaps %>% pull(sp), "_", " "))

# Loading predictors (Life-history, phylogeny, bioclimatic variables)
phyclimbine <- read.csv(".../data/all_predictors.csv")
dist_known_hosts <- read_csv(".../data/flvmean_phylo_dist_to_known_hosts.csv")

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
highevidence_set <- phyclimbine %>% filter(phylacine_binomial %in% flv_reservoirs) %>% 
  merge.data.frame(., dist_known_hosts, by.x = 1, by.y = 1) %>% 
  mutate(status = 1) %>% 
  dplyr::select(c(phylacine_binomial, iucn2020_binomial, order, family, 
                  adult_mass_g, max_longevity_d, gestation_length_d, litter_size_n, litters_per_year_n, 
                  interbirth_interval_d, weaning_age_d, trophic_level, area.1:mean_dist)) %>% 
  filter(complete.cases(.))

# Getting species-level research effort (based on reports recorded in VIRION)
species_level_re <- virion_m_clean %>% 
  filter(VirusGenus %in% c("flavivirus")) %>% 
  group_by(Host) %>%
  summarise(records = n())

# Getting family-level bcov research effort (based on reports recorded in VIRION)
family_level_re <- virion_m_clean %>% 
  filter(VirusGenus %in% c("flavivirus")) %>% 
  group_by(HostFamily) %>%
  summarise(records = n())

family_level_re$HostFamily <- str_to_sentence(family_level_re$HostFamily)

# Compute weights for high-evidence hosts (see Tonelli et al., 2025; doi.org/10.1111/2041-210X.14500)
highevidence_set_w <- highevidence_set %>%
  mutate(w = 1/length(highevidence_set$iucn2020_binomial),
         status = 1,
         real_status = "high-evidence host")

# Compute weights for susceptible species (see Tonelli et al., 2025; doi.org/10.1111/2041-210X.14500)
flv_lowevidence_set <- lowevidence_set %>%
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
flv_dataset <- rbind(highevidence_set_w[-35] %>% relocate(w, .after = last_col()),
                     flv_lowevidence_set,
                     pseudonegative_set %>% relocate(phylacine_binomial, iucn2020_binomial, order) %>% relocate(w, .after = last_col()))

# checking that sum(w of positive species) = sum(w of pseudo-negatives)
sum(flv_dataset$w[flv_dataset$status == "0"])
sum(flv_dataset$w[flv_dataset$status == "1"])

write.csv(flv_dataset, ".../data/flv_dataset_main.csv", row.names = F)