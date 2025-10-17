
######################## GOV2 Dataset Ecological zones inverse simpson added to our dataset

#load OTU table a
#Get from https://doi.org/10.1016/j.cell.2019.03.040
GOV2.tab = read.csv("Normalized_Viral_Abundances_ALL_95ANI_5kb_ECOLOGY.txt", sep = "\t", row.names=1) #GOV2 data
filtered_GOV2 <- GOV2.tab %>%
  select(ends_with("_SUR") | ends_with("_DCM"))

names(filtered_GOV2) <- gsub("^Station([0-9]+)_", "\\1_", names(filtered_GOV2))

ecological_zones <- read.csv("~/Documents/BATS R Studio/Outliers removed/Final_Scripts/import/PCoA1_ecological_zones copy.csv")

#From Fig1_script.R
inv_simp_whole_tib_All <- read.csv("~/Documents/BATS R Studio/Outliers removed/Final_Scripts/import/inv_simp_whole_tib_All.csv", row.names=1)
inv_simp_whole_tib_All$depth <- gsub("DCM", "BATS_DCM", inv_simp_whole_tib_All$depth)
inv_simp_whole_tib_All$depth <- gsub("SUR", "BATS_SUR", inv_simp_whole_tib_All$depth)



############### SUR - DCM inv. simpson ###############
library(vegan)
inv_simp_gov2 <- diversity(t(filtered_GOV2), index = "invsimpson")

inv_simp_df <- data.frame(
  sample = names(inv_simp_gov2),
  inv_simpson = inv_simp_gov2,
  stringsAsFactors = FALSE
)


# Extract sample ID (e.g., "11", "13") and layer (SUR or DCM)
inv_df <- inv_simp_df %>%
  mutate(
    id = sub("_(SUR|DCM)", "", sample),
    layer = sub(".*_(SUR|DCM)", "\\1", sample)
  )


wide_df <- inv_df %>%
  select(id, layer, inv_simpson) %>%
  pivot_wider(names_from = layer, values_from = inv_simpson)

wide_df <- wide_df %>%
  mutate(diff = SUR - DCM)

wide_df_filter <- wide_df %>%
  mutate(diff = SUR - DCM) %>%
  filter(!is.na(diff))


ecological_zone <- ecological_zones %>%
  mutate(Zone = gsub("_[^_]+$", "", Ecological_zone))

# First, extract the numeric ID from the ecological zone df
ecozone_df <- ecological_zone %>%
  mutate(ID = gsub("_.*", "", Sample))  # Gets the number before the underscore


# Join ecological zone information
wide_df_filter <- wide_df_filter %>%
  left_join(ecozone_df %>% select(ID, Zone) %>% distinct(), by = c("id" = "ID"))




inv_simp_whole_tib_All_gov2 <- inv_simp_whole_tib_All
inv_simp_whole_tib_All_gov2
inv_simp_whole_tib_All_gov2$depth_gov <- gsub("SUR", "BATS_SUR",
                                              gsub("DCM", "BATS_DCM",
                                                   inv_simp_whole_tib_All_gov2$depth))


# First, extract the shared site ID before the underscore from the `site` column
inv_simp_all_bats <- inv_simp_whole_tib_All_gov2 %>%
  mutate(site_id = sub("_[DS]$", "", site))

# Pivot wider so you get one row per site_id, with columns for SUR and DCM values
inv_simp_all_bats_wide <- inv_simp_all_bats %>%
  select(site_id, depth, inv_simp) %>%
  pivot_wider(names_from = depth, values_from = inv_simp)

# Now compute the difference: SUR - DCM
inv_simp_all_bats_wide <- inv_simp_all_bats_wide %>%
  mutate(diff = BATS_SUR - BATS_DCM) %>%
  filter(!is.na(diff))  # Remove rows with missing values

inv_simp_all_bats_wide$Zone <- "BATS"

colnames(inv_simp_all_bats_wide) <- c("id", "DCM", "SUR", "diff", "Zone")

gov2_invsimp_comparison <- rbind(inv_simp_all_bats_wide, wide_df_filter)
gov2_invsimp_comparison$abs <- abs(gov2_invsimp_comparison$diff)

pdf(file = "Supp_Fig2.pdf", width = 5, height = 5)
ggplot(gov2_invsimp_comparison, aes(x = "", y = diff)) +
  geom_boxplot(fill = "gray90", outlier.shape = NA) +
  geom_jitter(aes(color = Zone), width = 0.2, size = 4, alpha = 0.8) +
  facet_wrap(~Zone, ncol = 4) +
  theme_minimal() +
  labs(
    x = NULL,
    y = "Inverse Simpson Difference (SUR - DCM)",
    title = "Inverse Simpson Differences Across Ecological Zones"
  ) +
  theme(
    legend.title = element_blank()
  )
dev.off()
