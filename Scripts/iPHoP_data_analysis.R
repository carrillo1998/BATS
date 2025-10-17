######## Viral-Host Abundance ######



######### SURFACE ##############

library(BiocManager)
#BiocManager::install("microbiome")
library(vegan)
library(tidyverse)
library(apcluster)
library(corrplot)


#load OTU table a
meancovsur = read.csv("sur-meancov.txt", sep = "\t", row.names=1)


#rename column names to remove "_Trimmed.Mean"
colnames(meancovsur)
colnames(meancovsur)<-gsub("_R1_clean.fastq.gz.Trimmed.Mean","",colnames(meancovsur))
colnames(meancovsur)<-gsub("final.curated.subsampled.contigs.fa.self.blastn.clusters.fna.","",colnames(meancovsur))
columns_to_remove <- paste0("Contig.", 1:28)
columns_to_keep <- !(names(meancovsur) %in% columns_to_remove)
meancovsur <- meancovsur[, columns_to_keep]
meancovsur <- meancovsur[, c("T0_S", "T4_S", "T8_S", "T12_S", "T16_S", "T20_S", "T24_S", "T28_S", "T32_S",
                             "T36_S", "T40_S", "T44_S", "T48_S", "T52_S", "T56_S", "T60_S", "T64_S", "T68_S",
                             "T72_S", "T76_S", "T80_S", "T84_S", "T88_S", "T92_S", "T96_S", "T100_S", "T104_S", "T108_S", "T112_S")]
cleanedreads_sur <- scan("cleanedreads_sur.txt")
norm3 <- sweep(meancovsur, 2, cleanedreads_sur, FUN = '/')
vOTU_sur<- norm3*1000000000
vOTU_sur_no_outliers <- vOTU_sur[, !(colnames(vOTU_sur) %in% c("T52_S", "T88_S"))]
vOTU_sur <- vOTU_sur_no_outliers


library(dplyr)
combined_host_prediction_to_genome_m90 <- read.csv("~/Documents/BATS R Studio/Outliers removed/Final_Scripts/import/combined_host_prediction_to_genome_m90.csv") #iphop output
iphop_table <- combined_host_prediction_to_genome_m90

#Filter to only have 1 viral entry in dataframe that contains highest confidence.score instead of multiple entries for the same virus 
confidence_filtered_iphop <- iphop_table %>%
  group_by(Virus) %>%
  filter(Confidence.score == max(Confidence.score))

#make two dataframes, one with just genus info and one with just phyla info
iphop_p <- confidence_filtered_iphop %>%
  mutate(taxonomy_p = str_extract(Host.taxonomy, "(?<=p__)[^;]+")) %>%
  select(Virus, taxonomy_p)

iphop_g <- confidence_filtered_iphop %>%
  mutate(taxonomy_g = str_extract(Host.taxonomy, "(?<=g__)[^;]+")) %>%
  select(Virus, taxonomy_g)

vOTU_sur_iphop <- vOTU_sur

###### total host predictions ####
vOTU.tab_no_outlier <- read.csv("~/Documents/BATS R Studio/Outliers removed/Final_Scripts/import/vOTU.tab_no_outlier.csv", row.names=1) #from CoverM_All

vOTU_iphop <- vOTU.tab_no_outlier
# Find matching virus names
matching_names_all <- intersect(rownames(vOTU_iphop), iphop_p$Virus)

#Make df rownames into its own column
vOTU_iphop <- tibble::rownames_to_column(vOTU_iphop, var = "Virus")

vOTU_all_iphop <- left_join(vOTU_iphop, iphop_p, by = "Virus")

# Replace NA values in 'taxonomy_p' with "Unknown Host"
vOTU_all_iphop$taxonomy_p[is.na(vOTU_all_iphop$taxonomy_p)] <- "Unknown Host"

# Remove Unknown Host in 'taxonomy_p' "
vOTU_all_iphop <- vOTU_all_iphop[vOTU_all_iphop$taxonomy_p != "Unknown Host", ]

vOTU_all_iphop_unique <- vOTU_all_iphop[!duplicated(vOTU_all_iphop$Virus),]



#################################### Genus level analysis ##################
# Find matching virus names
matching_names_g <- intersect(rownames(vOTU_sur_iphop), iphop_g$Virus)

#Make df rownames into its own column
vOTU_sur_iphop_g <- vOTU_sur

vOTU_sur_iphop_g <- tibble::rownames_to_column(vOTU_sur_iphop_g, var = "Virus")

vOTU_sur_iphop_g <- left_join(vOTU_sur_iphop_g, iphop_g, by = "Virus")

# Replace NA values in 'taxonomy_p' with "Unknown Host"
vOTU_sur_iphop_g$taxonomy_g[is.na(vOTU_sur_iphop_g$taxonomy_g)] <- "Unknown Host"

# Remove Unknown Host in 'taxonomy_p' "
vOTU_sur_iphop_g <- vOTU_sur_iphop_g[vOTU_sur_iphop_g$taxonomy_g != "Unknown Host", ]

vOTU_sur_iphop_g <- vOTU_sur_iphop_g[!duplicated(vOTU_sur_iphop_g$Virus),]


# Remove the 'virus_name' column
vOTU_sur_iphop_g <- select(vOTU_sur_iphop_g, -Virus)

#Group by 'taxonomy_p' and summarize abundance values
aggregate_iphop_sur_g <- vOTU_sur_iphop_g %>%
  group_by(taxonomy_g) %>%
  summarize_all(sum)

filter_aggregate_sur_g <- aggregate_iphop_sur_g %>%
  filter(rowSums(select(., -taxonomy_g)) != 0)

filter_agg_sur_g <- column_to_rownames(filter_aggregate_sur_g, var = "taxonomy_g")


############################################## DCM ##############################



DCM_meancov <- read.csv("dcm-meancov.txt", sep = "\t", row.names = 1) #new file
#rename column names to remove "_Trimmed.Mean"
colnames(DCM_meancov)
colnames(DCM_meancov)<-gsub("_R1_clean.fastq.gz.Trimmed.Mean","",colnames(DCM_meancov))
colnames(DCM_meancov)<-gsub("vOTUs_5kb.fna.","",colnames(DCM_meancov))
colnames(DCM_meancov)<-gsub("final.curated.subsampled.contigs.fa.self.blastn.clusters.fna.","", colnames(DCM_meancov))
columns_to_remove <- paste0("Contig.", 1:8)
columns_to_keep <- !(names(DCM_meancov) %in% columns_to_remove)
DCM_meancov <- DCM_meancov[, columns_to_keep]
DCM_meancov <- DCM_meancov[, c("T4_D", "T16_D", "T40_D", "T52_D", "T64_D", "T76_D", "T88_D", "T100_D", "T112_D")]

#change these numbers based on the CoverM or FastQC results
cleanedreads_dcm <- scan("cleanedreads_dcm.txt")
norm4 <- sweep(DCM_meancov, 2, cleanedreads_dcm, FUN = '/')
vOTU_dcm<- norm4*1000000000

#remove outliers T4_D
vOTU_dcm_no_outlier <- vOTU_dcm[, !(colnames(vOTU_dcm) %in% c("T4_D"))]
vOTU_dcm <- vOTU_dcm_no_outlier
iphop_table <- combined_host_prediction_to_genome_m90

#Filter to only have 1 viral entry in dataframe that contains highest confidence.score instead of multiple entries for the same virus 
confidence_filtered_iphop <- iphop_table %>%
  group_by(Virus) %>%
  filter(Confidence.score == max(Confidence.score))

#make two dataframes, one with just genus info and one with just phyla info
iphop_p <- confidence_filtered_iphop %>%
  mutate(taxonomy_p = str_extract(Host.taxonomy, "(?<=p__)[^;]+")) %>%
  select(Virus, taxonomy_p)

iphop_g <- confidence_filtered_iphop %>%
  mutate(taxonomy_g = str_extract(Host.taxonomy, "(?<=g__)[^;]+")) %>%
  select(Virus, taxonomy_g)



vOTU_dcm_iphop <- vOTU_dcm


#################################### Genus level analysis ##################
# Find matching virus names
matching_names_g_dcm <- intersect(rownames(vOTU_dcm_iphop), iphop_g$Virus)

#Make df rownames into its own column
vOTU_dcm_iphop_g <- vOTU_dcm

vOTU_dcm_iphop_g <- tibble::rownames_to_column(vOTU_dcm_iphop_g, var = "Virus")

vOTU_dcm_iphop_g <- left_join(vOTU_dcm_iphop_g, iphop_g, by = "Virus")

# Replace NA values in 'taxonomy_g' with "Unknown Host"
vOTU_dcm_iphop_g$taxonomy_g[is.na(vOTU_dcm_iphop_g$taxonomy_g)] <- "Unknown Host"

# Remove Unknown Host in 'taxonomy_g' "
vOTU_dcm_iphop_g <- vOTU_dcm_iphop_g[vOTU_dcm_iphop_g$taxonomy_g != "Unknown Host", ]

vOTU_dcm_iphop_g <- vOTU_dcm_iphop_g[!duplicated(vOTU_dcm_iphop_g$Virus),]


# Remove the 'virus_name' column
vOTU_dcm_iphop_g <- select(vOTU_dcm_iphop_g, -Virus)

#Group by 'taxonomy_g' and summarize abundance values
aggregate_iphop_dcm_g <- vOTU_dcm_iphop_g %>%
  group_by(taxonomy_g) %>%
  summarize_all(sum)

filter_aggregate_dcm_g <- aggregate_iphop_dcm_g %>%
  filter(rowSums(select(., -taxonomy_g)) != 0)

filter_agg_dcm_g <- column_to_rownames(filter_aggregate_dcm_g, var = "taxonomy_g")



###################################################### Diel host predictions #####################################################
#### SURFACE
surface_signif_diel <- read.csv("~/Documents/BATS R Studio/Outliers removed/Final_Scripts/import/surface_signif_diel.csv")

# Find matching virus names
matching_names_diel_g_sur <- intersect(surface_signif_diel$vOTU_names, iphop_g$Virus)

diel_sur_vOTU_g <- vOTU_sur[rownames(vOTU_sur) %in% surface_signif_diel$vOTU_names, ]

#Make df rownames into its own column
diel_sur_iphop_g<- tibble::rownames_to_column(diel_sur_vOTU_g, var = "Virus")

diel_sur_iphop_g <- left_join(diel_sur_iphop_g, iphop_g, by = "Virus")

# Replace NA values in 'taxonomy_p' with "Unknown Host"
diel_sur_iphop_g$taxonomy_g[is.na(diel_sur_iphop_g$taxonomy_g)] <- "Unknown Host"


# Remove Unknown Host in 'taxonomy_p' "
diel_sur_iphop_g <- diel_sur_iphop_g[diel_sur_iphop_g$taxonomy_g != "Unknown Host", ]

diel_sur_iphop_g <- diel_sur_iphop_g[!duplicated(diel_sur_iphop_g$Virus),]

# Remove the 'virus_name' column
diel_sur_iphop_g <- select(diel_sur_iphop_g, -Virus)

#Group by 'taxonomy_p' and summarize abundance values
aggregate_iphop_sur_diel_g <- diel_sur_iphop_g %>%
  group_by(taxonomy_g) %>%
  summarize_all(sum)

filter_aggregate_sur_diel_g <- aggregate_iphop_sur_diel_g %>%
  filter(rowSums(select(., -taxonomy_g)) != 0)

filter_agg_sur_diel_g <- column_to_rownames(filter_aggregate_sur_diel_g, var = "taxonomy_g")


####################### Diel vs. Non-diel


# Remove the 'virus_name' column
diel_sur_iphop <- select(diel_sur_iphop, -Virus)

#Group by 'taxonomy_p' and summarize abundance values
aggregate_iphop_sur_diel <- diel_sur_iphop %>%
  group_by(taxonomy_p) %>%
  summarize_all(sum)

filter_aggregate_sur_diel <- aggregate_iphop_sur_diel %>%
  filter(rowSums(select(., -taxonomy_p)) != 0)

filter_agg_sur_diel <- column_to_rownames(filter_aggregate_sur_diel, var = "taxonomy_p")

col_sum_filter_agg_sur_diel <- rowSums(filter_agg_sur_diel)
sur_diel_overall <- as.data.frame(col_sum_filter_agg_sur_diel)
# Function to calculate proportion for each cell in a column
calculate_proportion <- function(col) {
  col / sum(col) # Divide each cell by the sum of its column excluding the current cell
}

# Apply the calculate_proportion function to each column of the dataframe
proportion_df_sur_diel_overall <- as.data.frame(apply(sur_diel_overall, 2, calculate_proportion))

colSums(proportion_df_sur_diel_overall)


# Reset rownames to a column
proportion_df_sur_diel_overall <- tibble::rownames_to_column(proportion_df_sur_diel_overall, var = "taxonomy_p")




################### Genus SUR Diel vs. Non-diel ########
diel_sur_vOTU_g <- vOTU_sur[rownames(vOTU_sur) %in% surface_signif_diel$vOTU_names, ]

#Make df rownames into its own column
diel_sur_iphop_g <- tibble::rownames_to_column(diel_sur_vOTU_g, var = "Virus")

diel_sur_iphop_g <- left_join(diel_sur_iphop_g, iphop_g, by = "Virus")

# Replace NA values in 'taxonomy_p' with "Unknown Host"
diel_sur_iphop_g$taxonomy_g[is.na(diel_sur_iphop_g$taxonomy_g)] <- "Unknown Host"


# Remove Unknown Host in 'taxonomy_p' "
diel_sur_iphop_g <- diel_sur_iphop_g[diel_sur_iphop_g$taxonomy_g != "Unknown Host", ]


# Remove the 'virus_name' column
diel_sur_iphop_g <- select(diel_sur_iphop_g, -Virus)

#Group by 'taxonomy_g' and summarize abundance values
aggregate_iphop_sur_diel_g <- diel_sur_iphop_g %>%
  group_by(taxonomy_g) %>%
  summarize_all(sum)

filter_aggregate_sur_diel_g <- aggregate_iphop_sur_diel_g %>%
  filter(rowSums(select(., -taxonomy_g)) != 0)

filter_agg_sur_diel_g <- column_to_rownames(filter_aggregate_sur_diel_g, var = "taxonomy_g")




diel_sur_virus_g <- tibble::rownames_to_column(diel_sur_vOTU_g, var = "Virus")


########### remake vOTU_sur_iphop_g with virus column

vOTU_sur_iphop <- vOTU_sur

#Make df rownames into its own column
vOTU_sur_iphop <- tibble::rownames_to_column(vOTU_sur_iphop, var = "Virus")

vOTU_sur_iphop_g <- left_join(vOTU_sur_iphop, iphop_g, by = "Virus")

# Remove duplicates
vOTU_sur_iphop_g <- vOTU_sur_iphop_g[!duplicated(vOTU_sur_iphop_g$Virus),]


# Replace NA values in 'taxonomy_p' with "Unknown Host"
vOTU_sur_iphop_g$taxonomy_g[is.na(vOTU_sur_iphop_g$taxonomy_g)] <- "Unknown Host"

# Remove Unknown Host in 'taxonomy_p' "
vOTU_sur_iphop_g <- vOTU_sur_iphop_g[vOTU_sur_iphop_g$taxonomy_g != "Unknown Host", ]

# Remove duplicates
vOTU_sur_iphop_g <- vOTU_sur_iphop_g[!duplicated(vOTU_sur_iphop_g$Virus),]


vOTU_sur_iphop_non_diel_g <- anti_join(vOTU_sur_iphop_g, diel_sur_virus_g, by = "Virus")

# Remove duplicates
vOTU_sur_iphop_non_diel_g <- vOTU_sur_iphop_non_diel_g[!duplicated(vOTU_sur_iphop_non_diel_g$Virus),]


# Replace NA values in 'taxonomy_p' with "Unknown Host"
vOTU_sur_iphop_non_diel_g$taxonomy_g[is.na(vOTU_sur_iphop_non_diel_g$taxonomy_g)] <- "Unknown Host"

# Remove Unknown Host in 'taxonomy_p' "
vOTU_sur_iphop_non_diel_g <- vOTU_sur_iphop_non_diel_g[vOTU_sur_iphop_non_diel_g$taxonomy_g != "Unknown Host", ]
vOTU_sur_iphop_non_diel_g_non0 <- vOTU_sur_iphop_non_diel_g %>%
  filter(rowSums(across(where(is.numeric))) != 0)

# Remove the 'virus_name' column
vOTU_sur_iphop_non_diel_g <- select(vOTU_sur_iphop_non_diel_g, -Virus)

#Group by 'taxonomy_p' and summarize abundance values
aggregate_iphop_sur_non_diel_g <- vOTU_sur_iphop_non_diel_g %>%
  group_by(taxonomy_g) %>%
  summarize_all(sum)

filter_aggregate_sur_non_diel_g <- aggregate_iphop_sur_non_diel_g %>%
  filter(rowSums(select(., -taxonomy_g)) != 0)

filter_agg_sur_non_diel_g <- column_to_rownames(filter_aggregate_sur_non_diel_g, var = "taxonomy_g")










