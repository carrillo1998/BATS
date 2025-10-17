


library(dplyr)
library(BiocManager)
#BiocManager::install("microbiome")
library(vegan)
library(tidyverse)
library(apcluster)
library(corrplot)
library(stringr)

####################### AMG Work #######################

all.annotations <- read.delim("~/Documents/BATS R Studio/Outliers removed/Final_Scripts/import/all-annotations.tsv") #annotate file
AMG_data <- all.annotations

amg_summary <- read.delim("~/Documents/BATS R Studio/Outliers removed/Final_Scripts/import/amg_summary.tsv") #distill summary file
AMG_annotations <- amg_summary


###### Generate permissive AMG catalog ######

permissive_catalog <- all.annotations %>%
  filter(auxiliary_score >= 1 & auxiliary_score <= 3,
         str_detect(amg_flags, "[MKE]"))

all.annotations
permissive_catalog


#### checkV viruses ####
combined.checkV.contam <- read.delim2("~/Documents/BATS R Studio/Outliers removed/Final_Scripts/import/combined-checkV-contam.tsv") #checkV file

viral_checkV <- combined.checkV.contam %>%
  filter(provirus == "Yes")

viral_checkV_text_edit <- viral_checkV %>%
  mutate(contig_id = str_extract(contig_id, "^[^|]+"))

permissive_text_edit <- permissive_catalog %>%
  mutate(scaffold = str_replace(scaffold, "_full.*|_partial.*", "")) %>%
  mutate(scaffold = str_replace_all(scaffold, "__", "_")) %>%
  mutate(scaffold = str_trim(scaffold)) 

permissive_text_edit <- permissive_text_edit %>%
  mutate(scaffold = str_remove(scaffold, "_$"))

permissive_checkV_filtered <- permissive_text_edit %>%
  filter(scaffold %in% viral_checkV_text_edit$contig_id)

permissive_checkV_filtered

all.annotations_text_edit <- all.annotations
all.annotations_text_edit$scaffold <- gsub("__full.*", "||full", all.annotations_text_edit$scaffold)
all.annotations_text_edit$scaffold <- gsub("__0_partial.*", "||0_partial", all.annotations_text_edit$scaffold)
all.annotations_text_edit$scaffold <- gsub("__k141", "_k141", all.annotations_text_edit$scaffold)


################## Permissive to conservative #################

############## checkV and VS2 ###########
combined.checkV.contam <- read.delim2("~/Documents/BATS R Studio/Outliers removed/Final_Scripts/import/combined-checkV-contam.tsv") #checkV contam file with data before subsample removal
final.viral.score <- read.delim("~/Documents/BATS R Studio/Outliers removed/Final_Scripts/import/final-viral-score.tsv") #final viral score after everything has been replaced with correct data



permissive_catalog #make it in same orientation as other files. No double __ and remove stuff after partial and full so the -cat_1 stuff in the scaffold column
# change _full to ||full


permissive_text_edit <- permissive_catalog %>%
  mutate(
    scaffold = str_replace(scaffold, "__", "_"),
    scaffold = str_remove(scaffold, "-cat.*"),
    scaffold = str_replace_all(scaffold, "__(full|partial)", "||\\1")
  )

permissive_text_edit$scaffold <- gsub("__0_partial.*", "||0_partial", permissive_text_edit$scaffold)
permissive_text_edit$scaffold <- gsub("full.*", "full", permissive_text_edit$scaffold)

#remove __ from final.viral.score
final.viral.score_edit <- final.viral.score %>%
  mutate(
    seqname = str_replace(seqname, "__", "_")
  )

############## filter using Virsorter2 data #############

final.viral.score_edit # VS2 output for viral score

VS2_output <- final.viral.score_edit %>%
  mutate(seqname = str_replace(seqname, "_full.*", "")) %>%
  mutate(seqname = str_replace_all(seqname, "__", "_")) %>%
  mutate(seqname = str_trim(seqname)) 

VS2_output <- VS2_output %>%
  mutate(seqname = str_replace(seqname, "_partial.*", "_partial")) %>%
  mutate(seqname = str_trim(seqname)) 

VS2_output$seqname <- gsub("full.*", "full", VS2_output$seqname)

VS2_filtered <- permissive_text_edit %>%
  semi_join(VS2_output %>% filter(max_score >= 0.95), 
            by = c("scaffold" = "seqname"))


######## VS2 and checkV passed AMGs ######
permissive_checkV_filtered$scaffold <- permissive_checkV_filtered$X

permissive_checkV_filtered <- permissive_checkV_filtered %>%
  mutate(
    scaffold = str_replace(scaffold, "__", "_"),
    scaffold = str_remove(scaffold, "-cat.*"),
    scaffold = str_replace_all(scaffold, "__(full|partial)", "||\\1")
  )


permissive_AMGs_step1 <- bind_rows(permissive_checkV_filtered, VS2_filtered) %>%
  distinct(X, X, .keep_all = TRUE) 

osc_permissive <- permissive_AMGs_step1
osc_permissive$names <- sub("-cat.*", "", osc_permissive$X)
osc_permissive$names <- gsub("__full", "||full", osc_permissive$names)
osc_permissive$names <- gsub("__0", "||0", osc_permissive$names)

osc_permissive_unique <- unique(osc_permissive$names)

write_lines(osc_permissive_unique, "AMG_scaffolds.txt")


###### tRNA-Scan filter ######

tRNA_scan_data <- read.delim("~/Documents/BATS R Studio/Outliers removed/Final_Scripts/import/tRNA_full_data (1).txt")
tRNA_scan_data_text_edit <- tRNA_scan_data %>%
  mutate(Sequence_Name = str_replace_all(Sequence_Name, "__", "_")) %>%
  mutate(Sequence_Name = str_trim(Sequence_Name)) 

tRNA_scan_data_text_edit$Sequence_Name <- gsub("full.*", "full", tRNA_scan_data_text_edit$Sequence_Name)
tRNA_scan_data_text_edit$Sequence_Name <- gsub("partial.*", "partial", tRNA_scan_data_text_edit$Sequence_Name)

permissive_AMGs_step1_tRNA_passed <- tRNA_scan_data_text_edit %>% #no tRNA
  filter(!Sequence_Name %in% permissive_AMGs_step1$scaffold)

tRNA_scan_matches <- tRNA_scan_data_text_edit %>% #tRNA hits
  filter(Sequence_Name %in% permissive_AMGs_step1$scaffold)

tRNA_scan_data_manual_check <- tRNA_scan_matches %>% #tRNA hit that was score >55 requires manual check
  filter(Cove_Score > 55)

tRNA_scan_data_manual_check_exempt <- tRNA_scan_matches %>% #tRNA hit was below 55 and doesn't require manual check 
  filter(!Cove_Score > 55)

add_to_permissive_tRNA <- tRNA_scan_data_text_edit %>% #AMGs to add back in
  filter(Sequence_Name %in% tRNA_scan_data_manual_check_exempt$Sequence_Name)




######## EMBOSS and tRNA check ######
conservative_emboss_summary <- read.csv("~/Documents/BATS R Studio/Outliers removed/Final_Scripts/import/conservative_emboss_summary.csv")
EMBOSS_output <- conservative_emboss_summary

EMBOSS_output$scaffold <- gsub("__k", "_k", EMBOSS_output$scaffold)
EMBOSS_output$scaffold <- gsub("full_", "full", EMBOSS_output$scaffold)

check2_match <- union(EMBOSS_output$scaffold, tRNA_scan_data_text_edit$Sequence_Name)

curation_pass2 <- permissive_AMGs_step1 %>% 
  filter(!scaffold %in% check2_match)
unique(curation_pass2$scaffold)
unique(curation_pass2$X)

#AMGs that have tRNA or repeat in scaffold
curation_check2 <- permissive_AMGs_step1 %>% 
  filter(scaffold %in% check2_match)
unique(curation_check2$scaffold)
unique(curation_check2$X)


############ tRNA and emboss for all AMGs, not just the ends ########
#AMGs that have tRNA or repeat in scaffold


# Step 1: Find positions for 1 gene left and 2 genes right
all.annotations_text_edit_search_range <- all.annotations_text_edit %>%
  arrange(scaffold, gene_position) %>%
  group_by(scaffold) %>%
  mutate(
    AMG_left1 = lag(gene_position, 1),   # 1 gene to the left
    AMG_left1_start = lag(start_position, 1),
    AMG_left1_end = lag(end_position, 1),
    
    AMG_left2 = lag(gene_position, 2),   # 1 gene to the left
    AMG_left2_start = lag(start_position, 2),
    AMG_left2_end = lag(end_position, 2),
    
    AMG_right2 = lead(gene_position, 2),  # 2 genes to the right
    AMG_right2_start = lead(start_position, 2),
    AMG_right2_end = lead(end_position, 2),
    
    AMG_right3 = lead(gene_position, 3),  # 2 genes to the right
    AMG_right3_start = lead(start_position, 3),
    AMG_right3_end = lead(end_position, 3)
  ) %>%
  ungroup()


# Step 2: Merge with curation to get target positions
curation_check2_extended <- curation_check2 %>%
  left_join(all.annotations_text_edit_search_range %>%
              select(scaffold, gene_position, AMG_left1_start, AMG_left1_end, AMG_left2_start, AMG_left2_end,
                     AMG_right2_start, AMG_right2_end, AMG_right3_start, AMG_right3_end),
            by = c("scaffold", "gene_position"))
colnames(curation_check2_extended)
tRNA_emboss_curation_all <- curation_check2_extended %>%
  select(X, fasta, scaffold, gene_position, start_position, end_position, auxiliary_score, AMG_left1_start, AMG_left1_end,
         AMG_left2_start, AMG_left2_end, AMG_right2_start, AMG_right2_end, AMG_right3_start, AMG_right3_end)
tRNA_emboss_curation_all$X 

tRNA_scan_data_text_edit <- tRNA_scan_data_text_edit %>%
  rename(scaffold = Sequence_Name)
EMBOSS_output <- EMBOSS_output %>%
  separate(positions, into = c("emboss_start", "emboss_end"), sep = "-")

# Step 3: Check for overlaps or adjacency within the range
tRNA_emboss_done_all <- curation_check2_extended %>%
  left_join(tRNA_scan_data_text_edit, by = "scaffold") %>%
  left_join(EMBOSS_output, by = "scaffold") %>%
  mutate(
    # tRNA or emboss hit overlaps within the left-to-right boundary
    tRNA_within = (!is.na(tRNA_Begin) & tRNA_Begin >= AMG_left1_start & tRNA_End <= AMG_right2_end),
    emboss_within = (!is.na(emboss_start) & emboss_start >= AMG_left1_start & emboss_end <= AMG_right2_end),
    
    # tRNA or emboss hit is adjacent or overlapping the left boundary
    tRNA_left_adjacent_or_overlap = (!is.na(tRNA_End) & (tRNA_End >= AMG_left2_start & tRNA_End <= AMG_left1_end & tRNA_Begin <= AMG_left1_end)),
    emboss_left_adjacent_or_overlap = (!is.na(emboss_end) & (emboss_end >= AMG_left2_start & emboss_end <= AMG_left1_end & emboss_start <= AMG_left1_end)),
    
    # tRNA or emboss hit is adjacent or overlapping the right boundary
    tRNA_right_adjacent_or_overlap = (!is.na(tRNA_Begin) & (tRNA_Begin <= AMG_right3_end & tRNA_Begin >= AMG_right2_start & tRNA_End >= AMG_right2_start)),
    emboss_right_adjacent_or_overlap = (!is.na(emboss_start) & (emboss_start <= AMG_right3_end & emboss_start >= AMG_right2_start & emboss_end >= AMG_right2_start)),
    
    # Final check flag
    check_flag = tRNA_within | emboss_within | tRNA_left_adjacent_or_overlap | emboss_left_adjacent_or_overlap | 
      tRNA_right_adjacent_or_overlap | emboss_right_adjacent_or_overlap
  )

colnames(tRNA_emboss_done_all)
tRNA_emboss_done_all <- tRNA_emboss_done_all %>%
  select(X, fasta, scaffold, gene_position, start_position, end_position, auxiliary_score, AMG_left1_start, AMG_left1_end,
         AMG_left2_start, AMG_left2_end, AMG_right2_start, AMG_right2_end, AMG_right3_start, AMG_right3_end, tRNA_Begin, 
         tRNA_End, emboss_start, emboss_end, tRNA_within, emboss_within, tRNA_left_adjacent_or_overlap,
         tRNA_right_adjacent_or_overlap, emboss_left_adjacent_or_overlap, emboss_right_adjacent_or_overlap, check_flag)

unique(tRNA_emboss_done_all$X)
# Step 4: Split into curation4_check and curation4_pass
curation_check2_2 <- tRNA_emboss_done_all %>% filter(check_flag == TRUE) %>% select(names(tRNA_emboss_done_all)) #546
curation_pass2_2 <- tRNA_emboss_done_all %>% filter(check_flag == FALSE) %>% select(names(tRNA_emboss_done_all)) #4280


unique(curation_check2_2$X)
unique(curation_pass2_2$X) 

intersect(curation_check2_2$X, curation_pass2_2$X)


###### Curation of those that have tRNA or emboss repeat near AMG #####
curation_failed <- curation_check2_2 %>% 
  filter(tRNA_within == TRUE | emboss_within == TRUE)
unique(curation_failed$X)

curation_check2_3 <- curation_check2_2 %>% 
  filter(tRNA_within == FALSE & emboss_within == FALSE)
unique(curation_check2_3$X) 


# Move rows where either tRNA_within or emboss_within is TRUE (or NA) to curation_failed
curation_failed <- curation_check2_2 %>% 
  filter(replace_na(tRNA_within, FALSE) == TRUE | replace_na(emboss_within, FALSE) == TRUE)
unique(curation_failed$X) 

# Move rows where both tRNA_within and emboss_within are FALSE (and not NA) to curation_check2_3
curation_check2_3 <- curation_check2_2 %>% 
  filter(replace_na(tRNA_within, FALSE) == FALSE & replace_na(emboss_within, FALSE) == FALSE)
unique(curation_check2_3$X) 


df1_unique <- curation_check2_2 %>%
  anti_join(curation_failed, by = colnames(curation_check2_2)) %>%  # Remove rows in df2
  anti_join(curation_check2_3, by = colnames(curation_check2_2)) 



###### tRNA and emboss length check #####
# If you have an AMG that is on the end of the contig, they will disappear in this step as they will not have AMG_right2. 
# This is not an issue because for it to be here means it has a tRNA or repeat region and it must have a repeat in the right to pass. 

curation_length_pass <- curation_check2_3 %>% 
  filter((AMG_right2_end - AMG_left1_start) <= 3000)
unique(curation_length_pass$X) 

curation_length_fail <- curation_check2_3 %>% 
  filter((AMG_right2_end - AMG_left1_start) > 3000)
unique(curation_length_fail$X) 


curation_length_pass <- curation_check2_3 %>% 
  filter(!is.na(AMG_left1_start) & !is.na(AMG_right2_end) & (AMG_right2_end - AMG_left1_start) <= 3000)
unique(curation_length_pass$X) 

curation_length_fail <- curation_check2_3 %>% 
  filter(is.na(AMG_left1_start) | is.na(AMG_right2_end) | (AMG_right2_end - AMG_left1_start) > 3000)
unique(curation_length_fail$X) 

df1_unique2 <- curation_check2_3 %>%
  anti_join(curation_length_pass, by = colnames(curation_check2_3)) %>%  # Remove rows in df2
  anti_join(curation_length_fail, by = colnames(curation_check2_3)) 


#### tRNA and repeat on rightside?
tRNA_and_repeat_on_right_pass <- curation_length_pass %>% 
  filter(tRNA_right_adjacent_or_overlap == TRUE & emboss_right_adjacent_or_overlap == TRUE)
unique(tRNA_and_repeat_on_right_pass$X)

tRNA_and_repeat_on_right_fail <- curation_length_pass %>% 
  filter(!(tRNA_right_adjacent_or_overlap == TRUE & emboss_right_adjacent_or_overlap == TRUE))
unique(tRNA_and_repeat_on_right_fail$X) 


# Identify the gene names (X) that have at least one row with both TRUE
AMGs_to_pass <- curation_length_pass %>%
  group_by(X) %>%
  filter(any(tRNA_right_adjacent_or_overlap == TRUE) & any(emboss_right_adjacent_or_overlap == TRUE)) %>%
  pull(X) %>%
  unique()

# Filter into pass/fail dataframes
position_pass <- curation_length_pass %>%
  filter(X %in% AMGs_to_pass)

position_fail <- curation_length_pass %>% 
  filter(!X %in% AMGs_to_pass)
unique(position_fail$X) 



############## Filter using output from steps ################

AMGs_to_remove <- unique(c(curation_failed$X, curation_length_fail$X, position_fail$X)) 

conservative_AMGs_step2 <- permissive_AMGs_step1 %>% 
  filter(!X %in% AMGs_to_remove)
unique(conservative_AMGs_step2$scaffold) 

gene_keywords_for_removal <- c("glycosyltransferase", "glycosyl transferase", "nucleotidyl transferase",
                               "carbohydrate kinase", "nucleotide sugar epimerase", "endonuclease",
                               "integrase", "plasmid stability")

# Create a filtered version of all.annotations that only contains relevant scaffolds
all_annotations_filtered_for_final_curation <- all.annotations_text_edit %>%
  filter(scaffold %in% conservative_AMGs_step2$scaffold)


# Identify scaffolds that contain any of the gene_keywords_for_removal in the specified columns
scaffolds_to_exclude <- all_annotations_filtered_for_final_curation %>%
  filter(str_detect(tolower(kegg_hit), str_c(gene_keywords_for_removal, collapse = "|")) |
           str_detect(tolower(pfam_hits), str_c(gene_keywords_for_removal, collapse = "|")) |
           str_detect(tolower(viral_hit), str_c(gene_keywords_for_removal, collapse = "|"))) %>%
  pull(scaffold) # Get the scaffolds that match the condition


# Separate into pass and fail based on scaffold matches
conservative_curation_fail_step3 <- conservative_AMGs_step2 %>% 
  filter(scaffold %in% scaffolds_to_exclude)

conservative_AMGs_final <- conservative_AMGs_step2 %>% 
  filter(!scaffold %in% scaffolds_to_exclude)
unique(conservative_AMGs_final$X)


conservative_AMGs_final_flag_check <- conservative_AMGs_final %>%
  filter(str_detect(amg_flags, "[MKE]"))
conservative_AMGs_final_flag_check <- conservative_AMGs_final_flag_check %>%
  filter(grepl("^[MKE]+$", amg_flags))

conservative_AMGs_final_flag_check_filter <- conservative_AMGs_final_flag_check %>%
  filter(!is.na(ko_id) & ko_id != "")

write.csv(conservative_AMGs_final_flag_check_filter, file = "Supp_table_6.csv", row.names = TRUE)






conservative_AMGs_final_flag_check_filter #annotations of all AMGs
vOTU.tab_no_outlier <- read.csv("~/Documents/BATS R Studio/Outliers removed/Final_Scripts/import/vOTU.tab_no_outlier.csv", row.names=1) #from CoverM_All



### Making table with AMG, vOTU, group, and abundance
amg_info_table <- select(conservative_AMGs_final_flag_check_filter, X, ko_id, scaffold)
amg_info_table$X <- gsub("__full", "||full", amg_info_table$X)
amg_info_table$X <- gsub("-cat.*", "", amg_info_table$X)
amg_info_table$X <- gsub("__0_partial", "||0_partial", amg_info_table$X)

amg_info_table <- select(amg_info_table, X, ko_id)
colnames(amg_info_table) <- c("vOTU", "ko_id")

vOTU.tab_no_outlier_edit <- rownames_to_column(vOTU.tab_no_outlier, var = "vOTU") #pulling abundance data

amg_info_table_abund <- semi_join(vOTU.tab_no_outlier_edit, amg_info_table, by = "vOTU")
amg_info_table_abund <- left_join(amg_info_table, vOTU.tab_no_outlier_edit, by = "vOTU")



######## Defining groups 
vOTU.tab_no_outlier

# Identify DCM and SUR columns
dcm_cols <- grep("_D$", names(vOTU.tab_no_outlier), value = TRUE)
sur_cols <- grep("_S$", names(vOTU.tab_no_outlier), value = TRUE)

# Calculate row sums
dcm_sum <- rowSums(vOTU.tab_no_outlier[, dcm_cols])
sur_sum <- rowSums(vOTU.tab_no_outlier[, sur_cols])

vOTU.tab_no_outlier_amg <- vOTU.tab_no_outlier #df to make groups
# Assign depth values
vOTU.tab_no_outlier_amg$depth <- ifelse(dcm_sum == 0 & sur_sum > 0, "SUR",
                    ifelse(sur_sum == 0 & dcm_sum > 0, "DCM",
                           ifelse(sur_sum > 0 & dcm_sum > 0, "both", NA)))
vOTU.tab_no_outlier_amg <- rownames_to_column(vOTU.tab_no_outlier_amg, var = "vOTU")

amg_info_table_abund <- left_join(amg_info_table_abund, vOTU.tab_no_outlier_amg[, c("vOTU", "depth")], by = "vOTU") #depth category added


#### adding diel or non-diel
surface_signif_diel <- read.csv("~/Documents/BATS R Studio/Outliers removed/Final_Scripts/import/surface_signif_diel.csv", row.names=1)
surface_signif_diel <- rownames_to_column(surface_signif_diel, var = "vOTU_names")
diel_vOTUs <- select(surface_signif_diel, vOTU_names)
diel_vOTUs$diel <- "Diel"
colnames(diel_vOTUs) <- c("vOTU", "diel")

amg_info_table_abund <- left_join(amg_info_table_abund, diel_vOTUs[, c("vOTU", "diel")], by = "vOTU") #depth category added

amg_info_table_abund$diel[is.na(amg_info_table_abund$diel)] <- "non-diel"


##### adding archetype from Fig6_script.R
SUR_som3_arche1_signif_subset <- read.csv("~/Documents/BATS R Studio/Outliers removed/Final_Scripts/import/SUR_som3_arche1_signif_subset.csv", row.names=1)
SUR_som3_arche2_signif_subset <- read.csv("~/Documents/BATS R Studio/Outliers removed/Final_Scripts/import/SUR_som3_arche2_signif_subset.csv", row.names=1)
SUR_som3_arche3_signif_subset <- read.csv("~/Documents/BATS R Studio/Outliers removed/Final_Scripts/import/SUR_som3_arche3_signif_subset.csv", row.names=1)

arche1_votu <- select(SUR_som3_arche1_signif_subset, virus)
arche1_votu$arche <- "Arche1"
colnames(arche1_votu) <- c("vOTU", "arche")

arche2_votu <- select(SUR_som3_arche2_signif_subset, virus)
arche2_votu$arche <- "Arche2"
colnames(arche2_votu) <- c("vOTU", "arche")

arche3_votu <- select(SUR_som3_arche3_signif_subset, virus)
arche3_votu$arche <- "Arche3"
colnames(arche3_votu) <- c("vOTU", "arche")

archetype_votus <- rbind(arche1_votu, arche2_votu, arche3_votu)

amg_info_table_abund <- left_join(amg_info_table_abund, archetype_votus[, c("vOTU", "arche")], by = "vOTU") #depth category added

amg_info_table_abund <- amg_info_table_abund[, c("vOTU", "ko_id", "depth", "diel", "arche", setdiff(names(amg_info_table_abund), c("vOTU", "ko_id", "depth", "diel", "arche")))]


amg_info_table_abund$diel[amg_info_table_abund$depth == "DCM"] <- NA

write.csv(amg_info_table_abund, file = "BATS_amg_abund.csv", row.names = FALSE)
write.csv(conservative_AMGs_final_flag_check_filter, file = "BATS_amg_annotation.csv", row.names = FALSE)


amg_arche_info <- select(amg_info_table_abund, vOTU, ko_id, depth, diel, arche)

arche_amg <- amg_arche_info[!is.na(amg_arche_info$arche), ]
arche_amg <- select(arche_amg, vOTU, ko_id, arche)


#### adding hosts from iPHoP_data_analysis.R script
iphop_g <- read.csv("~/Documents/BATS R Studio/Outliers removed/Final_Scripts/import/iphop_g.csv", row.names=1)

iphop_g_edit <- iphop_g
colnames(iphop_g_edit) <- c("vOTU", "taxonomy_g")

arche_amg <- left_join(arche_amg, iphop_g_edit[, c("vOTU", "taxonomy_g")], by = "vOTU") #depth category added



# Step 2: Count vOTUs per archetype per taxonomy_g and ko_id
counts <- arche_amg %>%
  group_by(taxonomy_g = taxonomy_g, ko_id, arche) %>%
  summarise(vOTU_count = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = arche,
    values_from = vOTU_count,
    names_prefix = "vOTUs_"
  )

# Step 3: Replace NAs with 0s for counts
counts <- counts %>%
  mutate(across(-c(taxonomy_g, ko_id), ~ replace_na(., 0)))

# Step 4: Determine Group classification
counts <- counts %>%
  rowwise() %>%
  mutate(Group = case_when(
    vOTUs_Arche1 > 0 & vOTUs_Arche2 == 0 & vOTUs_Arche3 == 0 ~ "Unique_Arche1",
    vOTUs_Arche2 > 0 & vOTUs_Arche1 == 0 & vOTUs_Arche3 == 0 ~ "Unique_Arche2",
    vOTUs_Arche3 > 0 & vOTUs_Arche1 == 0 & vOTUs_Arche2 == 0 ~ "Unique_Arche3",
    vOTUs_Arche1 > 0 & vOTUs_Arche2 > 0 & vOTUs_Arche3 > 0 ~ "shared_all",
    vOTUs_Arche1 > 0 & vOTUs_Arche2 > 0 & vOTUs_Arche3 == 0 ~ "shared_Arche1_Arche2",
    vOTUs_Arche1 > 0 & vOTUs_Arche3 > 0 & vOTUs_Arche2 == 0 ~ "shared_Arche1_Arche3",
    vOTUs_Arche2 > 0 & vOTUs_Arche3 > 0 & vOTUs_Arche1 == 0 ~ "shared_Arche2_Arche3",
    TRUE ~ "Other"
  )) %>%
  ungroup()

# Step 5: Rename and reorder columns to final output format
final_df <- counts %>%
  rename(
    host = taxonomy_g,
    AMG = ko_id
  ) %>%
  select(host, vOTUs_Arche1, vOTUs_Arche2, vOTUs_Arche3, AMG, Group)



arche_diel_amg <- amg_arche_info[!is.na(amg_arche_info$diel), ]
arche_diel_amg <- select(arche_diel_amg, vOTU, ko_id, diel, arche)

iphop_g_edit <- iphop_g
colnames(iphop_g_edit) <- c("vOTU", "taxonomy_g")

arche_diel_amg <- left_join(arche_diel_amg, iphop_g_edit[, c("vOTU", "taxonomy_g")], by = "vOTU") #depth category added

arche_diel_amg <- arche_diel_amg %>%
  distinct(vOTU, ko_id, diel, arche, taxonomy_g, .keep_all = TRUE)


# Step 2: Count vOTUs per archetype per taxonomy_g and ko_id
counts <- arche_diel_amg %>%
  group_by(taxonomy_g = taxonomy_g, ko_id, arche) %>%
  summarise(vOTU_count = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = arche,
    values_from = vOTU_count,
    names_prefix = "vOTUs_"
  )
colnames(counts) <- c("taxonomy_g", "ko_id", "vOTUs_non_diel", "vOTUs_Arche1", "vOTUs_Arche3", "vOTUs_Arche2")


# Step 3: Replace NAs with 0s for counts
counts <- counts %>%
  mutate(across(-c(taxonomy_g, ko_id), ~ replace_na(., 0)))

# Step 4: Determine Group classification
counts <- counts %>%
  rowwise() %>%
  mutate(Group = case_when(
    # Unique to Arche types
    vOTUs_Arche1 > 0 & vOTUs_Arche2 == 0 & vOTUs_Arche3 == 0 & vOTUs_non_diel == 0 ~ "Unique_Arche1",
    vOTUs_Arche2 > 0 & vOTUs_Arche1 == 0 & vOTUs_Arche3 == 0 & vOTUs_non_diel == 0 ~ "Unique_Arche2",
    vOTUs_Arche3 > 0 & vOTUs_Arche1 == 0 & vOTUs_Arche2 == 0 & vOTUs_non_diel == 0 ~ "Unique_Arche3",
    
    # Unique to non-diel
    vOTUs_non_diel > 0 & vOTUs_Arche1 == 0 & vOTUs_Arche2 == 0 & vOTUs_Arche3 == 0 ~ "Unique_non_diel",
    
    # Shared across archetypes
    vOTUs_Arche1 > 0 & vOTUs_Arche2 > 0 & vOTUs_Arche3 > 0 & vOTUs_non_diel == 0 ~ "shared_all",
    vOTUs_Arche1 > 0 & vOTUs_Arche2 > 0 & vOTUs_Arche3 == 0 & vOTUs_non_diel == 0 ~ "shared_Arche1_Arche2",
    vOTUs_Arche1 > 0 & vOTUs_Arche3 > 0 & vOTUs_Arche2 == 0 & vOTUs_non_diel == 0 ~ "shared_Arche1_Arche3",
    vOTUs_Arche2 > 0 & vOTUs_Arche3 > 0 & vOTUs_Arche1 == 0 & vOTUs_non_diel == 0 ~ "shared_Arche2_Arche3",
    
    # Shared with non-diel
    vOTUs_Arche1 > 0 & vOTUs_non_diel > 0 & vOTUs_Arche2 == 0 & vOTUs_Arche3 == 0 ~ "shared_Arche1_non_diel",
    vOTUs_Arche2 > 0 & vOTUs_non_diel > 0 & vOTUs_Arche1 == 0 & vOTUs_Arche3 == 0 ~ "shared_Arche2_non_diel",
    vOTUs_Arche3 > 0 & vOTUs_non_diel > 0 & vOTUs_Arche1 == 0 & vOTUs_Arche2 == 0 ~ "shared_Arche3_non_diel",
    vOTUs_Arche1 > 0 & vOTUs_Arche2 > 0 & vOTUs_non_diel > 0 & vOTUs_Arche3 == 0 ~ "shared_Arche1_Arche2_non_diel",
    vOTUs_Arche1 > 0 & vOTUs_Arche3 > 0 & vOTUs_non_diel > 0 & vOTUs_Arche2 == 0 ~ "shared_Arche1_Arche3_non_diel",
    vOTUs_Arche2 > 0 & vOTUs_Arche3 > 0 & vOTUs_non_diel > 0 & vOTUs_Arche1 == 0 ~ "shared_Arche2_Arche3_non_diel",
    vOTUs_Arche1 > 0 & vOTUs_Arche2 > 0 & vOTUs_Arche3 > 0 & vOTUs_non_diel > 0 ~ "shared_all_with_non_diel",
    
    TRUE ~ "Other"
  )) %>%
  ungroup()


# Step 5: Rename and reorder columns to final output format
final_df <- counts %>%
  rename(
    host = taxonomy_g,
    AMG = ko_id
  ) %>%
  select(host, vOTUs_non_diel, vOTUs_Arche1, vOTUs_Arche2, vOTUs_Arche3, AMG, Group)

write.csv(final_df, file = "AMG_with_host_and_arche.csv", row.names = FALSE)



depth_amg <- select(amg_arche_info, vOTU, ko_id, depth)

iphop_g_edit <- iphop_g
colnames(iphop_g_edit) <- c("vOTU", "taxonomy_g")

depth_amg <- left_join(depth_amg, iphop_g_edit[, c("vOTU", "taxonomy_g")], by = "vOTU") #depth category added

depth_amg <- depth_amg %>%
  distinct(vOTU, ko_id, depth, taxonomy_g, .keep_all = TRUE)


# Step 2: Count vOTUs per archetype per taxonomy_g and ko_id
depth_counts <- depth_amg %>%
  group_by(taxonomy_g = taxonomy_g, ko_id, depth) %>%
  summarise(vOTU_count = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = depth,
    values_from = vOTU_count,
    names_prefix = "vOTUs_"
  )


# Step 3: Replace NAs with 0s for counts
depth_counts <- depth_counts %>%
  mutate(across(-c(taxonomy_g, ko_id), ~ replace_na(., 0)))

# Step 4: Determine Group classification
depth_counts <- depth_counts %>%
  rowwise() %>%
  mutate(Group = case_when(
    vOTUs_SUR > 0 & vOTUs_DCM == 0 & vOTUs_both == 0 ~ "unique_SUR",
    vOTUs_DCM > 0 & vOTUs_SUR == 0 & vOTUs_both == 0 ~ "unique_DCM",
    vOTUs_SUR > 0 & vOTUs_DCM > 0 ~ "both",
    vOTUs_both > 0 ~ "both",
    TRUE ~ "Other"
  )) %>%
  ungroup()


# Step 5: Rename and reorder columns to final output format
depth_amg_df <- depth_counts %>%
  rename(
    host = taxonomy_g,
    AMG = ko_id
  ) %>%
  select(host, vOTUs_DCM, vOTUs_SUR, vOTUs_both, AMG, Group)

write.csv(depth_amg_df, file = "AMG_with_host_for_depth.csv", row.names = FALSE)









################ Making Supplemental Figures 11 and 12 ####################


install.packages("KEGGREST")
library(KEGGREST)
library(dplyr)
library(stringr)
library(purrr)

final_df
final_df_edit <- final_df %>%
  rename(ko_id = AMG)

final_df_edit_ann <- final_df_edit %>%
  mutate(
    ko_info = purrr::map(ko_id, function(id) {
      tryCatch({
        KEGGREST::keggGet(paste0("ko:", id))[[1]]
      }, error = function(e) NULL)
    }),
    ko_name = purrr::map_chr(ko_info, ~ if (!is.null(.x)) str_remove_all(.x$NAME, "\\s\\[.*\\]") else NA_character_),
    ko_symbol = purrr::map_chr(ko_info, ~ {
      if (!is.null(.x) && !is.null(.x$SYMBOL)) {
        str_split(.x$SYMBOL, " ", simplify = TRUE)[1]
      } else {
        NA_character_
      }
    }),
    pathway1 = purrr::map_chr(ko_info, ~ {
      if (!is.null(.x) && !is.null(.x$BRITE)) {
        brite_lines <- .x$BRITE
        if (length(brite_lines) >= 2) {
          pathway <- gsub("^[0-9]+\\s+", "", brite_lines[2])
          return(str_trim(pathway))
        }
      }
      return(NA_character_)
    }),
    pathway2 = purrr::map_chr(ko_info, ~ {
      if (!is.null(.x) && !is.null(.x$BRITE)) {
        brite_lines <- .x$BRITE
        if (length(brite_lines) >= 3) {
          pathway <- gsub("^[0-9]+\\s+", "", brite_lines[3])
          return(str_trim(pathway))
        }
      }
      return(NA_character_)
    }),
    pathway3 = purrr::map_chr(ko_info, ~ {
      if (!is.null(.x) && !is.null(.x$BRITE)) {
        brite_lines <- .x$BRITE
        if (length(brite_lines) >= 4) {
          pathway <- gsub("^[0-9]+\\s+", "", brite_lines[4])
          return(str_trim(pathway))
        }
      }
      return(NA_character_)
    })
  ) %>%
  select(-ko_info)


final_df_edit_ann

# remove digits from the pathway names
final_df_ann = final_df_edit_ann %>%
  mutate(across(c(pathway1, pathway2, pathway3), ~ str_remove(.x, "\\d+ ")))


write.csv(final_df_ann, file = "Supp_table12.csv", row.names = FALSE)









####### Depth 
depth_amg_df

depth_amg_file_edit <- depth_amg_df %>%
  rename(ko_id = AMG)

depth_amg_file_edit_ann <- depth_amg_file_edit %>%
  mutate(
    ko_info = purrr::map(ko_id, function(id) {
      tryCatch({
        KEGGREST::keggGet(paste0("ko:", id))[[1]]
      }, error = function(e) NULL)
    }),
    ko_name = purrr::map_chr(ko_info, ~ if (!is.null(.x)) str_remove_all(.x$NAME, "\\s\\[.*\\]") else NA_character_),
    ko_symbol = purrr::map_chr(ko_info, ~ {
      if (!is.null(.x) && !is.null(.x$SYMBOL)) {
        str_split(.x$SYMBOL, " ", simplify = TRUE)[1]
      } else {
        NA_character_
      }
    }),
    pathway1 = purrr::map_chr(ko_info, ~ {
      if (!is.null(.x) && !is.null(.x$BRITE)) {
        brite_lines <- .x$BRITE
        if (length(brite_lines) >= 2) {
          pathway <- gsub("^[0-9]+\\s+", "", brite_lines[2])
          return(str_trim(pathway))
        }
      }
      return(NA_character_)
    }),
    pathway2 = purrr::map_chr(ko_info, ~ {
      if (!is.null(.x) && !is.null(.x$BRITE)) {
        brite_lines <- .x$BRITE
        if (length(brite_lines) >= 3) {
          pathway <- gsub("^[0-9]+\\s+", "", brite_lines[3])
          return(str_trim(pathway))
        }
      }
      return(NA_character_)
    }),
    pathway3 = purrr::map_chr(ko_info, ~ {
      if (!is.null(.x) && !is.null(.x$BRITE)) {
        brite_lines <- .x$BRITE
        if (length(brite_lines) >= 4) {
          pathway <- gsub("^[0-9]+\\s+", "", brite_lines[4])
          return(str_trim(pathway))
        }
      }
      return(NA_character_)
    })
  ) %>%
  select(-ko_info)

depth_amg_file_edit_ann

# remove digits from the pathway names
depth_amg_file_edit_ann = depth_amg_file_edit_ann %>%
  mutate(across(c(pathway1, pathway2, pathway3), ~ str_remove(.x, "\\d+ ")))

write.csv(depth_amg_file_edit_ann, file = "Supp_table11.csv", row.names = FALSE)



