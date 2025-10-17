
library(dplyr)
library(BiocManager)
#BiocManager::install("microbiome")
library(vegan)
library(tidyverse)
library(apcluster)
library(corrplot)
library(stringr)
library(latticeExtra)
library(pracma)
library(rain)

dcm_pallete <- c("coral", "green","cadetblue1","burlywood","brown3","blue","aquamarine2","deeppink","darkslategrey",
                 "darksalmon", "darkseagreen4","magenta2","darkgoldenrod","gray90","gold","lightpink","lightcyan","khaki",
                 "lavender", "indianred","midnightblue","mediumspringgreen","mediumorchid","lightseagreen","olivedrab2",
                 "peachpuff","orangered1","orange","yellow","lightsalmon","purple","plum","snow4","steelblue","springgreen4",
                 "tan4","tomato2","violetred2")

############ Detrend and benjamini-hochberg FDR test on RAIN #############

####################################################### Surface Data generation ######

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
colnames(meancovsur)

cleanedreads_sur <- scan("cleanedreads_sur.txt")
norm3 <- sweep(meancovsur, 2, cleanedreads_sur, FUN = '/')
vOTU_sur<- norm3*1000000000

vOTU_sur_no_outliers <- vOTU_sur[, !(colnames(vOTU_sur) %in% c("T52_S", "T88_S"))]

vOTU_sur_detrend <- as.data.frame(t(apply(vOTU_sur_no_outliers, 1, detrend)))

vOTU_sur_colnames <- c("T0_S", "T4_S", "T8_S", "T12_S", "T16_S", "T20_S", "T24_S", "T28_S", "T32_S",
                       "T36_S", "T40_S", "T44_S", "T48_S", "T56_S", "T60_S", "T64_S", "T68_S",
                       "T72_S", "T76_S", "T80_S", "T84_S", "T92_S", "T96_S", "T100_S", "T104_S", "T108_S", "T112_S")
colnames(vOTU_sur_detrend) <- vOTU_sur_colnames

############################################################ DCM Data Generation #################################################

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
colnames(DCM_meancov)

#change these numbers based on the CoverM or FastQC results
cleanedreads_dcm <- scan("cleanedreads_dcm.txt")
norm4 <- sweep(DCM_meancov, 2, cleanedreads_dcm, FUN = '/')
vOTU_dcm<- norm4*1000000000



#remove outliers T4_D
vOTU_dcm_no_outlier <- vOTU_dcm[, !(colnames(vOTU_dcm) %in% c("T4_D"))]

vOTU_dcm <- vOTU_dcm_no_outlier
vOTU_dcm_detrend <- vOTU_dcm
vOTU_dcm_detrend <- as.data.frame(t(apply(vOTU_dcm, 1, detrend)))

vOTU_dcm_colnames <- c("T16_D", "T40_D", "T52_D", "T64_D", "T76_D", "T88_D", "T100_D", "T112_D")

colnames(vOTU_dcm_detrend) <- vOTU_dcm_colnames
####################### SUR Analysis ##################

### Running RAIN for rhythmic detection (Diel signals)
# measure.sequence used when samples are missing, in this case 0 indicates a removed sample (outliers)
# this will take a while
sur_results <- rain(t(vOTU_sur_detrend), deltat = 4, period = 24, peak.border = c(0.3, 0.7), measure.sequence = c(1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1),
                    verbose = FALSE); sur_results
sur_results$p.adj <- p.adjust(sur_results$pVal, method = "BH")


sur_best <- order(sur_results$p.adj)[1:10]; sur_best
times1 <- c(0:56) *4; times1

vOTU_sur[sur_best,]
vOTU_sur[sur_best, (0:27)]


custom_x_axis_4hour_intervals <- c(0:28 *4); custom_x_axis_4hour_intervals
custom_x_axis_12hour_intervals <- c(0:10 *12); custom_x_axis_12hour_intervals
custom_x_axis_24hour_intervals <- c(0:9 * 24); custom_x_axis_24hour_intervals


##### Surface  #####
SUR_best_data <- vOTU_sur[sur_best, (1:27)]


SUR_best_data$T52_S <- NA
SUR_best_data$T88_S <- NA
SUR_order <- c("T0_S", "T4_S", "T8_S", "T12_S", "T16_S", "T20_S", "T24_S", "T28_S", "T32_S", 
               "T36_S", "T40_S", "T44_S", "T48_S", "T52_S", "T56_S", "T60_S", "T64_S", "T68_S", 
               "T72_S", "T76_S", "T80_S", "T84_S", "T88_S", "T92_S", "T96_S", "T100_S", "T104_S", 
               "T108_S", "T112_S")


# Reshape your data into long format using pivot_longer
long_SUR_best_data <- SUR_best_data %>%
  rownames_to_column("Virus") %>%  # Add the row names (viruses) as a new column
  pivot_longer(cols = starts_with("T"),  # Selecting columns starting with "T" for time points
               names_to = "Time",  # New column for time points
               values_to = "Abundance")  # New column for abundance data


long_SUR_best_plot <- long_SUR_best_data
long_SUR_best_plot$Time <- as.numeric(sub("T([0-9]+)_S", "\\1", long_SUR_best_plot$Time))
long_SUR_best_plot$Time <- as.numeric(as.character(long_SUR_best_plot$Time))

long_SUR_best_plot$Time <- factor(long_SUR_best_plot$Time, levels = unique(long_SUR_best_plot$Time))
long_SUR_best_plot$Time <- as.numeric(as.character(long_SUR_best_plot$Time))

# Define a function to generate military time labels
military_time_labels <- function(x) {
  # Convert timepoints (e.g., 0, 4, 8, ...) to military time starting at 1600
  start_time <- 1600
  military_time <- (start_time + x * 100) %% 2400
  sprintf("%04d", military_time) # Ensure labels are zero-padded to 4 digits
}

ggplot(long_SUR_best_plot, aes(x = Time, y = Abundance, group = Virus, color = Virus)) +
  geom_line() +  # Create line graph
  geom_point() +  # Add points to the graph
  theme_minimal() +  # Minimal theme for the plot
  labs(title = "Abundance of Viruses Over Time", x = "Time Point", y = "Abundance") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels
  facet_wrap(~ Virus, scales = "free_y", ncol = 2)+ # Create 10 separate plots, 2 columns
  annotate("rect", xmin = c(2,26,50,74,98), xmax = c(14,38,62,86,110), ymin = 0.5, ymax = 80, #bars to show night time, doesn't work for multiple plots as height differs
           alpha = .5) +
  scale_x_continuous(breaks = seq(0, 112, by = 4), labels = military_time_labels, expand = c(0,0))




######################## DCM Attempt #################################################


dcm_time <- c(0,1,0,1,1,1,1,1,1,1)
dcm_results <- rain(t(vOTU_dcm_detrend), deltat = 12, period = 24, peak.border = c(0.3, 0.7), measure.sequence = c(0, 1, 0, 1, 1, 1, 1, 1, 1, 1), 
                    verbose = FALSE); dcm_results

dcm_results$p.adj <- p.adjust(dcm_results$pVal, method = "BH")
dcm_best <- order(dcm_results$pVal)[1:10]; dcm_best
DCM_best_data <- vOTU_dcm[dcm_best, (1:8)]

DCM_best_data$T4_D <- NA
DCM_best_data$T28_D <- NA
DCM_order <- c("T4_D", "T16_D", "T28_D", "T40_D", "T52_D", "T64_D", "T76_D", "T88_D", "T100_D", "T112_D")
DCM_best_data <- DCM_best_data[, DCM_order]


# Reshape your data into long format using pivot_longer
long_DCM_best_data <- DCM_best_data %>%
  rownames_to_column("Virus") %>%  # Add the row names (viruses) as a new column
  pivot_longer(cols = starts_with("T"),  # Selecting columns starting with "T" for time points
               names_to = "Time",  # New column for time points
               values_to = "Abundance")  # New column for abundance data


long_DCM_best_plot <- long_DCM_best_data
long_DCM_best_plot$Time <- as.numeric(sub("T([0-9]+)_D", "\\1", long_DCM_best_plot$Time))
long_DCM_best_plot$Time <- as.numeric(as.character(long_DCM_best_plot$Time))

long_DCM_best_plot$Time <- factor(long_DCM_best_plot$Time, levels = unique(long_DCM_best_plot$Time))
long_DCM_best_plot$Time <- as.numeric(as.character(long_DCM_best_plot$Time))

# Define a function to generate military time labels
military_time_labels_DCM <- function(x) {
  # Convert timepoints (e.g., 0, 4, 8, ...) to military time starting at 1600
  start_time <- 2000
  military_time <- (start_time + x * 100) %% 2400
  sprintf("%04d", military_time) # Ensure labels are zero-padded to 4 digits
}

ggplot(long_DCM_best_plot, aes(x = Time, y = Abundance, group = Virus, color = Virus)) +
  geom_line() +  # Create line graph
  geom_point() +  # Add points to the graph
  theme_minimal() +  # Minimal theme for the plot
  labs(title = "Abundance of Viruses Over Time", x = "Time Point", y = "Abundance") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels
  facet_wrap(~ Virus, scales = "free_y", ncol = 2)+ # Create 10 separate plots, 2 columns
  annotate("rect", xmin = c(2,26,50,74,98), xmax = c(14,38,62,86,110), ymin = 0.5, ymax = 80, 
           alpha = .5) +
  scale_x_continuous(breaks = seq(4, 112, by = 12), labels = military_time_labels, expand = c(0,0))





################### Separating sinusoidal based on day or night dominant ################
# Gather everything significant
surface_signif_diel <- subset(sur_results, p.adj <= 0.05);surface_signif_diel
dcm_signif_diel <- subset(dcm_results, p.adj <= 0.05); dcm_signif_diel



######################################### Day Dominant #####
# The way we are doing this is based on the phase information from the RAIN results. Phase tells us where the vOTU peaks each 24 time period.
# So if we tell it to give us things that peak during day hours we will get the day dominant diel cycling ones and if we do it night we should get the same. 
# For reference again, RAIN does phase based on its period information. For Surface we said the period is 24 hours at 4 hour intervals. Meaning each phase will
#   go up by 4. The first point will be 4 and the last point will be 24. Meaning although T0 is our first point and this is at 4pm, phase 4(the first phase option)
#   will be the phase for this. The next time point which is T4 will be the second in the model which will make it phase 8 and so on. 

Day_phase_numbers_sur <- c(4,20,24)

#### Surface

Day_dominant_diel_surface <- surface_signif_diel[surface_signif_diel$phase %in% c(4, 20, 24),]; Day_dominant_diel_surface


######################################### Night Dominant #####

Night_phase_numbers_sur <- c(8,12,16)

#### Surface

Night_dominant_diel_surface <- surface_signif_diel[surface_signif_diel$phase %in% c(8, 12, 16),]; Night_dominant_diel_surface

################### Plotting day or night dominant ##################

############################## Day
# Sort the DataFrame in order based on the column of interest
sur_day <- Day_dominant_diel_surface[order(Day_dominant_diel_surface$pVal), ]; sur_day


############################ Night

# Sort the DataFrame in order based on the column of interest
sur_night <- Night_dominant_diel_surface[order(Night_dominant_diel_surface$pVal), ]; sur_night


########################## Diel Names for vConTACT analysis #################

###### Surface Significant Diel
# All vOTUs Surface
surface_signif_diel <- surface_signif_diel %>%
  rownames_to_column(var = "vOTU_names"); surface_signif_diel
All_Surface_Diel_vOTUs <- surface_signif_diel$vOTU_names


# Day vOTUs
sur_day <- sur_day %>%
  rownames_to_column(var = "vOTU_names"); sur_day
Day_Surface_Diel_vOTUs <- sur_day$vOTU_names

# Night vOTUS
sur_night <- sur_night %>%
  rownames_to_column(var = "vOTU_names"); sur_night
Night_Surface_Diel_vOTUs <- sur_night$vOTU_names


###### Save and export ##########

# Make them into data_frames
All_Surface_vOTUs <- data.frame(All_Surface_Diel_vOTUs = All_Surface_Diel_vOTUs); All_Surface_Diel_vOTUs
Day_Surface_vOTUs <- data.frame(Day_Surface_Diel_vOTUs = Day_Surface_Diel_vOTUs)
Night_Surface_vOTUs <- data.frame(Night_Surface_vOTUs = Night_Surface_Diel_vOTUs)

# Save and export file
write.csv(All_Surface_vOTUs, file = "All_Surface_vOTUs.csv", row.names = FALSE)
write.csv(Day_Surface_vOTUs, file = "Day_Surface_vOTUs.csv", row.names = FALSE)
write.csv(Night_Surface_vOTUs, file = "Night_Surface_vOTUs.csv", row.names = FALSE)



######### SUR Diel Barplot #########
SUR_bar <- data.frame(Group = c("Diel", "Non-Diel"),
                      Viruses = c(3097, 25907))


SUR_bar <- SUR_bar %>%
  mutate(proportion = Viruses / sum(Viruses))


pdf(file = "Fig6a_1.pdf", width = 3, height = 5)
ggplot(SUR_bar, aes(x = "", y = proportion, fill = Group)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Proportional Stacked Barplot of DCM and SUR",
    x = "Depth",
    y = "Proportion of Diel Viral Population Dominance"
  ) +
  scale_fill_manual(values = c("Diel" = "maroon", "Non-Diel" = "black")) +
  theme_minimal()
dev.off()




######## Supplemental table 13
diel_vOTUs <- vOTU_sur[rownames(vOTU_sur) %in% surface_signif_diel$vOTU_names,]
write.csv(diel_vOTUs, file = "Supp_table_13.csv", row.names = TRUE)





