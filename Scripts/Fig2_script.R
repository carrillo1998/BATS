################ CoverM analysis ##################

#install and load all of our packages
#install.packages("tidyverse", dependencies = TRUE )
#install.packages('vegan')
#install.packages('pheatmap')
#install.packages('apcluster')
#install.packages('ggplot2')
#install.packages('corrplot')
library(BiocManager)
#BiocManager::install("microbiome")
library(vegan)
library(tidyverse)
library(apcluster)
library(corrplot)
library(ggplot2)
library(ggrepel)
library(sf)
library(maps)



############ Figure 2 scripts ############
#load OTU table a
OTU.tab = read.csv("all-meancov.txt", sep = "\t", row.names=1)

#rename column names to remove "_Trimmed.Mean"
colnames(OTU.tab)
colnames(OTU.tab)<-gsub("_R1_clean.fastq.gz.Trimmed.Mean.1","",colnames(OTU.tab))
colnames(OTU.tab)<-gsub("final.curated.subsampled.contigs.fa.self.blastn.clusters.fna.","",colnames(OTU.tab))
columns_to_remove <- paste0("Contig.", 1:75)
columns_to_keep <- !(names(OTU.tab) %in% columns_to_remove)
OTU.tab <- OTU.tab[, columns_to_keep]
OTU.tab <- OTU.tab[, c("T4_D", "T16_D", "T40_D", "T52_D", "T64_D", "T76_D", "T88_D", "T100_D", "T112_D", "T0_S", 
                       "T4_S", "T8_S", "T12_S", "T16_S", "T20_S", "T24_S", "T28_S", "T32_S",
                       "T36_S", "T40_S", "T44_S", "T48_S", "T52_S", "T56_S", "T60_S", "T64_S", "T68_S",
                       "T72_S", "T76_S", "T80_S", "T84_S", "T88_S", "T92_S", "T96_S", "T100_S", "T104_S", "T108_S", "T112_S")]

colnames(OTU.tab)


#### Normalization
cleanedreads_all <- scan("cleaned_reads.txt")
otu.tab2 <- sweep(OTU.tab, 2, cleanedreads_all, FUN = '/')
vOTU.tab <- otu.tab2 *1000000000

#remove outliers
vOTU.tab_no_outlier <- vOTU.tab[, !(colnames(vOTU.tab) %in% c("T4_D", "T52_S", "T88_S"))]
write.csv(vOTU.tab_no_outlier, file = "vOTU.tab_no_outlier.csv", row.names = TRUE)

#format abundance data
norm_abund_t_All <- t(vOTU.tab_no_outlier) %>%
  as_tibble(norm_abund_t_All, rownames= NA)
#log transformation
log_All <- log(norm_abund_t_All +1)


#bray curtis
bray_All <- vegdist(log_All, method = "bray", binary = FALSE, diag = FALSE, upper = FALSE, na.rm = FALSE)

#some data will need to be a matrix or dataframe later
bray_matrix_All <- as.matrix(bray_All)
write.table(bray_matrix_All, file = "All_bray_curtis.csv")


#my prefered way of generating figures - faster, works in base R
pdf(file = "All_bray.pdf", width = 10, height = 10);
heatmap(bray_matrix_All, main = "Bray Curtis")
dev.off()


#affinity propagation
All_ap_clust <- apcluster(negDistMat(r=3), bray_matrix_All, details = TRUE)

#r here is a power used to scale the data and allow for more easily recognizable separations
pdf(file = "All_ap_cluster.pdf", width = 10, height = 10);
heatmap(All_ap_clust)
dev.off()
All_ap_clust

color.palette2 = colorRampPalette(c("red4", "red3", "firebrick2", "indianred1", "salmon", "lightsalmon", "gold", "goldenrod1", "aquamarine2", "cyan3", "darkcyan", "darkslategrey"), space = "rgb")(100)

pdf(file = "Supp_Fig1.pdf", width = 10, height = 10);
heatmap(All_ap_clust, col = color.palette2)
dev.off()

#### Making a color_scale bar for figure
# Create a sequence of values from 0 to 1
values <- seq(0, 1, length.out = 100)
z_matrix <- matrix(values, nrow = 1, ncol = length(values))
# Plot the color scale as a rectangle using 'image'
par(mar = c(1, 2, 1, 2))  # Adjust margins to fit the color bar nicely
image(
  x = 1,                       # Single horizontal axis
  y = values,                   # The values corresponding to the colors
  z = z_matrix, # Create a matrix for image function
  col = color.palette2,          # Use the defined color palette
  xaxt = "n",                   # Remove x-axis (not needed)
  yaxt = "n",                   # Remove y-axis (will add a custom axis)
  bty = "n"                     # Remove box border
)
# Add the y-axis with labels 0 to 1 to represent the scale
axis(4, at = seq(0, 1, by = 0.2), labels = seq(0, 1, by = 0.2))





#writing to a tibble a list of samples and their respective cluster
#set the number of clusters generated in apcluster() here
All_clusters <- 5
#create blank tibble
All_types <- tibble()
types <- tibble()
#add each cluster to tibble
for (i in 1:All_clusters) {
  All_temp_types <- data.frame("site" = names(All_ap_clust[[i]]),
                               "type" = paste0("Cluster ", i))
  
  All_types <- bind_rows(All_types, All_temp_types)
}



#load the type data
All_type_no_change <- All_types
sites_All <- All_types$site
All_types <- as.data.frame(All_types)
All_types <- as.data.frame(All_types[, -1])
rownames(All_types) <- sites_All
names(All_types)[1] <- "type"
fac_All <-factor(All_types$type)
cols_All <- c("red", "blue", "orange", "green", "violet", "gray", "cyan")

#ordination
NMDS_All <- metaMDS(log_All, distance = "bray", k=2, trymax = 000, autotransform =TRUE, noshare = 0.1, expand = TRUE, trace = 1, plot = FALSE)
NMDS_All

# the stress value is printed here
#noshare means if there are fewer then this proportion of shared organisms, a stepacross (lowest possible similarity value) is prescribed
#k is the value of dimentions you want the data to be constrained to for a reported stress value
# we must reach convergence
pdf(file = "NMDS_stress_All.pdf", width = 12, height = 9);
stress_All <- stressplot(NMDS_All)
plot(stress_All)
dev.off()

#stress plot is just a relationship between ordination distance and dissimilarity
#goodness of fit values here, how well the visual representation of the data matches the dissimilarity matrix
# for stress < 0.1 good; 0.1<x<0.2 questionable; >0.2 bad

# now plot the ordination
pdf(file = "NMDS_All.pdf", width = 10, height = 10);
plot(NMDS_All, type="p", display="sites")
points(NMDS_All, display = "sites", col = cols_All[fac_All], pch = 19)
legend("topright", legend=levels(fac_All), col=cols_All, pch = 19)
# look at your stress data above and adjust this to be correct
text(x = -0.5, y = 0.4, "Stress = 9.538356e-05")
text(x = -0.5, y = 0.2, "Linear Fit = 1")
dev.off()


sur_dcm <- All_types
sur_dcm$type <- c(rep("DCM", 8), rep("SUR", 27))

sites_All <- sur_dcm$site
sur_dcm_types <- as.data.frame(sur_dcm)
sur_dcm_types <- as.data.frame(sur_dcm[, -1])
rownames(sur_dcm) <- sites_All
names(sur_dcm)[1] <- "types"
fac_sd <-factor(sur_dcm$type)
cols_sd <- c("salmon", "cyan3")

#ordination
NMDS_sd <- metaMDS(log_All, distance = "bray", k=2, trymax = 000, autotransform =TRUE, noshare = 0.1, expand = TRUE, trace = 1, plot = FALSE)
NMDS_sd


plot(NMDS_All, type="p", display="sites")
points(NMDS_All, display = "sites", col = cols_sd[fac_sd], pch = 19)
legend("topright", legend=levels(fac_sd), col=cols_sd, pch = 19, left)
# look at your stress data above and adjust this to be correct
text(x = -0.5, y = 0.4, "Stress = 9.538356e-05")
text(x = -1, y = 0.5, "P = 0.001")



#envfit, vector fitting for relating environmental data to ordination
env_data<-read.table("sample_data.txt", header=T, stringsAsFactors = T, sep = '\t') #read sample data table
attach(env_data)
env_data <- env_data[!(env_data$Sample_id %in% c("T4_D", "T52_S", "T88_S")), ]


# inverse simpson, simple alpha diversity
inv_simp_whole_All <- diversity(log_All, index = "invsimpson")


#tibble now with site, diversity and cluster
inv_simp_whole_tib_All <- tibble("site" = names(inv_simp_whole_All),
                                 "inv_simp" = inv_simp_whole_All,
                                 "cluster" = All_types)
inv_simp_whole_tib_All <- data.frame(site = names(inv_simp_whole_All),
                                     inv_simp = inv_simp_whole_All,
                                     cluster = All_types); inv_simp_whole_tib_All
colnames(inv_simp_whole_tib_All)[colnames(inv_simp_whole_tib_All) == 'type'] <- 'cluster'

write.table(inv_simp_whole_tib_All, file = "inv_simpson_t_All.csv")



alphaD_tib_All <- inner_join(x = All_type_no_change, y = inv_simp_whole_tib_All, by = "site")
pdf(file = "gginv_simp_All.pdf", width = 10, height = 10);
alphaD_plot_All <- ggplot(alphaD_tib_All, aes(x=cluster, y=inv_simp)) +
  geom_boxplot()
alphaD_plot_All
dev.off()


#kruskal wallis is there a signifigant difference between alpha diversities
# try playing with this to compare each group with eachother in a pairwise way - DEEP vs MEDITERR, DEEP vs INDIAN, INDIAN vs MEDITERR
kw_test_All <- kruskal.test(inv_simp_whole_All ~ cluster, data = inv_simp_whole_tib_All);kw_test_All
kw_test_All

########################### Shannon's H and making boxplots with kruskal wallis


# shannon's h, simple alpha diversity
shannons_whole_All <- diversity(log_All, index = "shannon")
#tibble now with site, diversity and cluster
shannons_whole_tib_All <- tibble("site" = names(shannons_whole_All),
                                 "shannon_h" = shannons_whole_All,
                                 "cluster" = All_types)
shannons_whole_tib_All <- data.frame(site = names(shannons_whole_All),
                                     shannon_h = shannons_whole_All,
                                     cluster = All_types); shannons_whole_tib_All
colnames(shannons_whole_tib_All)[colnames(shannons_whole_tib_All) == 'type'] <- 'cluster'
write.table(shannons_whole_tib_All, file = "shannons_t_All.csv")



shannon_alphaD_tib_All <- inner_join(x = All_type_no_change, y = shannons_whole_tib_All, by = "site")
pdf(file = "ggshannons_All.pdf", width = 10, height = 10);
shannon_alphaD_plot_All <- ggplot(shannon_alphaD_tib_All, aes(x=cluster, y=shannon_h)) +
  geom_boxplot()
shannon_alphaD_plot_All
dev.off()


#kruskal wallis is there a signifigant difference between alpha diversities
# try playing with this to compare each group with eachother in a pairwise way - DEEP vs MEDITERR, DEEP vs INDIAN, INDIAN vs MEDITERR
library("ggpubr")
library(ggplot2)
All_comparisons <- list(c("Cluster 1", "Cluster 2"), c("Cluster 1", "Cluster 3"), c("Cluster 1", "Cluster 4"), c("Cluster 1", "Cluster 5"),
                        c("Cluster 1", "Cluster 6"), c("Cluster 1", "Cluster 7"), c("Cluster 2", "Cluster 3"), c("Cluster 2", "Cluster 4"),
                        c("Cluster 2", "Cluster 5"), c("Cluster 2", "Cluster 6"), c("Cluster 2", "Cluster 7"), c("Cluster 3", "Cluster 4"),
                        c("Cluster 3", "Cluster 5"), c("Cluster 3", "Cluster 6"), c("Cluster 3", "Cluster 7"), c("Cluster 4", "Cluster 5"),
                        c("Cluster 4", "Cluster 6"), c("Cluster 4", "Cluster 7"), c("Cluster 5", "Cluster 6"), c("Cluster 5", "Cluster 7"),
                        c("Cluster 6", "Cluster 7"))

kw_shannons_all <- kruskal.test(shannons_whole_All ~ cluster, data = shannons_whole_tib_All)
p_value_shannons_all <- kw_shannons_all$p.value

Shannon_H_All_Boxplot <- ggplot(shannons_whole_tib_All, aes(x = cluster, y = shannon_h, fill = cluster)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               position = position_dodge(1)) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 8, outlier.size = 5) +
  ylab("Shannon's H") +
  xlab("Cluster")+
  theme_classic() +
  stat_compare_means(comparisons = All_comparisons) +
  stat_compare_means(label.y = 10.2)+
  theme(plot.title = element_text(color = "#0099f8", size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(face = "bold.italic", hjust = 0.5),
        plot.caption = element_text(face = "italic")) +
  ggtitle("Shannon's H All Data"); Shannon_H_All_Boxplot
ggsave("shannon's_H_all_comparisons.png", plot = Shannon_H_All_Boxplot, width = 6, height = 4)


########################## Inverse Simpson


Inv_Simp_All_Boxplot <- ggplot(inv_simp_whole_tib_All, aes(x = cluster, y = inv_simp, fill = cluster)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               position = position_dodge(1)) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 8, outlier.size = 5) +
  ylab("Inv_simp") +
  xlab("Cluster")+
  theme_classic() +
  stat_compare_means(comparisons = All_comparisons) +
  stat_compare_means(label.y = 27000, label.x = 1.25)+
  theme(plot.title = element_text(color = "#0099f8", size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(face = "bold.italic", hjust = 0.5),
        plot.caption = element_text(face = "italic")) +
  ggtitle("Inverse Simpson's All Data"); Inv_Simp_All_Boxplot
ggsave("inv_simp_all_comparisons.png", plot = Inv_Simp_All_Boxplot, width = 6, height = 4)


##############################################  MRPP test ############################################################
grouping_all <- c(rep("DCM", 8), rep("SUR", 27))
all_anosim_result <- anosim(bray_All, grouping_all)
all_nmds_p_value <- all_anosim_result$signif

mrpp_result <- mrpp(bray_All, grouping_all, permutations = 999, distance = "bray"); mrpp_result


###### Alpha depth ######
inv_simp_whole_tib_All$depth <- c(rep("DCM", 8), rep("SUR", 27))
kw_all_depth_inv <- kruskal.test(inv_simp_whole_All ~ depth, data = inv_simp_whole_tib_All)
p_value_all_depth <- kw_all_depth_inv$p.value

pdf(file = "Fig2_a.pdf", width = 6, height = 10)
all_depth_inv_simp <- ggplot(inv_simp_whole_tib_All, aes(x = depth, y = inv_simp, fill = depth)) +
  #  geom_dotplot(binaxis = 'y', stackdir = 'center',
  #              position = position_dodge(1)) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 8, outlier.size = 5, alpha = 0.5, color = c("salmon", "cyan3")) +
  ylab("Inverse Simpson's Index") +
  xlab("Depth")+
  theme_classic() +
  theme(text = element_text(size =18), plot.title = element_text(color = "#0099f8", size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(face = "bold.italic", hjust = 0.5),
        plot.caption = element_text(face = "italic")) +
  stat_compare_means(comparisons = list(c("DCM", "SUR")), 
                     method = "wilcox.test",
                     label = "p.signif") +
  scale_fill_manual(values = c("DCM" = "salmon", "SUR" = "cyan3")); all_depth_inv_simp
dev.off()
ggsave("Fig2a.png", plot = all_depth_inv_simp, width = 6, height = 4)



shannons_whole_tib_All$depth <- c(rep("DCM", 8), rep("SUR", 27))
kw_all_depth_shannons <- kruskal.test(shannons_whole_All ~ depth, data = shannons_whole_tib_All)
p_value_all_depth_shannon <- kw_all_depth_shannons$p.value

all_depth_shannon <- ggplot(shannons_whole_tib_All, aes(x = depth, y = shannon_h, fill = depth)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               position = position_dodge(1)) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 8, outlier.size = 5) +
  ylab("Shannon's H") +
  xlab("Depth")+
  geom_text(aes(label = sprintf("Kruskal-Wallis, p = %.3f", p_value_all_depth_shannon)), x = 1.2, y = max(shannons_whole_tib_All$shannon_h), hjust = 0.5) +
  theme_classic() +
  theme(plot.title = element_text(color = "#0099f8", size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(face = "bold.italic", hjust = 0.5),
        plot.caption = element_text(face = "italic")) +
  ggtitle("Shannon's H All"); all_depth_shannon

ggsave("shannon_h_all_depth.png", plot = all_depth_shannon, width = 6, height = 4)



Time_points <- c("T16","T40","T52","T64","T76","T88","T100","T112","T0","T4","T8","T12","T16","T20","T24","T28","T32","T36","T40","T44","T48",
                 "T56","T60","T64","T68","T72","T76","T80","T84","T92","T96","T100","T104","T108","T112")
inv_simp_whole_tib_All$Time_points <- Time_points


all_inv_simp <- ggplot(inv_simp_whole_tib_All, aes(x = Time_points, y = inv_simp, colour = depth)) +
  geom_point(size=3) +
  geom_line(aes(group=depth)) +
  scale_color_manual(values = c("salmon","cyan3")) +
  labs(color = "Depth") +  # Add this line to label the color legend
  theme(axis.text.x = element_text(angle = -90, hjust = 0)) +
  ggtitle("Inverse Simpson's All");all_inv_simp



# Reorder Time_points as a factor in the desired order
inv_simp_whole_tib_All$Time_points <- factor(inv_simp_whole_tib_All$Time_points, levels = c("T0","T4","T8","T12","T16","T20","T24","T28","T32",
                                                                                            "T36","T40","T44","T48","T52","T56","T60","T64",
                                                                                            "T68","T72","T76","T80","T84","T88","T92","T96",
                                                                                            "T100","T104","T108","T112"))

# Plot with reordered Time_points
all_inv_simp <- ggplot(inv_simp_whole_tib_All, aes(x = Time_points, y = inv_simp, colour = depth)) +
  geom_point(size=3) +
  geom_line(aes(group=depth)) +
  scale_color_manual(values = c("salmon","cyan3")) +
  labs(color = "Depth") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0)) +
  ggtitle("Inverse Simpson's All"); all_inv_simp





# Define breaks for the y-axis (intervals of 5,000)
y_breaks <- seq(22000, max(inv_simp_whole_tib_All$inv_simp), by = 2000)

# Create the plot with adjusted y-axis breaks
all_inv_simp <- ggplot(inv_simp_whole_tib_All, aes(x = Time_points, y = inv_simp, colour = depth)) +
  geom_point(size=3) +
  geom_line(aes(group=depth)) +
  scale_color_manual(values = c("salmon","cyan3")) +
  ylab("Inverse Simpsons") +
  xlab("Time Points") +
  labs(color = "Depth") +
  theme_minimal() +
  theme(text = element_text(size =16), axis.text.x = element_text(angle = -90, hjust = 0)) +
  scale_y_continuous(breaks = y_breaks) +  # Set y-axis breaks
  coord_cartesian(ylim = c(0, max(inv_simp_whole_tib_All$inv_simp) * 1.1)) ; all_inv_simp
ggsave("all_inv_simp_plot_zoomed_out.png", plot = all_inv_simp, width = 10, height = 6)



site_D <- c("T28_D", "T4_D", "T52_S", "T88_S")
side_d1 <- c("T28_D", "T4_D", "T52_S", "T88_S")
inv_simp_d <- c(NA)
cluster_d <- c("Cluster 0")
depth_d <- c("DCM", "DCM", "SUR", "SUR")
time_d <- c("T28","T4","T52", "T88")
removed <- data.frame(column = site_D, site = side_d1, inv_simp = inv_simp_d, cluster = cluster_d, depth = depth_d, Time_points = time_d)
row.names(removed) <- removed$column

# Remove the first column from the dataframe
removed <- removed[, -1]
inv_all_zoom_out <- inv_simp_whole_tib_All
inv_all_zoom_out <- rbind(removed, inv_all_zoom_out)
# Reorder Time_points as a factor in the desired order
inv_all_zoom_out$Time_points <- factor(inv_all_zoom_out$Time_points, levels = c("T0","T4","T8","T12","T16","T20","T24","T28","T32",
                                                                                "T36","T40","T44","T48","T52","T56","T60","T64",
                                                                                "T68","T72","T76","T80","T84","T88","T92","T96",
                                                                                "T100","T104","T108","T112"))
pdf(file = "Fig2a_2.pdf", width = 6, height = 6)
all_inv_simp2 <- ggplot(inv_all_zoom_out, aes(x = Time_points, y = inv_simp, colour = depth)) +
  geom_point(size=3) +
  geom_line(aes(group=depth)) +
  scale_color_manual(values = c("salmon","cyan3")) +
  ylab("Inverse Simpsons") +
  xlab("Time Points") +
  labs(color = "Depth") +
  theme(text = element_text(size =16), axis.text.x = element_text(angle = -90, hjust = 0)) +
  theme_minimal() +
  scale_y_continuous(breaks = y_breaks) +  # Set y-axis breaks
  coord_cartesian(ylim = c(0, max(inv_all_zoom_out$inv_simp) * 1.1)) ; all_inv_simp2
dev.off()
ggsave("Fig2a_2.png", plot = all_inv_simp, width = 10, height = 6)














################################ SUR Data Only ##################################
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

colnames(meancovsur)
cleanedreads_sur <- scan("cleanedreads_sur.txt")

norm3 <- sweep(meancovsur, 2, cleanedreads_sur, FUN = '/')
vOTU_sur<- norm3*1000000000
vOTU_sur_no_outliers <- vOTU_sur[, !(colnames(vOTU_sur) %in% c("T52_S", "T88_S"))]

#read in and format abundance data
norm_abund_t_SUR <- t(vOTU_sur_no_outliers) %>%
  as_tibble(norm_abund_t_SUR, rownames = NA)

#log transformation
#note this is ln
log_SUR <- log(norm_abund_t_SUR +1)

#bray curtis
bray_SUR <- vegdist(log_SUR, method = "bray", binary = FALSE, diag = FALSE, upper = FALSE, na.rm = FALSE)

#some data will need to be a matrix or dataframe later
bray_matrix_SUR <- as.matrix(bray_SUR)
write.table(bray_matrix_SUR, file = "SUR_bray_curtis.csv")

color.palette2 = colorRampPalette(c("red4", "red3", "firebrick2", "indianred1", "salmon", "lightsalmon", "gold", "goldenrod1", "aquamarine2", "cyan3", "darkcyan", "darkslategrey"), space = "rgb")(100)

#my prefered way of generating figures - faster, works in base R
pdf(file = "SUR_bray.pdf", width = 10, height = 10);
heatmap(bray_matrix_SUR, main = "Bray Curtis")
dev.off()


#affinity propagation
set.seed(111)
SUR_ap_clust <- apcluster(negDistMat(r=3), bray_matrix_SUR, details = TRUE)

#r here is a power used to scale the data and allow for more easily recognizable separations
pdf(file = "SUR_ap_cluster.pdf", width = 10, height = 10);
heatmap(SUR_ap_clust, col = color.palette2)
dev.off()
SUR_ap_clust

#writing to a tibble a list of samples and their respective cluster
#set the number of clusters generated in apcluster() here
SUR_clusters <- 7

#create blank tibble
SUR_types <- tibble()

types <- tibble()
#add each cluster to tibble
for (i in 1:SUR_clusters) {
  SUR_temp_types <- data.frame("site" = names(SUR_ap_clust[[i]]),
                               "type" = paste0("Cluster ", i))
  
  SUR_types <- bind_rows(SUR_types, SUR_temp_types)
}

#load the type data
SUR_type_no_change <- SUR_types
sites_SUR <- SUR_types$site
SUR_types <- as.data.frame(SUR_types)
SUR_types <- as.data.frame(SUR_types[, -1])
rownames(SUR_types) <- sites_SUR
names(SUR_types)[1] <- "type"
fac_SUR <-factor(SUR_types$type)
cols_SUR <- c("red", "blue", "orange", "green", "violet", "gray", "cyan", "black")

#ordination
NMDS_SUR <- metaMDS(log_SUR, distance = "bray", k=2, trymax = 000, autotransform =TRUE, noshare = 0.1, expand = TRUE, trace = 1, plot = FALSE)
NMDS_SUR

# the stress value is printed here
#noshare means if there are fewer then this proportion of shared organisms, a stepacross (lowest possible similarity value) is prescribed
#k is the value of dimentions you want the data to be constrained to for a reported stress value
# we must reach convergence
pdf(file = "NMDS_stress_SUR.pdf", width = 12, height = 9);
stress_SUR <- stressplot(NMDS_SUR)
plot(stress_SUR)
dev.off()

#stress plot is just a relationship between ordination distance and dissimilarity
#goodness of fit values here, how well the visual representation of the data matches the dissimilarity matrix
# for stress < 0.1 good; 0.1<x<0.2 questionable; >0.2 bad

# now plot the ordination
pdf(file = "NMDS_SUR.pdf", width = 10, height = 10);
plot(NMDS_SUR, type="p", display="sites")
points(NMDS_SUR, display = "sites", col = cols_SUR[fac_SUR], pch = 19)
legend("topright", legend=levels(fac_SUR), col=cols_SUR, pch = 19)
# look at your stress data above and adjust this to be correct
text(x = -0.1, y = 0.06, "Stress = 0.0652615")
text(x = -0.1, y = 0.05, "Linear Fit = 0.989")
dev.off()



#envfit, vector fitting for relating environmental data to ordination
Surface_env_table <-read.table("Surface_Data.txt", header=T, stringsAsFactors = T, sep = '\t') #read sample data table
Surface_env_table <- Surface_env_table[!(Surface_env_table$Sample_id %in% c("T52_S", "T88_S")), ]

attach(Surface_env_table)
NMDS_meta_SUR = envfit(NMDS_SUR~Time+Salinity+Temp, permutations=9999)
NMDS_meta_SUR
pdf(file = "envfit_SUR.pdf", width = 10, height = 10);
plot(NMDS_SUR, type="p", display="sites")
points(NMDS_SUR, display = "sites", col = cols_SUR[fac_SUR], pch = 19)
legend("topright", legend=levels(fac_SUR), col=cols_SUR, pch = 19)
text(x = -0.1, y = 0.06, "Stress = 0.07051535")
text(x = -0.1, y = 0.05, "Linear Fit = 0.988")
# please fill in the correct values above
plot(NMDS_meta_SUR)
dev.off()
#arrows represent strength of r2 values


# inverse simpson, simple alpha diversity
inv_simp_SUR <- diversity(log_SUR, index = "invsimpson", groups = SUR_types$type)
inv_simp_whole_SUR <- diversity(log_SUR, index = "invsimpson")

inv_simp_SUR
inv_simp_whole_SUR


#tibble now with site, diversity and cluster
inv_simp_whole_tib_SUR <- tibble("site" = names(inv_simp_whole_SUR),
                                 "inv_simp" = inv_simp_whole_SUR,
                                 "cluster" = SUR_types)
inv_simp_whole_tib_SUR <- data.frame(site = names(inv_simp_whole_SUR),
                                     inv_simp = inv_simp_whole_SUR,
                                     cluster = SUR_types); inv_simp_whole_tib_SUR
colnames(inv_simp_whole_tib_SUR)[colnames(inv_simp_whole_tib_SUR) == 'type'] <- 'cluster'

write.table(inv_simp_whole_tib_SUR, file = "inv_simpson_t_SUR.csv")



alphaD_tib_SUR <- inner_join(x = SUR_type_no_change, y = inv_simp_whole_tib_SUR, by = "site")
pdf(file = "gginv_simp_SUR.pdf", width = 10, height = 10);
alphaD_plot_SUR <- ggplot(alphaD_tib_SUR, aes(x=cluster, y=inv_simp)) +
  geom_boxplot()
alphaD_plot_SUR
dev.off()


#kruskal wallis is there a signifigant difference between alpha diversities
# try playing with this to compare each group with eachother in a pairwise way - DEEP vs MEDITERR, DEEP vs INDIAN, INDIAN vs MEDITERR
kw_test_SUR <- kruskal.test(inv_simp_whole_SUR ~ cluster, data = inv_simp_whole_tib_SUR);kw_test_SUR
kw_test_SUR

########################### Shannon's H and making boxplots with kruskal wallis


# shannon's h, simple alpha diversity
shannons_SUR <- diversity(log_SUR, index = "shannon", groups = SUR_types$type)
shannons_whole_SUR <- diversity(log_SUR, index = "shannon")
shannons_SUR
shannons_whole_SUR


#tibble now with site, diversity and cluster
shannons_whole_tib_SUR <- tibble("site" = names(shannons_whole_SUR),
                                 "shannon_h" = shannons_whole_SUR,
                                 "cluster" = SUR_types)
shannons_whole_tib_SUR <- data.frame(site = names(shannons_whole_SUR),
                                     shannon_h = shannons_whole_SUR,
                                     cluster = SUR_types); shannons_whole_tib_SUR
colnames(shannons_whole_tib_SUR)[colnames(shannons_whole_tib_SUR) == 'type'] <- 'cluster'

write.table(shannons_whole_tib_SUR, file = "shannons_t_SUR.csv")



shannon_alphaD_tib_SUR <- inner_join(x = SUR_type_no_change, y = shannons_whole_tib_SUR, by = "site")
pdf(file = "ggshannons_SUR.pdf", width = 10, height = 10);
shannon_alphaD_plot_SUR <- ggplot(shannon_alphaD_tib_SUR, aes(x=cluster, y=shannon_h)) +
  geom_boxplot()
shannon_alphaD_plot_SUR
dev.off()


#kruskal wallis is there a signifigant difference between alpha diversities
# try playing with this to compare each group with eachother in a pairwise way - DEEP vs MEDITERR, DEEP vs INDIAN, INDIAN vs MEDITERR
SUR_comparisons <- list(c("Cluster 1", "Cluster 2"), c("Cluster 1", "Cluster 3"), c("Cluster 1", "Cluster 4"), c("Cluster 1", "Cluster 5"),
                        c("Cluster 1", "Cluster 6"), c("Cluster 1", "Cluster 7"), c("Cluster 2", "Cluster 3"), 
                        c("Cluster 2", "Cluster 4"), c("Cluster 2", "Cluster 5"), c("Cluster 2", "Cluster 6"), c("Cluster 2", "Cluster 7"), 
                        c("Cluster 3", "Cluster 4"), c("Cluster 3", "Cluster 5"), c("Cluster 3", "Cluster 6"), 
                        c("Cluster 3", "Cluster 7"), c("Cluster 4", "Cluster 5"), c("Cluster 4", "Cluster 6"), 
                        c("Cluster 4", "Cluster 7"), c("Cluster 5", "Cluster 6"), c("Cluster 5", "Cluster 7"),
                        c("Cluster 6", "Cluster 7"))

kw_shannons_SUR <- kruskal.test(shannons_whole_SUR ~ cluster, data = shannons_whole_tib_SUR)
p_value_shannons_SUR <- kw_shannons_SUR$p.value

Shannon_H_SUR_Boxplot <- ggplot(shannons_whole_tib_SUR, aes(x = cluster, y = shannon_h, fill = cluster)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               position = position_dodge(1)) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 8, outlier.size = 5) +
  ylab("Shannon's H") +
  xlab("Cluster")+
  theme_classic() +
  stat_compare_means(comparisons = SUR_comparisons) +
  stat_compare_means(label.y = 10.6, label.x = 2)+
  theme(plot.title = element_text(color = "#0099f8", size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(face = "bold.italic", hjust = 0.5),
        plot.caption = element_text(face = "italic")) +
  ggtitle("Shannon's H Surface Data")
ggsave("shannon's_H_SUR_comparisons.png", plot = Shannon_H_SUR_Boxplot, width = 6, height = 4)



########################## Inverse Simpson


Inv_Simp_SUR_Boxplot <- ggplot(inv_simp_whole_tib_SUR, aes(x = cluster, y = inv_simp, fill = cluster)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               position = position_dodge(1)) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 8, outlier.size = 5) +
  ylab("Inv_simp") +
  xlab("Cluster")+
  theme_classic() +
  stat_compare_means(comparisons = SUR_comparisons) +
  stat_compare_means(label.y = 36000, label.x = 1.25)+
  theme(plot.title = element_text(color = "#0099f8", size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(face = "bold.italic", hjust = 0.5),
        plot.caption = element_text(face = "italic")) +
  ggtitle("Inverse Simpson's SUR Data"); Inv_Simp_SUR_Boxplot
ggsave("inv_simp_SUR_comparisons.png", plot = Inv_Simp_SUR_Boxplot, width = 6, height = 4)






##############################################  MRPP test ############################################################
bray_SUR
bray_matrix_SUR
inv_simp_whole_tib_SUR
colnames(vOTU_sur)
grouping_sur_day_night <- c("Day", "Night", "Night", "Night", "Day", "Day", "Day", "Night", "Night", "Night", 
                            "Day", "Day", "Day", "Night", "Night", "Day", "Day", "Day", "Night", "Night", 
                            "Night", "Day", "Day", "Night", "Night", "Night", "Day")

set.seed(111)
mrpp_result <- mrpp(bray_SUR, grouping_sur_day_night, permutations = 999, distance = "bray"); mrpp_result

cols_SUR_diel <- c("goldenrod3", "cyan")
SUR_types_diel <- SUR_types
SUR_types_diel$diel <- grouping_sur_day_night
fac_SUR_diel <- factor(SUR_types_diel$diel)

plot(NMDS_SUR, type="p", display="sites")
points(NMDS_SUR, display = "sites", col = cols_SUR_diel[fac_SUR_diel], pch = 19)
legend("topright", legend=levels(fac_SUR_diel), col=cols_SUR_diel, pch = 19)
text(x = -0.1, y = 0.08, "Stress = 0.07051535")
text(x = -0.1, y = 0.07, "Linear Fit = 0.988")
text(x = -0.1, y = 0.06, "MRPP = 0.068")

ordihull(NMDS_SUR, fac_SUR_diel, display = "sites", draw = c("polygon"), col = cols_SUR_diel, alpha = 100, label = FALSE)

Surface_env_table <- Surface_env_table[, -c(1,2,4,5)]
Surface_env_table$diel <- grouping_sur_day_night

attach(Surface_env_table)

NMDS_meta_SUR_diel = envfit(NMDS_SUR~Salinity+Temp+diel, permutations=9999)
NMDS_meta_SUR_diel

pdf(file = "envfit_SUR.pdf", width = 10, height = 10);
plot(NMDS_SUR, type="p", display="sites")
points(NMDS_SUR, display = "sites", col = cols_SUR_diel[fac_SUR_diel], pch = 19)
legend("topright", legend=levels(fac_SUR_diel), col=cols_SUR_diel, pch = 19)
text(x = -0.1, y = 0.08, "Stress = 0.07051535")
text(x = -0.1, y = 0.07, "Linear Fit = 0.989")
text(x = -0.1, y = 0.06, "MRPP = 0.068")

ordihull(NMDS_SUR, fac_SUR_diel, display = "sites", draw = c("polygon"), col = cols_SUR_diel, alpha = 100, label = FALSE)
# please fill in the correct values above
plot(NMDS_meta_SUR_diel)
dev.off()


########## Plotting the PCoA ##################


#read in and format abundance data
norm_abund_t_SUR <- t(vOTU_sur_no_outliers)
norm_abund_t_SUR <- as_tibble(norm_abund_t_SUR, rownames = NA)
#log transformation
#note this is ln
log_SUR <- log(norm_abund_t_SUR +1)
#bray curtis
bray_SUR <- vegdist(log_SUR, method = "bray", binary = FALSE, diag = FALSE, upper = FALSE, na.rm = FALSE)
#ordination
NMDS_SUR <- metaMDS(log_SUR, distance = "bray", k=2, trymax = 000, autotransform =TRUE, noshare = 0.1, expand = TRUE, trace = 1, plot = FALSE)
NMDS_SUR

pcoa_sur <- cmdscale(bray_SUR, k = 2, eig = T); pcoa_sur

# extract axis positions for each site from cmdscale object and create a dataframe for plotting
pcoa_plotting <- as.data.frame(pcoa_sur$points)
colnames(pcoa_plotting) <- c("axis_1", "axis_2")
pcoa_plotting$site <- rownames(pcoa_plotting)


# calculate the proportion of variance in the data which is explained by the first two PCoA axes
pcoa_sur$eig[1]/(sum(pcoa_sur$eig))
pcoa_sur$eig[2]/(sum(pcoa_sur$eig))





pcoa_plotting$group <- grouping_sur_day_night
pdf(file = "SUR_diel_PCoA.pdf", width = 10, height = 10)
pcoa <- ggplot(pcoa_plotting, aes(x = axis_1, y = axis_2, colour = group, fill = group)) +
  geom_point(alpha = 0.5, size = 3) +
  scale_color_manual(values = c("yellow3", "black")) +
  scale_fill_manual(values = c("yellow3", "black")) +
  theme(text = element_text(size = 16), axis.text = element_text(angle = -90, hjust = 0), legend.position = 'none') +
  xlab("PCoA 1 (30.89%)") +
  ylab("PCoA 2 (19.85%)") +
  theme_minimal() +
  annotate(geom = 'text', label = 'MRPP p = 0.068', x = 0.02, y = 0.045) +
  annotate(geom = 'text', label = 'Bray-Curtis', x = Inf, y = -Inf, hjust = 1.15, vjust = -1); pcoa
pcoa <- pcoa + theme(plot.title = element_text(hjust = 0.5)); pcoa
pcoa <- pcoa + theme(legend.position = 'right'); pcoa
ggsave("SUR_PCoA.png", plot = pcoa, width = 6, height = 4)


# Add ellipses
pcoa_ellipses <- pcoa +
  theme(text = element_text(size = 18)) +
  stat_ellipse(geom = "polygon", aes(group = group, fill = group), alpha = 0.2, colour = "black", linetype = "dashed"); pcoa_ellipses
ggsave("SUR_PCoA_ellipse.png", plot = pcoa_ellipses, width = 6, height = 4)
dev.off()



set.seed(111)
mrpp_result <- mrpp(bray_SUR, grouping_sur_day_night, permutations = 999, distance = "bray"); mrpp_result


ordination_scores_sur <- pcoa_plotting %>%
  select(axis_1, axis_2)  # Extract PCoA axes
sur_diel <- pcoa_plotting$group  # Group variable

# Use ordihull to calculate hull points
ordination_sur <- as.data.frame(ordination_scores_sur)
ordination_sur$diel <- sur_diel
hull_points_sur <- do.call(rbind, by(ordination_sur, ordination_sur$diel, function(df) {
  hull_idx <- chull(df$axis_1, df$axis_2)  # Find convex hull
  df[hull_idx, ]
}))
hull_points_sur$diel <- factor(hull_points_sur$diel)

SUR_permanova <- adonis2(bray_SUR ~ group,
                         data = pcoa_plotting,
                         by = "margin",
                         permutations = 999,
                         method = "bray"); SUR_permanova


pdf(file = "Fig2c.pdf", width = 10, height = 10)
pcoa <- ggplot(pcoa_plotting, aes(x = axis_1, y = axis_2, colour = group, fill = group)) +
  geom_point(alpha = 0.5, size = 6) +
  geom_polygon(data = hull_points_sur, aes(x = axis_1, y = axis_2, fill = diel, group = diel), 
               alpha = 0.2, colour = "black", linetype = "dotted") +  # Add convex hulls
  scale_color_manual(values = c("orange3", "black")) +
  scale_fill_manual(values = c("yellow2", "black")) +
  theme(text = element_text(size = 16), axis.text = element_text(angle = -90, hjust = 0), legend.position = 'none') +
  xlab("PCoA 1 (30.89%)") +
  ylab("PCoA 2 (19.85%)") +
  theme_minimal() +
  annotate(geom = 'text', label = 'MRPP p = 0.068', x = 0.02, y = 0.045) +
  annotate(geom = 'text', label = 'Permanova p = 0.074', x = 0.02, y = 0.045) +
  annotate(geom = 'text', label = 'Bray-Curtis', x = Inf, y = -Inf, hjust = 1.15, vjust = -1); pcoa
pcoa <- pcoa + theme(plot.title = element_text(hjust = 0.5)); pcoa
pcoa <- pcoa + theme(legend.position = 'right'); pcoa
ggsave("Fig2c.png", plot = pcoa, width = 6, height = 4)
dev.off()





inv_simp_whole_tib_SUR$diel <- grouping_sur_day_night

kw_sur_diel_inv <- kruskal.test(inv_simp_whole_SUR ~ diel, data = inv_simp_whole_tib_SUR)
p_value_sur_diel <- kw_sur_diel_inv$p.value

pdf(file = "Fig2b.pdf", width = 6, height = 10)
Inverse_simp_diel_sur <- ggplot(inv_simp_whole_tib_SUR, aes(x = diel, y = inv_simp, fill = diel)) +
  # geom_dotplot(binaxis = 'y', stackdir = 'center',
  #             position = position_dodge(1)) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 8, outlier.size = 5, alpha = 0.5, color = c("yellow3", "black")) +
  ylab("Inverse Simpsons") +
  xlab("Time")+
  theme_classic() +
  theme(text = element_text(size =18), plot.title = element_text(color = "#0099f8", size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(face = "bold.italic", hjust = 0.5),
        plot.caption = element_text(face = "italic")) +
  stat_compare_means(comparisons = list(c("Day", "Night")), 
                     method = "wilcox.test") +
  scale_fill_manual(values = c("Night" = "black", "Day" = "yellow")); Inverse_simp_diel_sur
ggsave("Fig2b.png", plot = Inverse_simp_diel_sur, width = 6, height = 4)
dev.off()









################################## DCM Data Only #######################################
DCM_meancov <- read.csv("dcm-meancov.txt", sep = "\t", row.names = 1) #new file
DCM <- DCM_meancov

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
cleanedreads_dcm <- scan("cleaned_reads.txt")
norm4 <- sweep(DCM_meancov, 2, cleanedreads_dcm, FUN = '/')
vOTU_dcm<- norm4*1000000000

#remove outliers T4_D
vOTU_dcm_no_outlier <- vOTU_dcm[, !(colnames(vOTU_dcm) %in% c("T4_D"))]


#read in and format abundance data
norm_abund_t_DCM <- t(vOTU_dcm_no_outlier)
norm_abund_t_DCM <- as_tibble(norm_abund_t_DCM, rownames = NA)

#log transformation
#note this is ln
log_DCM <- log(norm_abund_t_DCM +1)


#bray curtis
bray_DCM <- vegdist(log_DCM, method = "bray", binary = FALSE, diag = FALSE, upper = FALSE, na.rm = FALSE)

#some data will need to be a matrix or dataframe later
bray_matrix_DCM <- as.matrix(bray_DCM)
write.table(bray_matrix_DCM, file = "DCM_bray_curtis.csv")


#my prefered way of generating figures - faster, works in base R
pdf(file = "DCM_bray.pdf", width = 10, height = 10);
heatmap(bray_matrix_DCM, main = "Bray Curtis")
dev.off()


#affinity propagation
set.seed(111)
DCM_ap_clust <- apcluster(negDistMat(r=3), bray_matrix_DCM, details = TRUE)

#r here is a power used to scale the data and allow for more easily recognizable separations
pdf(file = "DCM_ap_cluster.pdf", width = 10, height = 10);
heatmap(DCM_ap_clust)
dev.off()
DCM_ap_clust

#writing to a tibble a list of samples and their respective cluster
#set the number of clusters generated in apcluster() here
DCM_clusters <- 3

#create blank tibble
DCM_types <- tibble()

types <- tibble()
#add each cluster to tibble
for (i in 1:DCM_clusters) {
  DCM_temp_types <- data.frame("site" = names(DCM_ap_clust[[i]]),
                               "type" = paste0("Cluster ", i))
  
  DCM_types <- bind_rows(DCM_types, DCM_temp_types)
}



#load the type data
DCM_type_no_change <- DCM_types
sites_DCM <- DCM_types$site
DCM_types <- as.data.frame(DCM_types)
DCM_types <- as.data.frame(DCM_types[, -1])
rownames(DCM_types) <- sites_DCM
names(DCM_types)[1] <- "type"
fac_DCM <-factor(DCM_types$type)
cols_DCM <- c("red", "blue", "green")

#ordination
NMDS_DCM <- metaMDS(log_DCM, distance = "bray", k=2, trymax = 000, autotransform =TRUE, noshare = 0.1, expand = TRUE, trace = 1, plot = FALSE)
NMDS_DCM

# the stress value is printed here
#noshare means if there are fewer then this proportion of shared organisms, a stepacross (lowest possible similarity value) is prescribed
#k is the value of dimentions you want the data to be constrained to for a reported stress value
# we must reach convergence
pdf(file = "NMDS_stress_DCM.pdf", width = 12, height = 9);
stress_DCM <- stressplot(NMDS_DCM)
plot(stress_DCM)
dev.off()

#stress plot is just a relationship between ordination distance and dissimilarity
#goodness of fit values here, how well the visual representation of the data matches the dissimilarity matrix
# for stress < 0.1 good; 0.1<x<0.2 questionable; >0.2 bad

# now plot the ordination
pdf(file = "NMDS_DCM.pdf", width = 10, height = 10);
plot(NMDS_DCM, type="p", display="sites")
points(NMDS_DCM, display = "sites", col = cols_DCM[fac_DCM], pch = 19)
legend("topright", legend=levels(fac_DCM), col=cols_DCM, pch = 19)
# look at your stress data above and adjust this to be correct
text(x = -0.1, y = 0.06, "Stress = 9.62045e-05  ")
text(x = -0.1, y = 0.05, "Linear Fit = 1")
dev.off()



#envfit, vector fitting for relating environmental data to ordination
DCM_env_table <-read.table("sampledata_dcm.txt", header=T, stringsAsFactors = T, sep = '\t') #read sample data table
DCM_env_table

DCM_env_table <- DCM_env_table[!(DCM_env_table$sample_id %in% c("T4_D")), ]; DCM_env_table
DCM_env_table$diel <- c("Day", "Day", "Night", "Day", "Night", "Day", "Night", "Day")
attach(DCM_env_table)

DCM_env_table

NMDS_meta_DCM = envfit(NMDS_DCM~sample_id+Date+Time+Lat+Long+Salinity+Temp+diel, permutations=9999)
NMDS_meta_DCM
pdf(file = "envfit_DCM.pdf", width = 10, height = 10);
plot(NMDS_DCM, type="p", display="sites")
points(NMDS_DCM, display = "sites", col = cols_DCM[fac_DCM], pch = 19)
legend("topright", legend=levels(fac_DCM), col=cols_DCM, pch = 19)
text(x = -0.1, y = 0.06, "Stress = 0.002970987")
text(x = -0.1, y = 0.05, "Linear Fit = 1")
# please fill in the correct values above
plot(NMDS_meta_DCM)
dev.off()
#arrows represent strength of r2 values

# inverse simpson, simple alpha diversity
inv_simp_DCM <- diversity(log_DCM, index = "invsimpson", groups = DCM_types$type)
inv_simp_whole_DCM <- diversity(log_DCM, index = "invsimpson")
inv_simp_DCM
inv_simp_whole_DCM


#tibble now with site, diversity and cluster
inv_simp_whole_tib_DCM <- tibble("site" = names(inv_simp_whole_DCM),
                                 "inv_simp" = inv_simp_whole_DCM,
                                 "cluster" = DCM_types)
inv_simp_whole_tib_DCM <- data.frame(site = names(inv_simp_whole_DCM),
                                     inv_simp = inv_simp_whole_DCM,
                                     cluster = DCM_types); inv_simp_whole_tib_DCM
colnames(inv_simp_whole_tib_DCM)[colnames(inv_simp_whole_tib_DCM) == 'type'] <- 'cluster'

write.table(inv_simp_whole_tib_DCM, file = "inv_simpson_t_DCM.csv")



alphaD_tib_DCM <- inner_join(x = DCM_type_no_change, y = inv_simp_whole_tib_DCM, by = "site")
pdf(file = "gginv_simp_DCM.pdf", width = 10, height = 10);
alphaD_plot_DCM <- ggplot(alphaD_tib_DCM, aes(x=cluster, y=inv_simp)) +
  geom_boxplot()
alphaD_plot_DCM
dev.off()


#kruskal wallis is there a signifigant difference between alpha diversities
# try playing with this to compare each group with eachother in a pairwise way - DEEP vs MEDITERR, DEEP vs INDIAN, INDIAN vs MEDITERR
kw_test_DCM <- kruskal.test(inv_simp_whole_DCM ~ cluster, data = inv_simp_whole_tib_DCM);kw_test_DCM
kw_test_DCM

########################### Shannon's H and making boxplots with kruskal wallis


# shannon's h, simple alpha diversity
shannons_DCM <- diversity(log_DCM, index = "shannon", groups = DCM_types$type)
shannons_whole_DCM <- diversity(log_DCM, index = "shannon")


shannons_DCM
shannons_whole_DCM


#tibble now with site, diversity and cluster
shannons_whole_tib_DCM <- tibble("site" = names(shannons_whole_DCM),
                                 "shannon_h" = shannons_whole_DCM,
                                 "cluster" = DCM_types)
shannons_whole_tib_DCM <- data.frame(site = names(shannons_whole_DCM),
                                     shannon_h = shannons_whole_DCM,
                                     cluster = DCM_types); shannons_whole_tib_DCM
colnames(shannons_whole_tib_DCM)[colnames(shannons_whole_tib_DCM) == 'type'] <- 'cluster'

write.table(shannons_whole_tib_DCM, file = "shannons_t_DCM.csv")



shannon_alphaD_tib_DCM <- inner_join(x = DCM_type_no_change, y = shannons_whole_tib_DCM, by = "site")
pdf(file = "ggshannons_DCM.pdf", width = 10, height = 10);
shannon_alphaD_plot_DCM <- ggplot(shannon_alphaD_tib_DCM, aes(x=cluster, y=shannon_h)) +
  geom_boxplot()
shannon_alphaD_plot_DCM
dev.off()


#kruskal wallis is there a signifigant difference between alpha diversities
# try playing with this to compare each group with eachother in a pairwise way - DEEP vs MEDITERR, DEEP vs INDIAN, INDIAN vs MEDITERR
DCM_comparisons <- list(c("Cluster 1", "Cluster 2"), c("Cluster 1", "Cluster 3"), c("Cluster 2", "Cluster 3"))

kw_shannons_DCM <- kruskal.test(shannons_whole_DCM ~ cluster, data = shannons_whole_tib_DCM)
p_value_shannons_DCM <- kw_shannons_DCM$p.value

Shannon_H_DCM_Boxplot <- ggplot(shannons_whole_tib_DCM, aes(x = cluster, y = shannon_h, fill = cluster)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               position = position_dodge(1)) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 8, outlier.size = 5) +
  ylab("Shannon's H") +
  xlab("Cluster")+
  theme_classic() +
  stat_compare_means(comparisons = DCM_comparisons) +
  stat_compare_means(label.y = 10.25, label.x = 2)+
  theme(plot.title = element_text(color = "#0099f8", size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(face = "bold.italic", hjust = 0.5),
        plot.caption = element_text(face = "italic")) +
  ggtitle("Shannon's H Deep Chlorophyll Maximum Data")
ggsave("shannon's_H_DCM_comparisons.png", plot = Shannon_H_DCM_Boxplot, width = 6, height = 4)



DCM_box_shannon <- ggplot(shannons_whole_tib_DCM, aes(x = cluster, y = shannon_h, fill = cluster)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               position = position_dodge(1)) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 8, outlier.size = 5) +
  ylab("Shannon's H") +
  xlab("Cluster")+
  geom_text(aes(label = sprintf("Kruskal-Wallis, p = %.3f", p_value_shannons_DCM)), x = 2.55, y = 10.175, hjust = 0.5) +
  theme_classic() +
  theme(plot.title = element_text(color = "#0099f8", size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(face = "bold.italic", hjust = 0.5),
        plot.caption = element_text(face = "italic")) +
  ggtitle("Shannon's H DCM"); DCM_box_shannon

ggsave("shannon_h_DCM.png", plot = DCM_box_shannon, width = 6, height = 4)




########################## Inverse Simpson


Inv_Simp_DCM_Boxplot <- ggplot(inv_simp_whole_tib_DCM, aes(x = cluster, y = inv_simp, fill = cluster)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               position = position_dodge(1)) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 8, outlier.size = 5) +
  ylab("Inv_simp") +
  xlab("Cluster")+
  theme_classic() +
  stat_compare_means(comparisons = DCM_comparisons) +
  stat_compare_means(label.y = 27000)+
  theme(plot.title = element_text(color = "#0099f8", size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(face = "bold.italic", hjust = 0.5),
        plot.caption = element_text(face = "italic")) +
  ggtitle("Inverse Simpson's DCM Data")
ggsave("inv_simp_DCM_comparisons.png", plot = Inv_Simp_DCM_Boxplot, width = 6, height = 4)


##############################################  MRPP test ############################################################
bray_DCM
bray_matrix_DCM
inv_simp_whole_tib_DCM
colnames(vOTU_dcm)
grouping <- inv_simp_whole_tib_DCM$cluster
mrpp_result <- mrpp(bray_DCM, grouping, permutations = 999, distance = "bray"); mrpp_result



inv_simp_whole_tib_DCM$diel <- c("Day", "Day", "Night", "Day", "Night", "Day", "Night", "Day")


kw_dcm_diel_inv <- kruskal.test(inv_simp_whole_DCM ~ diel, data = inv_simp_whole_tib_DCM)
p_value_dcm_diel <- kw_dcm_diel_inv$p.value

pdf(file = "DCM_invsimpson_diel_boxplot.pdf", width = 6, height = 10)
Inverse_simp_diel_dcm <- ggplot(inv_simp_whole_tib_DCM, aes(x = diel, y = inv_simp, fill = diel)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               position = position_dodge(1)) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 8, outlier.size = 5, alpha = 0, color = c("yellow3", "black")) +
  ylab("Inverse Simpsons") +
  xlab("Time")+
  theme_classic() +
  theme(text = element_text(size =18), plot.title = element_text(color = "#0099f8", size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(face = "bold.italic", hjust = 0.5),
        plot.caption = element_text(face = "italic")) +
  stat_compare_means(comparisons = list(c("Day", "Night")), 
                     method = "wilcox.test") +
  scale_fill_manual(values = c("Night" = "black", "Day" = "yellow")); Inverse_simp_diel_dcm
ggsave("inv_simp_dcm_diel.png", plot = Inverse_simp_diel_dcm, width = 6, height = 4)
dev.off()














