
library(cluster)
#install.packages("fpc")
library(fpc)
#install.packages("factoextra")
library(factoextra)
#install.packages("kohonen")
library(kohonen)
#install.packages("clustertend")
library(clustertend)
library(pracma)
library(gplots)
library(ggplot2)
library(tidyverse)
library(permute)
require(gridExtra)
library(grid)
library(fields)
#install.packages("ggforce")
library(ggforce)
library(vegan)


#From RAIN_detrend.R
surface_signif_diel <- read.csv("~/Documents/BATS R Studio/Outliers removed/Final_Scripts/import/surface_signif_diel.csv")
#From Fig2_script.R
vOTU_sur <- read.csv("~/Documents/BATS R Studio/Outliers removed/Final_Scripts/import/vOTU_sur.csv", row.names=1)


diel_sur_vOTU_g <- vOTU_sur[rownames(vOTU_sur) %in% surface_signif_diel$vOTU_names, ]
diel_sur_vOTU_g
diel_sur_vOTU_g <- diel_sur_vOTU_g[, !(colnames(diel_sur_vOTU_g) %in% c("T52_S", "T88_S"))]

# detrends and scales vOTUs abundance data that was considered to be diel
diel_vOTUs_sur_daniel_method <- t(scale(detrend(t(diel_sur_vOTU_g))))


####### diel vOTUs peak time #######
Surface_env_table <-read.table("Surface_Data.txt", header=T, stringsAsFactors = T, sep = '\t') #read sample data table
Surface_env_table <- Surface_env_table[!(Surface_env_table$Sample_id %in% c("T52_S", "T88_S")), ]
SUR_column_times <- Surface_env_table$Time

diel_sur_vOTU_g
diel_sur_vOTU_g_times <- diel_sur_vOTU_g
colnames(diel_sur_vOTU_g_times) <- SUR_column_times

time_map <- Surface_env_table$Time

time_labels <- paste0("T", seq(0, 112, 4))
time_labels <- time_labels[!time_labels %in% c("T52", "T88")]

time_mapping <- data.frame(Timepoint = time_labels, Time = time_map)

SUR_diel_peak_times <- diel_sur_vOTU_g

# Create a new dataframe with ranked values for each row
ranked_SUR_diel <- as.data.frame(t(apply(SUR_diel_peak_times, 1, function(x) rank(-x))))


# Assuming time_mapping is a dataframe with 'timepoint' and 'clock_time' columns
colnames(ranked_SUR_diel) <- Surface_env_table$Time

# Create a new dataframe to store the averaged values
averaged_SUR_diel <- as.data.frame(t(sapply(1:nrow(ranked_SUR_diel), function(i) {
  sapply(unique(colnames(ranked_SUR_diel)), function(time) {
    mean(as.numeric(ranked_SUR_diel[i, colnames(ranked_SUR_diel) == time]), na.rm = TRUE)
  })
})))

# Set the row names to be the same as the original dataframe
rownames(averaged_SUR_diel) <- rownames(ranked_SUR_diel)

# Identify rows with ties
# Create an empty list to store tie information
ties_list <- list()

# Loop through each row of the averaged dataframe
for (i in 1:nrow(averaged_SUR_diel)) {
  # Get the values from the current row
  row_values <- averaged_SUR_diel[i, ]
  
  # Find unique values and their frequencies
  freq_table <- table(row_values)
  
  # Check for ties (frequencies greater than 1)
  tied_values <- names(freq_table[freq_table > 1])
  
  if (length(tied_values) > 0) {
    # Find which columns have the tied values
    tied_columns <- colnames(averaged_SUR_diel)[row_values %in% tied_values]
    ties_list[[rownames(averaged_SUR_diel)[i]]] <- tied_columns
  }
}
print(ties_list)


# Create a vector to store the lowest average times
lowest_average_times_SUR_diel <- apply(averaged_SUR_diel, 1, function(x) {
  # Find the column name (time) corresponding to the minimum value
  colnames(averaged_SUR_diel)[which.min(x)]
})

# Add the result as a new column in the dataframe
averaged_SUR_diel$lowest_average_times_SUR_diel <- lowest_average_times_SUR_diel

table(averaged_SUR_diel$lowest_average_times_SUR_diel)






## Cluster validation statistics -- hopkins statistic identifies the degree to which
## modularity in the distance matrix exceeds that of uniformly distributed distances
## Values near 0.5 indicate uniform distribution, values near 1 indicate strong cluster tendency, values near 0 indicate regularly spaced (anti-cluster)
set.seed(123)
SUR_diel_h<-clustertend::hopkins(diel_vOTUs_sur_daniel_method,1000) ## clustertrend is 1-hopkins value. value is 0.7776
## Preliminarily calculate distance matrix to ease calculation
SUR_diel_dist<-dist(diel_vOTUs_sur_daniel_method)


## Calculate associated metrics for pam clustering algorithm
SUR_pam_asws<-c()
SUR_pam_chs<-c()
for(i in 1:8){
  SUR_pam_clust<-clara(diel_vOTUs_sur_daniel_method,i+1,metric="euclidean")
  SUR_pam_asws[i]<-summary(silhouette(SUR_pam_clust,full=TRUE))$avg.width
  SUR_pam_chs[i]<-calinhara(diel_vOTUs_sur_daniel_method,SUR_pam_clust$clustering)
  print(paste('completed cluster',i))
}

library(stats)
## Calculate associated metrics for hc clustering algorithm
SUR_ch_metrics<-c()
SUR_asw_metrics<-c()
for(i in 1:8){
  SUR_clust<-hcut(SUR_diel_dist,k=i+1,isdiss=TRUE)
  SUR_ch_metrics[i]<-calinhara(diel_vOTUs_sur_daniel_method,SUR_clust$cluster)
  SUR_asw_metrics[i]<-summary(silhouette(SUR_clust$cluster,SUR_diel_dist))$avg.width
  print(paste0('completed cluster ',i))
}

SUR_clust_3<-hcut(SUR_diel_dist,k=3,isdiss=TRUE, graph = TRUE)

## Calculate some SOMs and their associated metrics
SUR_x_grid_dim<-c(2,2,3,3,1,7,2,3,5,11,3,6,5,6,3,4)
SUR_y_grid_dim<-c(1,2,1,2,5,1,4,3,2,1,4,2,5,6,9,7)
SUR_ch_val<-c()
SUR_asw<-c()
for(i in 1:length(SUR_x_grid_dim)){
  set.seed(7583176)
  SUR_test_som<-som(as.matrix(diel_vOTUs_sur_daniel_method),grid=somgrid(SUR_x_grid_dim[i],SUR_y_grid_dim[i],'rectangular'))
  SUR_ch_val[i]<-calinhara(diel_vOTUs_sur_daniel_method,SUR_test_som$unit.classif)
  SUR_asw[i]<-summary(silhouette(SUR_test_som$unit.classif,SUR_diel_dist))$avg.width
  print(paste0('completed grid ',i))
}
SUR_som_nclust<-SUR_x_grid_dim*SUR_y_grid_dim
SUR_som_ch<-rep(NA,16)
SUR_som_asw<-rep(NA,16)
SUR_som_ch[SUR_som_nclust]<-SUR_ch_val
SUR_som_asw[SUR_som_nclust]<-SUR_asw

## Visualizing preliminary ODI for all data to get sense of data structure (WARNING: THIS IS VERY SLOW)
## ODI: Ordered dissimilarity image, representation of the distance matrix with rows/columns sorted to 
## maximize modularity
SUR_auto_odi<-get_clust_tendency(diel_vOTUs_sur_daniel_method,10,
                             gradient=list(low="Darkblue",mid="White",high="Gold"))


## Based off of manually inspecting these, it appears the data could roughly correspond to
## either three or four partially overlapping clusters. We can compare the metrics for these
## amongst algorithms

## So let's compare amongst algorithms + metrics
SUR_algs<-rep(c('pam','hc','som'),each=8)
SUR_chs<-c(SUR_pam_chs,SUR_ch_metrics[1:8],SUR_som_ch[2:9])
SUR_asws<-c(SUR_pam_asws,SUR_asw_metrics[1:8],SUR_som_asw[2:9])
SUR_metric_frame<-data.frame(alg=SUR_algs,ch=SUR_chs,asw=SUR_asws,nclust=rep(2:9,3))
SUR_ch_plot<-ggplot(SUR_metric_frame,aes(x=nclust,col=alg,y=ch))+
  geom_point(size=4)+
  xlab('Number Clusters')+ylab('C-H Value')+theme_bw()+
  scale_color_manual(name='Algorithm',labels=c('HC','PAM','SOM'),
                     values=c('gold','red','navy'))+
  theme(text=element_text(size=16,face='bold'))

ggsave(plot=ch_plot,filename='../figures/ch_metrics.pdf',device='pdf')


SUR_asw_plot<-ggplot(SUR_metric_frame,aes(x=nclust,col=alg,y=asw))+geom_point(size=4)+
  xlab('Number Clusters')+ylab('Avg Silhouette Width Value')+theme_bw()+
  scale_color_manual(name='Algorithm',labels=c('HC','PAM','SOM'),
                     values=c('gold','red','navy'))+
  theme(text=element_text(size=16,face='bold'))+
  ggrepel::geom_label_repel(data=data.frame(),aes(x=3,y=0.185,label='Elbow'),
                            col='navy',nudge_y=0.1,nudge_x=0.25)
ggsave('../figures/avg_silhouette_metric.pdf',plot=asw_plot,device='pdf')

## Under the decision to compare between 3 and 4, an SOM with 3 clusters appears to be the
## most supported

## Recapitulating 3 cluster SOM for SI
set.seed(1674)
SUR_som3<-som(as.matrix(diel_vOTUs_sur_daniel_method),grid=somgrid(3,1,'hexagonal'))


## Calculating individual silhouette profile
SUR_somsil3<-silhouette(dist=SUR_diel_dist,x=SUR_som3$unit.classif)
table(SUR_som3$unit.classif)
## Rearranging distance matrix for generating ODI
arche3<-which(SUR_som3$unit.classif==3)
arche2<-which(SUR_som3$unit.classif==2)
arche1<-which(SUR_som3$unit.classif==1)
SUR_rearranged_dmat3<-as.matrix(SUR_diel_dist)[c(arche1,arche2,arche3),
                                               c(arche1,arche2,arche3)]


################### Make figures after clustering ###############

## Making codebook figure for conceptual demonstration
set.seed(89781)
SUR_big_som<-som(as.matrix(diel_vOTUs_sur_daniel_method),grid=somgrid(6,4,'hexagonal'))
SUR_big_som_hc<-hcut((object.distances(SUR_big_som,'codes')),k=4)
png('SUR_codebook_figure.png',
    units='in',
    width=6,
    height=6,
    res=400)
plot(SUR_big_som,type='codes',shape='straight',main='')
add.cluster.boundaries(SUR_big_som,SUR_big_som_hc$cluster)
dev.off()


pdf(file = 'Fig_6a2.pdf',width=6,height=6)
png('Fig_6a2.png',width=6,height=6,units='in',res=275)
image(SUR_rearranged_dmat3,col=colorRampPalette(c('DarkViolet','White','Orange'))(n=100),
      axes=FALSE,
      main='Distance Matrix',
      cex.main=2)
abline(h=length(arche1)/length(SUR_som3$unit.classif),lwd=2,lty=5)
abline(h=(length(arche1)+length(arche2))/length(SUR_som3$unit.classif),lwd=2,lty=5)
abline(v=length(arche1)/length(SUR_som3$unit.classif),lwd=2,lty=5)
abline(v=(length(arche1)+length(arche2))/length(SUR_som3$unit.classif),lwd=2,lty=5)
text(x=c(length(arche1)/(2*length(SUR_som3$unit.classif)),
         (length(arche1)+(0.5*length(arche2)))/length(SUR_som3$unit.classif),
         (length(SUR_som3$unit.classif)-(0.5*length(arche3)))/length(SUR_som3$unit.classif)),
     y=c(length(arche1)/(2*length(SUR_som3$unit.classif)),
         (length(arche1)+(0.5*length(arche2)))/length(SUR_som3$unit.classif),
         (length(SUR_som3$unit.classif)-(0.5*length(arche3)))/length(SUR_som3$unit.classif)),
     labels=c('Archetype 1','Archetype 2','Archetype 3'),cex=1.5
)

image.plot(SUR_rearranged_dmat3[seq(1,nrow(SUR_rearranged_dmat3),by=100),
                                seq(1,nrow(SUR_rearranged_dmat3),by=100)],col=colorRampPalette(c('DarkViolet','White','Orange'))(n=100),
           axes=FALSE,
           horizontal=TRUE,
           main='Distance Matrix',
           cex=10,legend.only=TRUE,
           axis.args=list(cex.axis=2))
dev.off()





########### Labeling clusters #################
sur_sig_names <- surface_signif_diel
sur_sig_names <- sur_sig_names %>% rename(virus = vOTU_names)

SUR_som3$data
arche3<-which(SUR_som3$unit.classif==3)
SUR_som3_arche3_names <- rownames(diel_vOTUs_sur_daniel_method)[arche3]
SUR_som3_arche3_signif_subset <- sur_sig_names[sur_sig_names$virus %in% SUR_som3_arche3_names, ]

arche1<-which(SUR_som3$unit.classif==1)
SUR_som3_arche1_names <- rownames(diel_vOTUs_sur_daniel_method)[arche1]
SUR_som3_arche1_signif_subset <- sur_sig_names[sur_sig_names$virus %in% SUR_som3_arche1_names, ]


arche2<-which(SUR_som3$unit.classif==2)
SUR_som3_arche2_names <- rownames(diel_vOTUs_sur_daniel_method)[arche2]
SUR_som3_arche2_signif_subset <- sur_sig_names[sur_sig_names$virus %in% SUR_som3_arche2_names, ]

SUR_som3_arche3_signif_subset <- SUR_som3_arche3_signif_subset[!duplicated(SUR_som3_arche3_signif_subset$virus),]
SUR_som3_arche1_signif_subset <- SUR_som3_arche1_signif_subset[!duplicated(SUR_som3_arche1_signif_subset$virus),]
SUR_som3_arche2_signif_subset <- SUR_som3_arche2_signif_subset[!duplicated(SUR_som3_arche2_signif_subset$virus),]


table(SUR_som3_arche3_signif_subset$phase)
table(SUR_som3_arche1_signif_subset$phase)
table(SUR_som3_arche2_signif_subset$phase)

averaged_SUR_diel_names <- rownames_to_column(averaged_SUR_diel, var = "virus")
SUR_som3_signif_subset_arche1 <- averaged_SUR_diel_names[averaged_SUR_diel_names$virus %in% SUR_som3_arche1_names, ]
SUR_som3_signif_subset_arche2 <- averaged_SUR_diel_names[averaged_SUR_diel_names$virus %in% SUR_som3_arche2_names, ]
SUR_som3_signif_subset_arche3 <- averaged_SUR_diel_names[averaged_SUR_diel_names$virus %in% SUR_som3_arche3_names, ]

table(SUR_som3_signif_subset_arche1$lowest_average_times_SUR_diel)
table(SUR_som3_signif_subset_arche2$lowest_average_times_SUR_diel)
table(SUR_som3_signif_subset_arche3$lowest_average_times_SUR_diel)


################## Second part with line graph ################

SUR_samplers<-SUR_som3$unit.classif
set.seed(98014)
SUR_clust_ex_mat<-matrix(data=NA,ncol=27)
colnames(SUR_clust_ex_mat)<-colnames(diel_vOTUs_sur_daniel_method)
for(i in 1:3){
  SUR_clust_sample<-sample(which(SUR_samplers==i),140,replace=FALSE)
  SUR_clust_timeseries<-diel_vOTUs_sur_daniel_method[SUR_clust_sample,]
  SUR_clust_ex_mat<-rbind(SUR_clust_ex_mat,SUR_clust_timeseries)
}
SUR_full_code_mat<-rbind(SUR_clust_ex_mat[-1,],as.matrix(SUR_som3$codes[[1]]))
SUR_full_code_mat <- as.data.frame(SUR_full_code_mat)

SUR_full_code_mat$T52_S <- NA
SUR_full_code_mat$T88_S <- NA
SUR_order <- c("T0_S", "T4_S", "T8_S", "T12_S", "T16_S", "T20_S", "T24_S", "T28_S", "T32_S", 
               "T36_S", "T40_S", "T44_S", "T48_S", "T52_S", "T56_S", "T60_S", "T64_S", "T68_S", 
               "T72_S", "T76_S", "T80_S", "T84_S", "T88_S", "T92_S", "T96_S", "T100_S", "T104_S", 
               "T108_S", "T112_S")
SUR_full_code_mat <- SUR_full_code_mat[, SUR_order]


SUR_code_frame<-reshape(SUR_full_code_mat,direction='long',varying=list(1:29)) %>%
  mutate(time_label=factor(x=time,
                           levels=as.character(1:29),
                           labels=colnames(SUR_full_code_mat)),
         clust_label=factor(x=rep(c(rep(1:3,each=140),1,2,3),29),
                            levels=c(1,2,3),
                            labels=c('Archetype 1','Archetype 2', 'Archetype 3')),
         archetype=rep(c(rep('no',140*3),rep('yes',3)),29))
plot_labels_sur <- Surface_env_table$Time

pdf(file = "Fig_6b.pdf", width = 3, height = 4)
ggplot()+
  geom_rect(data=data.frame(xmin=c(1.5,7.5,13.5,19.5,25.5),xmax=c(4.5,10.5,16.5,22.5,28.5),ymin=rep(-Inf,5),ymax=rep(Inf,5)),
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.5,fill='grey')+
  geom_line(data=SUR_code_frame,aes(x=time,y=`T0_S`,group=id),alpha=0.2)+
  geom_line(data=filter(SUR_code_frame,archetype=='yes'),aes(x=time,y=`T0_S`,group=id,col=clust_label),lwd=1.5)+
  facet_wrap(~clust_label,ncol=1,strip.position='left')+
  scale_x_continuous(breaks = seq(1.5, 29, by = 3), 
                     labels = c("6PM", "6AM", "6PM", "6AM", "6PM", "6AM", "6PM", "6AM", "6PM", "6AM"), expand = c(0,0))+
  scale_color_manual(labels=c('Archetype 1','Archetype 2', 'Archetype 3'),
                     values=rev(c("midnightblue", "gold", "darkviolet")))+
  theme_bw()+
  theme(strip.background=element_blank(),
        text=element_text(size=16),
        axis.title=element_blank(),
        axis.text.x=element_text(angle=60,vjust=0.5),
        panel.border=element_blank(),
        strip.placement='outside',
        panel.spacing.y=unit(2,'lines'),
        panel.grid=element_blank(),
        legend.position='none',
        strip.text=element_blank())
dev.off()
ggsave(filename='som_archetypes_nolabels.png',device='png',width=6,units='in')



################## Bar graph figure ###############
diel_sur_vOTU_g <- vOTU_sur[rownames(vOTU_sur) %in% surface_signif_diel$vOTU_names, ]

#Make df rownames into its own column
diel_sur_iphop_g<- tibble::rownames_to_column(diel_sur_vOTU_g, var = "Virus")

#From iPHoP_data_analysis.R script
iphop_g <- read.csv("~/Documents/BATS R Studio/Outliers removed/Final_Scripts/import/iphop_g.csv", row.names=1) 
diel_sur_iphop_g <- left_join(diel_sur_iphop_g, iphop_g, by = "Virus")

# Replace NA values in 'taxonomy_p' with "Unknown Host"
diel_sur_iphop_g$taxonomy_g[is.na(diel_sur_iphop_g$taxonomy_g)] <- "Unknown Host"

# Remove Unknown Host in 'taxonomy_p' "
diel_sur_iphop_g <- diel_sur_iphop_g[diel_sur_iphop_g$taxonomy_g != "Unknown Host", ]
diel_sur_iphop_g <- diel_sur_iphop_g[!duplicated(diel_sur_iphop_g$Virus),]
diel_sur_iphop_g

unique(diel_sur_iphop_g$taxonomy_g)

#From iPHoP_data_analysis.R script
filter_agg_sur_diel_g <- read.csv("~/Documents/BATS R Studio/Outliers removed/Final_Scripts/import/filter_agg_sur_diel_g.csv", row.names=1)
SUR_sum_taxa_diel <- as.data.frame(rowSums(filter_agg_sur_diel_g))
SUR_sum_taxa_diel <- rownames_to_column(SUR_sum_taxa_diel, var = "Host")
colnames(SUR_sum_taxa_diel) <- c("Host", "sum_abundance")

top_20_SUR_abundance <- SUR_sum_taxa_diel %>%
  arrange(desc(sum_abundance)) %>%
  slice_head(n = 20)

prochloroB_syneE <- SUR_sum_taxa_diel %>%
  filter(Host %in% c("Prochlorococcus_B", "Synechococcus_E"))

# Combine the top 20 and specific rows
Hosts_for_figure <- bind_rows(top_20_SUR_abundance, prochloroB_syneE) %>%
  distinct()

# Create a dataframe with the remaining rows
remaining_df <- anti_join(SUR_sum_taxa_diel, Hosts_for_figure, by = colnames(SUR_sum_taxa_diel))

iphop_table_2 <- iphop_table
iphop_table_2$phylum <- sapply(strsplit(iphop_table_2$Host.taxonomy, ";"), function(x) {
  # Find the element that starts with 'p__'
  phylum_element <- grep("^p__", x, value = TRUE)
  # Remove the 'p__' prefix to get just the phylum name
  if (length(phylum_element) > 0) {
    return(sub("p__", "", phylum_element))
  } else {
    return(NA)  # Handle cases where no phylum is found
  }
})
iphop_table_2$genera <- sapply(strsplit(iphop_table_2$Host.taxonomy, ";"), function(x) {
  # Find the element that starts with 'g__'
  genera_element <- grep("^g__", x, value = TRUE)
  # Remove the 'g__' prefix to get just the phylum name
  if (length(genera_element) > 0) {
    return(sub("g__", "", genera_element))
  } else {
    return(NA)  # Handle cases where no phylum is found
  }
})

remaining_turn_to_phyla <- iphop_table_2[iphop_table_2$genera %in% remaining_df$Host, ]

remaining_turn_to_phyla <- remaining_turn_to_phyla %>% distinct(genera, .keep_all = TRUE)
colnames(remaining_df) <- c("genera", "sum_abundance")
remaining_phyla_names <- merge(remaining_df, remaining_turn_to_phyla[, c("genera", "phylum")], by = "genera", all.x = TRUE)

remaining_phyla_names <- remaining_phyla_names[, c("sum_abundance", "phylum")]
colnames(remaining_phyla_names) <- c("sum_abundance", "Host")

remaining_phyla_names_agg <- remaining_phyla_names %>%
  group_by(Host) %>%
  summarize(sum_abundance = sum(sum_abundance, na.rm = TRUE))




Hosts_for_figure <- bind_rows(Hosts_for_figure, remaining_phyla_names_agg)

arche1<-which(SUR_som3$unit.classif==1)

SOM3_1<-which(SUR_som3$unit.classif==1)
som3_1_names <- rownames(diel_vOTUs_sur_daniel_method)[SOM3_1]
som3_1_subset <- sur_sig_names[sur_sig_names$virus %in% som3_1_names, ]

SOM3_2<-which(SUR_som3$unit.classif==2)
som3_2_names <- rownames(diel_vOTUs_sur_daniel_method)[SOM3_2]
som3_2_subset <- sur_sig_names[sur_sig_names$virus %in% som3_2_names, ]

SOM3_3<-which(SUR_som3$unit.classif==3)
som3_3_names <- rownames(diel_vOTUs_sur_daniel_method)[SOM3_3]
som3_3_subset <- sur_sig_names[sur_sig_names$virus %in% som3_3_names, ]

combined_df <- bind_rows(
  som3_1_subset %>% select(virus) %>% mutate(cluster = 1),
  som3_2_subset %>% select(virus) %>% mutate(cluster = 2),
  som3_3_subset %>% select(virus) %>% mutate(cluster = 3)
)
colnames(combined_df) <- c("Virus", "cluster")

SOM_SUR_host <- left_join(combined_df, iphop_g, by = "Virus")
SOM_SUR_host <- left_join(SOM_SUR_host, iphop_p, by = "Virus")


# Replace NA values in 'taxonomy_p' with "Unknown Host"
SOM_SUR_host$taxonomy_g[is.na(SOM_SUR_host$taxonomy_g)] <- "Unknown Host"
SOM_SUR_host$taxonomy_p[is.na(SOM_SUR_host$taxonomy_p)] <- "Unknown Host"


# Remove Unknown Host in 'taxonomy_p' "
SOM_SUR_host <- SOM_SUR_host[SOM_SUR_host$taxonomy_g != "Unknown Host", ]
SOM_SUR_host <- SOM_SUR_host[SOM_SUR_host$taxonomy_p != "Unknown Host", ]


SOM_SUR_host <- SOM_SUR_host[!duplicated(SOM_SUR_host$Virus),]

SOM_SUR_host_g <- SOM_SUR_host %>%
  filter(taxonomy_g %in% Hosts_for_figure$Host)

SOM_SUR_host_p <- SOM_SUR_host %>%
  filter(!(taxonomy_g %in% Hosts_for_figure$Host))

SOM_SUR_host_g <- SOM_SUR_host_g[, c("Virus", "cluster", "taxonomy_g")]
colnames(SOM_SUR_host_g) <- c("Virus", "cluster", "host")
SOM_SUR_host_p <- SOM_SUR_host_p[, c("Virus", "cluster", "taxonomy_p")]
colnames(SOM_SUR_host_p) <- c("Virus", "cluster", "host")

SOM_SUR_host <- bind_rows(SOM_SUR_host_g, SOM_SUR_host_p)


diel_sur_vOTU_g <- vOTU_sur[rownames(vOTU_sur) %in% surface_signif_diel$vOTU_names, ]

#Make df rownames into its own column
diel_sur_iphop_g<- tibble::rownames_to_column(diel_sur_vOTU_g, var = "Virus")

# Add abundance data
SOM_SUR_host_abund <- left_join(SOM_SUR_host, diel_sur_iphop_g, by = "Virus")

# Sum of rows
SOM_SUR_host_abund$sum_abundance <- rowSums(SOM_SUR_host_abund[, !(colnames(SOM_SUR_host_abund) %in% c("Virus", "cluster", "host"))])
SOM_SUR_host_abund <- SOM_SUR_host_abund[, c("cluster", "host", "sum_abundance")]

aggregated_SUR_som_abund <- SOM_SUR_host_abund %>%
  group_by(cluster, host) %>%
  summarize(sum_abundance = sum(sum_abundance, na.rm = TRUE))

aggregated_SUR_som_abund <- aggregated_SUR_som_abund %>%
  group_by(host) %>%
  mutate(proportion = sum_abundance / sum(sum_abundance)) %>%
  ungroup()

#### Log count of taxa
log_count <- aggregated_SUR_som_abund

log_count <- log_count %>%
  group_by(host) %>%
  summarize(sum_all = sum(sum_abundance, na.rm = TRUE))

log_count$logcount <- log10(log_count$sum_all +1)


aggregated_SUR_som_abund$host<-factor(aggregated_SUR_som_abund$host,
                                            levels=log_count$host[order(log_count$logcount,decreasing=TRUE)])
log_count$host<-factor(log_count$host,levels=log_count$host[order(log_count$logcount,decreasing=TRUE)])
SUR_Palette <- c("1" = "midnightblue", 
                 "2" = "gold", 
                 "3" = "darkviolet")
aggregated_SUR_som_abund$cluster <- as.factor(aggregated_SUR_som_abund$cluster)

SUR_p1<-ggplot(aggregated_SUR_som_abund, aes(x=host, y=proportion, fill=cluster)) +
  geom_bar(stat="identity", colour="black", size=.3)+
  scale_fill_manual(values= SUR_Palette)+
  theme_classic() +
  theme(text = element_text(size=20),
        axis.title.y=element_text(hjust=0.66,size=10),
        axis.title=element_blank())+
  scale_y_continuous(expand = c(0,0))+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position='bottom',
        strip.placement='outside',
        strip.background.x=element_blank(),
        strip.text=element_text(angle=90))+
  theme(strip.text=element_blank(),
        axis.text.x=element_blank())+
  guides(fill=FALSE); SUR_p1

pdf(file = "test_all.pdf", width = 30, height = 10)
SUR_p1_labels<-SUR_p1 +
  theme(axis.text.x=element_text(angle=90,hjust=0.95,size=14),
        axis.text.y=element_text(size=12,angle=60,vjust=0.5))+
  guides(fill=guide_legend());SUR_p1_labels
ggsave(p1_labels,filename='../figures/fingerprint_with_labels.pdf',device='pdf')
dev.off()

SUR_som_p2<-ggplot(log_count,aes(x=host, y=logcount))+
  geom_bar(stat="identity", colour="black", size=.3)+
  theme_classic() +
  scale_y_continuous(expand = c(0,0))+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),strip.text.x=element_blank())+
  theme(text = element_text(size=15))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())+
  labs(y=expression(Log[10]~Count))+
  theme(#axis.text=element_blank(),
    axis.title.y=element_text(size=14));SUR_som_p2

SUR_plot_ob<-rbind(ggplotGrob(SUR_som_p2),ggplotGrob(SUR_p1_labels),size='first')
SUR_actual_plots<-unique(SUR_plot_ob$layout$t[grep('panel',SUR_plot_ob$layout$name)])
SUR_plot_ob$heights[SUR_actual_plots[1]]<-unit(1,'null')
SUR_plot_ob$heights[SUR_actual_plots[2]]<-unit(4,'null')

pdf(file = "cluster_taxonomy.pdf", height = 10, width = 10)
grid.newpage()
grid.draw(SUR_plot_ob)
dev.off()






############### Circular NMDS figure #################
diel_sur_vOTU_g <- diel_sur_vOTU_g[, !(colnames(diel_sur_vOTU_g) %in% c("T52_S", "T88_S"))]

# detrends and scales vOTUs abundance data that was considered to be diel
diel_vOTUs_sur_daniel_method <- t(scale(detrend(t(diel_sur_vOTU_g))))

# keep only rows that have host information
top_20_SUR_abundance <- SUR_sum_taxa_diel %>%
  arrange(desc(sum_abundance)) %>%
  slice_head(n = 20)

prochloroB_syneE <- SUR_sum_taxa_diel %>%
  filter(Host %in% c("Prochlorococcus_B", "Synechococcus_E"))

# Combine the top 20 and specific rows
nmds_to_keep <- bind_rows(top_20_SUR_abundance, prochloroB_syneE) %>%
  distinct()
nmds_host <- left_join(combined_df, iphop_g, by = "Virus")

# Replace NA values in 'taxonomy_p' with "Unknown Host"
nmds_host$taxonomy_g[is.na(nmds_host$taxonomy_g)] <- "Unknown Host"


# Remove Unknown Host in 'taxonomy_p' "
nmds_host <- nmds_host[nmds_host$taxonomy_g != "Unknown Host", ]


nmds_host <- nmds_host[!duplicated(nmds_host$Virus),]

nmds_hosts_filtered <- nmds_host %>%
  filter(taxonomy_g %in% nmds_to_keep$Host)

diel_vOTUs_sur_daniel_method_1 <- rownames_to_column(as.data.frame(diel_vOTUs_sur_daniel_method, var = "Host"))
diel_vOTUs_sur_daniel_method_1 <- diel_vOTUs_sur_daniel_method_1 %>%
  filter(rowname %in% nmds_hosts_filtered$Virus)
diel_vOTUs_sur_daniel_method_1 <- column_to_rownames(diel_vOTUs_sur_daniel_method_1, var = "rowname")

diel_vOTUs_sur_daniel_method_df <- as.data.frame(diel_vOTUs_sur_daniel_method)

## Calculating Euclidean distance matrix
SUR_diel_dist <-dist(diel_vOTUs_sur_daniel_method_1)
SUR_diel_dist_all <- dist(diel_vOTUs_sur_daniel_method_df)
SUR_diel_full_nmds<-metaMDS(SUR_diel_dist)
SUR_diel_full_all_nmds<-metaMDS(SUR_diel_dist_all)


## Additionally convert data to ranks to identify peak rank expression 
SUR_diel_rank_mat<-t(apply(diel_vOTUs_sur_daniel_method_1,1,rank))
SUR_diel_rank_all_mat<-t(apply(diel_vOTUs_sur_daniel_method_df,1,rank))


id_highest_mean_rank<-function(v){
  twelveam<-mean(v[c(3,9,14,20,25)],na.rm=TRUE)
  fouram<-mean(v[c(4,10,15,21,26)],na.rm=TRUE)
  eightam<-mean(v[c(5,11,16,27)],na.rm=TRUE)
  twelvepm<-mean(v[c(6,12,17,22)],na.rm=TRUE)
  fourpm<-mean(v[c(1,7,13,18,23)],na.rm=TRUE)
  eightpm<-mean(v[c(2,8,19,24)],na.rm=TRUE)
  meanranks<-c(twelveam,fouram,eightam,twelvepm,fourpm,eightpm)
  maxmeanrank<-which(meanranks==max(meanranks))
  if(1 %in% maxmeanrank & 6 %in% maxmeanrank){
    output<-6.5
  } else output<-mean(maxmeanrank)
  return(output)
}
SUR_mean_ranks<-apply(SUR_diel_rank_mat,1,id_highest_mean_rank)
SUR_mean_all_ranks<-apply(SUR_diel_rank_all_mat,1,id_highest_mean_rank)

SUR_diel_full_nmds_points <- SUR_diel_full_nmds$points
SUR_diel_full_all_nmds_points <- SUR_diel_full_all_nmds$points

SUR_taxa_names <- rownames(SUR_diel_full_nmds_points)
combined_df_text_edit <- combined_df
combined_df_text_edit$Virus <- gsub("__", "_", combined_df_text_edit$Virus) 

SUR_all_frame<-data.frame(x=SUR_diel_full_nmds$points[,1],
                      y=SUR_diel_full_nmds$points[,2],
                      Virus=rownames(SUR_diel_full_nmds$points),
                      taxa=SUR_taxa_names,time_rank=SUR_mean_ranks)
SUR_all_data_no_taxa_frame<-data.frame(x=SUR_diel_full_all_nmds$points[,1],
                          y=SUR_diel_full_all_nmds$points[,2],
                          Virus=rownames(SUR_diel_full_all_nmds$points),
                          time_rank=SUR_mean_all_ranks)
SUR_all_data_no_taxa_frame$Virus <- gsub("__", "_", SUR_all_data_no_taxa_frame$Virus) 

SUR_all_frame_hosts <- left_join(SUR_all_frame, nmds_hosts_filtered, by = "Virus")
SUR_all_data_no_taxa_frame <- left_join(SUR_all_data_no_taxa_frame, combined_df_text_edit, by = "Virus")

SUR_all_frame_hosts$distance <- sqrt(SUR_all_frame_hosts$x^2 + SUR_all_frame_hosts$y^2)
radius <- mean(SUR_all_frame_hosts$distance)
SUR_all_data_no_taxa_frame$distance <- sqrt(SUR_all_data_no_taxa_frame$x^2 + SUR_all_data_no_taxa_frame$y^2)
SUR_all_data_no_taxa_frame_radius <- mean(SUR_all_data_no_taxa_frame$distance)

ggplot(filter(SUR_all_frame_hosts, !taxonomy_g %in% c("Pelagibacter", "Flavobacterium", "Pelagibacter_A", 
                                                      "Winogradskyella")))+
  geom_point(aes(x=x,y=y,col=factor(floor(time_rank))),size=2.5)+
#  ggforce::geom_circle(aes(x0=0,y0=0,r=0.536))+
  coord_fixed()+
  theme_bw()+
  facet_wrap(~taxonomy_g,nrow=3)+
  scale_color_discrete(name='Peak Time',labels=c('0000','0400','0800','2000'))+
  theme(axis.title=element_blank(),
        legend.position='bottom',
        strip.background=element_blank())

pdf(file = "SUR_vOTU_all_taxonomy_nmds.pdf", width = 8, height = 5)
ggplot(SUR_all_frame_hosts)+
  geom_point(aes(x=x,y=y,col=factor(floor(time_rank))),size=2.5)+
#  ggforce::geom_circle(aes(x0=0,y0=0,r=0.536))+
  coord_fixed()+
  scale_color_manual(
    name = 'Peak Time',
    values = c("black", "gray", "goldenrod", "gold", "salmon", "gray45"),
    labels = c('0000', '0400', '0800', '1200', '1600', '2000'))+
  theme_bw()+
  facet_wrap(~taxonomy_g,nrow=5)+
 # scale_color_discrete(name='Peak Time',labels=c('0000','0400','0800','1200', '2000'))+
  theme(axis.title=element_blank(),
        legend.position='bottom',
        strip.background=element_blank())
dev.off()

pdf(file = "arche_layout_nmds.pdf", width = 8, height = 5)
ggplot(SUR_all_data_no_taxa_frame)+
  geom_point(aes(x=x,y=y,col=factor(floor(cluster))),size=2.5)+
  #  ggforce::geom_circle(aes(x0=0,y0=0,r=0.536))+
  coord_fixed()+
  theme_bw()+
#  facet_wrap(~taxonomy_g,nrow=5)+
  scale_color_discrete(name='Cluster',labels=c('1','2','3'))+
  theme(axis.title=element_blank(),
        legend.position='bottom',
        strip.background=element_blank())
dev.off()


SUR_all_data_no_taxa_frame_mrpp_time <- mrpp(SUR_diel_dist_all, SUR_all_data_no_taxa_frame$time_rank, permutations = 999, distance = "euclidean"); SUR_all_data_no_taxa_frame_mrpp_time
SUR_all_data_no_taxa_frame_mrpp_cluster <- mrpp(SUR_diel_dist_all, SUR_all_data_no_taxa_frame$cluster, permutations = 999, distance = "euclidean"); SUR_all_data_no_taxa_frame_mrpp_cluster

SUR_all_data_no_taxa_frame_permanova <- adonis2(SUR_diel_dist_all ~ time_rank + cluster,
                            data = SUR_all_data_no_taxa_frame,
                            by = "margin",
                            permutations = 999,
                            method = "euclidean"); SUR_all_data_no_taxa_frame_permanova


pdf(file = "Fig_6d.pdf", width = 8, height = 5)
ggplot(SUR_all_data_no_taxa_frame)+
  geom_point(aes(x=x,y=y,col=factor(floor(time_rank)), shape = as.factor(cluster)),size=2.5)+
  #  ggforce::geom_circle(aes(x0=0,y0=0,r=0.536))+
  coord_fixed()+
  theme_bw()+
  scale_shape_manual(
    name = 'cluster',
    values = c(16, 17, 18)  # You can customize the shape values as needed
  ) +
  scale_color_manual(
    name = 'Peak Time',
    values = c("black", "gray", "goldenrod", "gold", "salmon", "gray45"),
    labels = c('0000', '0400', '0800', '1200', '1600', '2000')
  ) +  #  facet_wrap(~taxonomy_g,nrow=5)+
 # scale_color_discrete(name='Peak Time',labels=c('0000','0400','0800','1200', '1600', '2000'))+
  theme(axis.title=element_blank(),
        legend.position='bottom',
        strip.background=element_blank())
dev.off()

ggsave(filename='Fig_6d.pdf',device='pdf',scale=2)

nmds_points <- full_nmds$points

##### Bar graph and pie chart of vOTU time data #####
SUR_vOTU_bar <- data.frame(Depth = c("SUR", "SUR"),
                      Time = c("Diel", "Non-diel"),
                      Viruses = c(3097, 25890))


SUR_vOTU_bar <- SUR_vOTU_bar %>%
  group_by(Depth) %>%
  mutate(Proportion = Viruses / sum(Viruses)) %>%
  ungroup()


SUR_vOTU_bar <- SUR_vOTU_bar %>%
  mutate(
    Depth = factor(Depth, levels = c("SUR")),
    Time = factor(Time, levels = c("Diel", "Non-diel"))
  )

pdf(file = "Fig_6a1.pdf", width = 3, height = 5)
ggplot(SUR_vOTU_bar, aes(x = Depth, y = Proportion, fill = Time)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Proportional Stacked Barplot of DCM and SUR",
    x = "Depth",
    y = "Proportion of Diel Viral Population Dominance"
  ) +
  scale_fill_manual(values = c("Diel" = "maroon", "Non-diel" = "black")) +
  theme_minimal()
dev.off()


########## pie chart
SUR_pie_diel <- data.frame(
  Time = c("0000", "0400", "0800", "1200", "1600", "2000"),
  Viruses = c(495, 1336, 78, 96, 9, 1083)
)

SUR_pie_diel <- SUR_pie_diel %>%
  arrange(desc(Time)) %>%
  mutate(prop = Viruses/ sum(SUR_pie_diel$Viruses) *100) %>%
  mutate(ypos = cumsum(prop) - 0.5*prop)


pdf(file = "SUR_pie_diel.pdf", width = 5, height = 5)
ggplot(SUR_pie_diel, aes(x = "", y = prop, fill = Time)) +
  geom_bar(stat = "identity", width = 1, color = "white")+
  coord_polar("y", start = 0) +
  theme_void()+
  theme(legend.position = "right") +
  geom_text(aes(y = ypos, label = Time), color = "white", size = 6)+
  scale_fill_manual(values = c("black", "gray", "goldenrod", "gold", "salmon", "gray45"))
dev.off()


write.csv(SUR_som3_arche1_signif_subset, file = "SUR_som3_arche1_signif_subset.csv", row.names = TRUE)
write.csv(SUR_som3_arche2_signif_subset, file = "SUR_som3_arche2_signif_subset.csv", row.names = TRUE)
write.csv(SUR_som3_arche3_signif_subset, file = "SUR_som3_arche3_signif_subset.csv", row.names = TRUE)

######## Archetype beta-dispersion ######

write.csv(SUR_all_data_no_taxa_frame, file = "Supp_table_14.csv", row.names = FALSE)
SUR_all_data_no_taxa_frame

# Ordered by time_rank
sur_data_time_order <- SUR_all_data_no_taxa_frame[order(SUR_all_data_no_taxa_frame$time_rank), ]

# Ordered by cluster
sur_data_time_cluster <- SUR_all_data_no_taxa_frame[order(SUR_all_data_no_taxa_frame$cluster), ]


#### Reorder df for later ease
#time
diel_vOTUs_sur_daniel_method_df_new <- as.data.frame(diel_vOTUs_sur_daniel_method_df)
diel_vOTUs_sur_daniel_method_df_new <- rownames_to_column(diel_vOTUs_sur_daniel_method_df_new)
diel_vOTUs_sur_daniel_method_df_new$rowname <- gsub("__", "_", diel_vOTUs_sur_daniel_method_df_new$rowname)
diel_vOTUs_sur_daniel_method_df_new <- column_to_rownames(diel_vOTUs_sur_daniel_method_df_new, var = "rowname")

diel_sur_reorder_time <- diel_vOTUs_sur_daniel_method_df_new[sur_data_time_order$Virus, , drop = FALSE]


#cluster
diel_sur_reorder_cluster <- diel_vOTUs_sur_daniel_method_df_new[sur_data_time_cluster$Virus, , drop = FALSE]

SUR_diel_dist_all_time <- dist(diel_sur_reorder_time)
SUR_diel_dist_all_cluster <- dist(diel_sur_reorder_cluster)

groups_time <- factor(c(rep(1,416), rep(1.5,25), rep(2,1217), rep(2.5,1), rep(3, 71), rep(3.5,1), rep(4, 108), rep(5, 7), rep(6,1241), rep(6.5, 10)), 
                      labels = c("0000","0200", "0400", "0600", "0800", "1000", "1200", "1600", "2000", "2200"))

groups_cluster <- factor(c(rep(1,943), rep(2,140), rep(3, 2014)), labels = c("Archetype 1","Archetype 2", "Archetype 3"))
groups_cluster <- factor(c(rep(1,1964), rep(2,148), rep(3, 985)), labels = c("Archetype 1","Archetype 2", "Archetype 3"))


## Calculate multivariate dispersions
SUR_time_betadisper <- betadisper(SUR_diel_dist_all_time, groups_time)
SUR_time_betadisper

SUR_cluster_betadisper <- betadisper(SUR_diel_dist_all_cluster, groups_cluster)
SUR_cluster_betadisper

## Perform test
anova(SUR_time_betadisper)

anova(SUR_cluster_betadisper)

## Permutation test for F
permutest(SUR_time_betadisper, pairwise = TRUE, permutations = 99)

permutest(SUR_cluster_betadisper, pairwise = TRUE, permutations = 99)

## Tukey's Honest Significant Differences
(SUR_time_betadisper.HSD <- TukeyHSD(SUR_time_betadisper))
plot(SUR_time_betadisper.HSD)


(SUR_cluster_betadisper.HSD <- TukeyHSD(SUR_cluster_betadisper))
plot(SUR_cluster_betadisper.HSD)
## Plot the groups and distances to centroids on the
## first two PCoA axes
plot(SUR_time_betadisper)

pdf(file = "Supp_Fig4.pdf", width =5, height =5)
plot(SUR_cluster_betadisper, col = c("midnightblue", "gold", "purple"))
dev.off()

## with data ellipses instead of hulls
plot(SUR_time_betadisper, ellipse = TRUE, hull = FALSE) # 1 sd data ellipse
plot(SUR_time_betadisper, ellipse = TRUE, hull = FALSE, conf = 0.90) # 90% data ellipse


plot(SUR_cluster_betadisper, ellipse = TRUE, hull = FALSE) # 1 sd data ellipse
plot(SUR_cluster_betadisper, ellipse = TRUE, hull = FALSE, conf = 0.90) # 90% data ellipse

# plot with manual colour specification
my_cols <- c("#1b9e77", "#7570b3")
plot(SUR_time_betadisper, col = my_cols, pch = c(16,17), cex = 1.1)

my_cols <- c("#1b9e77", "#7570b3")
plot(SUR_cluster_betadisper, col = my_cols, pch = c(16,17), cex = 1.1)

## Draw a boxplot of the distances to centroid for each group
boxplot(SUR_time_betadisper)


boxplot(SUR_cluster_betadisper)



############# Archetype figures ##########
times <- c("0000", "0400", "0800", "1200", "1600", "2000")
scale_fill_manual(values = c("black", "gray", "goldenrod", "gold", "salmon", "gray45"))

arche_data <- select(SUR_all_data_no_taxa_frame, Virus, cluster)
write.csv(arche_data, file = "arche_votu.csv", row.names = FALSE)


SUR_all_data_no_taxa_frame
arche1_time <- subset(SUR_all_data_no_taxa_frame, cluster == 1)
arche2_time <- subset(SUR_all_data_no_taxa_frame, cluster == 2)
arche3_time <- subset(SUR_all_data_no_taxa_frame, cluster == 3)

archetype_data <- data.frame(archetype = rep(c("Archetype 1", "Archetype 2", "Archetype 3"), each = 6),
                             time = rep(times, 3),
                             vOTUs = c(329, 647, 2, 6, 3, 977, 0, 18, 28, 98, 4, 0, 112, 553, 42, 4, 0, 274))


# Compute proportion of vOTUs per time within each archetype
prop_archetype <- archetype_data %>%
  group_by(archetype) %>%
  mutate(prop = vOTUs / sum(vOTUs)) %>%
  ungroup()

# Plot 1: stacked bar graph with proportions
pdf(file = "Fig_6c.pdf", width = 6, height = 5)
archetype_bar_stacked <- ggplot(prop_archetype, aes(x = archetype, y = prop, fill = time)) +
  geom_bar(stat = "identity") +
  labs(title = "Proportion of vOTUs by Time in Each Archetype",
       x = "Archetype", y = "Proportion") +
  scale_fill_manual(values = c("black", "gray", "goldenrod", "gold", "salmon", "gray45"))+
  theme_minimal(); archetype_bar_stacked
dev.off()

archetype_bar_unstacked <- ggplot(prop_archetype, aes(x = time, y = vOTUs, fill = time)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ archetype, ncol = 1) +
  labs(title = "Proportion of vOTUs by Time in Each Archetype",
       x = "Archetype", y = "Proportion") +
  scale_fill_manual(values = c("black", "gray", "goldenrod", "gold", "salmon", "gray45"))+
  theme_minimal(); archetype_bar_unstacked



