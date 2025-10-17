

################################## Pie Charts and Log2Fold Change Subfigures ############################
##### Use in conjunction with outputs from iPHoP_data_analysis.R

######## Figure 3b #####
# Vector of target rownames
target_rownames <- c("Pelagibacter", "AG-337-I02", "Pelagibacter_A", "Hyphobacterium", "Winogradskyella",
                     "Prochlorococcus_A", "Joostella", "ARS21", "Croceibacter", "GCA-2704625", "Actinomarina",
                     "TMED112", "MED-G52", "TMED25", "Aquimarina", "TMED108", "Flavobacterium", "Sediminibacterium",
                     "Prochlorococcus_B")
aggregate_iphop_sur_g #From iPHoP_data_analysis.R
sur_iphop <- column_to_rownames(aggregate_iphop_sur_g, var = "taxonomy_g")
sur_iphop <- as.data.frame(rowMeans(sur_iphop))
sur_iphop <- rownames_to_column(sur_iphop, var = "Host")
filtered_sur <- sur_iphop[sur_iphop$Host %in% target_rownames, ]
colnames(filtered_sur) <- c("Host", "Mean")

aggregate_iphop_dcm_g
dcm_iphop <- column_to_rownames(aggregate_iphop_dcm_g, var = "taxonomy_g")
dcm_iphop <- as.data.frame(rowMeans(dcm_iphop))
dcm_iphop <- rownames_to_column(dcm_iphop, var = "Host")
filtered_dcm <- dcm_iphop[dcm_iphop$Host %in% target_rownames, ]
colnames(filtered_dcm) <- c("Host", "Mean")


merged_depth <- merge(filtered_sur, filtered_dcm, by = "Host", suffixes = c("_sur", "_dcm"))
merged_depth$log2FC <- log2(merged_depth$Mean_sur / merged_depth$Mean_dcm)

pdf(file = "Fig3b.pdf", width = 6, height = 10);
ggplot(merged_depth, aes(x = reorder(Host, log2FC), y = log2FC)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    x = "Host",
    y = "Log2 Fold Change (SUR/DCM)",
    title = "Log2 Fold Change of Host Abundance (SUR vs DCM)"
  ) +
  theme_minimal()
dev.off()





########### Pie charts for Figures 3, 5, and 7 ######


####### SUR

SUR_hosts_pie <- data.frame(
  Group = c("Predictions", "No Predictions"),
  Viruses = c(7952, 21035)
)

SUR_hosts_pie <- SUR_hosts_pie %>%
  arrange(desc(Group)) %>%
  mutate(prop = Viruses/ sum(SUR_hosts_pie$Viruses) *100) %>%
  mutate(ypos = cumsum(prop) - 0.5*prop)

pdf(file = "SUR_hosts_pie.pdf", width = 5, height = 5)
ggplot(SUR_hosts_pie, aes(x = "", y = prop, fill = Group)) +
  geom_bar(stat = "identity", width = 1, color = "white")+
  coord_polar("y", start = 0) +
  theme_void()+
  theme(legend.position = "right") +
  geom_text(aes(y = ypos, label = Group), color = "black", size = 6)+
  scale_fill_manual(values = c("gray", "cyan3"))
dev.off()


###### DCM

DCM_hosts_pie <- data.frame(
  Group = c("Predictions", "No Predictions"),
  Viruses = c(11813, 16419)
)

DCM_hosts_pie <- DCM_hosts_pie %>%
  arrange(desc(Group)) %>%
  mutate(prop = Viruses/ sum(DCM_hosts_pie$Viruses) *100) %>%
  mutate(ypos = cumsum(prop) - 0.5*prop)

pdf(file = "DCM_hosts_pie.pdf", width = 5, height = 5)
ggplot(DCM_hosts_pie, aes(x = "", y = prop, fill = Group)) +
  geom_bar(stat = "identity", width = 1, color = "white")+
  coord_polar("y", start = 0) +
  theme_void()+
  theme(legend.position = "right") +
  geom_text(aes(y = ypos, label = Group), color = "black", size = 6)+
  scale_fill_manual(values = c("gray", "salmon"))
dev.off()


##### Diel SUR

SUR_diel_pie <- data.frame(
  Group = c("Predictions", "No Predictions"),
  Viruses = c(873, 2224)
)
SUR_nondiel_pie <- data.frame(
  Group = c("Predictions", "No Predictions"),
  Viruses = c(7079, 18811)
)

SUR_diel_pie <- SUR_diel_pie %>%
  arrange(desc(Group)) %>%
  mutate(prop = Viruses/ sum(SUR_diel_pie$Viruses) *100) %>%
  mutate(ypos = cumsum(prop) - 0.5*prop)
SUR_nondiel_pie <- SUR_nondiel_pie %>%
  arrange(desc(Group)) %>%
  mutate(prop = Viruses/ sum(SUR_nondiel_pie$Viruses) *100) %>%
  mutate(ypos = cumsum(prop) - 0.5*prop)

pdf(file = "SUR_diel_pie.pdf", width = 5, height = 5)
ggplot(SUR_diel_pie, aes(x = "", y = prop, fill = Group)) +
  geom_bar(stat = "identity", width = 1, color = "white")+
  coord_polar("y", start = 0) +
  theme_void()+
  theme(legend.position = "right") +
  geom_text(aes(y = ypos, label = Group), color = "black", size = 6)+
  scale_fill_manual(values = c("gray", "maroon"))
dev.off()

pdf(file = "SUR_nondiel_pie.pdf", width = 5, height = 5)
ggplot(SUR_nondiel_pie, aes(x = "", y = prop, fill = Group)) +
  geom_bar(stat = "identity", width = 1, color = "white")+
  coord_polar("y", start = 0) +
  theme_void()+
  theme(legend.position = "right") +
  geom_text(aes(y = ypos, label = Group), color = "black", size = 6)+
  scale_fill_manual(values = c("gray", "cyan3"))
dev.off()


###### arhcetypes

arche1_pie <- data.frame(
  Group = c("Predictions", "No Predictions"),
  Viruses = c(613, 1351)
)
arche2_pie <- data.frame(
  Group = c("Predictions", "No Predictions"),
  Viruses = c(22, 126)
)
arche3_pie <- data.frame(
  Group = c("Predictions", "No Predictions"),
  Viruses = c(238, 747)
)

arche1_pie <- arche1_pie %>%
  arrange(desc(Group)) %>%
  mutate(prop = Viruses/ sum(arche1_pie$Viruses) *100) %>%
  mutate(ypos = cumsum(prop) - 0.5*prop)
arche2_pie <- arche2_pie %>%
  arrange(desc(Group)) %>%
  mutate(prop = Viruses/ sum(arche2_pie$Viruses) *100) %>%
  mutate(ypos = cumsum(prop) - 0.5*prop)
arche3_pie <- arche3_pie %>%
  arrange(desc(Group)) %>%
  mutate(prop = Viruses/ sum(arche3_pie$Viruses) *100) %>%
  mutate(ypos = cumsum(prop) - 0.5*prop)

pdf(file = "arche1_pie.pdf", width = 5, height = 5)
ggplot(arche1_pie, aes(x = "", y = prop, fill = Group)) +
  geom_bar(stat = "identity", width = 1, color = "white")+
  coord_polar("y", start = 0) +
  theme_void()+
  theme(legend.position = "right") +
  geom_text(aes(y = ypos, label = Group), color = "black", size = 6)+
  scale_fill_manual(values = c("gray", "midnightblue"))
dev.off()

pdf(file = "arche2_pie.pdf", width = 5, height = 5)
ggplot(arche2_pie, aes(x = "", y = prop, fill = Group)) +
  geom_bar(stat = "identity", width = 1, color = "white")+
  coord_polar("y", start = 0) +
  theme_void()+
  theme(legend.position = "right") +
  geom_text(aes(y = ypos, label = Group), color = "black", size = 6)+
  scale_fill_manual(values = c("gray", "yellow"))
dev.off()

pdf(file = "arche3_pie.pdf", width = 5, height = 5)
ggplot(arche3_pie, aes(x = "", y = prop, fill = Group)) +
  geom_bar(stat = "identity", width = 1, color = "white")+
  coord_polar("y", start = 0) +
  theme_void()+
  theme(legend.position = "right") +
  geom_text(aes(y = ypos, label = Group), color = "black", size = 6)+
  scale_fill_manual(values = c("gray", "purple"))
dev.off()




















