library(tidyverse)
library(ggalluvial)


######################################
########### Archetypes ################
#########################################
votu_cluster = read_csv("./data/arche_votu.csv")
votu_cluster_host = votu_cluster %>% 
                          left_join(df %>% select(Virus, taxonomy_g)) %>%
                          drop_na(cluster) %>% 
                          distinct(Virus, cluster, taxonomy_g)

arch_cluster_host = df %>% 
                        distinct(Virus, .keep_all = TRUE) %>%
                        select(Virus, where(is.numeric)) %>%
                        left_join(votu_cluster_host) %>%
                        drop_na(cluster)

arch_colors = c("Archetype_1" = "#A63232", "Archetype_2" = "#F2D22E", "Archetype_3" = "#50BFBF")



# alluvial
color_gen = read_csv("./data/arche_genera_color.csv")

df_abund_arch = arch_cluster_host %>%
                pivot_longer(-c(Virus, cluster, taxonomy_g), names_to = "sample_name", values_to = "abundance") %>%
                rename(host = taxonomy_g) %>%
                group_by(host, cluster) %>%
                # summarise and create a column with the relative abundance
                dplyr::summarise(Abundance = sum(abundance)) %>%
                ungroup() %>%
                mutate(Archetype = paste0("Archetype_", cluster)) %>%
                pivot_wider(id_cols = host, names_from = Archetype, values_from = Abundance, values_fill = 0) %>%
                # reduce function combines (reduces) all of the elements into a single object `x` = sum (e.g. reduce(c(1, 2, 3), `+`))
                mutate(total = rowSums(across(where(is.numeric)))) %>% 
                arrange(desc(total)) %>%
                filter(host != "Unknown Host")

sum(df_abund_arch$Archetype_1)
sum(df_abund_arch$Archetype_2)
sum(df_abund_arch$Archetype_3)

host_abund_arch = df_abund_arch %>%
                        mutate(arch1_rel_ab = 100*(Archetype_1/sum(Archetype_1)),
                                arch2_rel_ab =  100*(Archetype_2/sum(Archetype_2)),
                                arch3_rel_ab =  100*(Archetype_3/sum(Archetype_3)),
                                total_rel = 100*(total/sum(total))) %>%
                        mutate(host = case_when(if_else(total_rel >=1, TRUE, FALSE) ~ host,
                                                TRUE ~ "others")) %>% 
                        select(host, ends_with("_rel_ab")) %>%
                        group_by(host) %>%
                        summarise(arch1_rel_ab = sum(arch1_rel_ab),
                                arch2_rel_ab = sum(arch2_rel_ab) ,
                                arch3_rel_ab =  sum(arch3_rel_ab))

col_pathway = host_abund_arch %>% 
    # mutate(family = replace_na(family, "others")) %>%
    left_join(color_gen, by = c("host" = "genus")) %>%
    distinct(host, color) %>% 
    deframe()
column_colors = c(col_pathway, arch_colors)


p_arch = host_abund_arch %>%
        pivot_longer(cols = c(arch1_rel_ab, arch2_rel_ab ,arch3_rel_ab), names_to = "Cluster", values_to = "rel_abund") %>%
        mutate(Cluster = case_when(Cluster == "arch1_rel_ab" ~ "Archetype 1",
                                     Cluster == "arch2_rel_ab" ~ "Archetype 2",
                                     TRUE ~ "Archetype 3"),
                Cluster = factor(Cluster, levels= c("Archetype 1", "Archetype 3", "Archetype 2"))) %>%
        mutate(host_id = paste0(Cluster,"_", host),) %>%
        group_by(Cluster) %>%
        arrange(desc(rel_abund), .by_group = TRUE) %>%
        mutate(host_id = factor(host_id, levels = unique(host_id))) %>%
        ggplot(aes(x = Cluster, stratum = host_id, alluvium = host, y = rel_abund, fill = host, label = host)) +
                              # geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "gray30") +
                              stat_alluvium() +
                              geom_stratum() +
                              geom_text(stat = "alluvium", size = 3, color = "black") +
                              scale_fill_manual(values = column_colors) +
                              # geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 5) +
                              scale_x_discrete(expand = c(.1,.1)) +
                              theme_void() + 
                              theme(
                                  axis.text = element_text(size = 10),
                                      axis.text.y = element_blank(),
                                      axis.ticks = element_blank(),
                                      panel.grid = element_blank(),
                                  # axis.text.y = element_text(size = 0),
                                  legend.position = "bottom"
                                  )
