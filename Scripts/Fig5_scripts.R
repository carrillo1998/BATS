library(tidyverse)
library(ggalluvial)


########### Fig 5 ##############


###############################
###### Diel vs no Diel ######
##############################
# only vOTU from the SUR since the Diel and non-diel were found only at the sruface
df_sur = read_csv("./data/SUR_vOTUs.csv") %>% rename(Virus = x)
# get the diel vOTUs
arch_votu = read_csv("./data/SUR_data_new.csv")
diel_votu = arch_votu %>% drop_na(archetype) %>% distinct(vOTU) %>% mutate(Group = "Diel")
no_diel_otu = df_sur %>% 
                mutate(Virus = str_replace(Virus, "\\|\\|", "_")) %>%
                filter(!Virus %in% diel_votu$vOTU) %>% distinct(Virus) %>%
                mutate(Group = "no Diel") %>%
                rename(vOTU = Virus)

df_ren = df %>% 
  mutate(Virus = str_replace(Virus, "\\|\\|", "_")) %>% 
  rename(vOTU = Virus)

df_diel = bind_rows(diel_votu, no_diel_otu) %>%
                    left_join(df_ren) %>%
                    distinct(vOTU, .keep_all = TRUE)


color_gen = read_csv("./arche_genera_color.csv")

df_abund_diel = df_diel %>%
                select(vOTU, Group, taxonomy_g, where(is.numeric)) %>%
                pivot_longer(-c(vOTU, Group, taxonomy_g), names_to = "sample_name", values_to = "abundance") %>%
                rename(host = taxonomy_g) %>%
                group_by(host, Group) %>%
                # summarise and create a column with the relative abundance
                dplyr::summarise(Abundance = sum(abundance)) %>%
                ungroup() %>%
                pivot_wider(id_cols = host, names_from = Group, values_from = Abundance, values_fill = 0) %>%
                # reduce function combines (reduces) all of the elements into a single object `x` = sum (e.g. reduce(c(1, 2, 3), `+`))
                mutate(total = rowSums(across(where(is.numeric)))) %>% 
                arrange(desc(total)) %>%
                filter(host != "Unknown Host")

sum(df_abund_diel$Diel)
sum(df_abund_diel$`no Diel`)

host_abund_diel = df_abund_diel %>%
                        mutate(diel_rel_ab = 100*(Diel/sum(Diel)),
                        no_diel_rel_ab =  100*(`no Diel`/sum(`no Diel`)),
                        total_rel = 100*(total/sum(total))) %>%
                        mutate(host = case_when(if_else(total_rel >=1, TRUE, FALSE) ~ host,
                                                TRUE ~ "others")) %>% 
                        select(host, diel_rel_ab, no_diel_rel_ab) %>%
                        group_by(host) %>%
                        summarise(diel_rel_ab  = sum(diel_rel_ab ),
                              no_diel_rel_ab  = sum(no_diel_rel_ab ))

col_pathway = host_abund_diel %>% 
    # mutate(family = replace_na(family, "others")) %>%
    left_join(color_gen, by = c("host" = "genus")) %>%
    distinct(host, color) %>% 
    deframe()
column_colors = c(col_pathway, Diel = "#b03061" , `no Diel` = "#00000052")

p_diel = host_abund_diel %>%
        pivot_longer(cols = c(diel_rel_ab, no_diel_rel_ab), names_to = "Group", values_to = "rel_abund") %>%
        mutate(Group = ifelse(Group == "diel_rel_ab", "Diel", "No Diel")) %>%
        mutate(host_id = paste0(Group,"_", host)) %>%
        group_by(Group) %>%
        arrange(desc(rel_abund), .by_group = TRUE) %>%
        mutate(host_id = factor(host_id, levels = unique(host_id))) %>%
        ggplot(aes(x = Group, stratum = host_id, alluvium = host, y = rel_abund, fill = host, label = host)) +
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
