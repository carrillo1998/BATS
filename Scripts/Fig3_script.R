library(tidyverse)
library(ggalluvial)

############ Figure 3 scripts ############
dcm_hosts = read_csv("../DCM_hosts.csv")
sur_hosts = read_csv("../SUR_hosts.csv")
df_host_ab = left_join(dcm_hosts, sur_hosts) %>% mutate(across(where(is.numeric), ~replace_na(., 0)))
color_gen = read_csv("../arche_genera_color.csv")

df_abund_depth = df_host_ab %>%
                    pivot_longer(-host, names_to = "sample_name", values_to = "abundance") %>%
                    mutate(depth = case_when(str_detect(sample_name, "_D") ~ "DCM",
                                                    TRUE ~ "SUR")) %>%
                    group_by(host, depth) %>%
                    # summarise and create a column with the relative abundance
                    dplyr::summarise(Abundance = sum(abundance)) %>%
                    ungroup() %>%
                    pivot_wider(id_cols = host, names_from = depth, values_from = Abundance, values_fill = 0) %>%
                    # reduce function combines (reduces) all of the elements into a single object `x` = sum (e.g. reduce(c(1, 2, 3), `+`))
                    mutate(total = rowSums(across(where(is.numeric)))) %>% 
                    arrange(desc(total))

host_abund_depth = df_abund_depth %>%
                        mutate(DCM_rel_ab = 100*(DCM/sum(DCM)),
                                SUR_rel_ab =  100*(SUR/sum(SUR)),
                                total_rel= 100*(total/sum(total))) %>%
                        arrange(desc(DCM_rel_ab)) %>%
                        mutate(host = case_when(if_else(total_rel >=1, TRUE, FALSE) ~ host,
                                                TRUE ~ "others")) 
sum(host_abund_depth$DCM)
sum(host_abund_depth$SUR)

host_abund_depth = host_abund_depth %>% 
                        select(host, DCM_rel_ab, SUR_rel_ab) %>%
                        group_by(host) %>%
                        summarise(DCM_rel_ab = sum(DCM_rel_ab),
                                SUR_rel_ab = sum(SUR_rel_ab))

col_pathway = host_abund_depth %>% 
    # mutate(family = replace_na(family, "others")) %>%
    left_join(color_gen, by = c("host" = "genus")) %>%
    distinct(host, color) %>% 
    deframe()
column_colors = c(col_pathway, depth_color)

p = host_abund_depth %>%
      pivot_longer(cols = c(DCM_rel_ab, SUR_rel_ab), names_to = "Depth", values_to = "rel_abund") %>%
        mutate(Depth = ifelse(Depth == "DCM_rel_ab", "DCM", "SUR")) %>%
        mutate(Depth = factor(Depth, levels = c("SUR", "DCM"))) %>%
        mutate(host_id = paste0(Depth,"_", host)) %>%
        group_by(Depth) %>%
        arrange(desc(rel_abund), .by_group = TRUE) %>%
        mutate(host_id = factor(host_id, levels = unique(host_id))) %>%
        ggplot(aes(x = Depth, stratum = host_id, alluvium = host, y = rel_abund, fill = host, label = host)) +
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
