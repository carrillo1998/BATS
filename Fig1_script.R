
library(tidyverse)
library(dplyr)
library(patchwork)
library(hrbrthemes)
library(vegan)
library(ggplot2)
library(pheatmap)
library(tidyr)
library(rstatix)
library(readxl)

############## Figure 1 Script ###########

#### Fig. 1a
### Latitude and Longitude Line graph 
# map of coordinates #
invirt_2019_ctd <- read.csv("~/Documents/BATS R Studio/Outliers removed/Scripts and background data/Daniel_Muratore_invirt_2019_ctd.csv")
env_map <- unique(Daniel_Muratore_invirt_2019_ctd[, c("cast", "binned_lat", "binned_lon")])

# Example: Add a row
C1 <- data.frame( #cast1 missing from dataset
  cast = 1,
  binned_lat = 31.68000,
  binned_lon = -64.17472
)

# Add to df
env_map <- rbind(env_map, C1)
env_map <- env_map[order(env_map$cast), ]


pdf(file = "Fig_1a.pdf", width = 5, height = 5)
ggplot(env_map, aes(x = binned_lon, y = binned_lat)) +
  geom_path(color = "gray50", size = 1) +   # Gray line
  geom_point(color = "black", size = 3) +   # Black points
  labs(x = "Longitude", y = "Latitude",
       title = "Path of Casts over Latitude and Longitude") +
  theme_minimal()
dev.off()





###### Fig. 1b

depth_table <- read_excel("depth_profile_data_continuous.xlsx")
depth_table <- depth_table %>% mutate(across(where(is.character), as.numeric))


p1 <- ggplot(depth_table, aes(x=Temperature_ITS_1, y=depth)) +
  geom_point(color="#69b3a2", size=2) + 
  ggtitle("Temperature") +
  theme_minimal()
#ggsave("plot1.png", plot = p1, width =6, height = 6)

p2 <- ggplot(depth_table, aes(x=Fluoresence_2, y=depth)) +
  geom_point(color="grey", linewidth=2) +
  ggtitle("Chlorophyll") +
  theme_minimal()

plot(p1 +p2)
plot(p1)

# Value used to transform the data
coeff <- 49.17

pdf(file = "Fig_1b.pdf", width = 10, height = 15);
ggplot(depth_table, aes(y = depth)) +
  geom_point(size = 4, aes(x = Temperature_ITS_1 / coeff), color = "darkblue") + 
  geom_point(size = 4, aes(x = Fluoresence_2), color = "green") +
  scale_fill_gradient(low = "blue", high = "lightblue", name = "Depth") +
  theme(text = element_text(size = 25)) +
  theme_minimal() +
  scale_x_continuous(
    name = "Chlorophyll(mmg/mm)",
    sec.axis = sec_axis(~.*coeff, name = "Temperature °C")
  ) +
  scale_y_reverse()
dev.off()













