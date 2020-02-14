#### code for extinction zones ####

# load libs
# data libs
library(data.table)
library(stringr)
library(dplyr)
library(purrr)
library(tidyr)

# plotting libs
library(ggplot2)
library(viridis)
library(ggthemes)

#### load data ####
data <- fread("data/extinction_data.csv")

# potential changes
tot_data <- expand(data, R, P, R_new, P_new) %>% 
  left_join(data)

# plot tiled facets
ggplot(tot_data)+
  geom_tile(aes(x = R_new, y = P_new, fill = gen_extinct),
            size = 0.1)+
  geom_tile(aes(x = R, y = P), fill = NA, col = "red", lwd = 1)+
  facet_grid(P~R, labeller = label_both, 
             switch = "both", as.table = FALSE)+
  scale_fill_viridis_c(option = "E",
                       na.value = viridis(1, begin = 1))+
  coord_cartesian(expand = F)+
  scale_y_continuous(breaks = c(0,0.5,1))+
  scale_x_continuous(breaks = c(0:5))+
  theme_few()+
  labs(x = NULL, y = NULL)+
  theme(strip.placement = "outside",
        legend.position = "none", 
        strip.background = element_rect(fill = "grey90", 
                                        colour = "black"))
  