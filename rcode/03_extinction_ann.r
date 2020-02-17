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

# plot tiled facets
ggplot(data)+
  geom_tile(aes(x = R_new, y = P_new, fill = gen_extinct),
            size = 0.1)+
  geom_tile(aes(x = R, y = P), fill = NA, col = "red", lwd = 1)+
  facet_grid(P~R, labeller = label_both, 
             switch = "both", as.table = FALSE)+
  scale_fill_viridis_c(option = "C")+
  coord_fixed(ratio = 5, expand = F)+
  scale_y_continuous(breaks = c(0,0.5,1))+
  scale_x_continuous(breaks = c(0:5))+
  theme_few()+
  labs(x = NULL, y = NULL)+
  theme(strip.placement = "outside",
        legend.position = "none",
        strip.background = element_rect(fill = "grey90", 
                                        colour = "black"))
  