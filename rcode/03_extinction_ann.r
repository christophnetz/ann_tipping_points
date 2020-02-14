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

#### load data ####
data <- fread("data/extinction_data.csv")

# plot tiled facets
ggplot(data)+
  geom_tile(aes(x = R_new, y = P_new, fill = gen_extinct))+
  facet_grid(P~R, labeller = label_both, 
             switch = "both", as.table = FALSE)+
  scale_fill_viridis_c()+
  coord_cartesian(expand = FALSE)+
  theme_bw()+
  theme(strip.placement = "outside")
