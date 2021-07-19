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
library(scico)
library(ggthemes)

setwd("C:/Users/P285100/Desktop/Botero/ann_tipping_points/ann_tipping_points/Run_9-7")

#### extinction landscape ####
data <- fread("extinction_data.csv")

# plot tiled facets
p <-ggplot(data)+
  geom_tile(aes(x = R_new, y = P_new, fill = extinct),
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


ggsave("extinctionplot.png", p, height = 3, width = 10, dpi = 300)

#### population trends ####
pop_data <- fread("data/pop_trend_data.csv")

# calculate if R and P are increased or not
pop_data[,`:=`(deltaP = P_new - P, deltaR = R_new - R)]
pop_data[,shape:=ifelse(deltaP > 0, 24, 25)]

# convert pop to proportion and add const val
pop_data[,prop:=(popsize/5e3)+(20*deltaP)]

ggplot(pop_data[deltaP %% 0.2 == 0,])+
  geom_line(aes(x = gen, y = popsize, col = deltaR,
                group = interaction(deltaR, deltaP)),
            size = 0.2, alpha = 0.9)+
  
  # geom_point(aes(x = gen, y = prop,
  #                fill = deltaR), alpha = 0.5, size = 0.1)+
  scale_size(range = c(0.1, 0.5))+
  scale_fill_scico(palette = "berlin")+
  scale_colour_scico(palette = "berlin")+
  facet_grid(R~P, labeller = label_both,
             switch = "both", as.table = FALSE)+
  theme_few()+
  coord_equal(ratio = 0.2)+
  theme(strip.placement = "outside",
        legend.position = "none", 
        strip.background = element_rect(fill = "grey90", 
                                        colour = "black"))