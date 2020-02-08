#### code to read preliminary output from botero sims ####

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
# list data files
data_files <- list.files("data/", full.names = TRUE)

# get params from filename
parameters <- str_extract_all(data_files, regex("[0-9].[0-9]")) %>% 
  map(function(vals){
    tibble(R = vals[1], P = vals[2])
  }) %>% 
  bind_rows()

# read data
data <- lapply(data_files, fread)

# attach parameters to data
data <- mutate(parameters, agents = data)

#### function to process cue ####
process_cue <- function(h, I01, I02, b1, b2, s, a){
  
  # work through parameters getting slope and intercept
  if(s <= 0.5){
    if(h <= 1){
      # name strategy
      data_rnorm <- data.table(intercept = c(I01*h, I02*(1-h)),
                               slope = 0, strategy = c("dbh1", "dbh2"))
      
    }else{
      data_rnorm <- data.table(intercept = I01, slope = 0, strategy = "cbh")
    }
  }else{
    if(h <= 1){
      # name strategy
      data_rnorm <- data.table(intercept = I02, slope = b2,
                               strategy = ifelse(a < 0, "dev_plast", "fen_plast"))
      
    }else{
      data_rnorm <- data.table(intercept = 0, slope = b1,
                               strategy = ifelse(a < 0, "dev_plast", "fen_plast"))
    }
  }
  
  return(data_rnorm)
  
}

#### get slope and intercept of evolved agents ####
data <- mutate(data, r_norms = map(agents, function(df){
  pmap_df(df[,2:ncol(df)], process_cue)
}))

# remove agents
data <- select(data, -agents) %>% 
  unnest(cols = "r_norms")

# plot figure
ggplot(data)+
  geom_abline(aes(intercept = intercept, slope = slope, col = strategy))+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  coord_cartesian(xlim=c(0,1), ylim=c(0, 1))+
  facet_grid(R~P, labeller = label_both)+
  theme_bw()
