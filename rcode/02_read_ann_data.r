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
{
  data_files <- list.files("data/", full.names = TRUE, pattern = "ann")
  
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
  
  # unnest data
  data <- unnest(data, cols = agents)}

#df <- count(data, R, P, cue, resp)

# plot by combo
ggplot(data)+
  # geom_line(aes(x = cue, y = resp - baseline, group = ind), col = "red", alpha=0.01)+
  #geom_hline(aes(yintercept = baseline), col = "blue", alpha=0.01)+
  geom_line(aes(x = cue, y = resp, group = ind), alpha=0.01)+
  facet_grid(R~P, labeller = label_both)
