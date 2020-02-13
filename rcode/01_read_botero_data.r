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
data_files <- list.files("data/botero/", full.names = TRUE)

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
    if(h < 1 && h > 0){
      # name strategy
      data_rnorm <- rbind(c(I01, 0, "dbh1"), 
                          c(I02, 0, "dbh2"))
      
    }else if (h >= 1 ){
      data_rnorm <- rbind(c(I01, 0, "cbh/adt"))
    }
    else if (h  <= 0){
      data_rnorm <- rbind(c(I02,0,"cbh/adt"))
    }
  }else{
    if(h < 1 && h > 0){
      # name strategy
      data_rnorm <- rbind(c(I01, b1, "pbh01"), c(I02, b2, "pbh02"))
      
    }else if (h >= 1 ){
      data_rnorm <- rbind(c(I01, b1, ifelse(a <= 0, "dev_plast", "fen_plast")))
    }
    else if (h <= 0 ){
      data_rnorm <- rbind(c(I02, b2, ifelse(a <= 0, "dev_plast", "fen_plast")))
    }
  }
  
  colnames(data_rnorm) <- c("intercept", "slope","strat")
  return(data_rnorm)
  
}

#### get slope and intercept of evolved agents ####
data <- mutate(data, r_norms = map(agents, function(df){
  temp_mat <- pmap(df[,2:ncol(df)], process_cue)
  temp_mat <- reduce(temp_mat, rbind)
  temp_mat <- as.data.table(temp_mat)
  temp_mat[,`:=`(slope = as.numeric(slope), intercept = as.numeric(intercept))]
}))

# remove agents
data2 <- select(data, -agents) %>% 
  unnest(cols = "r_norms")

# plot figure
ggplot(data2)+
  geom_abline(aes(intercept = intercept, slope = slope, col = strat), 
              alpha = 0.1)+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+
  coord_cartesian(xlim=c(-1,1), ylim=c(-1, 1))+
  scale_colour_manual(values = c("black", "green", "green", "blue", "red", "yellow", "yellow"))+
  facet_grid(P~R, labeller = label_both, as.table = FALSE)+
  theme(strip.placement = "inside")+
  theme_bw()
