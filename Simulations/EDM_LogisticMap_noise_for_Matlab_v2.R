rm(list=ls())

library(rEDM)
library(dplyr)
library(R.matlab)
library(ggplot2)
library(here)

data <-readMat(here::here("Simulations", 'LogisticMapTimeSeriesData.mat')) #Get time series data from Matlab
source(here::here("Simulations","nRMSE.R"))

#specify what portions of the data to use for constructing the simplex projection model and for testing the forecast skill.
lib <- c(1, 2800)
pred <- c(2801, 3000)

# select s-map tuning parameters
E_select = 2 #No need to use Simplex to determine the best embedding dimension with which to unfold the attactor; the embedding dimension for the logistic map is just 2. 
theta_select = 10
max_steps <- 1

# get growth rates
growth_rates <- data$rs
noise_level <- names(data)[1:8]

df_list <- vector("list", 100)
output_list <- vector("list")

for (k in 1:8){  #for each noise level 
  
  #X <- data$LM.Zpn
  X <- data[[noise_level[k]]]   #columns are time series for different values of the growth rate, r
  
  X.df <- as.data.frame(X)
  names(X.df) <- growth_rates
  
  
  # automatically select theta
  smap_output_list <- vector("list", 10)
  for (j in 1:ncol(X)) {

    ts <- X[,j]
    gr <- growth_rates[j]
    smap_output_list[[j]] <- lapply(1:max_steps, function(x) s_map(ts, lib, pred, E = E_select, theta=theta_select, tp=x, stats_only = FALSE))
    rmse_list <- lapply(1:max_steps, function(x) smap_output_list[[j]][[x]]["rmse"])
    df <- data.frame(noise_level = noise_level[k], growth_rate=rep(gr, max_steps), tp = 1:max_steps, rmse = unlist(rmse_list))
    df_list[[j]] <- df
    
  }
  
  output_list[[k]] <- bind_rows(df_list)
}

output <- bind_rows(output_list)

output$noise_level <- gsub("LM.", "", output$noise_level)

output$noise_type <- ifelse(grepl("on", output$noise_level), "obs_noise", "process_noise")
output$noise_level <- ifelse(grepl("H", output$noise_level), "3_high", 
                             ifelse(grepl("L", output$noise_level), "1_low", 
                                    ifelse(grepl("M", output$noise_level), "2_medium", 
                                           ifelse(grepl("Z", output$noise_level), "0_zero","4_super_high"))))

ggplot(data=output, aes(y=rmse, x=growth_rate, group=noise_level, colour=growth_rate)) + geom_point(alpha=.4) + 
  guides(colour=F) + scale_color_gradient(low="blue", high = "red") + facet_grid(noise_type~noise_level)

writeMat(here::here("Simulations","LogisticMapForecastError.mat"), output=output)

