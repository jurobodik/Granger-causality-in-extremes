# Meteorological application for the paper 'Granger causality in extremes'

# main function needed to be uploaded for this application is 'Extreme_causality_test' - this function can be found in the main file
source("./R/Main_functions.R")
source("./R/utils.R")

# Part 1: Data import. Note that the data can be accessed through hydrodaten.admin.ch and gate.meteoswiss.ch/idaweb after registration or by requesting the used data from the authors of Pasche et al.(2022). I do not have a permission to share them
# Part 2: Playground. You can try out uploading one station and conduct the test
# Part 3: Contains the test for each station while saving the Gamma values.

# The name of the meteostation is 'AIR', or 'Airolo'.
# Precipitation measured from 1883 to 2016
# Latitude: 46.5327
# Longitude: 8.6049
# Elevation	1139

library(readxl)
library(zoo)

data_folder <- "./data/Meteo_Application/"
results_folder <- "./Results/Meteo_Application/"
# precip_station <- "AIR" # the name of the station where precipitation was measured.
precip_station <- "AUB" # the name of the station where precipitation was measured.
# precip_station <- "BAU" # the name of the station where precipitation was measured.


## Part 1: Data import

discharges <- read.csv(file = paste0(data_folder, "Discharges_summer.csv")) # data about discharges at stations 1 to 68
precipitation <- read.csv(file = paste0(data_folder, "precip.csv")) # Precipitation data
meteo_stations_data <- read.csv(file = paste0(data_folder, "meteo_stations_data.csv")) # Other measurments such as temperature or humidity. These measurments were not from 'AIR' station but at the nearby 'LUZ' station, so there is a small data-collection discrepacy

names(precipitation)[1] <- "Date"
meteo_stations_data$Date <- as.Date(meteo_stations_data$Date)



read_station_precip_causes_discharge <- function(river_station, minimum_number_of_observation = 5000) {
  z1 <- discharges[, c(1, river_station + 1)]
  z1$Date <- as.Date(z1$Date)
  z2 <- precipitation[, c("Date", precip_station)]
  z2$Date <- as.Date(z2$Date)

  meteo_stations_data <- meteo_stations_data[, colSums(!is.na(meteo_stations_data)) > minimum_number_of_observation]
  meteo_stations_data$Date <- meteo_stations_data$Date + 2

  merged_df <- merge(z1, z2, by = "Date", all = TRUE)

  names(merged_df) <- c("Date", "discharge", "precip")
  merged_df <- merge(merged_df, meteo_stations_data, by = "Date", all = TRUE)
  # data <- subset(na.omit(merged_df), select = -c(Name, X)) # OP: what is X? it's not in my version of the dataset.
  data <- subset(na.omit(merged_df), select = -c(Name))

  for (j in 2:ncol(data)) {
    data[, j] <- as.numeric(data[, j])
  }
  return(data)
}





## Part 2: Playground

# There is 68 river discharge stations
# i <- sample(size = 1, 1:68) 
i <- 55
data <- read_station_precip_causes_discharge(i) # This is how you load the data for river station number i


x <- data$precip
y <- data$discharge
z1 <- data$`Air temperature 2 m above ground; daily maximum`
z2 <- data$`Relative air humidity 2 m above ground; daily maximum`
z <- data.frame(z1, z2) # You can add other covariates if you like. But we just really dont think they are significant confounders, so we ended up just using these two



Extreme_causality_test(x, y, z = z, p_value_computation = FALSE, both_tails = FALSE)
Extreme_causality_test(y, x, z = z, p_value_computation = FALSE, both_tails = FALSE)

# If you want to change lag or change the adjustment for confounding or consider also instantenius causality, change some hyperparameters
lag_future <- 1
lag_past <- 0
nu_x <- 0.3
q_y <- 0.2
q_z <- 0.1
instant <- FALSE

Extreme_causality_test(x, y, z = z, nu_x = nu_x, q_y = q_y, q_z = q_z, lag_future = lag_future, lag_past = lag_past, instant = instant, p_value_computation = FALSE, both_tails = FALSE)
Extreme_causality_test(y, x, z = z, nu_x = nu_x, q_y = q_y, q_z = q_z, lag_future = lag_future, lag_past = lag_past, instant = instant, p_value_computation = FALSE, both_tails = FALSE)

Extreme_causality_test(x, y, z = z, nu_x = nu_x, q_y = q_y, q_z = q_z, lag_future = lag_future, lag_past = lag_past, instant = instant, p_value_computation = TRUE, both_tails = FALSE)
Extreme_causality_test(y, x, z = z, nu_x = nu_x, q_y = q_y, q_z = q_z, lag_future = lag_future, lag_past = lag_past, instant = instant, p_value_computation = TRUE, both_tails = FALSE)








## Part 3: Test for each station

result_X_to_Y <- list()
result_Y_to_X <- list()

cat('Progress (/68): \n')
for (i in 1:68) {
  data <- read_station_precip_causes_discharge(i)

  x <- data$precip
  y <- data$discharge
  z1 <- data$`Air temperature 2 m above ground; daily maximum`
  z2 <- data$`Relative air humidity 2 m above ground; daily maximum`
  z <- data.frame(z1, z2)


  # result_X_to_Y <- append(result_X_to_Y, Extreme_causality_test(x, y, z, p_value_computation = TRUE, both_tails = FALSE)$p_value_tail)
  # result_Y_to_X <- append(result_Y_to_X, Extreme_causality_test(y, x, z, p_value_computation = TRUE, both_tails = FALSE)$p_value_tail)
  result_X_to_Y[[i]] <- Extreme_causality_test(x, y, z, p_value_computation = TRUE, both_tails = FALSE)
  result_Y_to_X[[i]] <- Extreme_causality_test(y, x, z, p_value_computation = TRUE, both_tails = FALSE)

  # cat(i, " = ", result_X_to_Y[[i]], "\n")
  cat(i, " = ", result_X_to_Y[[i]]$p_value_tail, "\n")
}

safe_save_rds(result_X_to_Y, file = paste0(results_folder, precip_station, "_result_X_to_Y.rds"))
safe_save_rds(result_Y_to_X, file = paste0(results_folder, precip_station, "_result_Y_to_X.rds"))

# Each element of the list is a list containing: output, p_value_tail, p_value_extreme, CTC, baseline
# Transform the list of lists into a data frame for easier analysis of all those variables:
result_X_to_Y_df <- do.call(rbind, lapply(result_X_to_Y, as.data.frame))
result_Y_to_X_df <- do.call(rbind, lapply(result_Y_to_X, as.data.frame))

write.csv(result_X_to_Y_df, file = paste0(results_folder, precip_station, "_result_X_to_Y_df.csv"))
write.csv(result_Y_to_X_df, file = paste0(results_folder, precip_station, "_result_Y_to_X_df.csv"))

# result_X_to_Y <- readRDS(paste0(results_folder, precip_station, "_result_X_to_Y.rds"))
# result_Y_to_X <- readRDS(paste0(results_folder, precip_station, "_result_Y_to_X.rds"))
# for (i in 1:68) {
#   cat(i, " = ", result_X_to_Y[[i]]$p_value_tail, "\n")
# }

# How many times we obtained the correct output?
# (sum(result_X_to_Y <= 0.05) + sum(result_Y_to_X >= 0.05)) / 136
(sum(result_X_to_Y$p_value_tail <= 0.05) + sum(result_Y_to_X$p_value_tail >= 0.05)) / 136

# Which stations gave the wrong result?
# c(which(result_X_to_Y <= 0.05), which(result_Y_to_X >= 0.05))
c(which(result_X_to_Y$p_value_tail <= 0.05), which(result_Y_to_X$p_value_tail >= 0.05))


