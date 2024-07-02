#Meteorological application for the paper 'Granger causality in extremes'

#main function needed to be uploaded for this application is 'Extreme_causality_test' - this function can be found in the main file


#First 30 lines is just data-upload. Note that the data can be accessed through hydrodaten.admin.ch and gate.meteoswiss.ch/idaweb after registration or by requesting the used data from the authors of Pasche et al.(2022). I do not have a permission to share them
#lines 50 - 100 are a playground - you can try out uploading one station and conduct the test
#lines 100-200 is just conducting the test for each station while saving the Gamma values. 
#lines 200-500 is the function 'Extreme_causality_test' that you need to upload

#The name of the meteostation is 'AIR', or 'Airolo'. 
#Precipitation measured from 1883 to 2016
#Latitude: 46.5327
#Longitude: 8.6049
#Elevation	1139

library("readxl")
library(zoo)

discharges = read.csv(file = 'Discharges_summer.csv')#data about discharges at stations 1 to 68
precipitation = read.csv(file = 'precip.csv') #Precipitation data
meteo_stations_data= read.csv(file = 'meteo_stations_data.txt') #Other measurments such as temperature or humidity. These measurments were not from 'AIR' station but at the nearby 'LUZ' station, so there is a small data-collection discrepacy

names(precipitation)[1] = 'Date'
meteo_stations_data$Date = as.Date(meteo_stations_data$Date)



read_station_precip_causes_discharge <- function(station, minimum_number_of_observation=5000){
  
  z1=discharges[,c(1,station+1)]; z1$Date = as.Date(z1$Date)
  z2=precipitation[, c(1, 3 )]; z2$Date = as.Date(z2$Date)
  
  meteo_stations_data = meteo_stations_data[,colSums(!is.na(meteo_stations_data))>minimum_number_of_observation]
  meteo_stations_data$Date = meteo_stations_data$Date +2
  
  merged_df <- merge(z1, z2, by = "Date", all = TRUE)
  
  names(merged_df) = c('Date', 'discharge', 'precip')
  merged_df <- merge(merged_df, meteo_stations_data, by = "Date", all = TRUE)
  data = subset(na.omit(merged_df), select = -c(Name, X))
  
  for (j in 2:ncol(data)) {data[,j] = as.numeric(data[,j])  } 
  return(data)
}








#Playground
i=sample(size =1, 1:68) #There is 68 river discharge stations
i=55
data = read_station_precip_causes_discharge(i) #This is how you upload the data from station number i


x=data$precip
y=data$discharge
z1=data$`Air temperature 2 m above ground; daily maximum`
z2=data$`Relative air humidity 2 m above ground; daily maximum`
z=data.frame(z1, z2) #You can add other covariates if you like. But we just really dont think they are significant confounders, so we ended up just using these two



Extreme_causality_test(x,y,z=z, p_value_computation = FALSE, both_tails = FALSE)
Extreme_causality_test(y,x,z=z, p_value_computation = FALSE, both_tails = FALSE)

#If you want to change lag or change the adjustment for confounding or consider also instantenius causality, change some hyperparameters
lag_future=1
lag_past = 0
nu_x = 0.3
q_y = 0.2
q_z = 0.1
instant=FALSE

Extreme_causality_test(x,y,z=z,  nu_x = nu_x,   q_y =q_y, q_z=q_z, lag_future = lag_future, lag_past = lag_past,instant=instant, p_value_computation = FALSE, both_tails = FALSE)
Extreme_causality_test(y,x,z=z,  nu_x = nu_x,   q_y =q_y, q_z=q_z, lag_future = lag_future, lag_past = lag_past,instant=instant, p_value_computation = FALSE, both_tails = FALSE)

Extreme_causality_test(x, y,z=z,  nu_x = nu_x,   q_y =q_y, q_z=q_z, lag_future = lag_future, lag_past = lag_past,instant=instant, p_value_computation = TRUE, both_tails = FALSE)
Extreme_causality_test(y, x,z=z,  nu_x = nu_x,   q_y =q_y, q_z=q_z, lag_future = lag_future, lag_past = lag_past,instant=instant, p_value_computation = TRUE, both_tails = FALSE)
















result_X_to_Y = list()
result_Y_to_X = list()

for (i in 1:68) {
  
  data = read_station_precip_causes_discharge(i) 
  
  x=data$precip
  y=data$discharge
  z1=data$`Air temperature 2 m above ground; daily maximum`
  z2=data$`Relative air humidity 2 m above ground; daily maximum`
  z=data.frame(z1, z2)
  
  
  result_X_to_Y = append(result_X_to_Y,Extreme_causality_test(x,y,z, p_value_computation = TRUE, both_tails = FALSE)$p_value_tail)
  result_Y_to_X = append(result_Y_to_X,Extreme_causality_test(y,x,z, p_value_computation = TRUE, both_tails = FALSE)$p_value_tail)
  
  cat(i, ' = ', result_X_to_Y[[i]], '\n')
}



#How many times we obtained the correct output?
(sum(result_X_to_Y<=0.05)  + sum(result_Y_to_X>=0.05))/136

#Which stations gave the wrong result?
c(which(result_X_to_Y<=0.05) , which(result_Y_to_X>=0.05))


