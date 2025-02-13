#### This file is for cleaning Daymet climate data for use in APSIM soybean yield gap analysis project.
#### The daymet files that are in use here are from Musgrave Research Farm NY, Cornell Agritech NY, and Arlington WI.

#Author: Sam Wallace
#Date Started: 7/29/2024
#Purpose: Clean data from Daymet into format for use in APSIM


####----Load in libraries

library(tidyverse)
library(janitor)
library(lubridate)
library(here)
library(daymetr)

####----Read in raw data

#Using here function, you can change where you read data in, as long as data folders are in folder where R project is
#Daymet files will need to be edited to remove information in rows 1-6 or until the column names to use in R

mus_daymet <- daymetr::download_daymet(site = "Musgrave", lat = 42.73, lon = -76.66, start = 2000, end = 2021)

mus_weather <- clean_names(mus_daymet$data)

mus_weather <- mus_weather %>%
  mutate(daily_rad_mj_m2 = (dayl_s*srad_w_m_2)/1e6) %>%
  select(year, yday, daily_rad_mj_m2, tmax_deg_c, tmin_deg_c, prcp_mm_day) %>%
  mutate(daily_rad_mj_m2 = round(daily_rad_mj_m2,2))

write_csv(mus_weather, "R_data/Weather_data/Daymet_clean/mus_long_term.csv")

musgrave_file <- here("R_data", "Weather_data", "Daymet_for_R", "musgrave_daymet_raw.csv")

#Use clean names to read in file for more R friendly naming in file
musgrave_raw <- clean_names(read_csv(musgrave_file))

#Calculate daily solar radiation based on equation provided from Daymet and select variables that are needed in APSIM
musgrave_update <- musgrave_raw %>%
  mutate(daily_rad_mj_m2 = (dayl_s*srad_w_m_2)/1000000) %>%
  mutate(date = as.Date(paste(year, yday), format = "%Y %j")) %>%
  select(year, yday, daily_rad_mj_m2, tmax_deg_c, tmin_deg_c, prcp_mm_day) %>%
  mutate(daily_rad_mj_m2 = round(daily_rad_mj_m2, 2)) 

musgrave_update <- musgrave_update %>%
#Change column names to the required column names in APSIM
colnames(musgrave_update) <- c('year', "day", 'radn', 'maxt', 'mint', 'rain')

#Check that you have no missing values in dataset; APSIM cannot handle missing weather data
sum(is.na(musgrave_update) == TRUE)

#Save file to your cleaned weather data folder
musgrave_file_clean <- here("R_data", "Weather_Data", "Daymet_clean","musgrave_daymet_clean.csv")
write_csv(musgrave_update, musgrave_file_clean)


#After saving data to folder, you will need to make constant file that at least includes latitude, but having
#average temperature, temperature amplitude and longitude can be helpful as well. The clean file can also be copied to
#the data folder outside of the R project for ease of access.


#### The above process can now be duplicated for the remaining sites in the study, in this case Geneva and Arlington

#File preparation for Geneva, NY

geneva_file <- here("R_data", "Weather_data", "Daymet_for_R", "geneva_daymet_raw.csv")

#Use clean names to read in file for more R friendly naming in file
geneva_raw <- clean_names(read_csv(geneva_file))

#Calculate daily solar radiation based on equation provided from Daymet and select variables that are needed in APSIM
geneva_update <- geneva_raw %>%
  mutate(daily_rad_mj_m2 = (dayl_s*srad_w_m_2)/1000000) %>%
  mutate(date = as.Date(paste(year, yday), format = "%Y %j")) %>%
  select(year, yday, daily_rad_mj_m2, tmax_deg_c, tmin_deg_c, prcp_mm_day) %>%
  mutate(daily_rad_mj_m2 = round(daily_rad_mj_m2, 2)) 

#Change column names to the required column names in APSIM
colnames(geneva_update) <- c('year', "day", 'radn', 'maxt', 'mint', 'rain')

#Check that you have no missing values in dataset; APSIM cannot handle missing weather data
sum(is.na(geneva_update) == TRUE)

#Save file to your cleaned weather data folder
geneva_file_clean <- here("R_data", "Weather_Data", "Daymet_clean","geneva_daymet_clean.csv")
write_csv(geneva_update, geneva_file_clean)

#File preparation for Arlington, WI

arlington_file <- here("R_data", "Weather_data", "Daymet_for_R", "arlington_daymet_raw.csv")

#Use clean names to read in file for more R friendly naming in file
arlington_raw <- clean_names(read_csv(arlington_file))

#Calculate daily solar radiation based on equation provided from Daymet and select variables that are needed in APSIM
arlington_update <- arlington_raw %>%
  mutate(daily_rad_mj_m2 = (dayl_s*srad_w_m_2)/1000000) %>%
  mutate(date = as.Date(paste(year, yday), format = "%Y %j")) %>%
  select(year, yday, daily_rad_mj_m2, tmax_deg_c, tmin_deg_c, prcp_mm_day) %>%
  mutate(daily_rad_mj_m2 = round(daily_rad_mj_m2, 2)) 

#Change column names to the required column names in APSIM
colnames(arlington_update) <- c('year', "day", 'radn', 'maxt', 'mint', 'rain')

#Check that you have no missing values in dataset; APSIM cannot handle missing weather data
sum(is.na(arlington_update) == TRUE)

#Save file to your cleaned weather data folder
arlington_file_clean <- here("R_data", "Weather_Data", "Daymet_clean","arlington_daymet_clean.csv")
write_csv(arlington_update, arlington_file_clean)
