#### This file is for finding the precipitation accumulation around a specific +/- day length around planting for each planting date in the SoYGAP project

library(tidyverse)
library(lubridate)
library(ggpubr)
library(ggpmisc)
library(knitr)
library(kableExtra)
library(readxl)

all_soy_rye <- read_csv("R_data/Soy_rye_data/NY_WI_NTvT_master.csv")

ny_soy_rye <- all_soy_rye %>%
  filter(State == "NY")

soy_planting_dates <- all_soy_rye %>%
  select(SoySeedingDate, ExperimentName, Year, State, Location, Rep) %>%
  distinct(SoySeedingDate, ExperimentName, Year, State, Location)


#### New York work ----

# Musgrave
mug_pl_dates <- soy_planting_dates %>%
  filter(State == "NY") %>% 
  filter(!is.na(ExperimentName)) %>%
  distinct(SoySeedingDate, Year, ExperimentName)

mus_met <- read_csv("~/Documents/clean_musgrave_prism.csv") 

mus_annual_rain <-mus_met %>%
  group_by(year_mug) %>%
  summarise(prcp_year = sum(prcp))

mean(mus_annual_rain$prcp_year)

mus_met <- mus_met %>%
  rename("Year" = year_mug)

mus_met_exp <- left_join(mus_met, mug_pl_dates, by = "Year")

base_temp_f <- 39.9
base_temp_c <- 4.4


mus_met_exp <- mus_met_exp%>%
  filter(!is.na(SoySeedingDate)) %>%
  mutate(SoySeedingDate = as.Date(SoySeedingDate, format = "%m/%d/%y")) %>%
  mutate(plant_day_of_year = yday(SoySeedingDate)) %>%
  group_by(SoySeedingDate) %>%
  mutate(plant_prcp = sum(prcp[day_of_year >= plant_day_of_year - 15 & day_of_year <= plant_day_of_year + 15], na.rm = TRUE)) %>%
  mutate(plant_before = sum(prcp[day_of_year >= plant_day_of_year - 15 & day_of_year <= plant_day_of_year], na.rm = TRUE)) %>%
  mutate(gdd_plant = sum((((maxt_c+mint_c)/2)-base_temp_c)[day_of_year >= 60 & day_of_year <= plant_day_of_year], na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(SoySeedingDate, Year, ExperimentName, plant_prcp, gdd_plant)

mus_soy_rye <- ny_soy_rye %>%
  filter(Location == "Musgrave") %>%
  group_by(Year, Till_NT) %>%
  mutate(yield_ratio_nt_t = mean(Yield_Mgha[Till_NT == "NT"]/Yield_Mgha[Till_NT == "Till"], na.rm = TRUE)) %>%
  mutate(yield_sd = sd(Yield_buac)) %>%
  mutate(mean_yield = mean(Yield_buac)) %>%
  distinct(ExperimentName, Year, Location, Yield_buac, Till_NT, yield_sd, mean_yield) %>%
  ungroup()

mus_met_yield <- left_join(mus_met_exp, mus_soy_rye)

# Geneva

geneva_planting_dates <- soy_planting_dates %>%
  filter(is.na(ExperimentName)) %>%
  distinct(Year, SoySeedingDate, Location)

gen_met <- read_csv("~/Documents/clean_geneva_prism.csv")
gen_met <- gen_met %>%
  rename("Year" = year)

gen_met_exp <- left_join(gen_met, geneva_planting_dates, by = "Year") 
gen_met_exp <- gen_met_exp %>%
  filter(!is.na(SoySeedingDate)) %>%
  mutate(SoySeedingDate = as.Date(SoySeedingDate, format = "%m/%d/%y")) %>%
  mutate(plant_day_of_year = yday(SoySeedingDate)) %>%
  group_by(SoySeedingDate) %>%
  mutate(plant_prcp = sum(prcp[day_of_year >= plant_day_of_year - 15 & day_of_year <= plant_day_of_year + 15], na.rm = TRUE)) %>%
  mutate(gdd_plant = sum((((maxt+mint)/2)-base_temp_f)[day_of_year >= plant_day_of_year - 10 & day_of_year <= plant_day_of_year + 10], na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(SoySeedingDate, Year, plant_prcp, gdd_plant)

gen_soy_rye <- ny_soy_rye %>%
  filter(Location == "Geneva") %>%
  group_by(Year, Till_NT) %>%
  mutate(yield_ratio_nt_t = mean(Yield_Mgha[Till_NT == "NT"]/Yield_Mgha[Till_NT == "Till"], na.rm = TRUE)) %>%
  mutate(yield_sd = sd(Yield_buac)) %>% 
  mutate(mean_yield = mean(Yield_buac)) %>%
  distinct(Year, Location, Yield_buac, Till_NT, yield_sd, mean_yield) %>%
  ungroup()

gen_met_yield <- left_join(gen_met_exp, gen_soy_rye)

gen_met_yield <- gen_met_yield %>%
  mutate(ExperimentName = rep(NA))

ny_met_yield <- rbind(mus_met_yield, gen_met_yield)

#### Wisconsin work

wi_pl_dates <- soy_planting_dates %>%
  filter(State == "WI") %>% 
  filter(!is.na(SoySeedingDate)) %>%
  distinct(ExperimentName, Year,SoySeedingDate) %>%
  filter(SoySeedingDate != "5/30/17") %>%
  mutate(Location = rep("AARS", 10))

wi_pl_dates %>%
  count(Year, Location) %>%
  filter(n > 1)

wi_met <- read_csv("~/Documents/clean_arlington_prism.csv")

wi_met <- wi_met %>%
  rename("Year" = year) %>%
  mutate(Location = rep("AARS", nrow(wi_met)))

wi_met_2021 <- wi_met %>%
  filter(Year == 2021)

wi_met_exp <- left_join(arlington_prism, wi_pl_dates, by = c("Year", "Location"))

wi_met_exp <- wi_met_exp %>%
  filter(!is.na(SoySeedingDate)) %>%
  mutate(SoySeedingDate = as.Date(SoySeedingDate, format = "%m/%d/%y")) %>%
  mutate(plant_day_of_year = yday(SoySeedingDate)) %>%
  group_by(SoySeedingDate) %>%
  mutate(plant_prcp = sum(prcp_cm[doy >= plant_day_of_year - 15 & doy <= plant_day_of_year + 15], na.rm = TRUE)) %>%
 # mutate(gdd_plant = sum((((maxt+mint)/2)-base_temp_f)[day_of_year >= plant_day_of_year - 10 & day_of_year <= plant_day_of_year + 10], na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(SoySeedingDate, Year, ExperimentName, plant_prcp)
head(wi_met_exp)

wi_soy_rye <- all_soy_rye %>%
  filter(State == "WI")

wi_soy_rye <- wi_soy_rye %>%
  group_by(ExperimentName, Year, Till_NT) %>%
  # Calculate the mean yield ratio of NT to Till
  mutate(yield_ratio_nt_t = mean(Yield_buac[Till_NT == "NT"], na.rm = TRUE)/ 
           mean(Yield_buac[Till_NT == "Till"], na.rm = TRUE)) %>%
  mutate(yield_sd = sd(Yield_buac)) %>%
  mutate(mean_yield = mean(Yield_buac)) %>%
  # Remove duplicates of yield ratios along with ExperimentName and Year
  distinct(ExperimentName, Year,Location, Yield_buac, Till_NT, yield_sd, mean_yield) %>%
  ungroup()

wi_met_yield <- left_join(wi_met_exp, wi_soy_rye)

wi_met_yield <- wi_met_yield %>%
  filter(SoySeedingDate != "2019-06-05")

wi_met_yield$Year[wi_met_yield$Year == 2021 & wi_met_yield$ExperimentName == "Seeding Depth Study"] <- "2021a" 
wi_met_yield$Year[wi_met_yield$Year == 2021 & wi_met_yield$ExperimentName == "Starter Fertilizer Study"] <- "2021b" 

wi_met_yield


comb_met_yield <- rbind(wi_met_yield, ny_met_yield)

comb_met_yield <- comb_met_yield %>%
  distinct(Year, Location, yield_sd, ExperimentName, plant_prcp, Till_NT, mean_yield, gdd_plant)

plant_prcp_exp <- comb_met_yield %>%
  distinct(Year, ExperimentName, Location, plant_prcp, gdd_plant) %>%
  mutate(prcp_cm = plant_prcp*2.54)

write_csv(plant_prcp_exp, "R_data/Weather_data/PRISM/plant_prcp_all.csv")


#### Make normal plots for each location

#Musgrave 
musgrave_normals <- read_xlsx("~/Documents/musgrave_climate_normals.xlsx")

musgrave_normals %>%
  rename("Month" = date, `Precipitation (in)` = pcpn_sum,`Average maximum temperature (F)` = maxt_mean,
         `Average minimum temperature (F)`=  mint_mean,`Average temperature (F)` = avgt_mean) %>%
  pivot_longer(c(`Average temperature (F)`, `Average maximum temperature (F)`, `Average minimum temperature (F)`),
               values_to = "temp_data", names_to = "measure") %>%
  mutate(Month = as.character(Month),  # Convert Month to character if it's numeric
         Month = case_when(
           Month == "1" ~ "Jan",
           Month == "2" ~ "Feb",
           Month == "3" ~ "Mar",
           Month == "4" ~ "Apr",
           Month == "5" ~ "May",
           Month == "6" ~ "Jun",
           Month == "7" ~ "Jul",
           Month == "8" ~ "Aug",
           Month == "9" ~ "Sep",
           Month == "10" ~ "Oct",
           Month == "11" ~ "Nov",
           Month == "12" ~ "Dec"
         )) %>%
  mutate(Month = factor(Month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul",
                                          "Aug", "Sep", "Oct", "Nov", "Dec")))->musgrave_normals

musgrave_normals%>%
  ggplot(aes(x = Month)) +
  geom_col(aes(y = `Precipitation (in)`*5),fill = "skyblue" ,alpha = 0.75) +
  geom_line(aes(y = temp_data, colour = measure, group = measure)) +
  geom_point(aes(y = temp_data, colour = measure, group = measure, shape = measure))+
  theme_bw() +
  scale_y_continuous(name = "Temperature (F)", sec.axis = sec_axis(~./15, name = "Precipitation (in)")) +
  labs(fill = "", colour = "", shape = "", title = "Aurora, NY Climate Normals 1991-2020") +
  theme(legend.position = c(.2,0.9),
        legend.background = element_rect(fill = NA, colour = NA),
        plot.title = element_text(hjust = 0.5))+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12))

sum(musgrave_normals$`Precipitation (in)`/3)

#Arlington
arlington_normals <- read_csv("~/Downloads/AARS 1991-2020-climate-normal.csv")

arlington_normals %>%
  pivot_longer(c(`Average temperature (F)`, `Average maximum temperature (F)`, `Average minimum temperature (F)`),
               values_to = "temp_data", names_to = "measure") %>%
  mutate(Month = factor(Month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul",
                                          "Aug", "Sep", "Oct", "Nov", "Dec"))) ->arlington_normals

arlington_normals%>%
  ggplot(aes(x = Month)) +
  geom_col(aes(y = `Precipitation (in)`*5),fill = "skyblue" ,alpha = 0.75) +
  geom_line(aes(y = temp_data, colour = measure, group = measure)) +
  geom_point(aes(y = temp_data, colour = measure, group = measure, shape = measure))+
  theme_bw() +
  scale_y_continuous(name = "Temperature (F)", sec.axis = sec_axis(~./15, name = "Precipitation (in)")) +
  labs(fill = "", colour = "", shape = "", title = "Arlington, WI Climate Normals 1991-2020") +
  theme(legend.position = c(.2,0.9),
        legend.background = element_rect(fill = NA, colour = NA),
        plot.title = element_text(hjust = 0.5))+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12))
#Geneva

geneva_normals <- read_xlsx("~/Documents/geneva_climate_normals.xlsx")

geneva_normals %>%
  rename("Month" = date, `Precipitation (in)` = pcpn_sum,`Average maximum temperature (F)` = maxt_mean,
         `Average minimum temperature (F)`=  mint_mean,`Average temperature (F)` = avgt_mean) %>%
pivot_longer(c(`Average temperature (F)`, `Average maximum temperature (F)`, `Average minimum temperature (F)`),
               values_to = "temp_data", names_to = "measure") %>%
  mutate(Month = as.character(Month),  # Convert Month to character if it's numeric
         Month = case_when(
           Month == "1" ~ "Jan",
           Month == "2" ~ "Feb",
           Month == "3" ~ "Mar",
           Month == "4" ~ "Apr",
           Month == "5" ~ "May",
           Month == "6" ~ "Jun",
           Month == "7" ~ "Jul",
           Month == "8" ~ "Aug",
           Month == "9" ~ "Sep",
           Month == "10" ~ "Oct",
           Month == "11" ~ "Nov",
           Month == "12" ~ "Dec"
         )) %>%
  mutate(Month = factor(Month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul",
                                          "Aug", "Sep", "Oct", "Nov", "Dec")))->geneva_normals

geneva_normals%>%
  ggplot(aes(x = Month)) +
  geom_col(aes(y = `Precipitation (in)`*5),fill = "skyblue" ,alpha = 0.75) +
  geom_line(aes(y = temp_data, colour = measure, group = measure)) +
  geom_point(aes(y = temp_data, colour = measure, group = measure, shape = measure))+
  theme_bw() +
  scale_y_continuous(name = "Temperature (F)", sec.axis = sec_axis(~./15, name = "Precipitation (in)")) +
  labs(fill = "", colour = "", shape = "", title = "Geneva, NY Climate Normals 1991-2020") +
  theme(legend.position = c(.2,0.9),
        legend.background = element_rect(fill = NA, colour = NA),
        plot.title = element_text(hjust = 0.5))+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12))

#Function to use if I need to convert these values into metric for weather
eng_to_metric <- eng_to_metric <- function(maxt_f,mint_f,avgt_f,prcp_in) {
  maxt_c <- (maxt_f-32)/(1.8)
  mint_c <- (mint_f-32)/(1.8)
  avgt_c <- (avgt_f-32)/(1.8)
  prcp_mm <- (prcp_in*25.4)
  
  return(data.frame(maxt_c, mint_c, avgt_c, prcp_mm))
}

geneva_normals %>%
  mutate(conversions = eng_to_metric(maxt_f = maxt_mean,mint_f = mint_mean,avgt_f = avgt_mean,prcp_in = pcpn_sum)) %>%
  unnest(conversions)

