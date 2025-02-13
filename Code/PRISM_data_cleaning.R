library(tidyverse)
library(janitor)
library(lubridate)
library(gghighlight)

getwd()
setwd("/Users/stw52/Documents")

#### Musgrave data ----
musgrave_prism_raw <- read_csv("musgrave_prism.csv")

colnames(musgrave_prism_raw) <- c("date_mug", "maxt", "mint", "prcp")

head(musgrave_prism_raw)

f_to_c <- 5/9

musgrave_prism_raw$date_mug <- as.Date(musgrave_prism_raw$date_mug, format = "%m/%d/%Y")


musgrave_prism_clean <- musgrave_prism_raw %>%
  mutate(date_mug = as.Date(date_mug, format = "%m/%d/%Y")) %>%
  mutate(maxt_c = (maxt-32)*f_to_c) %>%
  mutate(mint_c = (mint-32)*f_to_c) %>%
  mutate(prcp_mm = prcp*25.4) %>%
  mutate(year_mug = year(date_mug)) %>%
  mutate(day_of_year = yday(date_mug)) %>%
  group_by(year_mug) %>%
  mutate(plant_prcp = sum(prcp[day_of_year %in% 135:166])) %>%
  ungroup()

write_csv(musgrave_prism_clean, "clean_musgrave_prism.csv")

#### Geneva data ----

geneva_prism_raw <- read_csv("geneva_prism.csv")

colnames(geneva_prism_raw) <- c("date", "maxt", "mint", "prcp")

head(geneva_prism_raw)

f_to_c <- 5/9

geneva_prism_raw$date <- as.Date(geneva_prism_raw$date, format = "%m/%d/%Y")


geneva_prism_clean <- geneva_prism_raw %>%
  mutate(date = as.Date(date, format = "%m/%d/%Y")) %>%
  mutate(maxt_c = (maxt-32)*f_to_c) %>%
  mutate(mint_c = (mint-32)*f_to_c) %>%
  mutate(prcp_mm = prcp*25.4) %>%
  mutate(year = year(date)) %>%
  mutate(day_of_year = yday(date)) %>%
  group_by(year) %>%
  mutate(plant_prcp = sum(prcp[day_of_year %in% 135:166])) %>%
  ungroup()

write_csv(geneva_prism_clean, "clean_geneva_prism.csv")

#### Arlington data ---- 
arlington_prism_raw <- read_csv("arlington_prism.csv")

colnames(arlington_prism_raw) <- c("date", "maxt", "mint", "prcp")

head(arlington_prism_raw)

f_to_c <- 5/9

arlington_prism_raw$date <- as.Date(arlington_prism_raw$date, format = "%m/%d/%Y")


arlington_prism_clean <- arlington_prism_raw %>%
  mutate(date = as.Date(date, format = "%m/%d/%Y")) %>%
  mutate(maxt_c = (maxt-32)*f_to_c) %>%
  mutate(mint_c = (mint-32)*f_to_c) %>%
  mutate(prcp_mm = prcp*25.4) %>%
  mutate(year = year(date)) %>%
  mutate(day_of_year = yday(date)) %>%
  group_by(year) %>%
  mutate(plant_prcp = sum(prcp[day_of_year %in% 135:166])) %>%
  ungroup()

write_csv(arlington_prism_clean, "clean_arlington_prism.csv")


#### Plot of precipitation around planting since 2000 at Musgrave with years of study highlighted
musgrave_prism_clean %>% 
  distinct(plant_prcp, year_mug) %>%
  ggplot(aes(x = year_mug, y = plant_prcp)) +
  geom_col(fill = "darkgreen") +
  gghighlight(year_mug %in% c(2013,2014,2015,2016,2018,2019)) +
  theme_minimal() +
  labs(y = "Precipitation accum. May 15th to June 15th (in)", x = "Year")
              

t.test(musgrave_prism_clean$plant_prcp[musgrave_prism_clean$year_mug == 2016],
       musgrave_prism_clean$plant_prcp[musgrave_prism_clean$year_mug != 2016],
       alternative = "less")
