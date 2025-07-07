#### Making APSIM ready weather files with PRISM precip and temperature data supplemented with Daymet solar radiation
library(tidyverse)
library(lubridate)
library(daymetr)


mus_srad <- download_daymet(
  lat = 42.734, 
  lon = -76.657,
  start = 1981,
  end = 2024,
  internal = TRUE
)
mus_srad_df <- clean_names(mus_srad$data)

leap_years <- seq(1980, 2024, by = 4)

weather_df <- mus_srad_df

weather_df$date <- as.Date(weather_df$yday, origin = paste0(weather_df$year, "-01-01")) - 1

weather_df <- mus_srad_df %>%
  mutate(date = as.Date(yday - 1, origin = paste0(year, "-01-01")))

leap_years <- weather_df %>%
  filter(leap_year(year)) %>%
  pull(year) %>%
  unique()

dec_31_rows <- weather_df %>%
  filter(year %in% leap_years & yday == 364) %>%
  mutate(yday = 365,
         date = as.Date(paste0(year, "-12-31")))

df_extended <- bind_rows(weather_df, dec_31_rows) %>%
  arrange(year, yday)


prism_musgrave <- musgrave_prism %>%
  filter(Year < 2024)

apsim_musgrave <- left_join(prism_musgrave, df_extended) %>%
  dplyr::select(year, doy,srad_w_m_2,maxt_c, mint_c, prcp_cm, dayl_s) %>%
  mutate(radn = srad_w_m_2*dayl_s*(1/1e6),
         maxt = maxt_c,
         mint = mint_c,
         rain = prcp_cm*10,
         day = doy) %>%
  dplyr::select(year,day,radn,maxt,mint,rain)



write_csv(apsim_musgrave, file = "musgrave_apsim_weather.csv")

####
####
#### Same for Arlington
####
####

arl_srad <- download_daymet(
  lat = 43.30314599630097, 
  lon = -89.34653013390276,
  start = 1981,
  end = 2024,
  internal = TRUE
)
arl_srad_df <- clean_names(arl_srad$data)

leap_years <- seq(1980, 2024, by = 4)

weather_df <- arl_srad_df

weather_df$date <- as.Date(weather_df$yday, origin = paste0(weather_df$year, "-01-01")) - 1

weather_df <- arl_srad_df %>%
  mutate(date = as.Date(yday - 1, origin = paste0(year, "-01-01")))

leap_years <- weather_df %>%
  filter(leap_year(year)) %>%
  pull(year) %>%
  unique()

dec_31_rows <- weather_df %>%
  filter(year %in% leap_years & yday == 364) %>%
  mutate(yday = 365,
         date = as.Date(paste0(year, "-12-31")))

df_extended <- bind_rows(weather_df, dec_31_rows) %>%
  arrange(year, yday)


prism_arl <- arlington_prism %>%
  filter(Year < 2024)

apsim_arl <- left_join(prism_arl, df_extended) %>%
  dplyr::select(year, doy,srad_w_m_2,maxt_c, mint_c, prcp_cm, dayl_s) %>%
  mutate(radn = srad_w_m_2*dayl_s*(1/1e6),
         maxt = maxt_c,
         mint = mint_c,
         rain = prcp_cm*10,
         day = doy) %>%
  dplyr::select(year,day,radn,maxt,mint,rain)



write_csv(apsim_arl, file = "arlington_apsim_weather.csv")


