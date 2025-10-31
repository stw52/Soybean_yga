#### Making APSIM ready weather files with PRISM precip and temperature data supplemented with Daymet solar radiation
library(tidyverse)
library(lubridate)
library(daymetr)
library(janitor)
library(apsimx)
library(ggpubr)


#### Comparing Daymet and NASA Power to the years of observations

newa_mus <- clean_names(read_csv(c(
  "~/Downloads/mus_2022.csv",
  "~/Downloads/mus_2023.csv",
  "~/Downloads/mus_2024.csv"
)))

newa_metric <- newa_mus %>%
  mutate(date = as.Date(date, format = "%m/%d/%Y"),
         year = year(date),
         day = yday(date),
         radn = solar_radiation_langleys*0.041868,
         maxt = (max_air_temp_f-32)*5/9,
         mint = (min_air_temp_f-32)*5/9,
         rain = total_precipitation*25.4) %>%
  dplyr::select(year, day, radn, maxt, mint, rain) %>%
  mutate(source = rep("NEWA"))


daymet_comparison <- arranged_daymet %>%
  filter(year >= 2022) %>%
  dplyr::select(year, day, radn, maxt, mint, rain) %>%
  mutate(source = rep("Daymet"))
  

power_comparison <- power_mus %>%
  filter(year >= 2022) %>%
  dplyr::select(year, day, radn, maxt, mint, rain) %>%
  mutate(source = rep("NASA"))

prism_comparison <- prism_musgrave %>%
  filter(between(year, 2022, 2024)) %>%
  mutate(radn = rep(0),
         source = rep("PRISM"))

comp_met <- rbind(newa_metric, daymet_comparison, power_comparison, prism_comparison) 


pivoted <- comp_met %>%
  pivot_wider(
    id_cols = c(year, day),
    names_from = source,
    values_from = c(radn, maxt, mint, rain)
  )

# ---- 2. Function to calculate stats ----
error_stats <- function(sim, obs) {
  bias <- mean(sim - obs, na.rm = TRUE)
  mae  <- mean(abs(sim - obs), na.rm = TRUE)
  rmse <- sqrt(mean((sim - obs)^2, na.rm = TRUE))
  corr <- cor(sim, obs, use = "complete.obs")
  return(c(bias = bias, mae = mae, rmse = rmse, corr = corr))
}

# ---- 3. Apply function for each variable and source ----
results <- list(
  radn_power   = error_stats(pivoted$radn_NASA,   pivoted$radn_NEWA),
  radn_daymet  = error_stats(pivoted$radn_Daymet,  pivoted$radn_NEWA),
  maxt_power   = error_stats(pivoted$maxt_NASA,   pivoted$maxt_NEWA),
  maxt_daymet  = error_stats(pivoted$maxt_Daymet,  pivoted$maxt_NEWA),
  maxt_prism = error_stats(pivoted$maxt_PRISM, pivoted$maxt_NEWA),
  mint_power   = error_stats(pivoted$mint_NASA,   pivoted$mint_NEWA),
  mint_daymet  = error_stats(pivoted$mint_Daymet,  pivoted$mint_NEWA),
  mint_prism = error_stats(pivoted$mint_PRISM, pivoted$mint_NEWA),
  rain_power   = error_stats(pivoted$rain_NASA,   pivoted$rain_NEWA),
  rain_daymet  = error_stats(pivoted$rain_Daymet,  pivoted$rain_NEWA),
  rain_prism = error_stats(pivoted$rain_PRISM, pivoted$rain_NEWA)
)

# ---- 4. Convert to a clean summary table ----
results_df <- do.call(rbind, results) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("comparison") %>%
  separate(comparison, into = c("variable","dataset"), sep = "_")



qqplot(pivoted$rain_NEWA, pivoted$rain_NASA)

#### Arlington Nasa Power

power_aars <- get_power_apsim_met(lonlat = c( -89.34367303797328,43.301893058200726), dates = c("1984-01-01","2024-12-31"))

plot(power_aars, met.var = c("rain"), years = c(2013:2016, 2018, 2020, 2023), cumulative = T,
     climatology = T)
plot(power_aars, met.var = c("radn"), years = c(2013:2016, 2018, 2020, 2023), cumulative = T,
     climatology = T)


write_csv(power_aars, "power_aars.csv")

aars_meso <- clean_names(read_csv("~/Documents/aars_mesonet.csv"))

aars_24 <- aars_meso %>%
  mutate(date = mdy_hm(date_time_collected),
         year = year(date),
         day = yday(date),
         radn = daily_total_solar_radiation_mj,
         maxt = daily_maximum_air_temperature_c,
         mint = daily_minimum_air_temperature_c,
         rain = daily_total_rain_mm) %>%
  dplyr::select(year, day, radn, maxt, mint, rain) %>%
  filter(radn != "N/A") %>%
  mutate(radn = as.numeric(radn)) %>%
  filter(year == 2024) %>%
  mutate(source = rep("MESO"))

power_aars_24 <- power_aars %>%
  filter(year == 2024) %>%
  dplyr::select(year, day, radn, maxt, mint, rain) %>%
  filter(day %in% aars_24$day) %>%
  mutate(source = rep("NASA"))

daymet_aars <- get_daymet2_apsim_met(lonlat = c(-89.34367303797328,43.301893058200726), years = c(1984,2024))

daymet_aars <- daymet_aars %>%
  dplyr::select(year, day, radn, maxt, mint, rain)

for (y in leap_years) {
  has_365 <- nrow(filter(daymet_aars, year == y, day == 365)) > 0
  has_366 <- nrow(filter(daymet_aars, year == y, day == 366)) > 0
  if (has_365 & !has_366) {
    row_365 <- filter(daymet_aars, year == y, day == 365)
    new_row <- row_365 %>% mutate(day = 366)
    daymet_aars <- add_row(daymet_aars, !!!as.list(new_row[1, ]))
  }
}

arranged_daymet_aars <- arrange(daymet_aars, year, day) 
write_apsim_met(arranged_daymet_aars,filename =  "aars_apsim_daymet.met", wrt.dir = "~/Downloads")




#### Cumulative seasonal inputs per growing season period by each weather source


sim_years_ny <- c(2015,2016,2023,2024)
sim_years_wi <- c(2015,2017,2018,2019,2023)

power_sim_ny <- power_mus %>%
  filter(year %in% sim_years_ny,
         between(day, 150,300)) %>%
  group_by(year) %>%
  mutate(cum_rain = cumsum(rain),
            cum_radn = cumsum(radn),
         source = rep("NASA Power")) %>%
  dplyr::select(year, day, radn, maxt, mint, rain, cum_rain, cum_radn, source)

dm_sim_ny <- arranged_daymet %>%
  filter(year %in% sim_years_ny,
         between(day, 150,300)) %>%
  group_by(year) %>%
  mutate(cum_rain = cumsum(rain),
         cum_radn = cumsum(radn),
         source = rep("Daymet")) %>%
  dplyr::select(colnames(power_sim_ny))

np_dm_ny <- rbind(power_sim_ny, dm_sim_ny)

np_dm_ny_rad <- np_dm_ny %>%
  ggplot(aes(x = day, y = cum_radn)) +
  geom_line(aes(linetype = source)) +
  theme_pubr() +
  labs(y = expression(paste("Cumulative radiation " , "(Mj"," ",m^{-2},")")), x = "Day of year",
       color = "", linetype = "") +
  facet_wrap(~year) +
  theme(strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14))

ggsave(np_dm_ny_rad, filename = "Figures/np_dm_ny_radn.png",
       dpi = 300, width = 9, height = 6, units = "in")
  
np_dm_ny_rain <- np_dm_ny %>%
  ggplot(aes(x = day, y = cum_rain)) +
  geom_line(aes(linetype = source)) +
  theme_pubr() +
  labs(y = "Cumulative precipitation(mm)", x = "Day of year",
       color = "", linetype = "") +
  facet_wrap(~year) +
  theme(strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14))

ggsave(np_dm_ny_rain, filename = "Figures/np_dm_ny_rain.png",
       dpi = 300, width = 9, height = 6, units = "in")

power_sim_wi <- power_aars %>%
  filter(year %in% sim_years_wi,
         between(day, 150,300)) %>%
  group_by(year) %>%
  mutate(cum_rain = cumsum(rain),
         cum_radn = cumsum(radn),
         source = rep("NASA Power")) %>%
  dplyr::select(year, day, radn, maxt, mint, rain, cum_rain, cum_radn, source)

dm_sim_wi <- arranged_daymet_aars %>%
  filter(year %in% sim_years_wi,
         between(day, 150,300)) %>%
  group_by(year) %>%
  mutate(cum_rain = cumsum(rain),
         cum_radn = cumsum(radn),
         source = rep("Daymet")) %>%
  dplyr::select(colnames(power_sim_ny))

np_dm_wi <- rbind(power_sim_wi, dm_sim_wi)

np_dm_wi_rad <- np_dm_wi %>%
  ggplot(aes(x = day, y = cum_radn)) +
  geom_line(aes(linetype = source)) +
  theme_pubr() +
  labs(y = expression(paste("Cumulative radiation " , "(Mj"," ",m^{-2},")")), x = "Day of year",
       color = "", linetype = "") +
  facet_wrap(~year, nrow = 3) +
  theme(strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14))
ggsave(np_dm_wi_rad, filename = "Figures/np_dm_wi_radn.png",
       dpi = 300, width = 9, height = 6, units = "in")

np_dm_wi_rain <- np_dm_wi %>%
  ggplot(aes(x = day, y = cum_rain)) +
  geom_line(aes(linetype = source)) +
  theme_pubr() +
  labs(y = "Cumulative precipitation(mm)", x = "Day of year",
       color = "", linetype = "") +
  facet_wrap(~year, nrow = 3) +
  theme(strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14))

ggsave(np_dm_wi_rain, filename = "Figures/np_dm_wi_rain.png",
       dpi = 300, width = 9, height = 6, units = "in")


#### Bias comparison for temperature and precip for NASA Power, Daymet, and Climod

climod_mus <- read_csv("R_data/Weather_data/musgrave_climod.csv")

climod_mus_filtered <- climod_mus %>%
  filter(Precipitation != c("M", "T", "S"),
         !endsWith(Precipitation,"A"),
         MaxTemperature != "M",
         MinTemperature != "M") %>%
  mutate(rain = as.numeric(Precipitation)*25.4,
         maxt = (as.numeric(MaxTemperature)-32)*(5/9),
         mint = (as.numeric(MinTemperature)-32)*(5/9),
         date = as.Date(Date, format = "%m/%d/%y"),
         year = year(date),
         day = yday(date)) %>%
  dplyr::select(year, day, maxt, mint, rain) %>%
  mutate(year_day = paste0(year,"-",day)) %>%
  drop_na()
  
bias_power_mus <- power_mus %>%
  mutate(year_day = paste0(year,"-",day)) %>%
  filter(year_day %in% climod_mus_filtered$year_day) %>%
  dplyr::select(year, day, mint, maxt, rain, year_day) %>%
  mutate(source = rep("NASA"))

bias_daymet_mus <- arranged_daymet %>%
  mutate(year_day = paste0(year,"-",day)) %>%
  filter(year_day %in% climod_mus_filtered$year_day) %>%
  dplyr::select(year, day, mint, maxt, rain, year_day) %>%
  mutate(source = rep("Daymet"))

climod_mus_comp <- climod_mus_filtered %>%
  mutate(source = rep("Climod")) %>%
  arrange(year_day)

comp_met_climod <- rbind(bias_power_mus, bias_daymet_mus, climod_mus_comp, mus_iem_comp) 


pivoted <- comp_met_climod %>%
  pivot_wider(
    id_cols = c(year, day),
    names_from = source,
    values_from = c(maxt, mint, rain)
  )

# ---- 2. Function to calculate stats ----
error_stats <- function(sim, obs) {
  bias <- mean(sim - obs, na.rm = TRUE)
  mae  <- mean(abs(sim - obs), na.rm = TRUE)
  rmse <- sqrt(mean((sim - obs)^2, na.rm = TRUE))
  corr <- cor(sim, obs, use = "complete.obs")
  return(c(bias = bias, mae = mae, rmse = rmse, corr = corr))
}

# ---- 3. Apply function for each variable and source ----
results <- list(
  maxt_power   = error_stats(pivoted$maxt_NASA,   pivoted$maxt_Climod),
  maxt_daymet  = error_stats(pivoted$maxt_Daymet,  pivoted$maxt_Climod),
  maxt_iem     = error_stats(pivoted$maxt_IEM, pivoted$maxt_Climod),
  mint_power   = error_stats(pivoted$mint_NASA,   pivoted$mint_Climod),
  mint_daymet  = error_stats(pivoted$mint_Daymet,  pivoted$mint_Climod),
  mint_iem = error_stats(pivoted$mint_IEM, pivoted$mint_Climod),
  rain_power   = error_stats(pivoted$rain_NASA,   pivoted$rain_Climod),
  rain_daymet  = error_stats(pivoted$rain_Daymet,  pivoted$rain_Climod),
  rain_iem = error_stats(pivoted$rain_IEM, pivoted$rain_Climod)
)

# ---- 4. Convert to a clean summary table ----
results_df_mus <- do.call(rbind, results) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("comparison") %>%
  separate(comparison, into = c("variable","dataset"), sep = "_")

mus_table <- kable(results_df_mus) %>%
  kable_classic()
save_kable(aars_table, file = "Figures/mus_met_error.png", zoom = 2)


as.data.frame(comp_met_climod) %>%
  ggplot(aes(x = day, y = rain, colour = source)) +
  geom_point() +
  facet_wrap(~year)

power_mus_comp <- power_mus %>%
  select(colnames(mus_iem))
mus_iem <- get_iem_apsim_met(lonlat = c(-76.65,42.73), dates = c("1984-01-01","2024-12-31"), state = "NY")

write_apsim_met(mus_iem, wrt.dir = "~/Downloads", filename = "mus_iem.met")
write_csv(mus_iem, file = "R_data/Weather_data/mus_iem.csv")

mus_iem_comp <- mus_iem %>%
  mutate(year_day = paste0(year,"-",day)) %>%
  filter(year_day %in% climod_mus_filtered$year_day) %>%
  dplyr::select(year, day, mint, maxt, rain, year_day) %>%
  mutate(source = rep("IEM"))
  

summary(mus_iem)
iem_mus_rain <- plot(mus_iem, years = c(2015, 2016, 2023, 2024),cumulative = T, climatology = T, met.var = "rain")
plot(mus_iem, years = c(2015,2016,2023,2024), cumulative = T, climatology = T)

mus_iem_gdd <- mus_iem %>%
  filter(year %in% c(2015,2016,2023,2024)) %>%
  mutate(gdd = case_when(
    ((maxt+mint)/2)-10 > 0 ~ ((maxt+mint)/2)-10,
    .default = 0
  )) %>%
  group_by(year) %>%
  mutate(cum_gdd = cumsum(gdd))

mus_iem_gdd %>%
  ggplot(aes(x = day, y = cum_gdd)) +
  facet_wrap(~year) +
  geom_line() +
  ggpubr::theme_pubr()

iem_mus_rain %>%
  ggplot()

comp <- compare_apsim_met(mus_iem, power_mus_comp, met.var = "all", labels = c("iem","pwr"))

aars_iem <- get_iem_apsim_met(lonlat = c(-89.34367303797328,43.301893058200726), dates = c("1984-01-01","2024-12-31"), state = "WI")
write_apsim_met(aars_iem, wrt.dir = "~/Downloads", filename = "aars_iem.met")
write_csv(aars_iem, file = "R_data/Weather_data/aars_iem.csv")


plot(aars_iem, years = c(2015,2017,2018,2019,2020,2023), cumulative = T, climatology = T)

power_aars_comp <- power_aars %>%
  select(-c(rh, windspeed)) %>%
  mutate(source = rep("NASA"))

comp_dm_aars <- arranged_daymet_aars %>%
  mutate(source = rep("Daymet"))

comp_aars_iem <- aars_iem %>%
  mutate(source = rep("IEM"))


climod_aars <- read_csv("R_data/Weather_data/arlington_climod.csv")

climod_aars_filtered <- climod_aars %>%
  filter(Precipitation != c("M", "T", "S"),
         !endsWith(Precipitation,"A"),
         MaxTemperature != "M",
         MinTemperature != "M") %>%
  mutate(rain = as.numeric(Precipitation)*25.4,
         maxt = (as.numeric(MaxTemperature)-32)*(5/9),
         mint = (as.numeric(MinTemperature)-32)*(5/9),
         date = as.Date(Date, format = "%m/%d/%y"),
         year = year(date),
         day = yday(date)) %>%
  dplyr::select(year, day, maxt, mint, rain) %>%
  mutate(year_day = paste0(year,"-",day)) %>%
  drop_na()

bias_power_aars <- power_aars_comp %>%
  mutate(year_day = paste0(year,"-",day)) %>%
  filter(year_day %in% climod_aars_filtered$year_day) %>%
  dplyr::select(year, day, mint, maxt, rain, year_day) %>%
  mutate(source = rep("NASA"))

bias_daymet_aars <- comp_dm_aars %>%
  mutate(year_day = paste0(year,"-",day)) %>%
  filter(year_day %in% climod_aars_filtered$year_day) %>%
  dplyr::select(year, day, mint, maxt, rain, year_day) %>%
  mutate(source = rep("Daymet"))

climod_aars_comp <- climod_aars_filtered %>%
  mutate(source = rep("Climod")) %>%
  arrange(year_day)

aars_iem_comp <- comp_aars_iem %>%
  mutate(year_day = paste0(year,"-",day)) %>%
  filter(year_day %in% climod_aars_filtered$year_day) %>%
  dplyr::select(year, day, mint, maxt, rain, year_day) %>%
  mutate(source = rep("IEM"))

comp_met_aars <- rbind(bias_power_aars, aars_iem_comp, bias_daymet_aars, climod_aars_comp)



pivoted <- comp_met_aars %>%
  pivot_wider(
    id_cols = c(year, day),
    names_from = source,
    values_from = c(maxt, mint, rain)
  )

# ---- 2. Function to calculate stats ----
error_stats <- function(sim, obs) {
  bias <- mean(sim - obs, na.rm = TRUE)
  mae  <- mean(abs(sim - obs), na.rm = TRUE)
  rmse <- sqrt(mean((sim - obs)^2, na.rm = TRUE))
  corr <- cor(sim, obs, use = "complete.obs")
  return(c(bias = bias, mae = mae, rmse = rmse, corr = corr))
}

# ---- 3. Apply function for each variable and source ----
results <- list(
  maxt_power   = error_stats(pivoted$maxt_NASA,   pivoted$maxt_Climod),
  maxt_daymet  = error_stats(pivoted$maxt_Daymet,  pivoted$maxt_Climod),
  maxt_iem = error_stats(pivoted$maxt_IEM, pivoted$maxt_Climod),
  mint_power   = error_stats(pivoted$mint_NASA,   pivoted$mint_Climod),
  mint_daymet  = error_stats(pivoted$mint_Daymet,  pivoted$mint_Climod),
  mint_iem = error_stats(pivoted$mint_IEM, pivoted$mint_Climod),
  rain_power   = error_stats(pivoted$rain_NASA,   pivoted$rain_Climod),
  rain_daymet  = error_stats(pivoted$rain_Daymet,  pivoted$rain_Climod),
  rain_iem = error_stats(pivoted$rain_IEM, pivoted$rain_Climod)
)

# ---- 4. Convert to a clean summary table ----
results_df_aars <- do.call(rbind, results) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("comparison") %>%
  separate(comparison, into = c("variable","dataset"), sep = "_")

library(kableExtra)

aars_table <- kable(results_df_aars) %>%
  kable_classic()
save_kable(aars_table, file = "Figures/aars_met_error.png", zoom = 2)
