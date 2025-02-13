library(tidyverse)
library(emmeans)
library(lme4)
library(modelr)
library(ggpmisc)
library(lmerTest)
library(glmmTMB)
library(DHARMa)
library(shades)
library(meta)
library(metafor)
library(here)

# Bring in Soybean and rye data

soy_rye_raw <- read_csv(here("R_data/Soy_rye_data/NY_WI_NTvT_master.csv"))

soy_rye_raw %>%
  filter(Till_NT == "NT") %>%
  filter(RyeTerminationMethod != "Roller Crimped") %>%
  select(ExperimentName, Year, RyeTerminationMethod) #### 2022 includes mowing and rolling as methods of killing cereal rye

soy_rye_std <- soy_rye_raw %>%
  filter(Year != 2022) %>%
  filter(ExperimentName != "OCS")

all_soy_rye <- soy_rye_raw


# For Loop to change the number of days I am looking at for weather window


num_days <- seq(1,20,1)
planting_prcp_dataframes <- list()
for(i in 1:length(num_days)){
 
  mus_met <- read_csv("~/Documents/clean_musgrave_prism.csv")
  gen_met <- read_csv("~/Documents/clean_geneva_prism.csv")
  wi_met <- read_csv("~/Documents/clean_arlington_prism.csv")

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
    mutate(plant_prcp = sum(prcp[day_of_year >= plant_day_of_year - i & day_of_year <= plant_day_of_year + i], na.rm = TRUE)) %>%
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
  

  gen_met <- gen_met %>%
    rename("Year" = year)
  
  gen_met_exp <- left_join(gen_met, geneva_planting_dates, by = "Year") 
  gen_met_exp <- gen_met_exp %>%
    filter(!is.na(SoySeedingDate)) %>%
    mutate(SoySeedingDate = as.Date(SoySeedingDate, format = "%m/%d/%y")) %>%
    mutate(plant_day_of_year = yday(SoySeedingDate)) %>%
    group_by(SoySeedingDate) %>%
    mutate(plant_prcp = sum(prcp[day_of_year >= plant_day_of_year - i & day_of_year <= plant_day_of_year + i], na.rm = TRUE)) %>%
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
    filter(SoySeedingDate != "6/01/17")
  

  
  wi_met <- wi_met %>%
    rename("Year" = year)
  
  wi_met_2021 <- wi_met %>%
    filter(Year == 2021)
  
  wi_met_exp <- left_join(wi_met, wi_pl_dates, by = "Year")
  
  wi_met_exp <- wi_met_exp %>%
    filter(!is.na(SoySeedingDate)) %>%
    mutate(SoySeedingDate = as.Date(SoySeedingDate, format = "%m/%d/%y")) %>%
    mutate(plant_day_of_year = yday(SoySeedingDate)) %>%
    group_by(SoySeedingDate) %>%
    mutate(plant_prcp = sum(prcp[day_of_year >= plant_day_of_year - i & day_of_year <= plant_day_of_year + i], na.rm = TRUE)) %>%
    mutate(gdd_plant = sum((((maxt+mint)/2)-base_temp_f)[day_of_year >= plant_day_of_year - 10 & day_of_year <= plant_day_of_year + 10], na.rm = TRUE)) %>%
    ungroup() %>%
    distinct(SoySeedingDate, Year, ExperimentName, plant_prcp, gdd_plant)
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
    filter(SoySeedingDate != "2019-06-10")
  
  wi_met_yield$Year[wi_met_yield$Year == 2021 & wi_met_yield$ExperimentName == "Seeding Depth Study"] <- "2021a" 
  wi_met_yield$Year[wi_met_yield$Year == 2021 & wi_met_yield$ExperimentName == "Starter Fertilizer Study"] <- "2021b" 
  
  wi_met_yield
  
  
  comb_met_yield <- rbind(wi_met_yield, ny_met_yield)
  
  comb_met_yield <- comb_met_yield %>%
    distinct(Year, Location, yield_sd, ExperimentName, plant_prcp, Till_NT, mean_yield, gdd_plant)
  
  plant_prcp_exp <- comb_met_yield %>%
    distinct(Year, ExperimentName, Location, plant_prcp, gdd_plant) %>%
    mutate(prcp_cm = plant_prcp*2.54)
  
  planting_prcp_dataframes[[i]] <- plant_prcp_exp
  
}


walk2(planting_prcp_dataframes, seq_along(planting_prcp_dataframes), 
      ~write_csv(.x, str_c(here("R_data/Weather_data/PRISM/plant_prcp_"), .y, ".csv"), col_names = TRUE))


planting_prcp <- read_csv("R_data/Weather_data/PRISM/plant_prcp_all.csv")


rsq_weather_window <- data.frame(i = numeric(0), rsq_value = numeric(0))

for(i in 1:15){
  
  planting_prcp <- planting_prcp_dataframes[[i]]
  
  soy_rye_useful <- soy_rye_raw %>%
    select(Year, State, Location, ExperimentName, Rep, SoySeedingDate, Till_NT,
           RyeBiomass_kgha, SoyStand_ha, Per_Emergence, WeedBiomass_kgha, Yield_Mgha, Yield_buac, ExperimentalTrt) %>%
    filter(Year != 2022)
  
  wi_2014 <- soy_rye_useful %>%
    filter(State == "WI") %>%
    filter(Year == 2014) %>%
    mutate(plant_prcp = rep(NA)) %>%
    mutate(gdd_plant = rep(NA))
  
  soy_rye_useful$Year[soy_rye_useful$Year == 2021 & soy_rye_useful$ExperimentName == "Seeding Depth Study"] <- "2021a" 
  soy_rye_useful$Year[soy_rye_useful$Year == 2021 & soy_rye_useful$ExperimentName == "Starter Fertilizer Study"] <- "2021b" 
  
  soy_rye_useful <- soy_rye_useful %>%
    filter(SoySeedingDate != "2019-06-05")
  comb_soy_prcp <- left_join(soy_rye_useful, planting_prcp)
  
  comb_soy_prcp <- rbind(comb_soy_prcp) 
  
  comb_soy_prcp <- comb_soy_prcp %>%
    mutate(block = paste0(Rep,"-" ,State, "-", Location,"-",Year)) %>%
    mutate(ExperimentName = trimws(ExperimentName))
  
  comb_soy_prcp_revised <- comb_soy_prcp %>%
    filter(Location == "Musgrave" | Location == "AARS") %>% 
    mutate(site_year = trimws(paste0(State, "-", Year))) %>%
    filter(site_year != "WI-2014" & site_year != "WI-2022")
  
  barns <- comb_soy_prcp_revised %>%
    filter(ExperimentName == "BARNS") %>%
    mutate(block = Rep)
  
  barns_mod3 <- lmer(log(Yield_Mgha) ~ Till_NT*Year + (1|Year:block), data = barns) #Using this model for the table below as of now
  summary(barns_mod3)
  hist(resid(barns_mod3))
  plot(resid(barns_mod3)~predict(barns_mod3))
  barns_emm <- pairs(emmeans(barns_mod3, ~ Till_NT|Year))
  
  
  anova(barns_mod3, ddf = "Kenward-Roger")
  
  barns_model_data <- data.frame(
    count_nt = c(4, 4),
    count_till = c(4, 2),
    ExperimentName = c("BARNS", "BARNS"),
    state = c("NY","NY"),
    Year = c(2013, 2014),
    plant_prcp = c(planting_prcp$prcp_cm[planting_prcp$ExperimentName=="BARNS" & planting_prcp$Year==2013],
                   planting_prcp$prcp_cm[planting_prcp$ExperimentName=="BARNS" & planting_prcp$Year==2014]),
    mean_yield_gap = c(-0.0339, -0.2161),
    SE_TE = c(0.0681, 0.0900)
  )
  
  # NY2015 and NY2016
  moss <- comb_soy_prcp %>%
    filter(ExperimentName == "MOSS") %>%
    mutate(block = paste0(Rep,"-",Year))
  
  moss$Yield_Mgha
  
  moss_mod <- lmer(log(Yield_Mgha) ~ Till_NT*Year +  (1|block), data = moss)
  summary(moss_mod)
  hist(resid(moss_mod))
  plot(resid(moss_mod)~predict(moss_mod))
  moss_emmeans <- emmeans(moss_mod, ~Till_NT|Year)
  pairs(moss_emmeans)
  
  moss_model_data <- data.frame(
    count_nt = c(4, 4),
    count_till = c(4, 4),
    ExperimentName = c("MOSS", "MOSS"),
    state = c("NY","NY"),
    Year = c(2015, 2016),
    plant_prcp = c(planting_prcp$prcp_cm[planting_prcp$ExperimentName=="MOSS" & planting_prcp$Year==2015],
                   planting_prcp$prcp_cm[planting_prcp$ExperimentName=="MOSS" & planting_prcp$Year==2016]),
    mean_yield_gap = c(-0.0317, -0.6562),
    SE_TE = c(0.161, 0.161)
  )
  
  # WI2017, WI2018, WI2019, and WI2020
  rye_varieties_exp <- comb_soy_prcp %>%
    filter(ExperimentName == "CCBRT Soybean" | SoySeedingDate == "6/6/18" | Year == 2017) %>%
    mutate(block = paste0(Rep, "-", ExperimentName, "-", Year))
  
  rye_varieties_mod <- lmer(log(Yield_Mgha) ~ Till_NT*Year + 
                              (1|Year:ExperimentName:Rep), 
                            data = rye_varieties_exp)
  summary(rye_varieties_mod)
  hist(resid(rye_varieties_mod))
  plot(resid(rye_varieties_mod)~predict(rye_varieties_mod))
  rye_varieties_emmeans <- emmeans(rye_varieties_mod, ~Till_NT|Year)
  pairs(rye_varieties_emmeans)
  
  rye_var_model_data <- data.frame(
    count_nt = c(8, 24, 12, 16),
    count_till = c(4, 12, 4, 4),
    ExperimentName = c("NTSoy17", "NTSoy18", "CCRBT", "CCRBT"),
    state = c("WI","WI", "WI", "WI"),
    Year = c(2017, 2018, 2019, 2020),
    plant_prcp = c(planting_prcp$prcp_cm[planting_prcp$Location=="AARS" & planting_prcp$Year==2017],
                   planting_prcp$prcp_cm[planting_prcp$ExperimentName=="No-Till Soybeans" & planting_prcp$Year==2018],
                   planting_prcp$prcp_cm[planting_prcp$Location=="AARS" & planting_prcp$Year==2019], 
                   planting_prcp$prcp_cm[planting_prcp$ExperimentName=="CCBRT Soybean" & planting_prcp$Year==2020]),
    mean_yield_gap = c(-0.0835, -0.0397,-0.1341 ,-0.0336),
    SE_TE = c(0.0687, 0.0397, 0.0648, 0.0627)
  )
  
  # NY2018
  ocs <- comb_soy_prcp %>%
    filter(ExperimentName == "OCS") %>%
    mutate(block = Rep)
  
  ocs_mod <- lmer(log(Yield_Mgha) ~ Till_NT + (1|block:ExperimentalTrt), data = ocs)
  hist(resid(ocs_mod))
  plot(resid(ocs_mod)~predict(ocs_mod))
  ocs_emmeans <- emmeans(ocs_mod, ~Till_NT)
  pairs(ocs_emmeans)
  
  ocs_model_data <- data.frame(
    count_nt = c(32),
    count_till = c(32),
    ExperimentName = c("OCS"),
    state = c("NY"),
    Year = c(2018),
    plant_prcp = c(planting_prcp$prcp_cm[planting_prcp$ExperimentName=="OCS" & planting_prcp$Year==2018]),
    mean_yield_gap = c(-0.631),
    SE_TE = c(0.103)
  )
  
  # Mowtivation NY2023 and 2024
  mowtivation <- comb_soy_prcp %>%
    filter(ExperimentName == "Mowtivation") %>%
    mutate(block = Rep)
  
  mow_model <- lmer(log(Yield_Mgha) ~ Till_NT*Year + (1|Year:block:ExperimentalTrt), data = mowtivation)
  hist(resid(mow_model))
  plot(resid(mow_model)~predict(mow_model))
  mow_emmeans <- emmeans(mow_model, ~Till_NT|Year)
  pairs(mow_emmeans)
  
  mow_model_data <- data.frame(
    count_nt = c(18,18),
    count_till = c(12,12),
    ExperimentName = c("Mowtivation", "Mowtivation"),
    state = c("NY", "NY"),
    Year = c(2023,2024),
    plant_prcp = c(planting_prcp$prcp_cm[planting_prcp$ExperimentName=="Mowtivation" & planting_prcp$Year==2023],
                   planting_prcp$prcp_cm[planting_prcp$ExperimentName=="Mowtivation" & planting_prcp$Year==2024]),
    mean_yield_gap = c(0.04269, -0.00162),
    SE_TE = c(0.0488, 0.0348)
  )
  
  # WI2015
  wi2015 <- comb_soy_prcp %>%
    filter(Year == 2015 & State == "WI")%>%
    mutate(block = Rep)
  
  wi2015_mod <- lmer(log(Yield_Mgha) ~ Till_NT + (1|block), data = wi2015)
  hist(resid(wi2015_mod))
  plot(predict(wi2015_mod)~resid(wi2015_mod))
  wi2015_emmeans <- emmeans(wi2015_mod, ~Till_NT)
  pairs(wi2015_emmeans)
  
  wi2015_model_data <- data.frame(
    count_nt = c(9),
    count_till = c(3),
    ExperimentName = c("NTSoy2015"),
    state = c("WI"),
    Year = c(2015),
    plant_prcp = c(planting_prcp$prcp_cm[planting_prcp$ExperimentName=="No-Till Soybeans" & planting_prcp$Year==2015]),
    mean_yield_gap = c(-0.126),
    SE_TE = c(0.0279)
  )
  
  # WI2021b starter fert
  wi_sf <- comb_soy_prcp %>%
    filter(Year == "2021b" & State == "WI")%>%
    mutate(block = Rep)
  
  wi_sf_mod <- lmer(log(Yield_Mgha) ~ Till_NT + (1|block), data = wi_sf)
  hist(resid(wi_sf_mod))
  plot(predict(wi_sf_mod)~resid(wi_sf_mod))
  wi_sf_emmeans <- emmeans(wi_sf_mod, ~Till_NT)
  pairs(wi_sf_emmeans)
  
  wi_sf_model_data <- data.frame(
    count_nt = c(16),
    count_till = c(4),
    ExperimentName = c("SF"),
    state = c("WI"),
    Year = c(2021),
    plant_prcp = c(planting_prcp$prcp_cm[planting_prcp$Year=="2021b"]),
    mean_yield_gap = c(-0.555),
    SE_TE = c(0.0963)
  )
  # WI2021a seeding depth
  wi_seed_d <- comb_soy_prcp %>%
    filter((Year == "2021a"))%>%
    mutate(block = Rep)
  
  wi_seed_d_mod <- lmer(log(Yield_Mgha) ~ Till_NT + (1|block) , data = wi_seed_d)
  hist(resid(wi_seed_d_mod))
  plot(resid(wi_seed_d_mod)~predict(wi_seed_d_mod))
  wi_seed_d_emmeans <- emmeans(wi_seed_d_mod, ~Till_NT)
  pairs(wi_seed_d_emmeans)
  
  wi_sd_model_data <- data.frame(
    count_nt = c(20),
    count_till = c(4),
    ExperimentName = c("SeedD"),
    state = c("WI"),
    Year = c(2021),
    plant_prcp = c(planting_prcp$prcp_cm[planting_prcp$ExperimentName=="Seeding Depth Study" & planting_prcp$Year=="2021a"]),
    mean_yield_gap = c(-0.513),
    SE_TE = 0.0954
  )
  # WI2023
  wi_coulter <- comb_soy_prcp %>%
    filter(ExperimentName == "Coulter Study")%>%
    mutate(block = Rep)
  
  wi_coulter_mod <- lmer(log(Yield_Mgha) ~ Till_NT + (1|block), data = wi_coulter)
  hist(resid(wi_coulter_mod))
  plot(predict(wi_coulter_mod)~resid(wi_coulter_mod))
  wi_coulter_emmeans <- emmeans(wi_coulter_mod, ~Till_NT)
  pairs(wi_coulter_emmeans)
  
  coulter_model_data <- data.frame(
    count_nt = c(24),
    count_till = c(6),
    ExperimentName = c("Coulter"),
    state = c("WI"),
    Year = c(2023),
    plant_prcp = c(planting_prcp$prcp_cm[planting_prcp$ExperimentName=="Coulter Study" & planting_prcp$Year==2023]),
    mean_yield_gap = c(-0.39),
    SE_TE = 0.0533
  )
  # Hold out geneva for now ; data is super unclear and poorly labeled in sheet
  genevastudies <- comb_soy_prcp %>%
    filter(Location == "Geneva") %>%
    mutate(block = paste0(block, "-",Year))
  
  gen_mod <- lmer(log(Yield_Mgha) ~ Till_NT*Year + (1|block), data = genevastudies)
  hist(resid(gen_mod))
  plot(resid(gen_mod)~predict(gen_mod))
  gen_emmeans <- emmeans(gen_mod, ~Till_NT|Year)
  pairs(gen_emmeans)
  
  gen_model_data <- data.frame(
    count_nt = c(4,8),
    count_till = c(4,8),
    ExperimentName = c("Gen1", "Gen2"),
    state = c("NY", "NY"),
    Year = c(2017, 2019),
    plant_prcp = c(6.05, 9.58),
    SE_nt = c(NA, NA),
    SE_till = c(NA, NA),
    mean_yield_gap = c(-0.110,-0.437),
    nt_yield = c(2.25,NA ),
    till_yield = c(3.32,NA ),
    SE_TE = c(0.171, 0.123)
  )
  
  # May hold WI2014 until more information is found
  wi_2014 <- soy_rye_useful %>%
    filter(State == "WI") %>%
    filter(Year == 2014) %>%
    mutate(plant_prcp = rep(NA)) %>%
    mutate(gdd_plant = rep(NA))
  
  #### Initial meta-analysis
  
  meta_final_data <- rbind(barns_model_data, moss_model_data, ocs_model_data, wi2015_model_data, rye_var_model_data,
                           wi_sf_model_data, wi_sd_model_data, coulter_model_data)
  
  meta_final_data <- meta_final_data %>%
    mutate(site_year = paste0(state, "-", Year)) %>%
    mutate(var_te = SE_TE^2)
  
  metagen_analysis <- metagen(data = meta_final_data, TE = mean_yield_gap, seTE = SE_TE, studlab = site_year)
  summary(metagen_analysis)
  forest(metagen_analysis, )
  
  metagen_result <- data.frame(
    effect_size = metagen_analysis$TE,  # Treatment effect
    ci_lower = metagen_analysis$lower,  # Lower bound of confidence interval
    ci_upper = metagen_analysis$upper,  # Upper bound of confidence interval
    i2 = metagen_analysis$I2            # I² for heterogeneity
  )
  
  metagen_result_comb <- cbind(metagen_result, meta_final_data)
  
  metagen_result_final <- metagen_result_comb %>%
    mutate(yield_ratio = exp(mean_yield_gap),
           lowerci_trans = exp(ci_lower),
           upperci_trans = exp(ci_upper))
  
  rsq_value <- (cor(metagen_result_final$yield_ratio,metagen_result_final$plant_prcp))^2
  
  rsq_weather_window <- rbind(rsq_weather_window, data.frame(i = i, rsq_value = rsq_value))
}

write_csv(rsq_weather_window, "R_data/Soy_rye_data/weather_window_rsq.csv")



#### Maybe look at metafor package for meta-regression
library(metafor)

metafor_metareg <- rma.mv(yi = mean_yield_gap, V = var_te,
                          data = meta_final_data, slab = ExperimentName,
                          random = ~ 1 |ExperimentName/site_year,
                          test = "t", method = "REML",
                          mods = ~ plant_prcp)
summary(metafor_metareg)

exp(predict(metafor_metareg)$pred)

plot(exp(predict(metafor_metareg)$pred)~metafor_metareg$data$plant_prcp)

predicted_metafor <- data.frame(predict(metafor_metareg))

predicted_metafor_ratios <- predicted_metafor %>%
  mutate(yield_ratio = exp(pred),
         ci_lower_trans = exp(ci.lb),
         ci_upper_trans = exp(ci.ub)) 
predicted_metafor_ratios <- cbind(predicted_metafor_ratios, metafor_metareg$data$plant_prcp)

ggplot(predicted_metafor_ratios, aes(x = `metafor_metareg$data$plant_prcp`, y = yield_ratio))+
  geom_line() +
  geom_point() +
  geom_ribbon(aes(ymin = ci_lower_trans, ymax = ci_upper_trans, alpha = 0.5), show.legend = F) +
  theme_bw() +
  labs(y = "Yield ratio between no-till to till yield",
       x = "Precipitation accumulation +/- 15 days of planting (cm)",
       alpha = "") +
  stat_poly_eq(use_label(c("Eq","R2"))) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", alpha = 0.85)


print(predicted_metafor_ratios$pred)

class(metafor_metareg)

rainfall_sequence <- seq(2,15, 1)
coefs_metafor <- coef(metafor_metareg)  
predictions_metareg <- (coefs_metafor[2]*rainfall_sequence)+(coefs_metafor[1])
predictions_ratio <- exp(predictions_metareg)

predicted_df <- data.frame(predictions_ratio, rainfall_sequence)

ggplot(predicted_df, aes(x = rainfall_sequence, y = predictions_ratio)) +
  geom_point() +
  geom_line() +#### Line from predictions here; may want to add standard errors/CI
  theme_bw() +
  labs(y = "Predicted yield ratio", x = "Sequenced rainfall (cm)") 

mreg_df <- initial_mreg$data %>%
  select(ExperimentName, site_year, .TE, .seTE, plant_prcp) %>%
  mutate(yield_ratio = exp(.TE))

mreg_df %>%
  ggplot(aes(x = plant_prcp, y = yield_ratio, fill = site_year)) +
  geom_point()

metagen_result_final %>%
  filter(site_year == "NY-2015")
library(ggpubr)

ratio_yg_meta_check <- metagen_result_final %>%
  ggplot(aes(x = plant_prcp, y = yield_ratio))+
  geom_point(aes(color = state), size = 2) +
  geom_errorbar(aes(ymin = lowerci_trans, ymax = upperci_trans, colour = state), width = 0.2) + # Error bars for CI
  labs(x = "Precipitation accumulation +/- 15 days of planting (cm)",
       y = "Ratio of no-till yield to tilled yield",
       color = "") +
  geom_text(aes(label = site_year)) +
  stat_poly_eq(formula = y~x, use_label(c("R2", "eq")), size = 5)+
  theme_pubr() +
  scale_color_manual(values = c("black", "red")) +
  geom_smooth(method = "lm",formula = y ~ x,  se = FALSE, color = "gray", 
              linetype = "dashed")  +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))

ratio_yg_meta_english <- metagen_result_final %>%
  ggplot(aes(x = plant_prcp/2.54, y = yield_ratio))+
  geom_point(aes(color = state), size = 2) +
  geom_errorbar(aes(ymin = lowerci_trans, ymax = upperci_trans, colour = state), width = 0.1) + # Error bars for CI
  labs(x = "Precipitation accumulation +/- 15 days of planting (in)",
       y = "Ratio of no-till yield to tilled yield",
       color = "") +
  stat_poly_eq(formula = y~x, use_label(c("R2", "eq")), size = 5)+
  theme_pubr() +
  scale_color_manual(values = c("black", "red")) +
  geom_smooth(method = "lm",formula = y ~ x,  se = FALSE, color = "gray", 
              linetype = "dashed")  +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))

ggsave(ratio_yg_meta_check, filename = here("Figures/ratio_yg_metaplot_15day_metric.pdf"), dpi = 600, useDingbats = FALSE,
       height = 7, width = 9, units = "in")
ggsave(ratio_yg_meta_english, filename = here("Figures/ratio_yg_metaplot_15day_english.pdf"), dpi = 600, useDingbats = FALSE,
       height = 7, width = 9, units = "in")
