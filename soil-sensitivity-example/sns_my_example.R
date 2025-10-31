##### Doing sensitivity with my simulations based on the code Caio provided

## date: 2025-10-11
## apsim version: 2025.3.7681.0

library(apsimx)
library(ggplot2)
library(tidyverse)
library(ggpubr)


## running the simulation to make sure things work before
## we tweak them
lima_check <- get_ssurgo_soil_profile(lonlat = c(-76.65090, 42.73400))
plano_check <- get_ssurgo_soil_profile(lonlat = c(-89.34221, 43.30640), nlayers = 10)


lima_check[[1]]$soil$Soybean.XF <- c(1,1,1,1,0,0,0,0,0,0)
lima_check[[1]]$soil$Wheat.XF <- c(1,1,1,1,0,0,0,0,0,0)
lima_check[[1]]$soil$NO3N <- c(2,2,2,2,0,0,0,0,0,0)


plano_check[[1]]$soil$NO3N <- c(2,2,2,2,0.5,0.5,0.5,0.5,0.5,0.5)
plano_check[[1]]$soil$SAT <- c(0.544, 0.544, 0.544, 0.461, 0.461, 0.403, 0.403, 0.403, 0.389, 0.389)


lima_base <- lima_check[[1]]
plano_base <- plano_check[[1]]
sim0 <- apsimx('ny_soils_sensitivity_2.apsimx',
               cleanup = TRUE)
sim_name <- 'ny_soils_sensitivity_2-backup.apsimx'


check_apsimx_soil_profile(plano_base)
## it seems like things worked :)
head(sim0)

## creating different soil profiles


base.soil.phys <- extract_data_apsimx(sim_name, node = "Soil", root = "NY 2016", soil.child = "Physical")
base.soil.organic <- extract_data_apsimx(sim_name, node = "Soil", root = "NY 2016", soil.child = "Organic")
base.soil.no3 <- extract_data_apsimx(sim_name, node = "Soil", root = "NY 2016", soil.child = "NO3")
base.soil.chem <- extract_data_apsimx(sim_name, node = "Soil", root = "NY 2016", soil.child = "Chemical")
base.soil.meta <- extract_data_apsimx(sim_name, node = "Soil", root = "NY 2016", soil.child = "Metadata")
base.soil.soilwat <- base.soil.meta <- extract_data_apsimx(sim_name, node = "Soil", root = "NY 2016", soil.child = "SoilWater")
base.soil.initialwater <- base.soil.meta <- extract_data_apsimx(sim_name, node = "Soil", root = "NY 2016", soil.child = "InitialWater")



base.soil <- apsimx_soil_profile(nlayers = 7,
                    Thickness = base.soil.phys$soil.layers$Thickness,
                    ParticleSizeSand = base.soil.phys$soil.layers$ParticleSizeSand,
                    ParticleSizeSilt = base.soil.phys$soil.layers$ParticleSizeSilt,
                    ParticleSizeClay = base.soil.phys$soil.layers$ParticleSizeClay,
                    BD = base.soil.phys$soil.layers$BD,
                    AirDry = base.soil.phys$soil.layers$AirDry,
                    LL15 = base.soil.phys$soil.layers$LL15,
                    DUL = base.soil.phys$soil.layers$DUL,
                    SAT = base.soil.phys$soil.layers$SAT,
                    KS = base.soil.phys$soil.layers$KS,
                    crop.LL = base.soil.phys$crop$`Soybean LL`,
                    crop.XF = base.soil.phys$crop$`Soybean XF`,
                    Carbon = c(1.78,1.78,0.87,0.87,0.87,0.87,0.87),
                    SoilCNRatio = rep(8, 7),
                    FOM = c(100, 70, 60, 20, 10, 5, 2),
                    FOM.CN = 40,
                    FBiom = c(0.040, 0.031, 0.014, 0.010, 0.010, 0.010, 0.010),
                    FInert = c(0.193, 0.395, 0.667, 0.832, 0.990, 0.990, 0.990),
                    NO3N = c(17, 14, 12, 9, 6.5, 5, 3.5),
                    NH4N = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1),
                    PH = c(6.5,6.5,6.5,7,7,7,7),
                    soil.bottom = 200,
                    crops = c("Maize", "Wheat", "Soybean"),
                    metadata = base.soil.meta$soil)

lima_water_plot_df <- data.frame("depth" = base.soil$soil$Depth,
                                 "ll15" = base.soil$soil$LL15,
                                 "dul"  = base.soil$soil$DUL) 

lima_base_plot <- plot(base.soil, property = "water")

lima_base_plot +
  theme_bw() +
  labs(y = "Volumetric water content", x = "Soil depth (cm)") +
  geom_vline(xintercept = -86, color = "red", linetype = "dashed") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14))



## creating a function to make this less painful
## the trick when building a function like this is ensuring functional
## soil profiles. For example, by definition DUL < SAT, so when varying one, 
## it's important to check that it's not crossing any boundaries.
modify_soil_property <- function(soil.profile,
                                 soil.property = c('LL15', 'SAT', 'DUL', 'KS', "NO3N"), ## this can be expanded
                                 multiplier){
  soil.property <- match.arg(soil.property)
  new.soil <- soil.profile
  new.soil$soil[[soil.property]] <- new.soil$soil[[soil.property]] * multiplier 
  return(new.soil)
}



#### For DUL ####
sens.grid <- data.frame(DUL.Multiplier = c(0.8, 0.9, 1, 1.1))
sens.grid$soil.profile <- 1:nrow(sens.grid)

## creating different soil profiles for following the grid
soils.list <- lapply(sens.grid[['DUL.Multiplier']],
                     FUN = function(x){
                       modify_soil_property(base.soil,
                                            soil.property = "DUL",
                                            multiplier = x)
                     })

ending <- paste0("\\.", "apsimx","$")
sim_names <- list.files("MySims/NY", pattern = ending)
sim_names <- sim_names[sim_names != "ny2015_sim-backup.apsimx"]


sens <- list()
sens_df <- data.frame()


for (sim in sim_names) {
    sens[[sim]] <- sens_apsimx(
      file = sim,
      src.dir = "/Users/stw52/Downloads/soil-sensitivity-example/MySims/NY",
      grid = sens.grid[, 2, drop = FALSE],
      parm.paths = "soil.profile",
      summary = "none",
      soil.profiles = soils.list
    )
  sim_data <- sens[[sim]]$grid.sims
  
  sens_df <- rbind(sens_df, sim_data)  
  sens_merged <- merge(sens_df, sens.grid)
}

sens_output <- sens_merged %>%
  mutate(Year = year(Clock.Today),
         Month = month(Clock.Today),
         tillage = case_when(
           str_detect(Zone, "NT") ~ "NT",
           str_detect(Zone, "Till") ~ "Till"
         ))
sens_plot_dul <- sens_output %>%
  group_by(Zone) %>%
  mutate(Year_rank = dense_rank(Year)) %>%
  ungroup() %>%
  filter(Year_rank == 2, Month %in% 7:11) %>%
  ggplot(aes(x = Date,
         y = soy_kg_ha,
         colour = as.character(DUL.Multiplier)))+
  geom_line() +
  facet_grid(rows = vars(tillage), cols = vars(Year), scales = "free_x") +
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45)) +
  labs(color = "", y = "Soybean yield (kg/ha)") 


yield_summary <- sens_merged %>%
  mutate(Year = year(Clock.Today),
         tillage = case_when(
           str_detect(Zone, "NT") ~ "NT",
           str_detect(Zone, "Till") ~ "Till"
         )) %>%
  group_by(Zone, tillage, Year, DUL.Multiplier) %>%
  summarise(max_yield = max(soy_kg_ha, na.rm = TRUE), .groups = "drop") %>%
  filter(max_yield > 0)

sens_results <- yield_summary %>%
  group_by(tillage, DUL.Multiplier) %>%
  summarise(mean_yield = mean(max_yield, na.rm = TRUE),
            sd_yield = sd(max_yield, na.rm = TRUE),
            .groups = "drop")

ggplot(sens_results, aes(x = DUL.Multiplier, y = mean_yield, color = tillage)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean_yield - sd_yield, ymax = mean_yield + sd_yield), width = 0.02) +
  theme_minimal() +
  labs(x = "DUL multiplier", y = "Final soybean yield (kg/ha)")


sens_index <- yield_summary %>%
  arrange(tillage, DUL.Multiplier) %>%
  group_by(tillage) %>%
  mutate(
    delta_y = max_yield -  max_yield[DUL.Multiplier == 1],
    delta_x = DUL.Multiplier - 1,
    sens_index = (delta_y / max_yield[DUL.Multiplier == 1]) / (delta_x / 1)
  ) %>%
  drop_na()

sens_index_summary <- sens_index %>%
  group_by(tillage) %>%
  summarise(mean_S = mean(sens_index, na.rm = TRUE))

baseline_yield <- yield_summary %>%
  filter(DUL.Multiplier == 1) %>%
  select(tillage, base_yield = max_yield)

sens_from_base <- yield_summary %>%
  left_join(baseline_yield, by = "tillage") %>%
  mutate(
    delta_y = max_yield - base_yield,           # change from baseline
    rel_change = delta_y / base_yield           # proportional change
  )

ggplot(sens_from_base, aes(x = DUL.Multiplier, y = rel_change, color = tillage)) +
  #geom_line() +
  geom_point() +
  theme_minimal() +
  labs(
    x = "DUL multiplier",
    y = "Relative change in yield vs baseline (DUL=1)"
  )

#### For NO3 and Rye ####
sens.grid_no3 <- data.frame(NO3.Multiplier = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2))
sens.grid_no3$soil.profile <- 1:nrow(sens.grid_no3)

## creating different soil profiles for following the grid
soils.list_no3 <- lapply(sens.grid_no3[['NO3.Multiplier']],
                     FUN = function(x){
                       modify_soil_property(base.soil,
                                            soil.property = "NO3N",
                                            multiplier = x)
                     })


sens_no3 <- list()
sens_df_no3 <- data.frame()

for (sim in sim_names) {
  sens[[sim]] <- sens_apsimx(
    file = sim,
    src.dir = "/Users/stw52/Downloads/soil-sensitivity-example/MySims/NY",
    grid = sens.grid_no3[, 2, drop = FALSE],
    parm.paths = "soil.profile",
    summary = "none",
    soil.profiles = soils.list_no3
  )
  sim_data_no3 <- sens[[sim]]$grid.sims
  
  sens_df_no3 <- rbind(sens_df_no3, sim_data_no3)  
  sens_merged_no3 <- merge(sens_df_no3, sens.grid_no3)
}

okabe_ito <- c(
  "#999999",
  "#000000",
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#F0E442", # yellow
  "#0072B2", # blue
  "#D55E00", # vermilion
  "#CC79A7" # reddish purple
)

sens_plot_no3 <- sens_merged_no3 %>%
  filter(str_detect(Zone, "NT"),
         rye_bm_kg_ha > 0) %>%
  mutate(site_year = SimulationName) %>%
  ggplot(aes(x = Date,
             y = rye_bm_kg_ha,
             colour = as.character(NO3.Multiplier)))+
  geom_line() +
  scale_color_manual(values = okabe_ito) +
  facet_wrap(~site_year, scales = "free_x") +
  theme_classic() +
  theme(legend.position = "top",
        legend.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45),
        legend.box = "vertical",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.box.just = "right",
        legend.text = element_text(size = 12)) +
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 1)) + 
  labs(color = "Soil residual nitrate multiplier", y = "Cereal rye biomass (kg/ha)") 



rye_sens_output_ny <- sens_merged_no3 %>%
  filter(str_detect(Zone, "NT")) %>%
  group_by(Zone, NO3.Multiplier) %>%
  summarise(final_rye = max(rye_bm_kg_ha)) %>%
  separate_wider_delim(Zone, delim = " ", names = c("state", "tillage", "year", "soil_scenario", "plant_scenario")) %>%
  mutate(site_year = paste0(state, year))



#### Trying a modifier for XF

modify_soil_property_xf <- function(soil.profile,
                                    soil.properties = c("Soybean.XF", "Wheat.XF"),
                                    multiplier = NULL,
                                    additive_change = NULL,
                                    bottom_n = 2) {
  new.soil <- soil.profile
  
  # ensure soil exists and has the expected structure
  if (!"soil" %in% names(new.soil) || !is.data.frame(new.soil$soil)) {
    stop("The provided soil profile must contain a 'soil' data frame.")
  }
  
  n_rows <- nrow(new.soil$soil)
  idx <- (n_rows - bottom_n + 1):n_rows
  
  # apply to all provided soil property columns
  for (soil.property in soil.properties) {
    if (!soil.property %in% names(new.soil$soil)) {
      stop(paste("Soil property", soil.property, "not found in soil profile"))
    }
    
    # multiplicative change to entire column
    if (!is.null(multiplier)) {
      new.soil$soil[[soil.property]] <- new.soil$soil[[soil.property]] * multiplier
    }
    
    # additive change to bottom_n layers
    if (!is.null(additive_change)) {
      new.soil$soil[[soil.property]][idx] <-
        new.soil$soil[[soil.property]][idx] + additive_change
    }
    
    # optional: clip to [0, 1] range
    new.soil$soil[[soil.property]] <- pmax(pmin(new.soil$soil[[soil.property]], 1), 0)
  }
  
  return(new.soil)
}
sens.grid_xf <- data.frame(XF.Add = c(0, 0.25, 0.5, 0.75, 1))
sens.grid_xf$soil.profile <- 1:nrow(sens.grid_xf)

soils.list_xf <- lapply(sens.grid_xf[['XF.Add']], function(x) {
  modify_soil_property_xf(
    base.soil,
    soil.properties = c("Soybean.XF", "Wheat.XF"),
    additive_change = x
  )
})
ending <- paste0("\\.", "apsimx","$")
sim_names <- list.files("MySims/NY", pattern = ending)
sim_names <- sim_names[sim_names != "ny2015_sim-backup.apsimx"]


sens_xf <- list()
sens_df_xf <- data.frame()


for (sim in sim_names) {
  sens_xf[[sim]] <- sens_apsimx(
    file = sim,
    src.dir = "/Users/stw52/Downloads/soil-sensitivity-example/MySims/NY",
    grid = sens.grid_xf[, 2, drop = FALSE],
    parm.paths = "soil.profile",
    summary = "none",
    soil.profiles = soils.list_xf
  )
  sim_data_xf <- sens_xf[[sim]]$grid.sims
  
  sens_df_xf <- rbind(sens_df_xf, sim_data_xf)  
  sens_merged_xf <- merge(sens_df_xf, sens.grid_xf)
}

sens_output_xf <- sens_merged_xf %>%
  mutate(Year = year(Clock.Today),
         Month = month(Clock.Today),
         tillage = case_when(
           str_detect(Zone, "NT") ~ "NT",
           str_detect(Zone, "Till") ~ "Till"
         ))
sens_plot_xf <- sens_output_xf %>%
  group_by(Zone) %>%
  mutate(Year_rank = dense_rank(Year)) %>%
  ungroup() %>%
  filter(Year_rank == 2, Month %in% 7:11) %>%
  ggplot(aes(x = Date,
             y = soy_kg_ha,
             colour = as.character(XF.Add)))+
  geom_line() +
  facet_grid(rows = vars(tillage), cols = vars(Year), scales = "free_x") +
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45)) +
  labs(color = "", y = "Soybean yield (kg/ha)") 


####


##### Doing above for Wisconsin #####
#### First looking at NO3, then water params ####
sim_base_wi <- "WI2015_NT.apsimx"

base.soil.phys_wi <- extract_data_apsimx(sim_base_wi, node = "Soil", soil.child = "Physical", src.dir = "MySims/WI")
base.soil.organic_wi <- extract_data_apsimx(sim_base_wi, node = "Soil", soil.child = "Organic", src.dir = "MySims/WI")
base.soil.no3_wi <- extract_data_apsimx(sim_base_wi, node = "Soil", soil.child = "NO3", src.dir = "MySims/WI")
base.soil.chem_wi <- extract_data_apsimx(sim_base_wi, node = "Soil", soil.child = "Chemical", src.dir = "MySims/WI")
base.soil.meta_wi <- extract_data_apsimx(sim_base_wi, node = "Soil", soil.child = "Metadata", src.dir = "MySims/WI")
base.soil.soilwat_wi <- base.soil.meta <- extract_data_apsimx(sim_base_wi, node = "Soil", soil.child = "SoilWater", src.dir = "MySims/WI")
base.soil.initialwater_wi <- base.soil.meta <- extract_data_apsimx(sim_base_wi, node = "Soil", soil.child = "InitialWater", src.dir = "MySims/WI")



base.soil_wi <- apsimx_soil_profile(nlayers = 13,
                                 Thickness = base.soil.phys_wi$soil.layers$Thickness,
                                 ParticleSizeSand = base.soil.phys_wi$soil.layers$ParticleSizeSand,
                                 ParticleSizeSilt = base.soil.phys_wi$soil.layers$ParticleSizeSilt,
                                 ParticleSizeClay = base.soil.phys_wi$soil.layers$ParticleSizeClay,
                                 BD = base.soil.phys_wi$soil.layers$BD,
                                 AirDry = base.soil.phys_wi$soil.layers$AirDry,
                                 LL15 = base.soil.phys_wi$soil.layers$LL15,
                                 DUL = base.soil.phys_wi$soil.layers$DUL,
                                 SAT = base.soil.phys_wi$soil.layers$SAT,
                                 KS = plano_base$soil$KS,
                                 crop.LL = base.soil.phys_wi$crop$`soybean LL`,
                                 crop.XF = base.soil.phys_wi$crop$`soybean XF`,
                                 Carbon = base.soil.organic_wi$second$Carbon,
                                 SoilCNRatio = rep(8, 13),
                                 FOM = base.soil.organic_wi$second$FOM,
                                 FOM.CN = 40,
                                 FBiom = base.soil.organic_wi$second$FBiom,
                                 FInert = base.soil.organic_wi$second$FInert,
                                 #NO3N = c(10,10,10,10,10,10,10,10,10,10,10,10,10),
                                 NO3N = c(18,15,12,9,7,6,5,4,3.5,3,2.8,2.5,2.2),
                                 NH4N = rep(0, 13),
                                 PH = base.soil.chem_wi$second$PH,
                                 soil.bottom = 2510,
                                 crops = c("Maize", "Wheat", "Soybean"),
                                 metadata = base.soil.meta_wi$soil)




sens.grid <- data.frame(NO3.Multiplier = c(0 ,0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2))
sens.grid$soil.profile <- 1:nrow(sens.grid)

## creating different soil profiles for following the grid
soils.list <- lapply(sens.grid[['NO3.Multiplier']],
                     FUN = function(x){
                       modify_soil_property(base.soil_wi,
                                            soil.property = "NO3N",
                                            multiplier = x)
                     })

ending <- paste0("\\.", "apsimx","$")
sim_names <- list.files("MySims/WI", pattern = ending)
sim_names_wi <- sim_names[sim_names != c("WI2015_NT-backup-backup-backup-backup-backup-backup.apsimx", "WI2015_NT-backup-backup-backup-backup-backup.apsimx"       
                                      ,"WI2015_NT-backup-backup-backup-backup.apsimx", "WI2015_NT-backup-backup-backup.apsimx"                     
                                      ,"WI2015_NT-backup-backup.apsimx")]

sens_wi <- list()
sens_df_wi <- data.frame()

for (sim in sim_names_wi) {
  sens_wi[[sim]] <- sens_apsimx(
    file = sim,
    src.dir = "/Users/stw52/Downloads/soil-sensitivity-example/MySims/WI",
    grid = sens.grid[, 2, drop = FALSE],
    parm.paths = "soil.profile",
    summary = "none",
    soil.profiles = soils.list
  )
  sim_data_wi <- sens_wi[[sim]]$grid.sims
  
  sens_df_wi <- rbind(sens_df_wi, sim_data_wi)  
  sens_merged_wi <- merge(sens_df_wi, sens.grid)
}


sens_output_wi <- sens_merged_wi %>%
  mutate(Year = year(Clock.Today),
         Month = month(Clock.Today),
         tillage = case_when(
           str_detect(Zone, "NT") ~ "NT",
           str_detect(Zone, "Till") ~ "Till"
         ))



sens_plot_no3_wi <- sens_merged_wi %>%
  filter(str_detect(Zone, "NT"),
         rye_bm_kg_ha > 0) %>%
  mutate(site_year = SimulationName) %>%
  ggplot(aes(x = Date,
             y = rye_bm_kg_ha,
             colour = as.character(NO3.Multiplier)))+
  geom_line() +
  scale_color_manual(values = okabe_ito) +
  facet_wrap(~site_year, scales = "free_x") +
  theme_classic() +
  theme(legend.position = "top",
        legend.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45),
        legend.box = "vertical",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.box.just = "right",
        legend.text = element_text(size = 12)) +
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 1)) + 
  labs(color = "Soil residual nitrate multiplier", y = "Cereal rye biomass (kg/ha)") 

rye_sens_output_wi <- sens_output_wi %>%
  filter(str_detect(Zone, "NT")) %>%
  group_by(Zone, NO3.Multiplier) %>%
  summarise(final_rye = max(rye_bm_kg_ha)) %>%
  separate_wider_delim(Zone, delim = " ", names = c("state", "tillage", "year", "soil_scenario", "plant_scenario")) %>%
  mutate(site_year = paste0(state, year))


rye_sens_both <- rbind(rye_sens_output_ny, rye_sens_output_wi)

rye_obs <- read_csv("~/Documents/Soybean_yga/obs_rye_adj.csv")

mod_ef <- function(pred, obs){
  mse = sum((pred-obs)^2)
  mse_obs = sum((obs-mean(obs))^2)
  ef = 1-(mse/mse_obs)
}


sens_ny_no3_filter <- sens_merged_no3 %>%
  dplyr::select(Zone, Clock.Today, rye_bm_kg_ha, soy_kg_ha, Date, NO3.Multiplier, SimulationName)
  
sens_wi_no3_filter <- sens_merged_wi %>%
  dplyr::select(Zone, Clock.Today, rye_bm_kg_ha, soy_kg_ha, Date, NO3.Multiplier, SimulationName)

  
sens_no3_both <- rbind(sens_ny_no3_filter, sens_wi_no3_filter)


sens_no3_plot_all <-  sens_no3_both %>%
  filter(str_detect(Zone, "NT"),
         rye_bm_kg_ha > 0) %>%
  mutate(site_year = SimulationName) %>%
  ggplot(aes(x = Date,
             y = rye_bm_kg_ha,
             colour = as.character(NO3.Multiplier)))+
  geom_line() +
  scale_color_manual(values = okabe_ito) +
  facet_wrap(~site_year, scales = "free_x") +
  theme_classic() +
  theme(legend.position = "top",
        legend.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45),
        legend.box = "vertical",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.box.just = "right",
        legend.text = element_text(size = 12)) +
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 1)) + 
  labs(color = "Soil residual nitrate multiplier", y = "Cereal rye biomass (kg/ha)") 




rye_comparison <- left_join(rye_sens_both, rye_obs, by = "site_year") %>%
  group_by(NO3.Multiplier) %>%
  mutate(rmse_rye = Metrics::rmse(rye_bm_kg_ha, final_rye),
         rrmse_rye = rmse_rye/mean(rye_bm_kg_ha),
         ef_rye = mod_ef(pred = final_rye, obs = rye_bm_kg_ha),)


library(ggpmisc)

rye_comparison %>%
  ggplot(aes(x = rye_bm_kg_ha, y = final_rye)) +
  geom_point(aes(colour = state), size = 2.5) +
  facet_wrap(~NO3.Multiplier) + 
  theme_classic() +
  geom_abline(slope = 1) +
  scale_color_manual(values = c("red","black")) +
  coord_cartesian(ylim = c(0, 13000), xlim = c(0,13000)) +
  geom_text(aes(label = paste0("RMSE = ", round(rmse_rye,2)), x = 1000, y = 12000), size = 4,inherit.aes = F) +
  geom_text(aes(label = paste0("RRMSE = ", round(rrmse_rye*100,2),"%"), x = 1000, y = 11200), size = 4,inherit.aes = F) +
  geom_text(aes(label = paste0("EF = ", round(ef_rye,2)), x = 1000, y = 10500), size = 4, inherit.aes = F) +
  theme(legend.position = "top",
        legend.title = element_text(face = "bold"),
        legend.box = "vertical",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.box.just = "right",
        legend.text = element_text(size = 12)) +
  labs(x = "Observed cereal rye biomass (kg/ha)", y = "Predicted cereal rye biomass (kg/ha)", color = "") 

