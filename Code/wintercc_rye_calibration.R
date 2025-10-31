#### Using 2014-2015 data from winter cover crop database paper to help characterize cereal rye

library(tidyverse)
library(readxl)
library(ggpubr)
library(ggpmisc)
library(apsimx)

wintercc <- read_xlsx("~/Downloads/wintercc_musg_2014_2015.xlsx")

all_means <- wintercc %>%
  group_by(variety, kill_date) %>%
  summarise(mean_kg_ha = mean(biomass_kg_ha)) %>%
  mutate(kill_date = as.Date(kill_date))

write_csv(all_means, "~/Downloads/wintercc_all_vars.csv")

all_means %>%
  ggplot(aes(x = kill_date, y = mean_kg_ha, colour = variety)) +
  geom_point(size = 2) +
  geom_smooth(se = F, size = 0.5, method = "lm") +
  theme_pubr() +
  labs(color = "", y = "Cereal rye biomass (kg/ha)", x = "Date") +
  theme(legend.position = c(0.2,0.8))

aroos_means <- wintercc %>%
  filter(variety == "Aroostook") %>%
  group_by(kill_date) %>%
  summarise(mean_kg_ha = mean(biomass_kg_ha)) %>%
  mutate(date = as.Date(kill_date)) 

all_means_single <- all_means %>%
  group_by(kill_date) %>%
  summarise(mean_kg_ha = mean(mean_kg_ha))

write_csv(aroos_means, "~/Downloads/wintercc_aroostook.csv")

wintercc$zadoks <- as.numeric(gsub("Zadoks ", "", wintercc$zadoks))

mean_zadoks <- wintercc %>%
  group_by(kill_date) %>%
  summarise(zadoks = mean(zadoks)) %>%
  mutate(Date = as.Date(kill_date))

write_csv(mean_zadoks, "~/Downloads/mean_zadkos.csv")

#### Read in APSIM Simulation results

apsim_rye_2014 <- read_xlsx("~/Downloads/rye_model_experiments.xlsx", sheet = "Report")

kill_dates <- aroos_means$kill_date

apsim_rye_calib <- apsim_rye_2014 %>%
  filter(Clock.Today %in% kill_dates) %>%
  dplyr::select(Clock.Today, rye_bm_shoots)



apsim_plus_obs <- left_join(all_means_single, apsim_rye_calib, by = c("kill_date" = "Clock.Today"))

apsim_plus_obs %>%
  ggplot(aes(x = mean_kg_ha, y = rye_bm_shoots)) +
  geom_point() + 
  geom_abline(slope = 1) +
  theme_pubr() +
  coord_cartesian(xlim = c(0,7000), ylim = c(0,7000)) +
  stat_poly_line() +
  stat_poly_eq() +
  labs(x = "Observed cereal rye biomass (kg/ha)", y = "APSIM predicted cereal rye biomass (kg/ha)")
 
r_squared <- 1 - sum((apsim_plus_obs$mean_kg_ha - apsim_plus_obs$rye_bm_shoots)^2) / sum((apsim_plus_obs$mean_kg_ha - mean(apsim_plus_obs$mean_kg_ha))^2)
rmse <- sqrt(mean((apsim_plus_obs$rye_bm_shoots-apsim_plus_obs$mean_kg_ha)^2))


nse <- 1 - sum((apsim_plus_obs$rye_bm_shoots - apsim_plus_obs$mean_kg_ha)^2) / sum((apsim_plus_obs$mean_kg_ha - mean(apsim_plus_obs$mean_kg_ha))^2)


#### Running APSIM and trying Cultivar replacements in R ####

# Example
apsim_dir <- system.file("extdata", package = "apsimx")

inspect_apsimx("Wheat.apsimx", src.dir = apsim_dir, node = "Weather")

inspect_apsimx("Wheat.apsimx", src.dir = apsim_dir, 
               node = "Soil", soil.child = "Organic")

sim <- apsimx("Wheat.apsimx", src.dir = apsim_dir, value = "report")

summary(sim)

inspect_apsimx("Wheat.apsimx", src.dir = apsim_dir, node = "Crop")

inspect_apsimx("Wheat.apsimx", src.dir = apsim_dir, node = "Manager",
               parm = list("SowingRule1", NA))
inspect_apsimx_replacement("MaizeSoybean.apsimx", src.dir = apsim_dir,
                           node = "Maize", display.available = T)

inspect_apsimx_replacement("MaizeSoybean.apsimx", src.dir = apsim_dir,
                           node = "Maize", display.available = T, node.child = "Phenology",
                           node.subchild = "ThermalTime", node.subsubchild = "BaseThermalTime",
                           node.sub3child = "Response",
                           parm = "Y")

#### With rye experiment I set up
##apsimx_options(exe.path = "/Applications/APSIM2024.5.7502.0.app/Contents/MacOS/APSIM")


rye_dir <- "/Users/stw52/Downloads"

inspect_apsimx("rye_exp_197.apsimx", src.dir = rye_dir, node = "Crop")

sim_0 <- apsimx("rye_model_experiments.apsimx", src.dir = rye_dir, value = "report", simplify = F)

wheat <- get_apsimx_json(model = "Wheat")


insert_replacement_node("rye_exp_197.apsimx", src.dir = rye_dir, verbose = T,
                        rep.node = wheat)

inspect_apsimx_replacement("rye_exp_197.apsimx", src.dir = rye_dir,
                           node = "Wheat",display.available = T)

pp <- insert_replacement_node("rye_exp_197.apsimx", src.dir = rye_dir, wrt.dir = rye_dir, rep.node = wheat,
                              edit.tag = "-edited-STW",verbose = TRUE)

sim_1 <- apsimx("rye_exp_197-edited-SW.apsimx", src.dir = rye_dir, value = "report", simplify = F)
sim_1$Report$Clock.Today

pp1 <- inspect_apsimx_replacement("rye_exp_197-edited-SW.apsimx",
                                  src.dir = rye_dir,
                                  node = "Wheat",
                                  node.child = "Cultivars",
                                  node.subchild = "USA",
                                  node.subsubchild = "Yecora",
                                  display.available = TRUE,
                                  print.path = TRUE)

inspect_apsimx_replacement("rye_exp_197-edited-SW.apsimx",
                           src.dir = rye_dir,
                           node = "Wheat",
                           node.child = "Cultivars",
                           node.subchild = "USA",
                           node.subsubchild = "Yecora",
                           verbose = F)

param1 <- "Wheat.Leaf.Photosynthesis.RUE.FixedValue"
param2 <- "Wheat.Leaf.ExtinctionCoeff.VegetativePhase.FixedValue"

obs <- read_xlsx("~/Downloads/wintercc_all_vars.xlsx", sheet = "Sheet1")

sim_1_pred <- sim_1$Report
  
  
obs <- obs %>%
  mutate(Date = (Clock.Today)+dhours(12)) %>%
  dplyr::select(-Clock.Today)

rye_op <- optim_apsimx("rye_exp_197-edited-SW.apsimx", 
                    src.dir = rye_dir, 
                    parm.paths = c(param1, param2),
                    data = obs, 
                    weights = "mean",
                    replacement = c(TRUE, TRUE),
                    initial.values = c(1.2, 0.5))


#### Checking zadoks from simulation with variety observations

sim_zad <- read_xlsx("~/Downloads/rye_exp_197-edited-SW.xlsx", sheet = "Report")

kill

sims <- sim_zad %>%
  filter(Clock.Today %in% kill_dates) %>%
  mutate(variety = rep("pred"), Date = as.Date(Clock.Today)) %>%
  dplyr::select(Date, zadoks, variety)


obs_zadoks <- wintercc %>%
  group_by(variety, kill_date) %>%
  summarise(zadoks = mean(zadoks)) %>%
  mutate(Date = as.Date(kill_date)) %>%
  dplyr::select(Date, zadoks, variety)


sim_obs <- rbind(obs_zadoks, sims)

sim_obs %>%
  ggplot(aes(x = Date, y = zadoks, color = variety)) +
  geom_point() +
  geom_line() +
  theme_pubr() +
  labs(y = "Zadoks stage", x = "Date", color = "")

nse_zadok <- 1 - sum((sims$zadoks - mean_zadoks$zadoks)^2) / sum((mean_zadoks$zadoks - mean(mean_zadoks$zadoks))^2)

rmse <- sqrt(mean((sims$zadoks-mean_zadoks$zadoks)^2))

mean(mean_zadoks$zadoks)
