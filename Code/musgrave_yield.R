library(tidyverse)
library(here)


musgrave_soy <- here("R_data", "Soy_rye_data", "NY_WI_combined_STW_for_R.csv")

all_soy <- read_csv(musgrave_soy)

mug_soy <- mug_soy %>%
  filter(Location == "Musgrave")

mug_soy_yields <- mug_soy %>%
  filter(Till_NT == "NT") %>%
  group_by(Year) %>%
  summarise(annual_yield_bu_ac = max(Yield_buac))

emmerg_data <- all_soy %>%
  select(Till_NT, Per_Emergence, Location) %>%
  group_by(Till_NT, Location) %>%
  summarise(avg_emmerg = mean(!is.na(Per_Emergence)))
  
 


plantPop <- mug_soy %>%
  mutate(SoySeedingRate_ha = as.numeric(SoySeedingRate_ha)) %>%
  summarise(plant_pop = mean(SoySeedingRate_ha))

class(mug_soy$SoySeedingDate)

apsim_mug_input <- mug_soy %>%
  filter(Till_NT == "NT") %>%
  select(Year, Yield_buac) %>% 
  group_by(Year) %>%
  mutate(SimulationName = "Simulations") %>%
  mutate(Clock.Today = paste0(Year,"-10-17")) %>%
  mutate(Clock.Today = as.character(Clock.Today)) %>%
  mutate(soy_bu_ac = mean(Yield_buac)) %>%
  distinct(soy_bu_ac, .keep_all = TRUE) %>%
  ungroup(Year) %>%
  select(SimulationName, Clock.Today, soy_bu_ac)

write_csv(apsim_mug_input, "R_data/Soy_rye_data/apsim_musgrave.csv")


apsim_sims <- read_csv("R_data/Soy_rye_data/apsim_saved_data.csv")

apsim_sims$Clock.Today <- as.Date(apsim_sims$Clock.Today, format = "%m/%d/%y") 

apsim_sims %>%
  ggplot(aes(x = Clock.Today, y = soy_bu_ac, colour = SimulationName)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  labs(colour = "Simulation Name", y = "Soybean Yield (bu/ac)", x = "Year")
