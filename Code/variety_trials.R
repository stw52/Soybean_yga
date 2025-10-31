#### Aurora Variety Trials

library(tidyverse)
library(ggpubr)
library(ggpmisc)


var_trials <- read_csv("R_data/Soy_rye_data/aurora_variety_trials_2007_2013.csv")

trials_clean <- var_trials %>%
  mutate(sow_date = as.Date(sow_date, format = "%m/%d/%y"),
         stage_date = as.Date(stage_date, format = "%m/%d/%y"),
         harvest_date = as.Date(harvest_date, format = "%m/%d/%y")) %>%
  group_by(year, mg) %>%
  summarise(yield_avg = mean(yield_bu_ac),
            yield_std = sd(yield_bu_ac),
            sow_date = unique(sow_date),
            harvest_date = unique(harvest_date),
            stage_date = unique(stage_date),
            stage = unique(stage))
  
trials_clean %>%
  ggplot(aes(x = year, y = yield_avg, colour = mg)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  stat_poly_eq() +
  geom_errorbar(aes(ymin = yield_avg-yield_std, ymax = yield_avg+yield_std), width = 0.2) +
  theme_pubr() +
  labs(y = "Yield average (Bu/ac)", x = "Year", colour = "Maturity Group")
