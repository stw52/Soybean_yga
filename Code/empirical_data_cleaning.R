library(tidyverse)
library(ggpubr)
library(ggpmisc)
library(readxl)


setwd("~/")
all_soy_rye <- read_csv("Documents/NY_WI_NTvT_master_extra.csv")

ny_soy_rye <- all_soy_rye %>%
  filter(State == "NY")

wi_soy_rye <- all_soy_rye %>%
  filter(State == "WI")

soy_planting_dates <- all_soy_rye %>%
  select(SoySeedingDate, ExperimentName, Year, State, Location, Till_NT) %>%
  distinct(SoySeedingDate, ExperimentName, Year, State, Location, Till_NT)

all_yg <- all_soy_rye %>%
  select(Year,Till_NT, Yield_Mgha, Yield_buac, `Precip15May-15Jun_in`, ExperimentName, Rep, Location) %>%
  group_by(Year, ExperimentName, Rep) %>%
  mutate(yield_gap = mean(Yield_Mgha[Till_NT == "NT"]/Yield_Mgha[Till_NT == "Till"], na.rm = TRUE)) %>%
  distinct(Year, ExperimentName, Rep, yield_gap, `Precip15May-15Jun_in`, Location)

lm(all_yg$yield_gap~all_yg$`Precip15May-15Jun_in`)

#### Plots for Yield gap based on precipitation around planting
all_yg %>%
 ggplot(aes(x = `Precip15May-15Jun_in`, y = yield_gap)) +
 geom_point() +
 geom_smooth(method = "lm", color = "red") +
 theme_bw() +
 labs(x = "Precipitation May 15th to June 15th (in)", y = "Yield Gap between NT and Till", title = "No-Till vs Till Yield Gap related to Accumulated Precip. at Planting") +
 stat_poly_eq(use_label(c("eq","R2")))



ny_yg <- all_yg %>%
  filter(Location == "Musgrave")

ny_yg %>%
  ggplot(aes(x = `Precip15May-15Jun_in`, y = yield_gap, color = ExperimentName)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red") +
  theme_bw() +
  labs(x = "Precipitation May 15th to June 15th (in)", y = "Yield Gap between NT and Till", title = "No-Till vs Till Yield Gap related to Accumulated Precip. at Planting") +
  stat_poly_eq(use_label(c("eq","R2")))


ny_yg_distinct <- ny_yg %>%
  distinct(Year, yield_gap)


ny_yg_point <- ny_soy_rye %>% 
  select(Year,Till_NT, Yield_Mgha, Yield_buac, `Precip15May-15Jun_in`, ExperimentName) %>%
  group_by(ExperimentName) %>%
  mutate(yield_gap = mean(Yield_buac[Till_NT == "NT"]/Yield_buac[Till_NT == "Till"], na.rm = TRUE)) %>%
  ungroup()
  


#ny_yg_point %>%
#  ggplot(aes(x = `Precip15May-15Jun_in`, y = yield_gap)) +
 ## geom_point() +
#  theme_bw() +
 # labs(x = "Precipitation Accum. May 15th to June 15th (in)", y = "Soybean Yield Gap Ratio NT:Till") +
#  geom_smooth(method = "lm", formula = "y~x") +
 # stat_poly_eq(use_label(c("eq","R2")))

#ny_yg %>%
 # ggplot(aes(x = `Precip15May-15Jun_in`, y = yield_gap, label = Year)) +
  #geom_point() +
  #geom_smooth(method = "lm", color = "red") +
  #theme_bw() +
  #labs(x = "Precipitation May 15th to June 15th (in)", y = "Yield Gap between NT and Till", title = "New York No-Till vs Till Yield Gap related to Accumulated Precip. at Planting") +
  #stat_poly_eq(use_label(c("eq","R2"))) +
  #scale_color_viridis_c() +
  #geom_label()

wi_yg <- wi_soy_rye %>%
  select(Year,Till_NT, Yield_Mgha, Yield_buac, `Precip15May-15Jun_in`) %>%
  group_by(Year) %>%
  mutate(yield_gap = mean(Yield_Mgha[Till_NT == "NT"]/Yield_Mgha[Till_NT == "Till"], na.rm = TRUE))

#wi_yg %>%
#  ggplot(aes(x = `Precip15May-15Jun_in`, y = yield_gap)) +
 # geom_point() +
#  geom_smooth(method = "lm", color = "red") +
#  theme_bw() +
 # labs(x = "Precipitation May 15th to June 15th (in)", y = "Yield Gap between NT and Till", title = "Wisconsin No-Till vs Till Yield Gap related to Accumulated Precip. at Planting") +
 # stat_poly_eq(use_label(c("eq","R2")))


####Plots for Yield gap around emergence percentage

wi_emerg <- wi_soy_rye %>%
  select(Year, Per_Emergence, Yield_Mgha, Till_NT) %>%
  group_by(Year) %>%
  mutate(yield_gap = mean(Yield_Mgha[Till_NT == "NT"]/Yield_Mgha[Till_NT == "Till"], na.rm = TRUE)) %>%
  ungroup()

wi_emerg %>%
  ggplot(aes(x = Per_Emergence, y = yield_gap)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  stat_poly_eq(use_label(c("eq","R2")))


mean(wi_emerg$Per_Emergence[wi_emerg$Till_NT == "Till"], na.rm = TRUE)

sum(wi_emerg$Till_NT == "NT")

emerge_nt_wi <- wi_emerg$Per_Emergence[wi_emerg$Till_NT == "NT"]
emerge_till_wi <- wi_emerg$Per_Emergence[wi_emerg$Till_NT == "Till"]
t.test(emerge_nt_wi, emerge_till_wi, alternative = "less")


ny_emerg <- ny_soy_rye %>%
  select(Year, Per_Emergence, Yield_Mgha, Till_NT) %>%
  group_by(Year) %>%
  mutate(yield_gap = mean(Yield_Mgha[Till_NT == "NT"]/Yield_Mgha[Till_NT == "Till"], na.rm = TRUE)) %>%
  ungroup()

emerge_nt_ny <- ny_emerg$Per_Emergence[wi_emerg$Till_NT == "NT"]
emerge_till_ny <- ny_emerg$Per_Emergence[wi_emerg$Till_NT == "Till"]
#t.test(emerge_nt_ny, emerge_till_ny)
#emerge_till_ny


mean(emerge_nt_ny, na.rm = TRUE)

nt_yield <- all_soy_rye$Yield_buac[all_soy_rye$Till_NT == "NT"]
t_yield <- all_soy_rye$Yield_buac[all_soy_rye$Till_NT == "Till"]

t.test(nt_yield, t_yield, alternative = "less")



#### Read in new soybean and rye data from NY

new_projects <- read_xlsx(path = "R_data/Soy_rye_data/Additional Soybean Projects Spreadsheet.xlsx")

new_projects_keep <- new_projects %>%
  filter(ExperimentName %in% "Mowtivation") %>%
  mutate(SoySeedingDate = as.character(SoySeedingDate))


new_and_old_soy <- rbind(all_soy_rye, new_projects_keep)

write_csv(new_and_old_soy, "R_data/Soy_rye_data/all_soybean_data_combined.csv")
