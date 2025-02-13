library(tidyverse)
library(emmeans)
library(lme4)
library(modelr)
library(metafor)
library(ggpmisc)
library(lmerTest)
library(glmmTMB)
library(DHARMa)
library(shades)
##### Load in master spreadsheet for soy and rye data and the weather data

soy_rye_raw <- read_csv(here("R_data/Soy_rye_data/NY_WI_NTvT_master.csv"))

soy_rye_raw %>%
  filter(Till_NT == "NT") %>%
  filter(RyeTerminationMethod != "Roller Crimped") %>%
  select(ExperimentName, Year, RyeTerminationMethod) #### 2022 includes mowing and rolling as methods of killing cereal rye

soy_rye_std <- soy_rye_raw %>%
  filter(Year != 2022) %>%
  filter(ExperimentName != "OCS")

planting_prcp <- read_csv(here("R_data/Weather_data/PRISM/plant_prcp_all.csv"))


soy_rye_useful <- soy_rye_raw %>%
  select(Year, State, Location, ExperimentName, Rep, SoySeedingDate, Till_NT,
         RyeBiomass_kgha, SoyStand_ha, Per_Emergence, WeedBiomass_kgha, Yield_Mgha, Yield_buac, ExperimentalTrt) %>%
  filter(Year != 2022)

nt_t_counts <- soy_rye_useful %>%
  group_by(ExperimentName, Year) %>%
  mutate(nt_count = sum(Till_NT == "NT")) %>%
  mutate(till_count = sum(Till_NT == "Till")) %>%
  distinct(ExperimentName, till_count, nt_count)
  


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

#### Linear model to test overall difference in yield by tillage
lm_tillage_all <- lm(Yield_buac ~ Till_NT, data = comb_soy_prcp)
hist(residuals(lm_tillage_all))
plot(predict(lm_tillage_all) ~ residuals(lm_tillage_all))

emmeans(lm_tillage_all, pairwise~Till_NT) #### Indicates a significant difference between overall till and no-till yield values

lm_yield_mgha <- lm(Yield_Mgha ~ Till_NT, data = comb_soy_prcp_revised)

emmeans(lm_yield_mgha, pairwise ~ Till_NT)


#### Linear model to test difference between state first, then location

lm_state <- lm(Yield_buac ~ State, data = comb_soy_prcp)
emmeans(lm_state, pairwise~State) # Indicates a significant difference across states in Yield

lm_location <- lm(Yield_buac ~ Location, data = comb_soy_prcp)
emmeans(lm_location, pairwise~Location)

lm_loc_tillage <- lm(Yield_buac ~ Till_NT * Location, data = comb_soy_prcp) 
emmeans(lm_loc_tillage, pairwise~Till_NT|Location) # Indicates significant differences when we look across station location and tillage type in yield effect

emmip(lm_loc_tillage, Location ~ Till_NT , CI = TRUE) + # Plot shows yield difference based on tillage for each location
  theme_bw() +
  labs(y = expression(paste("Predicted Yield, " , "(Bu ", " ", ac^{-1}, ")")), x= "Tillage Type")


#### Plot observational data based on above anova
t_test_soy_yield <- comb_soy_prcp_revised %>%
  group_by(Year, Location) %>%
  summarize(t_test = list(t.test(Yield_Mgha ~ Till_NT))) %>%
  mutate(
    t_statistic = sapply(t_test, function(x) x$statistic),
    p_value = sapply(t_test, function(x) x$p.value),
    mean_tillage_1 = sapply(t_test, function(x) x$estimate[1]),
    mean_tillage_2 = sapply(t_test, function(x) x$estimate[2])
  )

comb_soy_prcp_revised <- comb_soy_prcp %>%
  filter(Location == "Musgrave" | Location == "AARS") %>% 
  mutate(site_year = trimws(paste0(State, "-", Year))) %>%
  filter(site_year != "WI-2014" & site_year != "WI-2022")

comb_soy_prcp_revised %>%
  group_by(Till_NT) %>%
  summarise(range_yield = range(Yield_Mgha),
            mean_yield = mean(Yield_Mgha))




soy_prcp_w_pvals <- left_join(comb_soy_prcp_revised, t_test_soy_yield)

soy_prcp_w_pvals <- soy_prcp_w_pvals %>%
  mutate(label = case_when(
      p_value < 0.05 & Till_NT == "Till"~ "*",
      p_value > 0.05 & Till_NT == "Till"~ "" ,
      !is.na(p_value) & Till_NT == "Till" ~ "*" 
    )
  ) %>%
  group_by(Year, Location, ExperimentName, Till_NT) %>%
  mutate(mean_yield = mean(Yield_Mgha)) %>%
  mutate(se_yield = sd(Yield_Mgha)/sqrt(n())) %>%
  ungroup() %>%
  distinct(Year, Location, ExperimentName, mean_yield, se_yield, p_value, label, Till_NT, site_year) %>%
  filter(ExperimentName != "Mowtivation")

soy_prcp_w_pvals$Year[soy_prcp_w_pvals$Year == 2021 & soy_prcp_w_pvals$ExperimentName == "Seeding Depth Study"] <- "2021a" 
soy_prcp_w_pvals$Year[soy_prcp_w_pvals$Year == 2021 & soy_prcp_w_pvals$ExperimentName == "Starter Fertilizer Study"] <- "2021b" 



comb_soy_prcp_revised %>%
  ggplot(aes(x = Year, y = Yield_Mgha, fill = Till_NT)) +
  stat_summary(geom = "col",fun = "mean", width = 0.6,  position = position_dodge(width = 0.6)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se",width = 0.2, position = position_dodge(width = 0.6)) +
  facet_wrap(~Location, nrow = 3, scales = "free_x", labeller = as_labeller(location_labels)) +
  scale_fill_manual(values = c("#3B528BFF", "#FDE779")) +
  theme_bw() +
  labs(y = expression(paste("Mean Yield, " , "(Mg ", " ", ha^{-1}, ")")),
       fill = "Tillage Type", 
       x = 'Site Year') +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) 


location_labels <- c("AARS" = "Arlington, WI", "Musgrave" = "Aurora, NY")

library(viridis)

soy_pval_graph <- soy_prcp_w_pvals %>%
  ggplot(aes(x = Year, y = mean_yield, fill = Till_NT)) +
  geom_col(width = 0.6,  position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = mean_yield-se_yield, ymax = mean_yield+se_yield), width = 0.2, position = position_dodge(width = 0.6)) +
  geom_text(aes(label = label, colour = "black"), size = 12, nudge_y = 0.02, show.legend = FALSE, color = "black") +
  facet_wrap(~Location, nrow = 3, scales = "free_x", labeller = as_labeller(location_labels)) +
  scale_fill_manual(values = c("#3B528BFF", "#FDE779"), 
                    labels = c("NT" = "No-till", "Till" = "Till")) +
  theme_bw() +
  labs(y = expression(paste("Mean yield, " , "(tons ", " ", ha^{-1}, ")")),
       fill = "", 
       x = 'Year') +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, vjust = 0.7),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 16)) +
  ylim(c(0,4.5))


ggsave(soy_pval_graph, filename = here("Figures/yield_by_tillage_pvals_metric.pdf"), dpi = 600, useDingbats = FALSE, 
       width = 8, height = 8, units = "in")

soy_pval_graph_english <- soy_prcp_w_pvals %>%
  mutate(yield_bu_ac = mean_yield*14.8696,
         se_yield_eng = se_yield*14.8696) %>%
  ggplot(aes(x = Year, y = yield_bu_ac, fill = Till_NT)) +
  geom_col(width = 0.6,  position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = yield_bu_ac-se_yield_eng, ymax = yield_bu_ac+se_yield_eng), width = 0.2, position = position_dodge(width = 0.6)) +
  geom_text(aes(label = label, colour = "black"), size = 12, nudge_y = 0.02, show.legend = FALSE, color = "black") +
  facet_wrap(~Location, nrow = 3, scales = "free_x", labeller = as_labeller(location_labels)) +
  scale_fill_manual(values = c("#3B528BFF", "#FDE779"), 
                    labels = c("NT" = "No-till", "Till" = "Till")) +
  theme_bw() +
  labs(y = expression(paste("Mean yield, " , "(Bu ", " ", ac^{-1}, ")")),
       fill = "", 
       x = 'Year') +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
  ylim(c(0,65))
ggsave(soy_pval_graph_english, filename = here("Figures/yield_by_tillage_pvals_english.pdf"), dpi = 600, useDingbats = FALSE, 
       width = 11, height = 6, units = "in")


comb_soy_prcp_revised %>%
  mutate(prcp_cm = plant_prcp*2.54) %>%
  ggplot(aes(x = Year, y = prcp_cm)) +
  stat_summary(geom = "col", fun = "mean") +
  facet_wrap(~Location, nrow = 1, scales = "free_x", labeller = as_labeller(location_labels)) +
  scale_fill_manual(values = c("#3B528BFF", "#FDE779")) +
  theme_bw() +
  labs(y = "Precipitation accumulation +/- 15 days of planting date (cm)",
       x = 'Year') +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) 



#### Splitting up dataset into experiments or groups based on similarities; then model for each group


# NY2013 and NY2014
barns <- comb_soy_prcp_revised %>%
  filter(ExperimentName == "BARNS") %>%
  mutate(block = Rep)

barns_mod <- lmer(log(Yield_Mgha) ~ Till_NT*Year + (1|Till_NT:Year), data = barns)
hist(resid(barns_mod))
plot(predict(barns_mod)~resid(barns_mod))
barns_emmeans <- emmeans(barns_mod, ~Till_NT|Year, type = "response")
pairs(barns_emmeans)

barns_emmeans

barns_mod2 <- glmmTMB(Yield_Mgha ~ Till_NT*Year + (1|block:Year), 
                      data = barns,
                      family = tweedie(link = "log"),
                      ziformula = ~1)
summary(barns_mod2)
simulateResiduals(barns_mod2, plot = TRUE)
pairs(emmeans(barns_mod2, ~Till_NT|Year, type = "response"))

barns_mod3 <- lmer(log(Yield_Mgha) ~ Till_NT*Year + (1|Year:block), data = barns) #Using this model for the table below as of now
summary(barns_mod3)
hist(resid(barns_mod3))
plot(resid(barns_mod3)~predict(barns_mod3))
pairs(emmeans(barns_mod3, ~ Till_NT|Year))

anova(barns_mod3, ddf = "Kenward-Roger")

barns_model_data <- data.frame(
  count_nt = c(4, 4),
  count_till = c(4, 2),
  ExperimentName = c("BARNS", "BARNS"),
  state = c("NY","NY"),
  Year = c(2013, 2014),
  plant_prcp = c(18.37, 6.78),
  SE_nt = c(0.145,0.133),
  SE_till = c(0.150,0.255),
  mean_yield_gap = c(-0.0339, -0.2161),
  nt_yield = c(2.89,2.65),
  till_yield = c(2.99,3.29),
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
  plant_prcp = c(16.89, 3.45),
  SE_nt = c(0.373,0.167),
  SE_till = c(0.385, 0.323),
  mean_yield_gap = c(-0.0317, -0.6562),
  nt_yield = c(2.63, 1.18),
  till_yield = c(2.71, 2.28),
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
  plant_prcp = c(8.17, 14.76, 8.2, 12.12),
  SE_nt = c(0.197, 0.459, 0.404, 0.512),
  SE_till = c(0.220, 0.503, 0.485, 0.560),
  mean_yield_gap = c(-0.0835, -0.0397,-0.1341 ,-0.0336),
  nt_yield = c(2.88, 3.38, 2.93, 3.74),
  till_yield = c(3.13, 3.65, 3.35, 3.87),
  SE_TE = c(0.0687, 0.0397, 0.0648, 0.0627)
)

# NY2018
ocs <- comb_soy_prcp %>%
  filter(ExperimentName == "OCS") %>%
  mutate(block = Rep)

ocs_mod <- lmer(log(Yield_Mgha) ~ Till_NT + (1|block), data = ocs)
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
  plant_prcp = c(4.45),
  SE_nt = c(0.113),
  SE_till = c(0.213),
  mean_yield_gap = c(-0.631),
  nt_yield = c(1.97),
  till_yield = c(3.7),
  SE_TE = c(0.103)
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
  plant_prcp = c(13.18),
  SE_nt = c(0.0671),
  SE_till = c(0.1119),
  mean_yield_gap = c(-0.126),
  nt_yield = c(3.67),
  till_yield = c(4.16),
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
  plant_prcp = c(7.44),
  SE_nt = c(0.125),
  SE_till = c(0.370),
  mean_yield_gap = c(-0.555),
  nt_yield = c(2.31),
  till_yield = c(4.01),
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
  plant_prcp = c(4.10),
  SE_nt = c(0.0927),
  SE_till = c(0.3459),
  mean_yield_gap = c(-0.513),
  nt_yield = c(2.38),
  till_yield = c(3.97),
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
  plant_prcp = c(1.12),
  SE_nt = c(0.062),
  SE_till = c(0.183),
  mean_yield_gap = c(-0.39),
  nt_yield = c(2.25),
  till_yield = c(3.32),
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
  mutate(sd_nt = SE_nt*sqrt(count_nt)) %>%
  mutate(sd_till = SE_till*sqrt(count_till)) %>%
  mutate(site_year = paste0(state, "-", Year)) %>%
  mutate(var_te = SE_TE^2)

library(meta)
library(metafor)

meta_analysis_yg <- metacont(data = meta_final_data, n.e = count_nt, n.c = count_till,
                             mean.e = nt_yield, mean.c = till_yield, sm = "SMD", sd.e = sd_nt, 
                             sd.c = sd_till, studlab = site_year)


summary(meta_analysis_yg)

forest(meta_analysis_yg, leftcols = "site_year", sortvar = TE)
funnel(meta_analysis_yg, studlab = TRUE)


meta_results <- data.frame(
  effect_size = meta_analysis_yg$TE,  # Treatment effect
  ci_lower = meta_analysis_yg$lower,  # Lower bound of confidence interval
  ci_upper = meta_analysis_yg$upper,  # Upper bound of confidence interval
  i2 = meta_analysis_yg$I2            # I² for heterogeneity
)

meta_result_full <- cbind(meta_final_data, meta_results)

metacont_plot <- meta_result_full %>%
  ggplot(aes(x = plant_prcp, y = effect_size))+
  geom_point(aes(color = site_year)) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper, colour = site_year), width = 0.2) + # Error bars for CI
  labs(x = "Precipitation Accumulation +/- 15 Days of Planting (cm)",
       y = "Standard Mean Difference in Yield",
       color = "Site-Year") +
  stat_poly_eq(use_label(c("Eq","R2")))+
  theme_bw() +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", alpha = 0.85)
  

ggsave(metacont_plot, filename = here("Figures/metacontplot.pdf"), dpi = 600, useDingbats = FALSE)

# Meta analysis different way
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

initial_mreg <- metareg(metagen_analysis, ~plant_prcp)
bubble(initial_mreg, studlab = TRUE)

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

ratio_yg_meta <- metagen_result_final %>%
  ggplot(aes(x = plant_prcp, y = yield_ratio))+
  geom_point(aes(color = state), size = 4) +
  geom_errorbar(aes(ymin = lowerci_trans, ymax = upperci_trans, colour = state), width = 0.4, linewidth = 2) + # Error bars for CI
  labs(x = "Precipitation accumulation +/- 15 days of planting (cm)",
       y = "Ratio of no-till yield to tilled yield",
       color = "") +
  stat_poly_eq(formula = y~x, use_label(c("R2", "eq")), size = 7)+
  theme_pubr() +
  scale_color_manual(values = c("black", "red")) +
  geom_smooth(method = "lm",formula = y ~ x,  se = FALSE, color = "gray", 
              linetype = "dashed", linewidth = 2)  +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))

ggsave(ratio_yg_meta, filename = here("Figures/ratio_yg_metaplot.pdf"), dpi = 600, useDingbats = FALSE,
       height = 7, width = 9, units = "in")

#### Random forest for predictor importance ; May not use for ASA but could be a cool dive into predictors for soybean yield gap
rf_soy_prcp_data<- comb_soy_prcp %>%
  group_by(Year, ExperimentName, Rep) %>%
  mutate(yield_gap = mean(Yield_buac[Till_NT == "NT"]/Yield_buac[Till_NT == "Till"], na.rm = TRUE)) %>%
  distinct(RyeBiomass_kgha, SoyStand_ha, Per_Emergence, WeedBiomass_kgha, yield_gap, plant_prcp, gdd_plant) %>%
  filter(!is.na(yield_gap)) %>%
  filter(RyeBiomass_kgha > 0)

rf_soy_prcp_data


rf_soy_prcp_data <- na.omit(rf_soy_prcp_data)

rf_soy_prcp_data <- rf_soy_prcp_data %>%
  mutate(WeedBiomass_kgha = as.numeric(WeedBiomass_kgha)) 

exclude_vars <- c("Year", "Rep", "ExperimentName")

rf_soy_prcp_subset <- rf_soy_prcp_data[, !names(rf_soy_prcp_data) %in% exclude_vars]

set.seed(222)
ind <- sample(2, nrow(rf_soy_prcp_subset), replace = TRUE, prob = c(0.7, 0.3))
train_soy_prcp <- rf_soy_prcp_subset[ind==1,]
test_soy_prcp <- rf_soy_prcp_subset[ind==2,]

library(caret)
library(modelr)
library(randomForest)
library(Boruta)
rf_soy_prcp <- randomForest(yield_gap ~.,
                     data = train_soy_prcp,
                     proximity = TRUE
)
varImpPlot(rf_soy_prcp)
rf_soy_prcp

plot(rf_soy_prcp)

feat_importance <- Boruta(yield_gap ~., data = rf_soy_prcp_subset, doTrace = 1)

feat_importance$ImpHistory

feat_stats <- attStats(feat_importance)
feat_stats$Feature <- rownames(feat_stats)

plot(feat_importance)
library(gghighlight)

feat_stats



feat_stats %>%
  mutate(Feature = c("Rye Biomass (kg/ha)", "Soy Population", "Percent Emergence", "Weed Biomass (kg/ha)", "Precip. at Planting","GDD at Planting")) %>%
  ggplot(aes(x = reorder(Feature, medianImp), y = medianImp)) +
  geom_col(fill = "darkgreen") +
  coord_flip() +
  gghighlight(meanImp > 20) +
  labs(x = "Predictor Variable for Soybean Yield Gap", y= "Mean Variable Importance") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))


#### Some PET Graphs from APSIM maybe?

apsim_musgrave_pet_nt <- read_xlsx("R_data/Soil_data/APSIM_No_till_2016_musgrave_waterbalance.xlsx")


mus_pet_nt <- apsim_musgrave_pet_nt %>%
  filter(Clock.Today > "2016-02-29" & Clock.Today < "2016-06-01") %>% 
  mutate(cum_total_PET = cumsum(ET_potential),
         cum_soil_ev = cumsum(soil_evap_nt),
         cum_rye_et = cumsum(cc_trans),
         cum_soy_et = cumsum(Soy_ET)) %>%
  pivot_longer(names_to = "measure",
               values_to = "values",
               cols = c(cum_total_PET, cum_soil_ev, cum_rye_et, cum_soy_et))
  
mus_pet_nt %>%
  ggplot(aes(x = Clock.Today, values, color = measure)) +
  geom_line() +
  theme_bw() +
  scale_color_discrete(labels = c("cum_rye_et" = "Cereal Rye Transpiration",
                                  "cum_soil_ev" = "Soil Evaporation",
                                  "cum_soy_et" = "Soybean Transpiration",
                                  "cum_total_PET" = "Cummulative Total PET")) +
  labs(y = "Cummulative Potential Evapotranspiration (mm)", 
       x = "Month",
       color = "PET Measure")


mus_pet_plant <- left_join(mus_pet_nt, mus_pet_nt)

apsim_musgrave_pet_till <- read_xlsx("R_data/Soil_data/till_water_balance_musg_2016_apsim.xlsx")


mus_pet_till <- apsim_musgrave_pet_till %>%
  filter(Clock.Today < "2016-06-01") %>%
  rename() %>%
  mutate(cum_total_PET = cumsum(ET_potential_till),
         cum_soil_ev = cumsum(soil_evap_till),
         cum_soy_et = cumsum(soy_et_till)) %>%
  pivot_longer(names_to = "measure",
               values_to = "values",
               cols = c(cum_total_PET, cum_soil_ev, cum_soy_et))

mus_pet_till %>%
  ggplot(aes(x = Clock.Today, values, color = measure)) +
  geom_line() +
  theme_bw() +
  scale_color_discrete(labels = c("cum_soil_ev" = "Soil Evaporation",
                                  "cum_soy_et" = "Soybean Transpiration",
                                  "cum_total_PET" = "Cummulative Total PET")) +
  labs(y = "Cummulative Potential Evapotranspiration (mm)", 
       x = "Month",
       color = "PET Measure")


mus_water_def_2016_nt <- read_csv("~/Downloads/water-deficit-for-2016.csv")
mus_water_def_2016_nt <- mus_water_def_2016_nt %>%
  mutate(DateTime = as.Date(DateTime)) %>%
  mutate(Deficit_metric = Deficit*0.833333) %>%
  mutate(tillage = rep("No-Till"))

mus_water_def_2016_nt %>%
  ggplot(aes(x = DateTime, y = Deficit_metric)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = as.Date("2016-05-20"), color = "red") +
  theme_bw() +
  labs(x = "Month",
       y = "Soil water deficit (cm of water per m of soil)")

mus_water_def_2016_till <- read_csv("~/Downloads/water-deficit-for-2016_till_mus.csv")

mus_water_def_2016_till <- mus_water_def_2016_till %>%
  mutate(DateTime = as.Date(DateTime)) %>%
  mutate(Deficit_metric = Deficit*0.833333) %>%
  mutate(tillage = rep("Till"))

mus_water_def_2016_till %>%
  ggplot(aes(x = DateTime, y = Deficit_metric)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = as.Date("2016-05-20"), color = "red") +
  theme_bw() +
  labs(x = "Month",
       y = "Soil water deficit (cm of water per m of soil)")

total_water_def_2016 <- rbind(mus_water_def_2016_nt, mus_water_def_2016_till)


cer_rye_2016_evap <- total_water_def_2016 %>%
  filter(DateTime < "2016-06-02") %>%
  ggplot(aes(x = DateTime, y = Deficit_metric)) +
  geom_point(aes(colour = tillage), size = 2) +
  geom_line(aes(colour = tillage), linewidth = 1.5) +
  scale_color_manual(values = c("#3B528BFF", "#FDE779")) +
 # geom_vline(xintercept = as.Date("2016-05-20"), color = "red") +
  theme_bw() +
  labs(x = "Month",
       y = "Soil water deficit (mm of water per cm of soil)",
       color = "Tillage") +
  theme(axis.title = element_text(size = 18),
                                 axis.text = element_text(size = 16),
                                 legend.text = element_text(size = 16),
                                 legend.title = element_text(size = 16))

ggsave(cer_rye_2016_evap, filename = "Figures/spring_evap_2016_mus.pdf", dpi = 600,
       useDingbats = F)

# Model decision diagram
library(DiagrammeR)

grViz("digraph {
  graph [layout = dot, rankdir = TB]
  
  node [shape = rectangle]        
  rec1 [label = 'Filter data to standard scenarios']
  # edge definitions with the node IDs
  rec1
  }",
      height = 400, width = 650)

grViz("digraph {
  graph [layout = dot, rankdir = TB]
  
  node [shape = rectangle]        
  rec1 [label = 'Filter data to standard scenarios']
  rec2 [label = 'Create linear mixed effects models similar experiments']
  
  # edge definitions with the node IDs
  rec1 -> rec2
  }",
      height = 400, width = 650)

grViz("digraph {
  graph [layout = dot, rankdir = TB]
  
  node [shape = rectangle]        
  rec1 [label = 'Filter data to standard scenarios']
  rec2 [label = 'Create linear mixed effects models similar experiments']
  rec3 [label =  'Extract treatment effects and standard errors from each model']

  # edge definitions with the node IDs
  rec1 -> rec2 -> rec3 
  }",
      height = 400, width = 650)


grViz("digraph {
  graph [layout = dot, rankdir = TB]
  
  node [shape = rectangle]        
  rec1 [label = 'Filter data to standard till and no-till scenarios']
  rec2 [label = 'Collect gridded precipitation data from PRISM']
  rec3 [label = 'Create linear mixed effects models for similar experiments']
  rec4 [label =  'Extract treatment effects and standard errors from each model']
  rec5 [label = 'Utilize meta-analysis techniques to explore yield gap relation to planting precipitation']
  
  # edge definitions with the node IDs
  rec1 -> rec2 -> rec3 -> rec4 -> rec5
  }",
      height = 400, width = 650)


# Experiments timeline 

library(timevis)
data <- data.frame(
  id = 1:3,
  content = c("Event 1", "Event 2", "Event 3"),
  start = c("2023-01-01", "2023-02-15", "2023-03-30"),
  end = c(NA, "2023-03-01", NA)
)
timevis(data, options = list(zoomable = F))


