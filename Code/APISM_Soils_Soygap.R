library(tidyverse)
library(soilDB)
library(apsimx)
library(sf)

mus_soil_wss <- get_ssurgo_soil_profile(lonlat = c(-76.65055, 42.73376))

mus_soil_apsim <- mus_soil_wss[[1]]$soil

mus_soil_tex <- mus_soil_apsim %>%
  dplyr::select(ParticleSizeClay, ParticleSizeSilt, ParticleSizeSand, BD, Depth, LL15, DUL, SAT, KS) %>%
  rename("sand" = ParticleSizeSand, "silt" = ParticleSizeSilt, "clay" = ParticleSizeClay)


mus_key <- SDA_spatialQuery(st_sfc(st_point(c(-76.65055, 42.73376)), crs = 4236), what = "mukey")

mus_soildb <- fetchSDA("mukey = ('289444')", duplicates = T)


lima_cokey <- mus_soildb@site$cokey[mus_soildb@site$compname == "Lima"]

mus_horizons <- mus_soildb@horizons 

lima_horizons <- mus_horizons %>%
  filter(cokey == lima_cokey) %>%
  dplyr::select(hzname, hzdept_r, hzdepb_r, fragvol_r, sandtotal_r, silttotal_r, claytotal_r,
                om_r, dbthirdbar_r, ksat_r, awc_r, ph1to1h2o_h) %>%
  unite(col = "depth", c("hzdept_r", "hzdepb_r"), sep = "-") %>%
  rename("rock" = fragvol_r, "sand" = sandtotal_r, "silt" = silttotal_r, "clay" = claytotal_r,
         "om" = om_r, "bd" = dbthirdbar_r, "ksat" = ksat_r, "awc" = awc_r, "ph" = ph1to1h2o_h)
lima_text <- lima_horizons %>%
  dplyr::select(sand, silt, clay, bd)

lima_rosetta <-ROSETTA(lima_text, vars = names(lima_text), v = "3", include.sd = T)

kpa_33 <- 33*10.197
kpa_1500 <- 1500*10.197

lima_rosetta_vwc <- lima_rosetta %>%
  mutate(scale_a = 10^alpha,
         scale_n = 10^npar,
         scale_ksat = (10^ksat)*10,
         theta_33 = theta_r+(theta_s-theta_r)/(1+(scale_a*kpa_33)^scale_n)^(1-(1/scale_n)),
         theta_1500 = theta_r+(theta_s-theta_r)/(1+(scale_a*kpa_1500)^scale_n)^(1-(1/scale_n)),
         awc = theta_33-theta_1500)

lima_s <- c(41.4, 41.4, 41.8, 42.1, 38.8, 38.6)
lima_si <- c(41.2, 41.2, 38.7, 34, 47.6, 48)
lima_c <- c(17.4, 17.4, 19.5, 23.9, 13.6, 13.4)
lima_bd <- c(1.38,1.38,1.46,1.55,1.55,1.73)

lima_lab_data <- data.frame(sand = lima_s,
                            silt = lima_si, 
                            clay = lima_c,
                            bd = lima_bd)

lima_lab_ros <- ROSETTA(lima_lab_data, vars = names(lima_lab_data), v = "3", include.sd = T)
lima_rosetta_vwc_lab <- lima_lab_ros %>%
  mutate(scale_a = 10^alpha,
         scale_a_low = 10^(alpha-sd_alpha),
         scale_a_high = 10^(alpha+sd_alpha),
         scale_n = 10^npar,
         scale_n_low = 10^(npar-sd_npar),
         scale_n_high = 10^(npar+sd_npar),
         scale_ksat = (10^ksat)*10,
         scale_ksat_low = (10^(ksat-sd_ksat))*10,
         scale_ksat_low = (10^(ksat+sd_ksat))*10,
         theta_r_low = theta_r-sd_theta_r,
         theta_r_high = theta_r+sd_theta_r,
         theta_s_low = theta_s-sd_theta_s,
         theta_s_high = theta_s+sd_theta_s,
         theta_33 = theta_r+(theta_s-theta_r)/(1+(scale_a*kpa_33)^scale_n)^(1-(1/scale_n)),
         theta_33_low = theta_r_low+(theta_s_low-theta_r_low)/(1+(scale_a_low*kpa_33)^scale_n_low)^(1-(1/scale_n_low)),
         theta_33_high = theta_r_high+(theta_s_high-theta_r_high)/(1+(scale_a_high*kpa_33)^scale_n_high)^(1-(1/scale_n_high)),
         theta_1500 = theta_r+(theta_s-theta_r)/(1+(scale_a*kpa_1500)^scale_n)^(1-(1/scale_n)),
         awc = theta_33-theta_1500,
         depth = c(150,100,50,160,400,180),
         taw = awc*depth)

#### First chatgpt rec 
vg_theta <- function(h_kpa, thetas, thetar, alpha, n_par) {
  h_cm <- h_kpa * 10.197        # IMPORTANT: convert kPa -> cm
  m <- 1 - 1 / n_par
  Se <- (1 + (alpha * h_cm)^n_par)^(-m)
  thetar + (thetas - thetar) * Se
}

summarize <- function(x) {
  c(mean = mean(x), low95 = quantile(x, 0.025),
    high95 = quantile(x, 0.975), min = min(x), max = max(x))
}

mc_soil_layer <- function(mu, sd, N=20000, seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  # thetas, thetar are expected as linear (mu_thetas, mu_thetar)
  thetas <- rnorm(N, mu$thetas, sd$thetas)
  thetar <- rnorm(N, mu$thetar, sd$thetar)
  
  # alpha: accept either mu_log10_alpha or mu_alpha
  if ("log10_alpha" %in% names(mu)) {
    log10_alpha <- rnorm(N, mu$log10_alpha, sd$log10_alpha)
    alpha <- 10^log10_alpha
  } else {
    alpha <- rnorm(N, mu$alpha, sd$alpha)
  }
  
  # n: accept either log10(n) or linear n
  if ("log10_n" %in% names(mu)) {
    log10_n <- rnorm(N, mu$log10_n, sd$log10_n)
    n_par <- 10^log10_n
  } else {
    n_par <- rnorm(N, mu$n, sd$n)
  }
  
  # ks: accept either log10_ks or ks (positive)
  if ("log10_ks" %in% names(mu)) {
    log10_ks <- rnorm(N, mu$log10_ks, sd$log10_ks)
    ks <- 10^log10_ks
  } else {
    ks <- rnorm(N, mu$ks, sd$ks)
    ks <- pmax(ks, 1e-12)
  }
  
  # physical constraints
  thetas <- pmin(pmax(thetas, 0), 1)
  thetar <- pmin(pmax(thetar, 0), thetas - 1e-8)
  n_par  <- pmax(n_par, 1.01)
  ks     <- pmax(ks, 1e-12)
  
  theta33   <- vg_theta(33, thetas, thetar, alpha, n_par)
  theta1500 <- vg_theta(1500, thetas, thetar, alpha, n_par)
  paw       <- theta33 - theta1500
  
  list(thetas = summarize(thetas),
       theta33 = summarize(theta33),
       theta1500 = summarize(theta1500),
       ks = summarize(ks),
       paw = summarize(paw))
}

# wrapper over a dataframe with rows = layers (ordered top->bottom)
mc_soil_profile <- function(df, N=20000, seed=123) {
  results <- vector("list", nrow(df))
  for (i in seq_len(nrow(df))) {
    mu <- as.list(df[i, grep("^mu_", names(df))])
    names(mu) <- sub("^mu_", "", names(mu))
    sd <- as.list(df[i, grep("^sd_", names(df))])
    names(sd) <- sub("^sd_", "", names(sd))
    results[[i]] <- mc_soil_layer(mu, sd, N, seed+i)
  }
  results
}

soil_df <- data.frame(
  mu_thetas = lima_lab_ros$theta_s,
  mu_thetar = lima_lab_ros$theta_r,
  mu_log10_alpha = lima_lab_ros$alpha,
  mu_log10_n = lima_lab_ros$npar,
  mu_log10_ks = lima_lab_ros$ksat,
  sd_thetas = lima_lab_ros$sd_theta_s,
  sd_thetar = lima_lab_ros$sd_theta_r,
  sd_log10_alpha = lima_lab_ros$sd_alpha,
  sd_log10_n = lima_lab_ros$sd_npar,
  sd_log10_ks = lima_lab_ros$sd_ksat
)

out_nonmc <- mc_soil_profile(soil_df)

# safe vg_theta (kPa -> cm)
vg_theta <- function(h_kpa, thetas, thetar, alpha, n_par) {
  h_cm <- h_kpa * 10.197
  m <- 1 - 1 / n_par
  Se <- (1 + (alpha * h_cm)^n_par)^(-m)
  thetar + (thetas - thetar) * Se
}

# one-layer sampler + summarizer (robust)
mc_soil_profile_df2 <- function(df, N = 20000, seed = 123) {
  set.seed(seed)
  nlay <- nrow(df)
  out_rows <- vector("list", nlay)
  
  # helper to compute stats safely
  safe_stats <- function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (length(x) == 0) return(rep(NA_real_, 5))
    q <- as.numeric(quantile(x, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE))
    c(
      mean = mean(x, na.rm = TRUE),
      low95 = q[1],
      high95 = q[2],
      min = min(x, na.rm = TRUE),
      max = max(x, na.rm = TRUE)
    )
  }
  
  for (i in seq_len(nlay)) {
    # extract mu and sd columns for this row (strip prefixes)
    mu_row <- df[i, grep("^mu_", names(df)), drop = FALSE]
    sd_row <- df[i, grep("^sd_", names(df)), drop = FALSE]
    names(mu_row) <- sub("^mu_", "", names(mu_row))
    names(sd_row) <- sub("^sd_", "", names(sd_row))
    mu <- as.list(mu_row)
    sd <- as.list(sd_row)
    
    # required checks (stop early with informative error if missing)
    need_mu <- c("thetas", "thetar")
    if (!all(need_mu %in% names(mu))) {
      stop("Missing mu_thetas or mu_thetar in row ", i)
    }
    need_sd <- c("thetas", "thetar")
    if (!all(need_sd %in% names(sd))) {
      stop("Missing sd_thetas or sd_thetar in row ", i)
    }
    
    # sample linear-theta params
    thetas <- rnorm(N, mu$thetas, sd$thetas)
    thetar <- rnorm(N, mu$thetar, sd$thetar)
    
    # handle alpha (either log10_alpha or alpha)
    if ("log10_alpha" %in% names(mu)) {
      if (!("log10_alpha" %in% names(sd))) stop("sd_log10_alpha missing for row ", i)
      log10_alpha <- rnorm(N, mu$log10_alpha, sd$log10_alpha)
      alpha <- 10^log10_alpha
    } else if ("alpha" %in% names(mu)) {
      if (!("alpha" %in% names(sd))) stop("sd_alpha missing for row ", i)
      alpha <- rnorm(N, mu$alpha, sd$alpha)
    } else {
      stop("No alpha or log10_alpha in mu for row ", i)
    }
    
    # handle n (either log10_n or n)
    if ("log10_n" %in% names(mu)) {
      if (!("log10_n" %in% names(sd))) stop("sd_log10_n missing for row ", i)
      log10_n <- rnorm(N, mu$log10_n, sd$log10_n)
      n_par <- 10^log10_n
    } else if ("n" %in% names(mu)) {
      if (!("n" %in% names(sd))) stop("sd_n missing for row ", i)
      n_par <- rnorm(N, mu$n, sd$n)
    } else {
      stop("No n or log10_n in mu for row ", i)
    }
    
    # handle ks (optional)
    if ("log10_ks" %in% names(mu)) {
      if (!("log10_ks" %in% names(sd))) stop("sd_log10_ks missing for row ", i)
      log10_ks <- rnorm(N, mu$log10_ks, sd$log10_ks)
      ks <- 10^log10_ks
    } else if ("ks" %in% names(mu)) {
      if (!("ks" %in% names(sd))) stop("sd_ks missing for row ", i)
      ks <- rnorm(N, mu$ks, sd$ks)
    } else {
      # if ks missing entirely, fill with NA
      ks <- rep(NA_real_, N)
    }
    
    # enforce physical constraints
    thetas <- pmin(pmax(thetas, 0), 1)
    thetar <- pmin(pmax(thetar, 0), thetas - 1e-8)
    n_par  <- pmax(n_par, 1.01)
    ks     <- ifelse(is.finite(ks), pmax(ks, 1e-12), NA_real_)
    
    # compute retention points and PAW
    theta33   <- vg_theta(33, thetas, thetar, alpha, n_par)
    theta1500 <- vg_theta(1500, thetas, thetar, alpha, n_par)
    paw       <- theta33 - theta1500
    
    # build row of stats
    thetas_stats   <- safe_stats(thetas)
    theta33_stats  <- safe_stats(theta33)
    theta1500_stats<- safe_stats(theta1500)
    ks_stats       <- safe_stats(ks)
    paw_stats      <- safe_stats(paw)
    
    out_rows[[i]] <- data.frame(
      layer_index = i,
      
      thetas_mean   = thetas_stats["mean"],
      thetas_low95  = thetas_stats["low95"],
      thetas_high95 = thetas_stats["high95"],
      thetas_min    = thetas_stats["min"],
      thetas_max    = thetas_stats["max"],
      
      theta33_mean   = theta33_stats["mean"],
      theta33_low95  = theta33_stats["low95"],
      theta33_high95 = theta33_stats["high95"],
      theta33_min    = theta33_stats["min"],
      theta33_max    = theta33_stats["max"],
      
      theta1500_mean   = theta1500_stats["mean"],
      theta1500_low95  = theta1500_stats["low95"],
      theta1500_high95 = theta1500_stats["high95"],
      theta1500_min    = theta1500_stats["min"],
      theta1500_max    = theta1500_stats["max"],
      
      ks_mean   = ks_stats["mean"],
      ks_low95  = ks_stats["low95"],
      ks_high95 = ks_stats["high95"],
      ks_min    = ks_stats["min"],
      ks_max    = ks_stats["max"],
      
      paw_mean   = paw_stats["mean"],
      paw_low95  = paw_stats["low95"],
      paw_high95 = paw_stats["high95"],
      paw_min    = paw_stats["min"],
      paw_max    = paw_stats["max"],
      stringsAsFactors = FALSE
    )
  } # end layers loop
  
  out_df <- do.call(rbind, out_rows)
  rownames(out_df) <- NULL
  out_df
}
# run it (example N smaller for speed while testing)
out_df <- mc_soil_profile_df2(soil_df, N = 10000, seed = 42)
head(out_df)

out_ny_low <- out_df %>%
  dplyr::select(layer_index, thetas_mean, theta33_low95, theta1500_high95, ks_mean) 

out_ny_high <- out_df %>%
  dplyr::select(layer_index, thetas_mean, theta33_high95, theta1500_low95, ks_mean) 

out_ny_mean <- out_df %>%
  dplyr::select(layer_index, thetas_mean, theta33_mean, theta1500_mean, ks_mean) 

out_df_wi <- mc_soil_profile_df2(soil_df_plano, N = 10000, seed = 42)

out_wi_low <- out_df_wi %>%
  dplyr::select(layer_index, thetas_mean, theta33_low95, theta1500_high95, ks_mean) 

out_wi_high <- out_df_wi %>%
  dplyr::select(layer_index, thetas_mean, theta33_high95, theta1500_low95, ks_mean) 

out_wi_mean <- out_df_wi %>%
  dplyr::select(layer_index, thetas_mean, theta33_mean, theta1500_mean, ks_mean) 


#### ChatGPT recmmendation to get high and low values: 

vg_theta <- function(h_kpa, thetas, thetar, alpha, n_par) {
  h_cm <- h_kpa * 10.197        
  m <- 1 - 1 / n_par
  Se <- (1 + (alpha * h_cm)^n_par)^(-m)
  thetar + (thetas - thetar) * Se
}


summarize <- function(x) {
  c(
    mean   = mean(x),
    low95  = quantile(x, 0.025),
    high95 = quantile(x, 0.975),
    min    = min(x),
    max    = max(x)
  )
}


mc_soil_layer <- function(mu, sd, N=20000, seed=NULL, keep_draws=FALSE) {
  if (!is.null(seed)) set.seed(seed)
  
  thetas <- rnorm(N, mu$thetas, sd$thetas)
  thetar <- rnorm(N, mu$thetar, sd$thetar)
  
  log10_alpha <- rnorm(N, mu$log10_alpha, sd$log10_alpha)
  alpha <- 10^log10_alpha
  log10_n <- rnorm(N, mu$log10_n, sd$log10_n)
  n_par <- 10^log10_n
  log10_ks <- rnorm(N, mu$log10_ks, sd$log10_ks)
  ks <- 10^log10_ks
  
  thetas <- pmin(pmax(thetas, 0), 1)
  thetar <- pmin(pmax(thetar, 0), thetas - 1e-8)
  n_par  <- pmax(n_par, 1.01)
  
  theta33   <- vg_theta(33, thetas, thetar, alpha, n_par)
  theta1500 <- vg_theta(1500, thetas, thetar, alpha, n_par)
  paw <- theta33 - theta1500
  
  summaries <- list(
    thetas   = summarize(thetas),
    theta33  = summarize(theta33),
    theta1500= summarize(theta1500),
    ks       = summarize(ks),
    paw      = summarize(paw)
  )
  
  if (keep_draws) {
    return(list(summary = summaries,
                draws = data.frame(thetas, thetar, alpha, n_par, ks,
                                   theta33, theta1500, paw)))
  } else {
    return(summaries)
  }
}

mc_soil_profile <- function(df, N=20000, seed=123, keep_draws=FALSE) {
  results <- vector("list", nrow(df))
  for (i in seq_len(nrow(df))) {
    mu <- as.list(df[i, grep("^mu_", names(df))])
    names(mu) <- sub("^mu_", "", names(mu))
    sd <- as.list(df[i, grep("^sd_", names(df))])
    names(sd) <- sub("^sd_", "", names(sd))
    results[[i]] <- mc_soil_layer(mu, sd, N, seed+i, keep_draws=keep_draws)
  }
  results
}



make_high_low_profiles <- function(out_list, probs = c(0.025, 0.5, 0.975)) {
  nlay <- length(out_list)
  
  high_list <- vector("list", nlay)
  mean_list <- vector("list", nlay)
  low_list  <- vector("list", nlay)
  
  for (i in seq_len(nlay)) {
    draws <- out_list[[i]]$draws
    
    # Identify the draw closest to the upper and lower quantiles of PAW
    q_low  <- quantile(draws$paw, probs[1])
    q_mean <- quantile(draws$paw, probs[2])
    q_high <- quantile(draws$paw, probs[3])
    
    # which draw is closest to those quantiles?
    idx_low  <- which.min(abs(draws$paw - q_low))
    idx_mean <- which.min(abs(draws$paw - q_mean))
    idx_high <- which.min(abs(draws$paw - q_high))
    
    low_list[[i]]  <- draws[idx_low, ]
    mean_list[[i]] <- draws[idx_mean, ]
    high_list[[i]] <- draws[idx_high, ]
  }
  
  # stack into data.frames
  low_df  <- do.call(rbind, low_list)
  mean_df <- do.call(rbind, mean_list)
  high_df <- do.call(rbind, high_list)
  
  list(low = low_df, high = high_df, mean = mean_df)
}

out <- mc_soil_profile(soil_df, N=20000, seed=123, keep_draws=TRUE)

profiles_lima <- make_high_low_profiles(out)
high_soil_lima <- profiles_lima$high %>%
  mutate(layer = as.factor(1:6), 
         scenario = rep("High"))
low_soil_lima  <- profiles_lima$low %>%
  mutate(layer = as.factor(1:6),
         scenario = rep("Low"))
mean_soil_lima  <- profiles_lima$mean %>%
  mutate(layer = as.factor(1:6),
         scenario = rep("Mean"))


lima_high_low_mean <- rbind(high_soil_lima, low_soil_lima, mean_soil_lima) %>%
  group_by(scenario) %>%
  mutate(across(thetas:paw,
                ~ ifelse(layer == 2, .[layer == 1], .)))

lima_bias <- lima_high_low_mean %>%
  group_by(layer) %>%
  summarise(bias_thetas_high = thetas[scenario == "High"]-thetas[scenario == "Mean"],
            bias_thetas_low = thetas[scenario == "Low"]-thetas[scenario == "Mean"],
            bias_theta33_low = theta33[scenario == "Low"]-theta33[scenario == "Mean"],
            bias_theta33_high = theta33[scenario == "High"]-theta33[scenario == "Mean"],
            bias_theta1500_low = theta1500[scenario == "Low"]-theta1500[scenario == "Mean"],
            bias_theta1500_high = theta1500[scenario == "High"]-theta1500[scenario == "Mean"])

lima_high_low_curves <- lima_high_low_mean %>%
  group_by(scenario, layer) %>%
  summarise(matric = seq(1,2000,1),
            vwc = vg_theta(matric, thetas, thetar, alpha, n_par))

layer_labels <- c("1" = "0-15 cm",
                  "2" = "15-25 cm",
                  "3" = "25-30 cm",
                  "4" = "30-46 cm",
                  "5" = "46-86 cm",
                  "6" = "86-104 cm")

lima_vang <- lima_high_low_curves %>%
  ggplot(aes(x = matric, y = vwc, color = scenario)) +
  geom_line() +
  geom_vline(xintercept = c(33,1500), linetype = "dashed", color = "gray", alpha = 0.6) +
  facet_wrap(~layer, labeller = as_labeller(layer_labels)) +
  theme_pubr() +
  scale_x_log10() +
  scale_color_manual(values = c("black","red", "blue")) +
  labs(y = "Volumetric water content",
       x = "Matric potential (kPa)",
       color = "")
ggsave(lima_vang, filename = "Figures/lima_vang.png",
       dpi = 300, height = 6, width = 6, units = "in")

#### Creating a soil water retention curve example 

mc_retention_curve <- function(mu, sd, h_seq = exp(seq(log(1), log(15000), length.out=200)),
                               N=20000, seed=123) {
  set.seed(seed)
  
  # Linear params
  thetas <- rnorm(N, mu$thetas, sd$thetas)
  thetar <- rnorm(N, mu$thetar, sd$thetar)
  
  # Log-space params
  log10_alpha <- rnorm(N, mu$log10_alpha, sd$log10_alpha)
  alpha <- 10^log10_alpha
  log10_n <- rnorm(N, mu$log10_n, sd$log10_n)
  n_par <- 10^log10_n
  
  # Constraints
  thetas <- pmin(pmax(thetas, 0), 1)
  thetar <- pmin(pmax(thetar, 0), thetas - 1e-8)
  n_par  <- pmax(n_par, 1.01)
  
  # Storage for results
  res <- data.frame(h_kpa = h_seq,
                    mean = NA, low95 = NA, high95 = NA)
  
  # Loop over each suction value
  for (j in seq_along(h_seq)) {
    theta_vals <- vg_theta(h_seq[j], thetas, thetar, alpha, n_par)
    res$mean[j]   <- mean(theta_vals)
    res$low95[j]  <- quantile(theta_vals, 0.025)
    res$high95[j] <- quantile(theta_vals, 0.975)
  }
  
  res
}

mu <- list(thetas = lima_lab_ros$theta_s[1], thetar = lima_lab_ros$theta_r[1],
           log10_alpha = lima_lab_ros$alpha[1], log10_n = lima_lab_ros$npar[1])
sd <- list(thetas = lima_lab_ros$sd_theta_s[1], thetar = lima_lab_ros$sd_theta_r[1],
           log10_alpha = lima_lab_ros$sd_alpha[1], log10_n = lima_lab_ros$sd_npar[1])

curve_out <- mc_retention_curve(mu, sd, N=5000)
ros_curv_ex <- ggplot(curve_out, aes(x=h_kpa)) +
  geom_ribbon(aes(ymin=low95, ymax=high95),color = "black",linetype = "dashed", alpha = 0) +
  geom_line(aes(y=mean), size=1) +
  scale_x_log10() +
  labs(x="Matric potential (kPa, log scale)",
       y="Volumetric water content (cm³/cm³)") +
  theme_pubr() +
  coord_cartesian(xlim = c(1, 2000)) +
  theme(axis.title = element_text(size = 14))
ggsave(ros_curv_ex, filename = "Figures/lima_top_ros_curve.png",
       dpi = 300, width = 7, height = 6, units = "in")

#### Now for Wisconsin
arl_soil <- get_ssurgo_soil_profile(lonlat = c(-89.33146, 43.30458))
get_ssurgo_soil_profile(lonlat = c(-89.33146, 43.30458), nlayers = 7)
get_ss
arl_key <- SDA_spatialQuery(st_sfc(st_point(c(-89.33146, 43.30458)), crs = 4236), what = "mukey")

arl_soildb <- fetchSDA("mukey = ('423333')", duplicates = T)

plano_osd <- fetchOSD("Plano")

plano_osd_horizons <- plano_osd@site

plano_cokey <- arl_soildb@site$cokey[arl_soildb@site$compname == "Plano"]

arl_horizons <- arl_soildb@horizons 

plano_horizons <- arl_horizons %>%
  filter(cokey == plano_cokey) %>%
  dplyr::select(hzname, hzdept_r, hzdepb_r, fragvol_r, sandtotal_r, silttotal_r, claytotal_r,
                om_r, dbthirdbar_r, ksat_r, awc_r, ph1to1h2o_r) %>%
  unite(col = "depth", c("hzdept_r", "hzdepb_r"), sep = "-") %>%
  rename("rock" = fragvol_r, "sand" = sandtotal_r, "silt" = silttotal_r, "clay" = claytotal_r,
         "om" = om_r, "bd" = dbthirdbar_r, "ksat" = ksat_r, "awc" = awc_r, "ph" = ph1to1h2o_r)

plano_text <- data.frame(
  sand = c(3.9, 3.4, 3.8, 3.6, 3.2, 3.5, 12.6, 50.4, 62.4, 72.4, 72.9, 72.4),
  silt = c(71.7, 73.6, 77, 70.6, 65.9, 68.4, 64.6, 35.1, 24.6, 20.7, 22.2, 23.1),
  clay = c(24.4, 23, 19.2, 25.8, 30.9, 28.1, 22.8, 14.5, 13, 6.9, 4.9, 4.5),
  bd = c(1.34, 1.34, 1.42, 1.42, 1.42, 1.42, 1.6, 1.6, 1.5, 1.5, 1.5, 1.5)
)


plano_rosetta <-ROSETTA(plano_text, vars = names(plano_text), v = "3", include.sd = T)

soil_df_plano <- data.frame(
  mu_thetas = plano_rosetta$theta_s,
  mu_thetar = plano_rosetta$theta_r,
  mu_log10_alpha = plano_rosetta$alpha,
  mu_log10_n = plano_rosetta$npar,
  mu_log10_ks = plano_rosetta$ksat,
  sd_thetas = plano_rosetta$sd_theta_s,
  sd_thetar = plano_rosetta$sd_theta_r,
  sd_log10_alpha = plano_rosetta$sd_alpha,
  sd_log10_n = plano_rosetta$sd_npar,
  sd_log10_ks = plano_rosetta$sd_ksat
)

out_plano <- mc_soil_profile(soil_df_plano, N=20000, seed=123, keep_draws=TRUE)

profiles_plano <- make_high_low_profiles(out_plano)
high_soil_plano <- profiles_plano$high %>%
  mutate(layer = as.factor(1:12), 
         scenario = rep("High"))
low_soil_plano  <- profiles_plano$low %>%
  mutate(layer = as.factor(1:12),
         scenario = rep("Low"))
mean_soil_plano  <- profiles_plano$mean %>%
  mutate(layer = as.factor(1:12),
         scenario = rep("Mean"))


plano_high_low_mean <- rbind(high_soil_plano, low_soil_plano, mean_soil_plano) %>%
  group_by(scenario) %>%
  mutate(across(thetas:paw,
                ~ ifelse(layer == 2, .[layer == 1], .)))
plano_bias <- plano_high_low_mean %>%
  group_by(layer) %>%
  summarise(bias_thetas_high = thetas[scenario == "High"]-thetas[scenario == "Mean"],
            bias_thetas_low = thetas[scenario == "Low"]-thetas[scenario == "Mean"],
            bias_theta33_low = theta33[scenario == "Low"]-theta33[scenario == "Mean"],
            bias_theta33_high = theta33[scenario == "High"]-theta33[scenario == "Mean"],
            bias_theta1500_low = theta1500[scenario == "Low"]-theta1500[scenario == "Mean"],
            bias_theta1500_high = theta1500[scenario == "High"]-theta1500[scenario == "Mean"])


plano_high_low_curves <- plano_high_low_mean %>%
  group_by(scenario, layer) %>%
  summarise(matric = seq(1,2000,1),
            vwc = vg_theta(matric, thetas, thetar, alpha, n_par))

layer_labels <- c("1" = "0-20 cm",
                  "2" = "20-38 cm",
                  "3" = "38-51 cm",
                  "4" = "51-66 cm",
                  "5" = "66-94 cm",
                  "6" = "94-112 cm",
                  "7" = "112-150 cm",
                  "8" = "150-163 cm",
                  "9" = "163-178 cm",
                  "10" = "178-190 cm",
                  "11" = "190-221 cm",
                  "12" = "221-251 cm")

plano_vang <- plano_high_low_curves %>%
  ggplot(aes(x = matric, y = vwc, color = scenario)) +
  geom_line() +
  geom_vline(xintercept = c(33,1500), linetype = "dashed", color = "gray", alpha = 0.6) +
  facet_wrap(~layer, labeller = as_labeller(layer_labels)) +
  theme_pubr() +
  scale_x_log10() +
  scale_color_manual(values = c("black","red", "blue")) +
  labs(y = "Volumetric water content",
       x = "Matric potential (kPa)",
       color = "")
ggsave(plano_vang, filename = "Figures/plano_vang.png",
       dpi = 300, height = 6, width = 6, units = "in")

#### Following examples in miguez video

mus_sps <- get_ssurgo_tables(lonlat = c(-76.65055, 42.73376), shift = 300)
mus_sps$mapunit.shp$mukey <- as.factor(mus_sps$mapunit.shp$mukey)
plot(mus_sps$mapunit.shp[, "mukey"], key.pos = 1)


mus_lima  <- get_ssurgo_soil_profile(lonlat = c(-76.6507612677327, 42.734075011700575))
mus_lima_apsimx <- mus_lima[[1]]$soil

muslima_text <- mus_lima_apsimx %>%
  dplyr::select(ParticleSizeSand, ParticleSizeSilt, ParticleSizeClay, BD) %>%
  rename("sand" = ParticleSizeSand, "silt" = ParticleSizeSilt, "clay" = ParticleSizeClay,"bd" = BD)
muslima_ros <- ROSETTA(muslima_text, vars = names(muslima_text), v = "3", include.sd = T)

soil_df_ssurgo <- data.frame(
  mu_thetas = muslima_ros$theta_s,
  mu_thetar = muslima_ros$theta_r,
  mu_log10_alpha = muslima_ros$alpha,
  mu_log10_n = muslima_ros$npar,
  mu_log10_ks = muslima_ros$ksat,
  sd_thetas = muslima_ros$sd_theta_s,
  sd_thetar = muslima_ros$sd_theta_r,
  sd_log10_alpha = muslima_ros$sd_alpha,
  sd_log10_n = muslima_ros$sd_npar,
  sd_log10_ks = muslima_ros$sd_ksat
)
out_ssurgo <- mc_soil_profile(soil_df_ssurgo, N=20000, seed=123, keep_draws=TRUE)

profiles_lima_apsimx <- make_high_low_profiles(out_ssurgo)


#### Getting other soil profiles for APSIM use from top soil series in NY and WI

ontario_apsimx <- get_ssurgo_soil_profile(lonlat = c(-77.28298, 43.03276))

hilton_apsimx <- get_ssurgo_soil_profile(lonlat = c(-77.0841, 43.0152))

honeyoye_apsimx <- get_ssurgo_soil_profile(lonlat = c(-76.6126, 42.6985))

appleton_apsimx <- get_ssurgo_soil_profile(lonlat = c(-76.6126, 42.6985), nsoil = 4)

lima_apsimx_representative <- get_ssurgo_soil_profile(lonlat = c(-76.6126, 42.6985), nsoil = 2)

southern_tier <- get_ssurgo_soil_profile(lonlat = c(-76.3941, 42.1728), 2)
