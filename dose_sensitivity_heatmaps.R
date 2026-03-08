library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)

this_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(this_dir)
source(file.path(this_dir, "insulinRscript.R"))
source(file.path(this_dir, "Simulate.R"))

# ============================================================================
# Dose Sensitivity Heatmaps: Dose x Injection Site for Cmax, Tmax, AUC
# ============================================================================

sweep_doses <- c(5, 7, 10, 15)

dose_sweep_grid <- expand.grid(
  site_idx = 1:nrow(site_profiles),
  dose_val = sweep_doses,
  stringsAsFactors = FALSE
)

dose_sweep_grid$site <- site_profiles$site[dose_sweep_grid$site_idx]
dose_sweep_grid$Cmax <- NA_real_
dose_sweep_grid$Tmax <- NA_real_
dose_sweep_grid$AUC  <- NA_real_

for (r in seq_len(nrow(dose_sweep_grid))) {
  i <- dose_sweep_grid$site_idx[r]
  d <- dose_sweep_grid$dose_val[r]

  local_dose_sizes <- rep(d, length(dose_times))

  local_ode <- function(t, state, params) {
    with(as.list(c(state, params)), {
      t_day <- t %% 24

      u_input <- 0
      for (j in seq_along(dose_times)) {
        if (t_day >= dose_times[j] && t_day < dose_times[j] + inj_dur) {
          u_input <- u_input + local_dose_sizes[j] / inj_dur
        }
      }

      G_in <- 0
      for (j in seq_along(meal_times)) {
        t_rel <- t_day - meal_times[j]
        if (t_rel >= 0 && t_rel <= meal_durs[j]) {
          x <- t_rel / meal_durs[j]
          G_in <- G_in + meal_sizes[j] * 4 * x * (1 - x)
        }
      }

      S_pos <- max(S, 0); P_pos <- max(P, 0); G_pos <- max(G, 0)

      kabs <- Vmax * S_pos / (Km + S_pos)
      E_G   <- Emax * G_pos^hillN / (EC50^hillN + G_pos^hillN)
      C_t   <- 1 / (1 + P_pos / IC50)
      sens  <- sens_func(t_day)
      U_PG  <- Umax * P_pos * G_pos / ((KP + P_pos) * (KG + G_pos))
      H     <- kH * (G_pos - G_set)

      dS <- u_input - kabs - kdeg * S_pos + kre * P_pos
      dP <- Vc * (kabs + E_G * C_t + I_basal) - kclr * P_pos - kenz * P_pos - kre * P_pos
      dG <- Gb + G_in - sens * U_PG - k0 * G_pos - H

      return(list(c(dS, dP, dG)))
    })
  }

  params_site <- modifyList(base_params, list(
    Vmax = site_profiles$Vmax[i],
    Km   = site_profiles$Km[i],
    day_length = day_length
  ))

  out <- ode(y = state0, times = times, func = local_ode,
             parms = params_site, method = "lsoda")
  out_df_local <- as.data.frame(out)

  dose_sweep_grid$Cmax[r] <- max(out_df_local$P)
  dose_sweep_grid$Tmax[r] <- out_df_local$time[which.max(out_df_local$P)]
  dose_sweep_grid$AUC[r]  <- sum(out_df_local$P) * (diff(out_df_local$time)[1])
}

dose_sweep_grid <- dose_sweep_grid %>%
  mutate(Dose = factor(dose_val, levels = sweep_doses))

# -----------------------------
# Shared theme
heatmap_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey50", size = 11),
    panel.grid    = element_blank()
  )

# -----------------------------
# 1. Cmax Heatmap
p_cmax <- ggplot(dose_sweep_grid, aes(x = site, y = Dose, fill = Cmax)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = round(Cmax, 1)), size = 4) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = median(dose_sweep_grid$Cmax),
                       name = "Cmax") +
  labs(title = "Dose Sensitivity: Cmax (Peak Plasma Insulin)",
       subtitle = "Uniform per-meal dose swept across all injection sites",
       x = "Injection Site", y = "Dose (Units)") +
  heatmap_theme

p_cmax

# -----------------------------
# 2. Tmax Heatmap
p_tmax <- ggplot(dose_sweep_grid, aes(x = site, y = Dose, fill = Tmax)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = round(Tmax, 2)), size = 4) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = median(dose_sweep_grid$Tmax),
                       name = "Tmax (hr)") +
  labs(title = "Dose Sensitivity: Tmax (Time to Peak)",
       subtitle = "Uniform per-meal dose swept across all injection sites",
       x = "Injection Site", y = "Dose (Units)") +
  heatmap_theme

p_tmax

# -----------------------------
# 3. AUC Heatmap
p_auc <- ggplot(dose_sweep_grid, aes(x = site, y = Dose, fill = AUC)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = round(AUC, 0)), size = 4) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = median(dose_sweep_grid$AUC),
                       name = "AUC") +
  labs(title = "Dose Sensitivity: AUC (Total Insulin Exposure)",
       subtitle = "Uniform per-meal dose swept across all injection sites",
       x = "Injection Site", y = "Dose (Units)") +
  heatmap_theme

p_auc
