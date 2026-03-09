library(deSolve)
library(dplyr)
library(ggplot2)
library(patchwork)
library(shadowtext)

# get system ODE
this_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(this_dir)
source(file.path(this_dir, "insulinRscript.R"))

# Time grid
times <- seq(0, 24, by = 0.01)

# Initial conditions
state_init <- c(
  S = 0,
  P = 10,
  G = base_params$G_set
)

# Parameter ranges
Vmax_vals <- seq(5, 20)
Km_vals   <- seq(25, 55, length.out = length(Vmax_vals))

results <- data.frame()

# Run scan
for (v in Vmax_vals) {
  for (k in Km_vals) {
    
    params <- base_params
    params$Vmax <- v
    params$Km   <- k
    
    out <- ode(
      y = state_init,
      times = times,
      func = insulin_ode,
      parms = params
    )
    
    df <- as.data.frame(out)
    
    # PK metrics
    Cmax <- max(df$P)
    Tmax <- df$time[which.max(df$P)]
    
    # trapezoidal AUC
    AUC <- sum(diff(df$time) * 
                 (head(df$P,-1) + tail(df$P,-1)) / 2)
    
    results <- rbind(
      results,
      data.frame(
        Vmax = v,
        Km = k,
        Cmax = Cmax,
        Tmax = Tmax,
        AUC = AUC
      )
    )
  }
}


# Plots (stolen from aditya)
# site parameters:
site_params <- data.frame(
  site  = c("Abdomen", "Upper Arm", "Buttock", "Thigh"),
  Vmax  = c(18, 13, 10, 8),
  Km    = c(30, 40, 42, 50),
  color = c("#c53030", "#d69e2e", "#2b6cb0", "#718096"),
  stringsAsFactors = FALSE
)


heatmap_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey50", size = 11),
    panel.grid    = element_blank()
  )

# Cmax plot
p_cmax <- ggplot(results, aes(x = Vmax, y = Km, fill = Cmax)) +
  geom_tile(color = "white", linewidth = 0.8) +
  scale_fill_viridis_c(option = "viridis", name = "Cmax") +
  labs(
    title = "Absorption Parameter Sensitivity: Cmax",
    subtitle = "Peak plasma insulin across Vmax–Km parameter space",
    x = "Vmax (Absorption Capacity)",
    y = "Km (Diffusion Resistance)"
  ) +
  heatmap_theme + geom_point(
    data = site_params,
    aes(x = Vmax, y = Km),
    inherit.aes = FALSE,
    size = 3,
    shape = 24,
    fill = "black",
    color = "white",
    stroke = 1.2
  ) + shadowtext::geom_shadowtext(
    data = site_params,
    aes(x = Vmax, y = Km, label = site),
    inherit.aes = FALSE,
    size = 3.5,
    vjust = 1.5,
    bg.colour = "black",
    bg.r = 0.1
  )

p_cmax

# tmax plot
p_tmax <- ggplot(results, aes(x = Vmax, y = Km, fill = Tmax)) +
  geom_tile(color = "white", linewidth = 0.8) +
  scale_fill_viridis_c(option = "viridis", name = "Tmax") +
  labs(
    title = "Absorption Parameter Sensitivity: Tmax",
    subtitle = "Time to peak plasma insulin",
    x = "Vmax",
    y = "Km"
  ) +
  heatmap_theme + geom_point(
    data = site_params,
    aes(x = Vmax, y = Km),
    inherit.aes = FALSE,
    size = 3,
    shape = 24,
    fill = "black",
    color = "white",
    stroke = 1.2
  ) + shadowtext::geom_shadowtext(
    data = site_params,
    aes(x = Vmax, y = Km, label = site),
    inherit.aes = FALSE,
    size = 3.5,
    vjust = 1.5,
    bg.colour = "black",
    bg.r = 0.1
  )

p_tmax

# AUC plot
p_auc <- ggplot(results, aes(x = Vmax, y = Km, fill = AUC)) +
  geom_tile(color = "white", linewidth = 0.8) +
  scale_fill_viridis_c(option = "viridis", name = "AUC") +
  labs(
    title = "Absorption Parameter Sensitivity: AUC",
    subtitle = "Total systemic insulin exposure",
    x = "Vmax",
    y = "Km"
  ) +
  heatmap_theme + geom_point(
    data = site_params,
    aes(x = Vmax, y = Km),
    inherit.aes = FALSE,
    size = 3,
    shape = 24,
    fill = "black",
    color = "white",
    stroke = 1.2
  ) + shadowtext::geom_shadowtext(
    data = site_params,
    aes(x = Vmax, y = Km, label = site),
    inherit.aes = FALSE,
    size = 3.5,
    vjust = 1.5,
    bg.colour = "black",
    bg.r = 0.1
  )

p_auc

# Plot all
p_cmax | p_tmax | p_auc
