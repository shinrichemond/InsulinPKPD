library(deSolve)
library(dplyr)
library(ggplot2)
library(patchwork)
library(shadowtext)
library(metR)

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
Vmax_vals <- seq(5, 20, 0.25)
Km_vals   <- seq(25, 55, length.out = length(Vmax_vals))

# Create progress bar
n_iter <- length(Vmax_vals) * length(Km_vals)
pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
iter <- 0

results <- data.frame()

# Run scan
for (v in Vmax_vals) {
  for (k in Km_vals) {

    iter <- iter + 1
    setTxtProgressBar(pb, iter)
    
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
    
    last_dose <- max(dose_times)
    window <- 4   # hours after injection to search
    
    df_window <- df %>%
      filter(time >= last_dose & time <= last_dose + window)
    
    if (nrow(df_window) > 0) {
      
      # Cmax
      idx_max <- which.max(df_window$P)
      Cmax <- df_window$P[idx_max]
      
      # Tmax (relative to injection)
      Tmax <- df_window$time[idx_max] - last_dose
      
      # AUC (trapezoidal)
      AUC <- sum(diff(df_window$time) *
                   (head(df_window$P, -1) + tail(df_window$P, -1)) / 2)
      
    } else {
      Cmax <- NA
      Tmax <- NA
      AUC  <- NA
    }
    
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

close(pb)  # Close progress bar

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
breaks <- seq(0, max(results$Cmax, na.rm = TRUE), by = 2)

p_cmax <- ggplot(results, aes(x = Vmax, y = Km)) +
  geom_tile(aes(fill = Cmax)) + 
  geom_contour(aes(z = Cmax), breaks = breaks, color = "white") +
  metR::geom_text_contour(aes(z = Cmax), breaks=breaks, stroke = 0.08) +
  scale_fill_viridis_c(option = "viridis", name = "Cmax") +
  labs(
    title = "Cmax",
    #subtitle = "Peak plasma insulin across Vmax–Km parameter space",
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

p_cmax

# tmax plot
breaks <- seq(0, max(results$Tmax, na.rm = TRUE), by = 0.2)

p_tmax <- ggplot(results, aes(x = Vmax, y = Km)) +
  geom_tile(aes(fill = Tmax)) + 
  geom_contour(aes(z = Tmax), breaks = breaks, color = "white") +
  metR::geom_text_contour(aes(z = Tmax), breaks=breaks, stroke = 0.08) +  
  scale_fill_viridis_c(option = "viridis", name = "Tmax") +
  labs(
    title = "Tmax",
    #subtitle = "Time to peak plasma insulin",
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
    vjust = -1.5,
    bg.colour = "black",
    bg.r = 0.1
  )

p_tmax

# AUC plot
breaks <- seq(0, max(results$AUC, na.rm = TRUE), by = 10)

p_auc <- ggplot(results, aes(x = Vmax, y = Km)) +
  geom_tile(aes(fill = AUC)) +
  geom_contour(aes(z = AUC), breaks = breaks, color = "white") +
  metR::geom_text_contour(aes(z = AUC), breaks=breaks, stroke = 0.08) +  
  scale_fill_viridis_c(option = "viridis", name = "AUC") +
  labs(
    title = "AUC",
    #subtitle = "Total systemic insulin exposure",
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

