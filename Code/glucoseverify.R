library(deSolve)
library(dplyr)
library(ggplot2)
library(tidyr)

Vmax_vals <- seq(5, 100, length.out = 10)
Km_vals   <- seq(5, 250, length.out = 10)
param_grid <- expand.grid(Vmax = Vmax_vals, Km = Km_vals)

# store results
all_results <- list()

for(i in seq_len(nrow(param_grid))){
  params <- base_params
  params$Vmax <- param_grid$Vmax[i]
  params$Km   <- param_grid$Km[i]
  
  out <- ode(
    y = c(S = 0, P = 0, G = 90), 
    times = seq(0, 24, by = 0.05),  # 24 hours with 3 min resolution
    func = insulin_ode,
    parms = params
  )
  
  out_df <- as.data.frame(out)
  out_df$Vmax <- param_grid$Vmax[i]
  out_df$Km   <- param_grid$Km[i]
  
  all_results[[i]] <- out_df
}

results_df <- bind_rows(all_results)

results_df %>%
  group_by(time) %>%
  summarize(
    G_mean = mean(G),
    G_min  = min(G),
    G_max  = max(G)
  ) %>%
  ggplot(aes(x = time)) +
  geom_ribbon(aes(ymin = G_min, ymax = G_max), fill = "steelblue", alpha = 0.3) +
  geom_line(aes(y = G_mean), color = "steelblue", linewidth = 1) +
  labs(
    title = "Plasma Glucose Robustness Across Absorption Parameters",
    subtitle = "Shaded area shows min–max glucose for all Vmax/Km combinations",
    x = "Time (hours)", y = "Glucose (mg/dL)"
  ) +
  theme_minimal(base_size = 12)