library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(plotly)
library(knitr)
library(kableExtra)

this_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(this_dir)
source(file.path(this_dir, "insulinRscript.R"))
source(file.path(this_dir, "Simulate.R"))

out_df <- readRDS(file.path(this_dir, "data.rds"))

# ============================================================================
# TIME SERIES: S(t), P(t), G(t) FOR ALL INJECTION SITES

site_colors <- setNames(site_profiles$color, site_profiles$site)

p_S <- ggplot(out_df, aes(x = time, y = S, color = site)) +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  scale_color_manual(values = site_colors) +
  labs(title="Subcutaneous Depot S(t)", y="S", color="Site") +
  theme_minimal()

p_P <- ggplot(out_df, aes(x = time, y = P, color = site)) +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  scale_color_manual(values = site_colors) +
  labs(title="Plasma Insulin P(t)", y="P", color="Site") +
  theme_minimal()

p_G <- ggplot(out_df, aes(x = time, y = G, color = site)) +
  geom_hline(yintercept=c(70,180), linetype="dashed", color="grey60") +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=70, ymax=180,
           fill="green", alpha=0.04) +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  scale_color_manual(values = site_colors) +
  labs(title="Blood Glucose G(t)", y="G", x="Time (hours)", color="Site") +
  theme_minimal()

p_S / p_P / p_G +
  plot_annotation(
    title = "24-Hour Dynamics by Injection Site",
    subtitle = "Same meals, same dosing — only absorption kinetics vary",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(color = "grey40", size = 12)
    )
  )

# =====================
# Glucose Drift Check

# Define physiological range
fasting_min <- 70
fasting_max <- 180

# Glucose Drift Plot
p_drift <- ggplot(out_df, aes(x = time/24, y = G)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = fasting_min, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = fasting_max, linetype = "dashed", color = "grey50") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = fasting_min, ymax = fasting_max,
           fill = "green", alpha = 0.05) +
  labs(
    title = "Glucose Drift Check",
    subtitle = "Green band = physiological normal range (70–180 mg/dL)",
    x = "Time (days)", y = "Glucose (mg/dL)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey50", size = 11),
    legend.position = "bottom"
  )

p_drift

# -----------------------------
# Phase Plane: G vs P
phase_plane <- ggplot(out_df, aes(x = P, y = G, color = site)) +
  geom_path(linewidth = 0.8, alpha = 0.8) +
  labs(title = "Phase Plane: Glucose vs Plasma Insulin",
       x = "Plasma Insulin (P)", y = "Glucose (G, mg/dL)") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")
phase_plane

# -----------------------------
# Cmax / Tmax / AUC per site
insulin_metrics <- out_df %>%
  group_by(site) %>%
  summarise(
    Cmax = max(P),
    Tmax = time[which.max(P)],
    AUC  = sum(P) * (diff(time)[1]) # approximate integral via sum*dt
  ) %>%
  pivot_longer(cols = c(Cmax, Tmax, AUC), names_to = "Metric", values_to = "Value")

cmax_plot <- ggplot(insulin_metrics %>% filter(Metric=="Cmax"), 
                    aes(x = site, y = Value, fill = site)) +
  geom_col() + labs(title="Cmax by Site", y="Peak Plasma Insulin") + theme_minimal()

tmax_plot <- ggplot(insulin_metrics %>% filter(Metric=="Tmax"), 
                    aes(x = site, y = Value, fill = site)) +
  geom_col() + labs(title="Tmax by Site", y="Time to Peak (hours)") + theme_minimal()

auc_plot <- ggplot(insulin_metrics %>% filter(Metric=="AUC"), 
                   aes(x = site, y = Value, fill = site)) +
  geom_col() + labs(title="AUC by Site", y="Total Exposure") + theme_minimal()

metrics_plot <- cmax_plot / tmax_plot / auc_plot
metrics_plot

# -----------------------------
# Glucose Time-in-Range Heatmap
# Define TIR (70–180 mg/dL)
tir_summary <- out_df %>%
  mutate(in_range = G >= 70 & G <= 180) %>%
  group_by(site) %>%
  summarise(TIR = mean(in_range)*100) # percent of time in range

tir_plot <- ggplot(tir_summary, aes(x = site, y = 1, fill = TIR)) +
  geom_tile() +
  scale_fill_gradient(low="red", high="green") +
  labs(title="Glucose Time-in-Range (%) by Site", x="Injection Site", y=NULL, fill="TIR %") +
  theme_minimal() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

tir_plot

# -----------------------------
# Histogram of Glucose Distribution

glucose_hist <- ggplot(out_df, aes(x = G, fill = site)) +
  geom_histogram(alpha = 0.6, position="identity", bins = 50) +
  labs(title="Distribution of Glucose Levels", x="G (mg/dL)", y="Count") +
  theme_minimal()

glucose_hist

# -----------------------------
# Depot vs Plasma Insulin Scatter / kabs overlay
depot_scatter <- ggplot(out_df, aes(x = S, y = P, color = site)) +
  geom_path(size = 0.7, alpha=0.8) +
  geom_line(aes(y = kabs), linetype="dashed") +
  labs(title="Subcutaneous Depot vs Plasma Insulin", x="S (Depot Units)", y="P / kabs") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")
depot_scatter


