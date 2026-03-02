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