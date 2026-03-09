library(deSolve)
library(dplyr)
library(ggplot2)
library(patchwork)

# get system ODE
this_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(this_dir)
source(file.path(this_dir, "insulinRscript.R"))

# finding peaks
find_peaks <- function(time, P) {
  # find local maxima
  idx <- which(diff(sign(diff(P))) == -2) + 1
  
  data.frame(
    t_peak = time[idx],
    P_peak = P[idx]
  )
}

# Time grid
times <- seq(0, 24, by = 0.01)

# Initial conditions
state_init <- c(
  S = 0,
  P = 10,
  G = base_params$G_set
)

# Parameter ranges
Vmax_vals <- seq(8, 30, length.out = 12)
Km_vals   <- seq(10, 80, length.out = 12)

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
    
    peaks <- find_peaks(df$time, df$P)
    
    # if no peaks detected (rare but possible)
    if(nrow(peaks) == 0){
      Cmax <- NA
      Tmax <- NA
      peak_id <- NA
    } else {
      Cmax <- max(peaks$P_peak)
      Tmax <- peaks$t_peak[which.max(peaks$P_peak)]
      
      # which peak is dominant
      peak_id <- which.max(peaks$P_peak)
    }
    
    # trapezoidal AUC
    AUC <- sum(diff(df$time) *
                 (head(df$P,-1) + tail(df$P,-1))/2)
    
    results <- rbind(results,
                     data.frame(
                       Vmax = v,
                       Km = k,
                       Cmax = Cmax,
                       Tmax = Tmax,
                       AUC = AUC,
                       peak_id = peak_id
                     )
    )
  }
}

ggplot(results, aes(Vmax, Km, fill = factor(peak_id))) +
  geom_tile(color="white") +
  labs(
    title = "Dominant Peak Index",
    fill = "Peak Number"
  ) +
  heatmap_theme

slice <- results %>% filter(Km == 10)

ggplot(slice, aes(Vmax, Tmax)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Tmax vs Vmax (Km fixed)"
  )

ggplot(df, aes(time, P)) +
  geom_line() +
  labs(title="Plasma insulin time series")
