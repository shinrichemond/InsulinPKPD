# import libraries
library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(plotly)
library(knitr)
library(kableExtra)

# "Base" Parameters
base_params <- list(
  Vmax = 18,    # Parameter for Insulin Absorption; Will be overridden per site
  Km   = 30,    # Parameter for Insulin Absorption; Will be overridden per site
  
  kdeg = 0.05,  # Insulin SC degradation
  kre = 0.01,   # Insulin SC absorption from central Plasma
  
  kclr = 0.15,  # Insulin Plasma clearence as waste from kidney, liver, etc...
  kenz = 0.08,  # Insulin Plasma enzymatic degredation
  
  Vc = 5.0,     # Volume of central plasma
  
  Emax = 3.0,   # maximum endogenous glucose production
  EC50 = 120,   # half saturation term for glucose production
  hillN = 2.5,  # saturation power term
  
  IC50 = 25,    # half saturation term for glucose secretion inhibition
  
  Gb = 1.8,     # Basal glucose production
  Umax = 5.5,   # Maximal glucose uptake
  KP = 15,      # Half saturation for insulin in utilization
  KG = 80       # Hald saturation for glucose in utilization
)

# Site Profiles
site_profiles <- data.frame(
  site  = c("Abdomen", "Upper Arm", "Buttock", "Thigh"),
  Vmax  = c(18, 13, 10, 8),
  Km    = c(30, 40, 42, 50),
  color = c("#c53030", "#d69e2e", "#2b6cb0", "#718096"),
  stringsAsFactors = FALSE
)

# Exogenous Input parameters (Glucose + Insulin)
meal_times  <- c(8, 13, 19)        # breakfast, lunch, dinner (hours in day)
meal_sizes  <- c(50, 40, 55)       # glucose load scale
meal_durs   <- c(0.75, 0.75, 1.0)  # absorption windows

dose_times  <- c(7.75, 12.75, 18.75)  # pre-meal injections
dose_sizes  <- c(10, 8, 10)
inj_dur <- 0.05   # 3 min

# Main System ODE
insulin_ode <- function(t, state, params) {
  with(as.list(c(state, params)), {
    
    # Set Time (mod 24)
    t_day <- t %% 24
    
    # Injection Site Input
    u_input <- 0
    for (i in seq_along(dose_times)) {
      t_dose <- dose_times[i]
      dose   <- dose_sizes[i]
      if (t_day >= t_dose && t_day < t_dose + inj_dur) {
        u_input <- u_input + dose / inj_dur
      }
    }
    
    # Meal Glucose Input
    G_in <- 0
    for (i in seq_along(meal_times)) {
      t_meal  <- meal_times[i]
      meal    <- meal_sizes[i]
      dur     <- meal_durs[i]
      t_rel   <- t_day - t_meal
      
      if (t_rel >= 0 && t_rel < dur) {
        t_peak  <- dur * 0.3
        t_decay <- dur * 0.4
        G_in <- G_in + meal * (t_rel / t_peak) * exp(-(t_rel - t_peak) / t_decay)
      }
    }
    
    # Bounding SVs to 0
    S_pos <- max(S, 0)
    P_pos <- max(P, 0)
    G_pos <- max(G, 0)
    
    # Auxiliary Functions
    kabs <- Vmax * S_pos / (Km + S_pos)                                # Insulin Absorption
    E_G   <- Emax * G_pos^hillN   / (EC50^hillN + G_pos^hillN )        # Endogenous glucose production
    C_t <- 1 / (1 + P_pos / IC50)                                      # glucose production inhibitory term
    U_PG <- Umax * P_pos * G_pos / ((KP + P_pos) * (KG + G_pos))       # glucose utilization 
    
    # Main ODEs
    dS <- u_input - kabs - kdeg * S_pos + kre * P_pos
    dP <- Vc * (kabs + E_G * C_t) - kclr * P_pos - kenz * P_pos - kre * P_pos
    dG <- Gb + G_in - U_PG
    
    # Return
    return(list(
      c(dS, dP, dG),
      kabs    = kabs,
      E_G     = E_G,
      C_t     = C_t,
      U_PG    = U_PG,
      G_in    = G_in,
      u_input = u_input
    ))
  })
}

# Initial Condition
state0 <- c(
  S = 0,   # subcutaneous insulin depot
  P = 0,   # plasma insulin
  G = 90   # baseline glucose
)

# Time God
day_length <- 24
end_day <- 7
t_end <- day_length * end_day     # run for however many days
dt    <- 1/60                     # set time resolution
times <- seq(0, t_end, by = dt)

# Set up parameters
params <- c(
  base_params,
  list(day_length = day_length)
)

out <- ode(
  y     = state0,
  times = times,
  func  = insulin_ode,
  parms = params,
  method = "lsoda"
)

out_df <- as.data.frame(out)
