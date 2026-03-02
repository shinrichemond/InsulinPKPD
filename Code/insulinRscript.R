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
  I_basal = 0.5,# Insulin-Glucose feedback term
  
  Vc = 5.0,     # Volume of central plasma
  
  Emax = 3.0,   # maximum endogenous glucose production
  EC50 = 90,    # half saturation term for glucose production
  hillN = 2.5,  # saturation power term
  
  IC50 = 25,    # half saturation term for glucose secretion inhibition
  
  Gb = 0.3,     # Basal glucose production
  Umax = 9,     # Maximal glucose uptake
  KP = 10,      # Half saturation for insulin in utilization
  KG = 80       # Half saturation for glucose in utilization
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
meal_sizes  <- c(25, 35, 30)       # glucose load scale
meal_durs   <- c(1.0, 1.0, 1.0)    # absorption windows

dose_times  <- c(7.75, 12.75, 18.75)  # pre-meal injections
dose_sizes  <- c(7, 8, 7)
inj_dur <- 0.05   # 3 min

# Insulin Sensitivity functions
sens_linear <- function(t_day, wake=6, sleep=22, m=0.02){ # sensitivity decline
  if(t_day < wake) return(1 - m*(sleep-wake))  # night baseline
  if(t_day > sleep) return(1 - m*(sleep-wake))
  
  x <- t_day - wake
  y <- 1 - m*x   # y = -m x + 1
  max(y, 0.5)    # physiological floor (avoid zero sensitivity)
}
sens_dawn <- function(t_day, mu=6, sigma=1.5, depth=0.35){ # dawn phenomenon
  1 - depth * exp(-((t_day - mu)^2)/(2*sigma^2))
}

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
    kabs <- Vmax * S_pos / (Km + S_pos)                           # Insulin Absorption
    E_G   <- Emax * G_pos^hillN   / (EC50^hillN + G_pos^hillN )   # Endogenous glucose production
    C_t <- 1 / (1 + P_pos / IC50)                                 # glucose production inhibitory term
    sens <- sens_linear(t_day) * sens_dawn(t_day)                 # time-dependent insulin sensitivity
    U_PG <- Umax * P_pos * G_pos / ((KP + P_pos) * (KG + G_pos))  # glucose utilization 
    
    # Main ODEs
    dS <- u_input - kabs - kdeg * S_pos + kre * P_pos
    dP <- Vc * (kabs + E_G * C_t + I_basal) - kclr * P_pos - kenz * P_pos - kre * P_pos
    dG <- Gb + G_in - sens * U_PG
    
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
end_day <- 14
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

# Graphs
ggplot(out_df, aes(x = time/24, y = G)) +
  geom_line() +
  theme_minimal()

ggplot(out_df, aes(x = time/24, y = G)) +
  geom_line() +
  geom_hline(yintercept = c(70, 100), linetype="dashed") +
  labs(title="Fasting glucose range (70–100 mg/dL)") +
  theme_minimal()

ggplot(out_df, aes(x = time/24, y = G)) +
  geom_line() +
  geom_hline(yintercept = c(120,160), linetype="dotted") +
  labs(title="Post-prandial glucose excursion window") +
  theme_minimal()

ggplot(out_df, aes(x = time/24, y = P)) +
  geom_line() +
  labs(title="Plasma insulin concentration dynamics") +
  theme_minimal()

ggplot(out_df, aes(x = time/24, y = S)) +
  geom_line() +
  labs(title="Subcutaneous insulin depot") +
  theme_minimal()
sens_df <- data.frame(
  time = out_df$time,
  t_day = out_df$time %% 24,
  sens = sapply(out_df$time %% 24, function(td){
    sens_linear(td) * sens_dawn(td)
  })
)

sens_df <- data.frame(
  time = out_df$time,
  t_day = out_df$time %% 24,
  sens = sapply(out_df$time %% 24, function(td){
    sens_linear(td) * sens_dawn(td)
  })
)
ggplot(sens_df, aes(x = time/24, y = sens)) +
  geom_line() +
  labs(title="Time-varying insulin sensitivity") +
  theme_minimal()

out_df$U_eff <- with(out_df,
                     sapply(time %% 24, function(td){
                       sens_linear(td) * sens_dawn(td)
                     }) * base_params$Umax * P * G / ((base_params$KP + P)*(base_params$KG + G))
)
ggplot(out_df, aes(x = time/24, y = U_eff)) +
  geom_line() +
  labs(title="Effective glucose utilization") +
  theme_minimal()

out_df$E_prod <- with(out_df, base_params$Emax * G^base_params$hillN / (base_params$EC50^base_params$hillN + G^base_params$hillN))
ggplot(out_df, aes(x = time/24, y = E_prod)) +
  geom_line() +
  labs(title="Endogenous glucose production") +
  theme_minimal()

out_df$net_G <- with(out_df, base_params$Gb + G_in - U_eff)
ggplot(out_df, aes(x = time/24, y = net_G)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype="dashed") +
  labs(title="Glucose mass balance") +
  theme_minimal()

ggplot(out_df, aes(x = P, y = G)) +
  geom_path(alpha=0.5) +
  labs(title="Phase plane: Plasma insulin vs Glucose") +
  theme_minimal()

ggplot(out_df, aes(x = time/24, y = G)) +
  geom_line() +
  geom_hline(yintercept = c(70,140), linetype="dashed") +
  labs(title="Glucose time-in-range (70–140)") +
  theme_minimal()
