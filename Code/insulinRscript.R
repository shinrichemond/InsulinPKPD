# import libraries
library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(plotly)
library(knitr)
library(kableExtra)
library(crayon)

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
  hillN = 2.5,    # saturation power term
  
  IC50 = 10,    # half saturation term for glucose secretion inhibition
  
  Gb = 0.5,     # Basal glucose production
  k0 = 0.01,    # Basal glucose uptake
  Umax = 6,     # Maximal glucose uptake
  KP = 15,      # Half saturation for insulin in utilization
  KG = 40,       # Half saturation for glucose in utilization
  
  # Circadian rhythm stuff
  m = 0.02,      # slope for linear decrease
  wake=6,         # wake for linear decrease
  depth=0.25,     # stuff for dawn effect
  mu=6,
  sigma=2,
  
  # Hepatic buffering
  G_set = 102,      # setpoint
  kH = .2      # buffering strength
)

# Exogenous Input parameters (Glucose + Insulin)
meal_times  <- c(8, 13, 19)        # breakfast, lunch, dinner (hours in day)
meal_sizes  <- c(35, 50, 50)       # glucose load scale
meal_durs   <- c(1.0, 1.0, 1.0)    # absorption windows

dose_times  <- c(7.75, 12.75, 18.75)  # pre-meal injections
dose_sizes  <- c(7, 8, 7)
inj_dur <- 0.05   # 3 min

# Circadian rhythm function
sens_func <- function(t_day,
                      wake=6, sleep=22,
                      m=0.015,
                      mu=6, sigma=1.5, depth=0.25,
                      floor=0.6){
  
  # linear daily resistance accumulation
  if(t_day < wake){
    lin_loss <- m*(sleep-wake)
  } else if(t_day > sleep){
    lin_loss <- m*(sleep-wake)
  } else {
    lin_loss <- m*(t_day - wake)
  }
  
  # dawn phenomenon (resistance bump)
  dawn_loss <- depth * exp(-((t_day - mu)^2)/(2*sigma^2))
  
  # total sensitivity
  sens <- 1 - (lin_loss + dawn_loss)
  
  # physiological bounds
  sens <- max(sens, floor)
  sens <- min(sens, 1.1)
  
  return(sens)
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
      
      if (t_rel >= 0 && t_rel <= dur) {
        x <- t_rel / dur              # normalize to [0,1]
        
        # parabolic dome: 4x(1-x)
        shape <- 4 * x * (1 - x)
        
        G_in <- G_in + meal * shape
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
    sens <- sens_func(t_day)                                      # time-dependent insulin sensitivity
    U_PG <- Umax * P_pos * G_pos / ((KP + P_pos) * (KG + G_pos))  # glucose utilization 
    H <- kH * (G_pos - G_set)                                     # Hepatic Buffering

    # Main ODEs
    dS <- u_input - kabs - kdeg * S_pos + kre * P_pos
    dP <- Vc * (kabs + E_G * C_t + I_basal) - kclr * P_pos - kenz * P_pos - kre * P_pos
    dG <- Gb + G_in - sens*U_PG - k0*G_pos - H 
    
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