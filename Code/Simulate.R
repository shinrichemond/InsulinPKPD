library(deSolve)

this_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(this_dir)
source(file.path(this_dir, "insulinRscript.R"))

# Site Profiles
site_profiles <- data.frame(
  site  = c("Abdomen", "Upper Arm", "Buttock", "Thigh"),
  Vmax  = c(18, 10, 2, 5),
  Km    = c(30, 50, 60, 90),
  color = c("#c53030", "#d69e2e", "#2b6cb0", "#718096"),
  stringsAsFactors = FALSE
)

# Initial Condition
state0 <- c(
  S = 0,   # subcutaneous insulin depot
  P = 0,   # plasma insulin
  G = 90   # baseline glucose
)

# Time God
day_length <- 24
end_day <- 1
t_end <- day_length * end_day     # run for however many days
dt    <- 1/60                     # set time resolution
times <- seq(0, t_end, by = dt)

# Set up parameters
params <- c(
  base_params,
  list(day_length = day_length)
)

# Storage
results <- lapply(1:nrow(site_profiles), function(i){
  
  site <- site_profiles$site[i]
  params_site <- modifyList(
    base_params,
    list(
      Vmax = site_profiles$Vmax[i],
      Km   = site_profiles$Km[i],
      day_length = day_length
    )
  )
  
  out <- ode(
    y = state0,
    times = times,
    func = insulin_ode,
    parms = params_site,
    method = "lsoda"
  )
  
  df <- as.data.frame(out)
  df$site <- site
  df
})

# Combined dataframe (all sites)
all_sites_df <- dplyr::bind_rows(results)

# Save result
saveRDS(all_sites_df, file.path(this_dir, "data.rds"))