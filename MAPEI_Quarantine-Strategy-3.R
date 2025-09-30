library(odin)
library(tidyverse)
library(plotly)

# 1. Define model

seird_quarantine <- odin::odin({
  
  # 1.1 Compartment dynamics
  deriv(S)  <- -beta * S * I / N
  deriv(E)  <-  beta * S * I / N - lambda * E
  
  #Making q, the probability of catching a patient before turns infectious
  #At the same time, modelling in the day that quarantine starts
  e_Q <- (t >= Q_day) * e_Quara
  deriv(QR) <- e_Q * (1 - p_D) * lambda * E - gamma_R * QR   # Recovery path (Qua)
  deriv(QD) <- e_Q * p_D * lambda * E       - gamma_D * QD   # Death path (Qua)
  deriv(IR) <- (1-e_Q) * (1 - p_D) * lambda * E - gamma_R * IR   # Recovery path
  deriv(ID) <- (1-e_Q) * p_D * lambda * E       - gamma_D * ID   # Death path
  deriv(R)  <- gamma_R * IR
  deriv(D)  <- gamma_D * ID
  
  # 1.2 Aggregated totals 
  I <- IR + ID                        # Total infectious
  Q <- QR + QD                        # Total Quarantine
  N <- S + E + IR + ID + R            # Living population (excluding D)
  
  # 1.3 Initial conditions 
  initial(S)  <- N_ini - I_ini
  initial(E)  <- 0
  initial(IR) <- I_ini * (1 - p_D)
  initial(ID) <- I_ini * p_D
  initial(QR) <- 0
  initial(QD) <- 0
  initial(R)  <- 0
  initial(D)  <- 0
  
  # 1.4 User-defined parameters
  N_ini   <- user(600000)
  I_ini   <- user(1)
  beta    <- user(0.2323)           # Transmission rate
  lambda  <- user(1 / 9.399205)     # Incubation rate
  gamma_R <- user(1 / 16.09345)     # Recovery rate
  gamma_D <- user(1 / 7.49748)      # Death rate
  p_D     <- user(0.564)            # Case fatality probability
  Q_day <- user(9999999)            # Default not starting
  e_Quara <- user(0)                # Default quarantine completely ineff
  
  # 1.5 Derived quantities
  T_inf_eff <- (1 - p_D) / gamma_R + p_D / gamma_D
  gamma_eff <- 1 / T_inf_eff
  R0        <- user(2.2708)
})


# 2. Scenario Creation
# a. For simple 3x3
scenarios <- expand.grid(
  Q_day = c(30, 60, 90),
  e_Quara = c(0.2, 0.4, 0.6)
)
scenarios <- bind_rows(
  data.frame(Q_day = 600, e_Quara = 0), # Intervention starts after simulation ends
  scenarios
)

# b. For fancy 3D plot
scenarios <- expand.grid(
  Q_day = 1:300,
  e_Quara = seq(0.0, 0.8, by = 0.01)
)

# 3. Create function for running each scenario 600 days
run_scenario <- function(day, q_val) {
  # 3.1 initialize mode
  seird_mod <- seird_quarantine$new(     p_D     = 0.564,
                                         N_ini   = 600000,
                                         I_ini   = 1,
                                         beta    = 0.2323,
                                         lambda  = 1 / 9.399205,
                                         gamma_R = 1 / 16.09345,
                                         gamma_D = 1 / 7.49748, 
                                         Q_day = day, 
                                         e_Quara = q_val)
  
  # 3.2 Duration to run
  times <- seq(1,600)
  
  # 3.3 Run the model
  seird_run <- seird_mod$run(times)
  
  # 3.4 Convert results to a data frame
  results_df <- as.data.frame(seird_run)
  
  # 3.5 Extract key metrics
  final_deaths <- tail(results_df$D, 1)
  peak_infections <- max(results_df$IR+results_df$ID)
  peak_quarantine <- max(results_df$QR+results_df$QD)
  peak_quarantine_day <-results_df$t[which.max(results_df$QR+results_df$QD)]
  peak_infections_day <- results_df$t[which.max(results_df$IR+results_df$ID)]
  
  # 3.6 Return a single row data frame with the results
  return(
    data.frame(
      quarantine_day = day,
      q_effectiveness = q_val,
      peak_infections = peak_infections,
      peak_infections_day = peak_infections_day,
      peak_quarantine = peak_quarantine,
      peak_quarantine_day = peak_quarantine_day,
      final_deaths_at_600_days = final_deaths
    )
  )
}

# 4. Run scenario to get results
all_results <- purrr::map2_dfr(
  scenarios$Q_day,
  scenarios$e_Quara,
  run_scenario
)


# 5. Plot compartments
plot_ly(data = all_results, 
        x = ~quarantine_day, 
        y = ~q_effectiveness, 
        z = ~final_deaths_at_600_days,
        color = ~peak_infections,
        type = 'scatter3d', 
        mode = 'markers',
        # Set the color bar title here
        marker = list(colorbar = list(title = "Peak Quarantine"))) %>%
  layout(title = "Final death count by varying start and effectiveness",
         # For 3D plots, axis titles go inside 'scene'
         scene = list(
           xaxis = list(title = "Quarantine Start"),
           yaxis = list(title = "Quarantine Effectiveness"),
           zaxis = list(title = "Final Deaths")
         ))
