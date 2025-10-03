# check distancing 
# obtain the parameters from gt_ver 2
source("gt_main_parm_vaccination_plot.R")

seird_user_ode_distancing <- odin::odin({
  deriv(S) <- -beta_eff * S * (I_R + I_D) / N
  deriv(E) <-  beta_eff * S * (I_R + I_D) / N - lambda * E
  deriv(I_R) <-  lambda * E * proportion_recovery - recovery_obtained * I_R
  deriv(I_D) <-  lambda * E * proportion_death - death_obtained * I_D
  deriv(R) <-  recovery_obtained * I_R 
  deriv(D) <- death_obtained * I_D
  
  deriv(C) <- beta_eff * S * (I_R + I_D) / N    # cumulative infections
  initial(C) <- 0
  N <- S + E + I_R + I_D + R 
  
  initial(S) <- N_ini - I_ini
  initial(E) <- 0
  initial(I_R) <- proportion_recovery * I_ini 
  initial(I_D) <- proportion_death    * I_ini
  initial(R) <- 0
  initial(D) <- 0
  
  proportion_recovery <- 0.436
  proportion_death <- 0.564
  recovery_obtained <- user()
  death_obtained <- user()
  
  # user parameters
  N_ini <- user(6e5)
  I_ini <- user(1)
  lambda <- user()
  beta <- user()
  
  intervention_day <- user(90)
  # changed parameter
  beta_eff <- if (t < intervention_day) beta else 0.5 * beta
  
  
}, verbose = FALSE)

# initialize model
seird_mod <- seird_user_ode_distancing$new(beta = beta0, lambda = lambda_obtained, recovery_obtained = recovery_obtained,death_obtained = death_obtained )
# how long to run
times <- seq(1,1500*1)
# run the model
seird_mod <- seird_mod$run(times)

p_seird <- seird_mod %>%
  ggplot(aes(x = t)) +
  geom_line(aes(y = S, color = "Susceptible")) +
  geom_line(aes(y = E, color = "Exposed")) +
  geom_line(aes(y = I_R+I_D, color = "Infected")) +
  geom_line(aes(y = R, color = "Removed")) +
  geom_line(aes(y = D, color = "Death")) +
  xlab("Time (days)") +
  ylab("Number of Individuals") +
  scale_color_manual(
    values = c(
      "Susceptible" = "#0072B2",  # blue
      "Exposed"    = "grey", #
      "Infected" = "#D55E00",   # red-orange
      "Removed"    = "#009E73",    # green
      "Death"    = "black"    # black
    ),
    breaks = c("Susceptible", "Exposed", "Infected", "Removed", "Death")  # set order here
  ) +
  theme_classic(base_size = 16) +
  labs(color = "Compartment") +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.7)   # adjust as needed
  )

#p_seird


# 
result_df_d <- as.data.frame(seird_mod)
result_df_d$R[1500]
result_df_d$D[1500]
result_df_d$S[1500]
max(result_df_d$I_R + result_df_d$I_D)
which((result_df_d$I_R + result_df_d$I_D) == max(result_df_d$I_R + result_df_d$I_D))
