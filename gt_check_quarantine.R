# this is for checking the quarantine
# obtain the parameters from gt_ver 2
source("gt_vaccination_plot.R")
### then is the quarantine
# assume we capture and quarantine(securely) x % Exposed individuals
###### SEIRDQ
seirdq_user_ode <- odin::odin({
  deriv(S) <- -beta * S * (I_R + I_D) / N 
  deriv(E) <-  beta * S * (I_R + I_D) / N  - lambda * E * (1-quarantine_eff) - lambda * quarantine_eff * E
  deriv(I_R) <-  lambda * E *  (1-quarantine_eff) * proportion_recovery - recovery_obtained * I_R 
  deriv(I_D) <-  lambda * E *  (1-quarantine_eff) * proportion_death - death_obtained * I_D 
  deriv(R) <-  recovery_obtained * I_R + recovery_obtained * Q_R
  deriv(D) <- death_obtained * I_D + death_obtained * Q_D
  deriv(Q_R) <- lambda * quarantine_eff * E * proportion_recovery - recovery_obtained * Q_R
  deriv(Q_D) <- lambda * quarantine_eff * E * proportion_death - death_obtained * Q_D
  
  deriv(C) <- beta * S * (I_R + I_D) / N      # cumulative infections
  initial(C) <- 0
  
  N <-  S + E + I_R + I_D + R #+ D + Q_R + Q_D
  
  initial(S) <- N_ini - I_ini
  initial(E) <- 0
  initial(I_R) <- proportion_recovery * I_ini 
  initial(I_D) <- proportion_death    * I_ini
  initial(R) <- 0
  initial(D) <- 0
  initial(Q_R) <- 0
  initial(Q_D) <- 0
  
  proportion_recovery <- 0.436
  proportion_death <- 0.564
  
  # user parameters
  N_ini <- user(6e5)
  I_ini <- user(1)
  lambda <- user()
  beta <- user()
  quarantine <- user(0.6)
  recovery_obtained <- user()
  death_obtained <- user()
  
  ## vaccination schedule & number (constant number per day)
  quarantine_start     <- user(30)   # day intervention begins
  
  quarantine_eff <- if (t < quarantine_start) 0  else quarantine
  
  
}, verbose = FALSE)


# initialize model
seirdq_mod <- seirdq_user_ode$new(beta = beta0, lambda = lambda_obtained, recovery_obtained = recovery_obtained,death_obtained = death_obtained,
                                  quarantine = 0.4, quarantine_start = 60)
# run the model
seirdq_mod <- seirdq_mod$run(times)

result_df <- as.data.frame(seirdq_mod)
result_df[450,]$D
max(result_df$I_R + result_df$I_D)
which((result_df$I_R + result_df$I_D) == max(result_df$I_R + result_df$I_D))



########### if we want to stop the quarantine after 30 days with daily cases < 0.3


seirdq_user_ode <- odin::odin({
  ## -------- Core dynamics --------
  deriv(S)   <- -beta * S * (I_R + I_D) / N
  deriv(E)   <-  beta * S * (I_R + I_D) / N  - lambda * E * (1 - quarantine_eff) - lambda * quarantine_eff * E
  deriv(I_R) <-  lambda * E * (1 - quarantine_eff) * proportion_recovery - recovery_obtained * I_R
  deriv(I_D) <-  lambda * E * (1 - quarantine_eff) * proportion_death    - death_obtained    * I_D
  deriv(R)   <-  recovery_obtained * I_R + recovery_obtained * Q_R
  deriv(D)   <-  death_obtained    * I_D + death_obtained    * Q_D
  deriv(Q_R) <-  lambda * quarantine_eff * E * proportion_recovery - recovery_obtained * Q_R
  deriv(Q_D) <-  lambda * quarantine_eff * E * proportion_death    - death_obtained    * Q_D
  
  ## Cumulative infections
  deriv(C) <- beta * S * (I_R + I_D) / N
  initial(C) <- 0
  
  ## Population
  N <- S + E + I_R + I_D + R # + Q_R + Q_D 
  
  ## -------- Initial conditions --------
  initial(S)   <- N_ini - I_ini
  initial(E)   <- 0
  initial(I_R) <- proportion_recovery * I_ini
  initial(I_D) <- proportion_death    * I_ini
  initial(R)   <- 0
  initial(D)   <- 0
  initial(Q_R) <- 0
  initial(Q_D) <- 0
  
  proportion_recovery <- 0.436
  proportion_death    <- 0.564
  
  ## -------- User parameters --------
  N_ini <- user(6e5)
  I_ini <- user(1)
  lambda <- user()
  beta <- user()
  recovery_obtained <- user()
  death_obtained <- user()
  
  quarantine       <- user(0.6)   # strength when ON
  quarantine_start <- user(30)
  
  ## Auto-stop rule
  case_thresh <- user(0.3)
  window_days <- user(30)
  
  ## -------- Indicators implemented as algebra --------
  ## (lambda*E < case_thresh) → 1 if true, 0 if false
  below_cond <- 1.0 * (lambda * E < case_thresh)
  
  ## (t >= quarantine_start) → 1 if true, 0 if false
  after_start <- 1.0 * (t >= quarantine_start)
  
  ## Combined condition
  cond_all <- after_start * below_cond * (1 - q_off)
  
  ## -------- Counter for consecutive days --------
  ## Integrates +1 when condition true, resets to 0 when false
  deriv(below_days) <- cond_all - (1 - cond_all) * below_days
  initial(below_days) <- 0
  
  ## -------- Latch to turn quarantine off --------
  off_trigger <- 1.0 * (below_days >= window_days)
  deriv(q_off) <- off_trigger * (1 - q_off)
  initial(q_off) <- 0
  
  ## -------- Effective quarantine --------
  quarantine_eff <- quarantine * after_start * (1 - q_off)
  
}, verbose = FALSE)


# initialize model
seirdq_mod <- seirdq_user_ode$new(beta = beta0, lambda = lambda_obtained, recovery_obtained = recovery_obtained,death_obtained = death_obtained,
                                  quarantine = 0.2, quarantine_start = 30)
# run the model
seirdq_mod <- seirdq_mod$run(times)

result_df <- as.data.frame(seirdq_mod)
result_df[600,]$D
max(result_df$I_R + result_df$I_D)
which((result_df$I_R + result_df$I_D) == max(result_df$I_R + result_df$I_D))

# while compared with previous one this one basically yields very similar result
