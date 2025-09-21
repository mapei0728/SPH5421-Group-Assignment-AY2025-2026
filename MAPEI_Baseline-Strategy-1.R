# SEIRD Model with Split Infection Paths (Recovered vs Deaths)
# Population Size: 600,000
# Model Output: Baseline Compartment Model + Key Summary Statistics

library(odin)
library(tidyverse)


# 1. Define the SEIRD model with split IR (recover) / ID (death)

seird_split <- odin::odin({
  
  # 1.1 Compartment dynamics
  deriv(S)  <- -beta * S * I / N
  deriv(E)  <-  beta * S * I / N - lambda * E
  deriv(IR) <- (1 - p_D) * lambda * E - gamma_R * IR   # Recovery path
  deriv(ID) <- p_D * lambda * E       - gamma_D * ID   # Death path
  deriv(R)  <- gamma_R * IR
  deriv(D)  <- gamma_D * ID
  
  # 1.2 Aggregated totals 
  I <- IR + ID                        # Total infectious
  N <- S + E + IR + ID + R            # Living population (excluding D)
  
  # 1.3 Initial conditions 
  initial(S)  <- N_ini - I_ini
  initial(E)  <- 0
  initial(IR) <- I_ini * (1 - p_D)
  initial(ID) <- I_ini * p_D
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
  
  # 1.5 Derived quantities
  T_inf_eff <- (1 - p_D) / gamma_R + p_D / gamma_D
  gamma_eff <- 1 / T_inf_eff
  R0        <- user(2.2708)
})


# 2. Initialize model
seird_mod <- seird_split$new(
  p_D     = 0.564,
  N_ini   = 600000,
  I_ini   = 1,
  beta    = 0.2323,
  lambda  = 1 / 9.399205,
  gamma_R = 1 / 16.09345,
  gamma_D = 1 / 7.49748
)


# 3. Run simulation over 600 days
times <- seq(0, 600, by = 1)
seird_out <- seird_mod$run(times)

df_seird_out <- as.data.frame(seird_out) %>%
  mutate(I = IR + ID)   # Add total infectious for plotting


# 4. Plot compartments
p_seird <- df_seird_out %>%
  ggplot(aes(x = t)) +
  geom_line(aes(y = S, color = "S")) +
  geom_line(aes(y = E, color = "E")) +
  geom_line(aes(y = I, color = "I")) +
  geom_line(aes(y = R, color = "R")) +
  geom_line(aes(y = D, color = "D"), linetype = "dashed") +
  scale_colour_manual(values = c(
    "S" = "steelblue", "E" = "orange", "I" = "red",
    "R" = "green", "D" = "black"
  )) +
  labs(
    title = "SEIRD Model with Split Infectious Paths",
    x = "Time (days)",
    y = "Number of individuals",
    color = "Compartment"
  ) +
  theme_classic(base_size = 18) +
  theme(legend.position = c(0.9, 0.75))

print(p_seird)


# 5. Extract Key Summary Statistics
N0         <- with(df_seird_out[1,], S + E + IR + ID + R)
peak_I     <- max(df_seird_out$I)
day_peak_I <- df_seird_out$t[which.max(df_seird_out$I)][1]

peak_E     <- max(df_seird_out$E)
day_peak_E <- df_seird_out$t[which.max(df_seird_out$E)][1]

final_R <- tail(df_seird_out$R, 1)
final_D <- tail(df_seird_out$D, 1)
final_S <- tail(df_seird_out$S, 1)

attack_rate <- (final_R + final_D) / N0        # Infection attack rate
CFR         <- final_D / (final_R + final_D)   # Case fatality ratio


# 6. Print Summary
cat("=== SEIRD Summary ===\n")
cat("Initial population (alive)     :", N0, "\n")
cat("Peak Infections (I)            :", peak_I, "on day", day_peak_I, "\n")
cat("Peak Exposed (E)               :", peak_E, "on day", day_peak_E, "\n")
cat("Final Recovered                :", final_R, "\n")
cat("Final Deaths                   :", final_D, "\n")
cat("Attack Rate (R + D) / N0       :", attack_rate, "\n")
cat("Case Fatality Ratio (D / R+D)  :", CFR, "\n\n")


# 7. Mass Balance Check 
mass_check <- with(tail(df_seird_out, 1), S + E + IR + ID + R + D)
cat("Mass balance check (final total):", mass_check,
    " | Should ≈", N0 + final_D, "\n\n")


# 8. Print Final Compartment Values
final_vals <- tail(df_seird_out, 1)

cat("=== Final Compartment Values ===\n")
cat("S  (Susceptible)           :", final_vals$S,  "\n")
cat("E  (Exposed)               :", final_vals$E,  "\n")
cat("IR (Infectious - Recover) :", final_vals$IR, "\n")
cat("ID (Infectious - Death)   :", final_vals$ID, "\n")
cat("R  (Recovered)             :", final_vals$R,  "\n")
cat("D  (Deceased)              :", final_vals$D,  "\n")
