# SEIRD Model: Time-Varying Beta (Intervention at t = 90)
# Baseline vs Distancing Comparison

library(odin)
library(dplyr)
library(ggplot2)


# 1. ODIN Model Definition with Time-Varying Beta
seird_split_intervene <- odin::odin({
  
  deriv(S)  <- -beta * S * I / N
  deriv(E)  <-  beta * S * I / N - lambda * E
  deriv(IR) <- (1 - p_D) * lambda * E - gamma_R * IR
  deriv(ID) <- p_D * lambda * E       - gamma_D * ID
  deriv(R)  <- gamma_R * IR
  deriv(D)  <- gamma_D * ID
  
  I <- IR + ID
  N <- S + E + IR + ID + R + D
  
  initial(S)  <- N_ini - I_ini
  initial(E)  <- 0
  initial(IR) <- I_ini * (1 - p_D)
  initial(ID) <- I_ini * p_D
  initial(R)  <- 0
  initial(D)  <- 0
  
  beta <- if (t < t_intervene) beta0 else beta1
  
  N_ini       <- user()
  I_ini       <- user()
  beta0       <- user()
  beta1       <- user()
  t_intervene <- user()
  lambda      <- user()
  gamma_R     <- user()
  gamma_D     <- user()
  p_D         <- user()
})



# 2. Common Parameters
pars <- list(
  N_ini        = 600000,
  I_ini        = 1,
  beta0        = 0.2323,
  beta1        = 0.2323 * 0.5,
  t_intervene  = 90,
  lambda       = 1 / 9.399205,
  gamma_R      = 1 / 16.09345,
  gamma_D      = 1 / 7.49748,
  p_D          = 0.564
)
times <- seq(0, 1200, by = 1)


# 3. Baseline Model (No Intervention)
mod_base <- seird_split_intervene$new(
  N_ini = pars$N_ini,
  I_ini = pars$I_ini,
  beta0 = pars$beta0,
  beta1 = pars$beta0,    # No change in beta
  t_intervene = 9999,    # Never intervene
  lambda = pars$lambda,
  gamma_R = pars$gamma_R,
  gamma_D = pars$gamma_D,
  p_D = pars$p_D
)
df_base <- as.data.frame(mod_base$run(times)) %>%
  mutate(I = IR + ID)


# 4. Distancing Model (Intervention at t = 90)
mod_dist <- seird_split_intervene$new(
  N_ini = pars$N_ini,
  I_ini = pars$I_ini,
  beta0 = pars$beta0,
  beta1 = pars$beta1,
  t_intervene = pars$t_intervene,
  lambda = pars$lambda,
  gamma_R = pars$gamma_R,
  gamma_D = pars$gamma_D,
  p_D = pars$p_D
)
df_dist <- as.data.frame(mod_dist$run(times)) %>%
  mutate(I = IR + ID)


# 5. Plot Distancing Scenario
p_seird_dist <- df_dist %>%
  ggplot(aes(x = t)) +
  geom_line(aes(y = S, color = "S")) +
  geom_line(aes(y = E * 10, color = "E*10")) +
  geom_line(aes(y = I * 10, color = "I*10")) +
  geom_line(aes(y = R, color = "R")) +
  geom_line(aes(y = D, color = "D"), linetype = "dashed") +
  scale_colour_manual(
    values = c(
      "S" = "steelblue", "E*10" = "orange",
      "I*10" = "black", "R" = "green4", "D" = "gray30"
    )
  ) +
  labs(
    title = "SEIRD Compartments Over Time (With Distancing)",
    x = "Time (days)",
    y = "Number of Individuals",
    color = "Compartment"
  ) +
  theme_classic(base_size = 16) +
  theme(legend.position = c(0.9, 0.75))

print(p_seird_dist)
ggsave("seird_distancing_plot.png", p_seird_dist, width = 15, height = 8, dpi = 300)


# 6. Summary: Distancing Scenario
peak_I     <- max(df_dist$I)
day_peak_I <- df_dist$t[which.max(df_dist$I)]
peak_E     <- max(df_dist$E)
day_peak_E <- df_dist$t[which.max(df_dist$E)]
final_R    <- tail(df_dist$R, 1)
final_D    <- tail(df_dist$D, 1)
final_S    <- tail(df_dist$S, 1)

CFR    <- final_D / (final_R + final_D)
attack <- (pars$N_ini - final_S) / pars$N_ini

cat("=== SEIRD Distancing Summary ===\n")
cat("Peak infections (I):", peak_I, "on day", day_peak_I, "\n")
cat("Peak exposed (E)   :", peak_E, "on day", day_peak_E, "\n")
cat("Final recovered     :", final_R, "\n")
cat("Final deaths        :", final_D, "\n")
cat("Case Fatality Ratio :", CFR, "\n")
cat("Attack rate         :", attack, "\n\n")


# 7. Summary: Baseline Scenario
peak_I_base     <- max(df_base$I)
day_peak_I_base <- df_base$t[which.max(df_base$I)]
peak_E_base     <- max(df_base$E)
day_peak_E_base <- df_base$t[which.max(df_base$E)]
final_R_base    <- tail(df_base$R, 1)
final_D_base    <- tail(df_base$D, 1)
final_S_base    <- tail(df_base$S, 1)

CFR_base    <- final_D_base / (final_R_base + final_D_base)
attack_base <- (pars$N_ini - final_S_base) / pars$N_ini

cat("=== SEIRD Baseline Summary ===\n")
cat("Peak infections (I):", peak_I_base, "on day", day_peak_I_base, "\n")
cat("Peak exposed (E)   :", peak_E_base, "on day", day_peak_E_base, "\n")
cat("Final recovered     :", final_R_base, "\n")
cat("Final deaths        :", final_D_base, "\n")
cat("Case Fatality Ratio :", CFR_base, "\n")
cat("Attack rate         :", attack_base, "\n")
cat("Susceptible (End)   :", final_S_base, "\n")
cat("Susceptible (%)     :", 100 * final_S_base / pars$N_ini, "%\n\n")


# 8. Effectiveness Comparison
eff_attack <- 1 - (attack / attack_base)
eff_deaths <- 1 - (final_D / final_D_base)
eff_peakI  <- 1 - (peak_I / peak_I_base)
delay_peak_day <- day_peak_I - day_peak_I_base

cat("=== Intervention Effectiveness ===\n")
cat("↓ Attack Rate Reduction:", 100 * eff_attack, "%\n")
cat("↓ Deaths Reduction     :", 100 * eff_deaths, "%\n")
cat("↓ Peak Infections Drop :", 100 * eff_peakI, "%\n")
cat("→ Peak Day Delayed by  :", delay_peak_day, "days\n")
