# SEIRD Model with Vaccination Strategies
# Baseline vs 2.5% Daily Vaccination vs Fixed 3000/day

library(odin)
library(tidyverse)
library(patchwork)


# 1. ODIN Model Definition: SEIRD + Vaccination 
seird_vaccine <- odin::odin({
  
  ## SEIRD + Vaccination Compartments
  deriv(S)  <- -beta * S * I / N - vax_effective
  deriv(E)  <-  beta * S * I / N - lambda * E
  deriv(IR) <- (1 - p_D) * lambda * E - gamma_R * IR
  deriv(ID) <- p_D * lambda * E - gamma_D * ID
  deriv(R)  <- gamma_R * IR                     # Recovered from infection
  deriv(D)  <- gamma_D * ID
  deriv(V)  <- vax_effective                    # Immune via vaccination
  
  I <- IR + ID
  N <- S + E + IR + ID + R + V
  
  ## Vaccination Strategy
  vax_1 <- if (t >= t_vax_start && vax_type == 1) vax_rate_prop * S else 0
  vax_2 <- if (t >= t_vax_start && vax_type == 2) vax_rate_fixed     else 0
  vax_admin  <- vax_1 + vax_2
  vax_demand <- if (vax_admin < S) vax_admin else S
  vax_effective <- vax_demand * vax_eff
  
  output(vax_daily_admin)     <- vax_demand
  output(vax_daily_effective) <- vax_effective
  
  ## Initial Conditions
  initial(S) <- N_ini - I_ini
  initial(E) <- 0
  initial(IR) <- I_ini * (1 - p_D)
  initial(ID) <- I_ini * p_D
  initial(R)  <- 0
  initial(D)  <- 0
  initial(V)  <- 0
  
  ## User Parameters
  N_ini          <- user()
  I_ini          <- user()
  beta           <- user()
  lambda         <- user()
  gamma_R        <- user()
  gamma_D        <- user()
  p_D            <- user()
  vax_type       <- user()
  vax_rate_prop  <- user()
  vax_rate_fixed <- user()
  vax_eff        <- user()
  t_vax_start    <- user()
})


# 2. Parameter Setup & Model Execution
pars <- list(
  N_ini          = 600000,
  I_ini          = 1,
  beta           = 0.2323,
  lambda         = 1 / 9.399205,
  gamma_R        = 1 / 16.09345,
  gamma_D        = 1 / 7.49748,
  p_D            = 0.564,
  vax_eff        = 0.84,
  vax_rate_prop  = 0.025,
  vax_rate_fixed = 3000,
  t_vax_start    = 90
)
times <- seq(0, 600, 1)

# Function to run a scenario
run_model <- function(vax_type) {
  mod <- seird_vaccine$new(
    N_ini=pars$N_ini, I_ini=pars$I_ini,
    beta=pars$beta, lambda=pars$lambda,
    gamma_R=pars$gamma_R, gamma_D=pars$gamma_D,
    p_D=pars$p_D, vax_type=vax_type,
    vax_eff=pars$vax_eff,
    vax_rate_prop=pars$vax_rate_prop,
    vax_rate_fixed=pars$vax_rate_fixed,
    t_vax_start=pars$t_vax_start
  )
  as.data.frame(mod$run(times)) %>%
    mutate(
      I = IR + ID,
      scenario = c("Baseline", "Vax_2.5%", "Vax_3000")[vax_type + 1]
    ) %>%
    rename(R_recov = R)
}

# Run three scenarios
df_baseline <- run_model(0)
df_vax25    <- run_model(1)
df_vax3000  <- run_model(2)
df_all      <- bind_rows(df_baseline, df_vax25, df_vax3000)


# 3. Scenario Summary Tables
# Parameter table
param_table <- tibble(
  Scenario = c("Baseline", "Vax_2.5%", "Vax_3000"),
  N_ini    = pars$N_ini,
  I_ini    = pars$I_ini,
  beta     = pars$beta,
  lambda   = pars$lambda,
  gamma_R  = pars$gamma_R,
  gamma_D  = pars$gamma_D,
  p_D      = pars$p_D,
  vax_policy = c(
    "None",
    sprintf("Prop: %.3f×S/day; eff=%.2f; start=%d", pars$vax_rate_prop,  pars$vax_eff, pars$t_vax_start),
    sprintf("Fixed: %d/day;    eff=%.2f; start=%d", pars$vax_rate_fixed, pars$vax_eff, pars$t_vax_start)
  )
)

# Output metrics
summarise_metrics <- function(df) {
  peak_idx <- which.max(df$I)
  tibble(
    Scenario        = unique(df$scenario),
    Recovered_cases = tail(df$R_recov, 1),
    Deaths_cases    = tail(df$D, 1),
    Peak_Infections = df$I[peak_idx],
    Peak_Day        = df$t[peak_idx]
  )
}
results_table <- bind_rows(
  summarise_metrics(df_baseline),
  summarise_metrics(df_vax25),
  summarise_metrics(df_vax3000)
) %>%
  mutate(across(where(is.numeric), round, 0))

cat("Model Parameters by Scenario \n")
print(param_table)

cat("\n Baseline Outcomes by Scenario \n")
print(results_table)


# 4. Intervention Effectiveness Comparison
metrics <- function(df, N_ini) {
  final_Rrecov <- tail(df$R_recov, 1)
  final_D      <- tail(df$D, 1)
  list(
    AR      = (final_Rrecov + final_D) / N_ini,
    D       = final_D,
    peakI   = max(df$I),
    peakDay = df$t[which.max(df$I)]
  )
}
m_base   <- metrics(df_baseline, pars$N_ini)
m_v25    <- metrics(df_vax25,    pars$N_ini)
m_v3000  <- metrics(df_vax3000,  pars$N_ini)

eff_rel <- function(m, m_base) {
  list(
    eff_attack = 1 - m$AR     / m_base$AR,
    eff_deaths = 1 - m$D      / m_base$D,
    eff_peakI  = 1 - m$peakI  / m_base$peakI,
    delay_peak =    m$peakDay - m_base$peakDay
  )
}
e25    <- eff_rel(m_v25,   m_base)
e3000  <- eff_rel(m_v3000, m_base)

cat("\n Effectiveness: Vax 2.5% \n")
cat(" Attack Rate:", sprintf("%.2f%%", 100 * e25$eff_attack), "\n")
cat(" Deaths     :", sprintf("%.2f%%", 100 * e25$eff_deaths), "\n")
cat(" Peak I     :", sprintf("%.2f%%", 100 * e25$eff_peakI), "\n")
cat(" Delay Peak :", e25$delay_peak, "days\n")

cat("\n Effectiveness: Vax 3000/day \n")
cat(" Attack Rate:", sprintf("%.2f%%", 100 * e3000$eff_attack), "\n")
cat(" Deaths     :", sprintf("%.2f%%", 100 * e3000$eff_deaths), "\n")
cat(" Peak I     :", sprintf("%.2f%%", 100 * e3000$eff_peakI), "\n")
cat(" Delay Peak :", e3000$delay_peak, "days\n")


# 5. Infection Curves Over Time
ggplot(df_all, aes(x = t, y = I, color = scenario)) +
  geom_line(size = 1) +
  labs(
    title = "Infected Over Time (I)",
    x = "Time (days)", y = "Number of Infectious",
    color = "Scenario"
  ) +
  theme_minimal(base_size = 14)


# 6. Detailed Compartment Plots for Each Scenario
plot_seird_compartments <- function(df, title_text) {
  legend_levels <- c(
    "Susceptible", "Exposed", "Infectious",
    "Recovered (infection)", "Vaccinated (immune)", "Deaths"
  )
  legend_colors <- c(
    "Susceptible"           = "steelblue",
    "Exposed"              = "orange",
    "Infectious"            = "red",
    "Recovered (infection)" = "green4",
    "Vaccinated (immune)"   = "purple",
    "Deaths"                = "black"
  )
  ggplot(df, aes(x = t)) +
    geom_line(aes(y = S,       color = "Susceptible")) +
    geom_line(aes(y = E,       color = "Exposed")) +
    geom_line(aes(y = I,       color = "Infectious")) +
    geom_line(aes(y = R_recov*10, color = "Recovered (infection)")) +
    geom_line(aes(y = V,       color = "Vaccinated (immune)")) +
    geom_line(aes(y = D,       color = "Deaths"), linetype = "dashed") +
    scale_colour_manual(
      values = legend_colors, breaks = legend_levels, limits = legend_levels,
      name = "Compartment"
    ) +
    labs(title = title_text, x = "Time (days)", y = "Individuals") +
    theme_classic(base_size = 18) +
    theme(legend.position = "right")
}

p_vax25   <- plot_seird_compartments(df_vax25,   "Vaccination: 2.5% Daily")
p_vax3000 <- plot_seird_compartments(df_vax3000, "Vaccination: 3000/day")
(p_vax25 + p_vax3000) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")


# 7. Faceted Comparison Grid by Compartment

# Define desired compartment order and label mapping
compartment_levels <- c(
  "Susceptible", "Exposed", "Infectious",
  "Recovered (infection)", "Vaccinated (immune)", "Deaths"
)

compartment_labels <- c(
  S = "Susceptible",
  E = "Exposed",
  I = "Infectious",
  R_recov = "Recovered (infection)",
  V = "Vaccinated (immune)",
  D = "Deaths"
)

# Combine three scenarios and reshape to long format
df_comp <- bind_rows(df_baseline, df_vax25, df_vax3000) %>%
  select(t, scenario, S, E, I, R_recov, V, D) %>%
  pivot_longer(
    cols = c(S, E, I, R_recov, V, D),
    names_to = "compartment_raw",
    values_to = "value"
  ) %>%
  mutate(
    # Recode to human-readable labels
    compartment = recode(compartment_raw, !!!compartment_labels),
    # Convert to factor with defined order
    compartment = factor(compartment, levels = compartment_levels),
    # Ensure consistent scenario factor levels
    scenario = factor(scenario, levels = c("Baseline", "Vax_2.5%", "Vax_3000"))
  )


# Define custom color palette for scenarios
scenario_colors <- c(
  "Baseline"  = "#6c757d",   # grey
  "Vax_2.5%"  = "#1f77b4",   # blue
  "Vax_3000"  = "#d62728"    # red
)

# Plot with dashed lines for Deaths only
ggplot(df_comp, aes(x = t, y = value, color = scenario)) +
  geom_line(aes(linetype = compartment == "Deaths"), linewidth = 0.9) +
  scale_linetype_manual(
    values = c("FALSE" = "solid", "TRUE" = "dashed"),
    guide = "none"
  ) +
  scale_color_manual(values = scenario_colors, name = "Scenario") +
  facet_wrap(~ compartment, ncol = 3, scales = "free_y") +
  labs(
    title = "SEIRD — Vaccination Scenarios per Compartment",
    x = "Time (days)",
    y = "Number of individuals"
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11)
  )
