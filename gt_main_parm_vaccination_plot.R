# install and load necessary packages
#install.packages("outbreaks")
library(outbreaks)
library(dplyr)
library(ggplot2)
library(tidyr)
library(odin)
library(scales)
library(patchwork)

# read in data
ebola <- outbreaks::ebola_sim
ebola_clean <- ebola_sim$linelist
ebola_clean <- ebola_clean[, -c(8:11)]

# clean the data
ebola_clean_clean <- ebola_clean %>%
  filter((as.numeric(as.Date(date_of_onset) - as.Date(date_of_infection))) >=0) %>%
  filter((as.numeric(as.Date(date_of_outcome) - as.Date(date_of_onset))) >=0)

#  not all rows with NA removed, just rows with NA in desired columns
ebola_clean_clean <- ebola_clean %>%
  filter(
    is.na(date_of_onset) | is.na(date_of_infection) |
      (as.numeric(as.Date(date_of_onset) - as.Date(date_of_infection))) >= 0
  ) %>%
  filter(
    is.na(date_of_outcome) | is.na(date_of_onset) |
      (as.numeric(as.Date(date_of_outcome) - as.Date(date_of_onset))) >= 0
  )

nrow(ebola_clean_clean)

# prepare the observed cases

all_days <- seq(min(as.Date(ebola_clean_clean$date_of_onset), na.rm = TRUE),
                max(as.Date(ebola_clean_clean$date_of_onset), na.rm = TRUE),
                by = "day")
daily_new_cases_full <- data.frame(
  date  = all_days,
  cases = as.integer(table(factor(as.Date(ebola_clean_clean$date_of_onset), levels = all_days)))
)

# plot to checl
# plot(daily_new_cases_full$date, daily_new_cases_full$cases,type = "p",
#      xlab = "Days",
#      ylab = "Daily Cases",
#      cex.lab = 1.5,
#      cex.axis = 1.2)

# prepare the cumulative cases
cumulative_cases_full <- within(daily_new_cases_full, cum_cases <- cumsum(cases))
cumulative_df <- subset(cumulative_cases_full, cum_cases > 0)
cumulative_df$day_index <- seq_len(nrow(cumulative_df))

# plot to check
# plot(cumulative_cases_full$date, log(cumulative_cases_full$cum_cases),type = "p",
#      xlab = "Days",
#      ylab = "Daily Cases",
#      cex.lab = 1.5,
#      cex.axis = 1.2)

# fit and summary
cutoff <- as.Date("2014-05-27") # day 51
fit_df <- subset(cumulative_df, as.Date(date) < cutoff)
fit_df$day_index <- seq_len(nrow(fit_df))
fit_df <- fit_df[which(fit_df$cum_cases >= 5),]

# the glm model to obtain r
glm_model <- glm(log(cum_cases) ~ day_index, data = fit_df)
summary(glm_model)

# plot to check
ggplot(cumulative_df, aes(day_index, log(cum_cases))) +
  geom_point(size = 1) +
  geom_smooth(data = fit_df, aes(day_index, log(cum_cases)),
              method = "glm", formula = y ~ x, se = FALSE, color = "red") +
  #geom_vline(xintercept = nrow(fit_df), linetype = "dashed") +
  xlab("Days since first case") +
  ylab("Log cumulative infections") +
  theme_classic(base_size = 16)

# estimate the generation time
# from one being infected to the other being infected
# serial interval is symptom onset

# contacts
ct <- ebola$contacts

# lookup: case_id -> infection date
ebola_clean_clean$date_of_infection <- as.Date(ebola_clean_clean$date_of_infection)
inf_date <- setNames(ebola_clean_clean$date_of_infection, ebola_clean_clean$case_id)

# generation times (days) for pairs with both dates available
gt_days <- with(ct, as.numeric(inf_date[case_id] - inf_date[infector]))
gt_days <- gt_days[!is.na(gt_days)]

# summary statsics
summary(gt_days)

# remove outlier
gt_days_cutoff <- quantile(gt_days, 0.75) + 1.5*IQR(gt_days)
gt_days <- gt_days[gt_days <= gt_days_cutoff]
mean(gt_days)

R0_estimated <- exp(mean(gt_days)*coef(glm_model)[2])
cat("Estimated R0 is", R0_estimated, "\n")

# incubation period estimation
ebola_period_incubation <- ebola_clean_clean[, c(3,4)]
ebola_period_incubation <- drop_na(ebola_period_incubation)
incubation_period <- difftime(as.POSIXct(ebola_period_incubation$date_of_onset), as.POSIXct(ebola_period_incubation$date_of_infection), units = "days")
incubation_period <- incubation_period[which(incubation_period>0)]
outlier_point <-  quantile(incubation_period, 0.75) + 1.5* IQR(incubation_period)
incubation_period <- incubation_period[which(incubation_period<=outlier_point)]

#summary
summary(incubation_period)
mean(incubation_period)

# infectious period estimation
ebola_period_infectious<- ebola_clean_clean[, c(4,6,7)]
#ebola_period_infectious<-drop_na(ebola_period_infectious)
infectious_period<-difftime(as.POSIXct(ebola_period_infectious$date_of_outcome), as.POSIXct(ebola_period_infectious$date_of_onset), units = "days")
infectious_period<-infectious_period[infectious_period>=0]
mean(infectious_period, na.rm = T)

outlier_point_infectious <- quantile(infectious_period, 0.75, na.rm = T) + 1.5* IQR(infectious_period, na.rm = T)
infectious_period <- infectious_period[which(infectious_period<=outlier_point_infectious)]
mean(infectious_period, na.rm = T)

# infection - hospital period estimation, not used here
# ebola_period_hospital<-ebola_clean_clean[, c(4,5)]
# ebola_period_hospital<-drop_na(ebola_period_hospital)
# hospital_period<-difftime(as.POSIXct(ebola_period_hospital$date_of_hospitalisation), as.POSIXct(ebola_period_hospital$date_of_onset), units = "days")
# hospital_period<-hospital_period[hospital_period>=0]
# summary(hospital_period)

# different events

recovery_events <- ebola_period_infectious %>%
  filter(outcome == "Recover") %>%
  #filter(!is.na(date_of_onset) & !is.na(date_of_outcome)) %>%
  mutate(days_diff = as.numeric(as.Date(date_of_outcome) - as.Date(date_of_onset))) #%>%
  #filter(days_diff >= 0)
#summary(recovery_events$days_diff)

outlier_point_infectious_recovery <- quantile(recovery_events$days_diff, 0.75, na.rm = T) + 1.5* IQR(recovery_events$days_diff, na.rm = T)
infectious_period_recovery <- recovery_events$days_diff[which(recovery_events$days_diff<=outlier_point_infectious_recovery)]
mean(infectious_period_recovery)

death_events <- ebola_period_infectious %>%
  filter(outcome == "Death") %>%
  #filter(!is.na(date_of_onset) & !is.na(date_of_outcome)) %>%
  mutate(days_diff = as.numeric(as.Date(date_of_outcome) - as.Date(date_of_onset))) #%>%
  #filter(days_diff > 0)
summary(death_events$days_diff)

outlier_point_infectious_death <- quantile(death_events$days_diff, 0.75, na.rm = T) + 1.5* IQR(death_events$days_diff, na.rm = T)
infectious_period_death <- death_events$days_diff[which(death_events$days_diff<=outlier_point_infectious_death)]
mean(infectious_period_death)

# test if recovery and death difference in days are significantly different
# two-sample t-test (default two-sided)
t.test(recovery_events$days_diff, death_events$days_diff, var.equal = FALSE)

# proportion
proportion_recovery <- length(recovery_events$days_diff)/(length(recovery_events$days_diff) + length(death_events$days_diff))
proportion_death <- length(death_events$days_diff)/(length(recovery_events$days_diff) + length(death_events$days_diff))
cat("Sum is : ", proportion_death + proportion_recovery )

# from obtained parameters
################################
lambda_obtained <- 1/mean(as.numeric(incubation_period))
gamma_obtained <-  proportion_death* 1/mean(infectious_period_death) + proportion_recovery* 1/mean(infectious_period_recovery)  
#1/mean(as.numeric(infectious_period))
recovery_obtained <- 1/mean(infectious_period_recovery)
death_obtained <- 1/mean(infectious_period_death)

N_ini_val <- 6e5
beta0 <- 0.2323 # as.numeric(R0_estimated) * gamma_obtained
as.numeric(R0_estimated) * gamma_obtained

# the model baseline
###### SEIRD
seird_user_ode <- odin::odin({
  deriv(S) <- -beta * S * (I_R + I_D) / N
  deriv(E) <-  beta * S * (I_R + I_D) / N - lambda * E
  deriv(I_R) <-  lambda * E * proportion_recovery - recovery_obtained * I_R
  deriv(I_D) <-  lambda * E * proportion_death - death_obtained * I_D
  deriv(R) <-  recovery_obtained * I_R 
  deriv(D) <- death_obtained * I_D
  
  deriv(C) <- beta * S * (I_R + I_D) / N    # cumulative infections
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


}, verbose = FALSE)

# initialize model
seird_mod <- seird_user_ode$new(beta = beta0, lambda = lambda_obtained, recovery_obtained = recovery_obtained,death_obtained = death_obtained )
# how long to run
times <- seq(1,600*1)
# run the model
seird_mod <- seird_mod$run(times)

######
result_df_s1 <- as.data.frame(seird_mod)
N_ini_val - result_df_s1[600,]$S
result_df_s1[600,]$D
max(result_df_s1$I_R + result_df_s1$I_D)
which((result_df_s1$I_R + result_df_s1$I_D) == max(result_df_s1$I_R + result_df_s1$I_D))

# generate the vaccination intervention:
###### SEIRDV
seirdv_user_ode <- odin::odin({
  deriv(S) <- -beta * S * (I_R + I_D)/ N  - vacc_flow * S
  deriv(V) <- vacc_flow * S
  deriv(E) <-  beta * S * (I_R + I_D)/ N  - lambda * E
  deriv(I_R) <-  lambda * E * proportion_recovery - recovery_obtained * I_R
  deriv(I_D) <-  lambda * E * proportion_death - death_obtained * I_D
  deriv(R) <-  recovery_obtained * I_R 
  deriv(D) <- death_obtained * I_D
  
  deriv(C) <- beta * S * (I_R + I_D) / N      # cumulative infections
  initial(C) <- 0
  N <- S + E + I_R + I_D + R + V
  
  initial(S) <- N_ini - I_ini
  initial(E) <- 0
  initial(I_R) <- proportion_recovery * I_ini 
  initial(I_D) <- proportion_death    * I_ini
  initial(R) <- 0
  initial(D) <- 0
  initial(V) <- 0

  
  proportion_recovery <- 0.436
  proportion_death <- 0.564
  
  # user parameters
  N_ini <- user(6e5)
  I_ini <- user(1)
  lambda <- user()
  beta <- user()
  recovery_obtained <- user()
  death_obtained <- user()
  vacc_start <- user(120)
  vacc_rate <- user(0.025)
  vacc_eff <- user(0.84)
  
  # changed parameter
  vacc_flow <- if (t < vacc_start) 0 else vacc_rate * vacc_eff
  
}, verbose = FALSE)

# initialize model
seirdv_mod <- seirdv_user_ode$new(beta = beta0, lambda = lambda_obtained, recovery_obtained = recovery_obtained, death_obtained = death_obtained, 
                                  vacc_start = 60 )
# run the model
seirdv_mod <- seirdv_mod$run(times)
result_df_v1_60 <- as.data.frame(seirdv_mod)
result_df_v1_60$R[600]
result_df_v1_60$D[600]
min(result_df_v1_60$S)
max(result_df_v1_60$I_R + result_df_v1_60$I_D)
which((result_df_v1_60$I_R + result_df_v1_60$I_D) == max(result_df_v1_60$I_R + result_df_v1_60$I_D))

# initialize model
seirdv_mod <- seirdv_user_ode$new(beta = beta0, lambda = lambda_obtained, recovery_obtained = recovery_obtained, death_obtained = death_obtained, 
                                  vacc_start = 90 )
# run the model
seirdv_mod <- seirdv_mod$run(times)
result_df_v1_90 <- as.data.frame(seirdv_mod)
result_df_v1_90$R[600]
result_df_v1_90$D[600]
min(result_df_v1_90$S)
max(result_df_v1_90$I_R + result_df_v1_90$I_D)
which((result_df_v1_90$I_R + result_df_v1_90$I_D) == max(result_df_v1_90$I_R + result_df_v1_90$I_D))


##### if we want to use constant vaccination number
### assume each day vaccinate susceptible population 
# assume 84% effectiveness
seirdv_user_ode_v2 <- odin::odin({
  ## Flows
  deriv(S)   <- -beta * S * (I_R + I_D)/ N  - vacc_capped * vacc_eff
  deriv(V)   <-  vacc_capped * vacc_eff
  deriv(E)   <-  beta * S * (I_R + I_D)/ N  - lambda * E
  deriv(I_R) <-  lambda * E * proportion_recovery - recovery_obtained * I_R
  deriv(I_D) <-  lambda * E * proportion_death   - death_obtained    * I_D
  deriv(R)   <-  recovery_obtained * I_R
  deriv(D)   <-  death_obtained    * I_D
  
  ## Cumulative infections
  deriv(C) <- beta * S * (I_R + I_D)/ N 
  initial(C) <- 0
  
  N <- S + E + I_R + I_D + R + V
  
  initial(S) <- N_ini - I_ini
  initial(E) <- 0
  initial(I_R) <- proportion_recovery * I_ini 
  initial(I_D) <- proportion_death    * I_ini
  initial(R) <- 0
  initial(D) <- 0
  initial(V) <- 0
  
  ## Fixed fractions/outcomes
  proportion_recovery <- 0.436
  proportion_death    <- 0.564
  
  ## User parameters
  N_ini <- user(6e5)
  I_ini <- user(1)
  lambda <- user()
  beta   <- user()
  recovery_obtained <- user()
  death_obtained    <- user()
  vacc_eff <- user(0.84)
  
  # vaccination schedule & number (constant number per day)
  vacc_start     <- user(120)   # day intervention begins
  vacc_per_day   <- user(3000)   # constant number of vaccinations per day
  
  # daily vaccination flow turns on at vacc_start
  vacc_flow <- if (t < vacc_start) 0 else vacc_per_day
  
  # we can't vaccinate more than S remaining (S cannot < 0)
  vacc_capped <- if (S > vacc_flow) vacc_flow else S
  
  # outputs in case
  output(daily_vaccinations) <- vacc_capped
}, verbose = FALSE)


# initialize model
seirdv_mod <- seirdv_user_ode_v2$new(beta = beta0, lambda = lambda_obtained, recovery_obtained = recovery_obtained, death_obtained = death_obtained, 
                                  vacc_start = 60 )
# run the model
seirdv_mod <- seirdv_mod$run(times)
result_df_v2_60 <- as.data.frame(seirdv_mod)
result_df_v2_60$R[600]
result_df_v2_60$D[600]
max(result_df_v2_60$I_R + result_df_v2_60$I_D)
which((result_df_v2_60$I_R + result_df_v2_60$I_D) == max(result_df_v2_60$I_R + result_df_v2_60$I_D))

# 90 days
# initialize model
seirdv_mod <- seirdv_user_ode$new(beta = beta0, lambda = lambda_obtained, recovery_obtained = recovery_obtained, death_obtained = death_obtained, 
                                  vacc_start = 90 )
# run the model
seirdv_mod <- seirdv_mod$run(times)
result_df_v2_90 <- as.data.frame(seirdv_mod)
result_df_v2_90$R[600]
result_df_v2_90$D[600]
max(result_df_v2_90$I_R + result_df_v2_90$I_D)
which((result_df_v2_90$I_R + result_df_v2_90$I_D) == max(result_df_v2_90$I_R + result_df_v2_90$I_D))


# these are the plot settings
df_I <- bind_rows(
  result_df_s1      %>% transmute(t, value = I_R + I_D, scenario = "S1",    variant = "baseline / day 60"),
  result_df_v1_60      %>% transmute(t, value = I_R + I_D, scenario = "S4(a)", variant = "baseline / day 60"),
  result_df_v2_60      %>% transmute(t, value = I_R + I_D, scenario = "S4(b)", variant = "baseline / day 60"),
  result_df_v1_90   %>% transmute(t, value = I_R + I_D, scenario = "S4(a)", variant = "day 90"),
  result_df_v2_90   %>% transmute(t, value = I_R + I_D, scenario = "S4(b)", variant = "day 90")
)

df_V <- bind_rows(
  result_df_v1_60      %>% transmute(t, value = V, scenario = "S4(a)", variant = "baseline / day 60"),
  result_df_v2_60      %>% transmute(t, value = V, scenario = "S4(b)", variant = "baseline / day 60"),
  result_df_v1_90   %>% transmute(t, value = V, scenario = "S4(a)", variant = "day 90"),
  result_df_v2_90   %>% transmute(t, value = V, scenario = "S4(b)", variant = "day 90")
)

df_C <- bind_rows(
  result_df_s1      %>% transmute(t, value = C, scenario = "S1",    variant = "baseline / day 60"),
  result_df_v1_60      %>% transmute(t, value = C, scenario = "S4(a)", variant = "baseline / day 60"),
  result_df_v2_60      %>% transmute(t, value = C, scenario = "S4(b)", variant = "baseline / day 60"),
  result_df_v1_90   %>% transmute(t, value = C, scenario = "S4(a)", variant = "day 90"),
  result_df_v2_90   %>% transmute(t, value = C, scenario = "S4(b)", variant = "day 90")
)

# colors and the transparency etc.
cols  <- c("S1"="#1f77b4","S4(a)"="#E69F00","S4(b)"="#2ca02c")
alph  <- c("S1"=0.55, "S4(a)"=1, "S4(b)"=1)
fmt_k <- label_number(scale = 1/1000, big.mark = ",")
lts   <- c("baseline / day 60"="solid", "day 90"="dotdash")

base_theme <- theme_classic(base_size = 12) +
  theme(legend.position = "bottom")

# use panels to show
make_panel <- function(dat, ylab, tag_letter, sqrt_scale = FALSE,
                       show_legend = FALSE, custom_breaks = NULL) {
  p <- ggplot(dat, aes(t, value, color = scenario, alpha = scenario, linetype = variant)) +
    geom_line(linewidth = 1.1) +
    scale_color_manual(
      values = cols, name = "Scenario",
      limits = c("S1","S4(a)","S4(b)"), drop = FALSE
    ) +
    scale_linetype_manual(
      values = lts,  name = "Start",
      limits = c("baseline / day 60","day 90"), drop = FALSE
    ) +
    scale_alpha_manual(values = alph, guide = "none") +
    labs(x = "Time (days)", y = ylab) +
    base_theme +
    annotate("text", x = -Inf, y = Inf, label = paste0(tag_letter, ")"),
             hjust = -0.4, vjust = 1.0, fontface = "bold", size = 4)
  
  if (sqrt_scale) {
    p <- p + scale_y_sqrt(breaks = custom_breaks, labels = fmt_k,
                          expand = expansion(mult = c(0, 0.02)))
  } else {
    p <- p + scale_y_continuous(breaks = custom_breaks, labels = fmt_k,
                                expand = expansion(mult = c(0, 0.02)))
  }
  
  if (!show_legend) {
    p <- p + guides(color = "none", linetype = "none") + theme(legend.position = "none")
  } else {
    # --- dummy layer to force all legend entries, even if not in data for this panel ---
    dummy <- expand.grid(
      t = NA_real_, value = NA_real_,
      scenario = c("S1","S4(a)","S4(b)"),
      variant  = c("baseline / day 60","day 90")
    )
    p <- p +
      geom_line(data = dummy,
                aes(t, value, color = scenario, linetype = variant),
                inherit.aes = FALSE, linewidth = 1.1, show.legend = TRUE, na.rm = TRUE) +
      guides(color = guide_legend(nrow = 1, order = 1),
             linetype = guide_legend(nrow = 1, order = 2))
  }
  p
}


# determined fixed breaks
p_a <- make_panel(
  df_I, "Infected individuals (k)", "a",
  sqrt_scale = T, show_legend = F,
  custom_breaks = c( 2000, 5000, 10000, 20000, 40000, 60000, 80000)
)

p_b <- make_panel(
  df_V, "Vaccinated individuals (k)", "b",
  sqrt_scale = F, show_legend = T,
  custom_breaks = seq(100000, 600000, by = 100000)
)

p_c <- make_panel(
  df_C, "Cumulative infections (k)", "c",
  sqrt_scale = T, show_legend = F,
  custom_breaks = c(10000, 50000, 100000, 200000, 300000, 500000)
)

# plot it out tgt
(p_a | p_b | p_c) +
  plot_layout(widths = c(1,1,1)) &
  theme(legend.position = "bottom")





