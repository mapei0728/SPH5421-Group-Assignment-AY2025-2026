# Ebola Cleaned Data Analysis: Delays, Generation Time, Growth Rate, R0, β

library(tidyverse)
library(gt)
library(scales)

# 0. Load Cleaned Data
ll_clean  <- read_csv("linelist_clean.csv",  show_col_types = FALSE)
contacts  <- read_csv("contacts_clean.csv", show_col_types = FALSE)


# 1. Incubation Delay (Infection to Onset)
incubation <- ll_clean %>%
  filter(!is.na(delay_incubation), delay_incubation > 0) %>%
  pull(delay_incubation)

cutoff_incub <- quantile(incubation, 0.75) + 1.5 * IQR(incubation)
mean_raw_incub    <- mean(incubation)
mean_cutoff_incub <- mean(incubation[incubation <= cutoff_incub])

cat("Delay = Incubation (infection → onset):\n")
cat("Raw mean               :", round(mean_raw_incub, 6), "\n")
cat("Upper IQR cutoff       :", round(cutoff_incub, 6), "\n")
cat("Mean after IQR cutoff  :", round(mean_cutoff_incub, 6), "\n\n")


# 2. Onset to Outcome Delay (All Outcomes)
all_outcomes <- ll_clean %>%
  pull(delay_onset_to_outcome) %>%
  na.omit()

cutoff_all <- quantile(all_outcomes, 0.75) + 1.5 * IQR(all_outcomes)
mean_raw_all    <- mean(all_outcomes)
mean_cutoff_all <- mean(all_outcomes[all_outcomes <= cutoff_all])

cat("All outcomes combined:\n")
cat("Raw mean               :", round(mean_raw_all, 6), "\n")
cat("Upper IQR cutoff       :", round(cutoff_all, 6), "\n")
cat("Mean after IQR cutoff  :", round(mean_cutoff_all, 6), "\n\n")


# 3. Onset to Death Delay
death <- ll_clean %>%
  filter(outcome == "Death") %>%
  pull(delay_onset_to_outcome) %>%
  na.omit()

cutoff_death <- quantile(death, 0.75) + 1.5 * IQR(death)
mean_raw_death    <- mean(death)
mean_cutoff_death <- mean(death[death <= cutoff_death])

cat("Outcome = Death:\n")
cat("Raw mean               :", round(mean_raw_death, 6), "\n")
cat("Upper IQR cutoff       :", round(cutoff_death, 6), "\n")
cat("Mean after IQR cutoff  :", round(mean_cutoff_death, 6), "\n\n")


# 4. Onset to Recover Delay
recover <- ll_clean %>%
  filter(outcome == "Recover") %>%
  pull(delay_onset_to_outcome) %>%
  na.omit()

cutoff_recover <- quantile(recover, 0.75) + 1.5 * IQR(recover)
mean_raw_recover    <- mean(recover)
mean_cutoff_recover <- mean(recover[recover <= cutoff_recover])

cat("Outcome = Recover:\n")
cat("Raw mean               :", round(mean_raw_recover, 6), "\n")
cat("Upper IQR cutoff       :", round(cutoff_recover, 6), "\n")
cat("Mean after IQR cutoff  :", round(mean_cutoff_recover, 6), "\n\n")

# 6.Outcome Distribution Summary
ll_clean %>%
  filter(!is.na(outcome)) %>%
  count(outcome) %>%
  mutate(ratio = round(n / sum(n), 6))


# 7. Generation Time Estimation (Infector to Infectee)
gt_data <- contacts %>%
  left_join(ll_clean %>% select(case_id, date_of_infection), by = c("infector" = "case_id")) %>%
  rename(infector_infection = date_of_infection) %>%
  left_join(ll_clean %>% select(case_id, date_of_infection), by = c("infectee" = "case_id")) %>%
  rename(infectee_infection = date_of_infection) %>%
  mutate(gt = as.numeric(infectee_infection - infector_infection)) %>%
  filter(!is.na(gt) & gt >= 0)

cat("Raw GT distribution:\n")
print(summary(gt_data$gt))

# Cutoff = 30 days
gt_cutoff <- 30
gt_data_clean <- gt_data %>% filter(gt <= gt_cutoff)

cat("\nRecords after filtering (gt ≤", gt_cutoff, "):", nrow(gt_data_clean), "\n")

# Function to summarise GT
make_gt_summary <- function(df, label) {
  df %>%
    summarise(
      n = n(),
      mean = mean(gt),
      sd = sd(gt),
      median = median(gt),
      q25 = quantile(gt, 0.25),
      q75 = quantile(gt, 0.75),
      min = min(gt),
      max = max(gt)
    ) %>%
    mutate(version = label)
}

summary_raw   <- make_gt_summary(gt_data, "Raw (no cutoff)")
summary_clean <- make_gt_summary(gt_data_clean, paste0("Cutoff ≤ ", gt_cutoff, " days"))

# Combine and display summary
gt_summary_compare <- bind_rows(summary_raw, summary_clean) %>%
  relocate(version, .before = n)

gt_summary_compare %>%
  gt(groupname_col = "version") %>%
  tab_header(title = paste0("Generation Time Summary (Raw vs Cutoff ≤ ", gt_cutoff, " days)")) %>%
  cols_label(
    version = "Data Version",
    n = "N",
    mean = "Mean",
    sd = "SD",
    median = "Median",
    q25 = "Q25",
    q75 = "Q75",
    min = "Min",
    max = "Max"
  )


# 8. Epidemic Growth Rate (r)
cutoff_days <- 49     # First 50 days
min_cum     <- 5      # Minimum cumulative count to avoid log(0)

# Prepare daily case count
cases <- ll_clean %>%
  filter(!is.na(date_of_onset)) %>%
  mutate(day = as.numeric(date_of_onset - min(date_of_onset)))

times <- tibble(day = 0:cutoff_days)

daily_cases <- cases %>%
  count(day, name = "cases") %>%
  right_join(times, by = "day") %>%
  arrange(day) %>%
  mutate(
    cases = replace_na(cases, 0),
    cumulative = cumsum(cases),
    date_of_onset = min(ll_clean$date_of_onset) + day
  )

growth_data <- daily_cases %>%
  filter(cumulative >= min_cum)

# Fit exponential growth (log-linear)
fit  <- glm(log(cumulative) ~ day, data = growth_data)
smry <- summary(fit)

r_est <- coef(fit)[["day"]]
se_r  <- smry$coefficients["day", "Std. Error"]
z_val <- qnorm(0.975)
r_CI  <- c(r_est - z_val * se_r, r_est + z_val * se_r)

# Plot r fit
label_text <- sprintf("r = %.4f (95%% CI: %.4f – %.4f)", r_est, r_CI[1], r_CI[2])

p_r <- ggplot(growth_data, aes(x = day, y = log(cumulative))) +
  geom_point(color = "steelblue", alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  annotate("text", x = max(growth_data$day) * 0.6,
           y = max(log(growth_data$cumulative)) * 0.9,
           label = label_text, hjust = 0, size = 4.5) +
  labs(
    title = "Exponential Growth Fit (Activity 2 Method)",
    subtitle = sprintf("Window: First %d days; Cumulative ≥ %d", cutoff_days, min_cum),
    x = "Days Since First Onset",
    y = "log(Cumulative Cases)"
  ) +
  theme_minimal()

print(p_r)


# 9. Estimate R0 and β
r     <- r_est
tg    <- 10.29650 
R0    <- exp(r * tg)  

cat("\nEstimated R0 =", R0, "\n")

# Infectious periods
T_IR <- 16.09345    # Mean time to recovery
T_ID <- 7.49748     # Mean time to death

# Recovery and death rates
gamma_R <- 1 / T_IR
gamma_D <- 1 / T_ID
gamma   <- 0.436 * gamma_R + 0.564 * gamma_D   # Weighted average

# Transmission rate β
beta <- R0 * gamma

# Output
cat("gamma_R =", round(gamma_R, 6), "\n")
cat("gamma_D =", round(gamma_D, 6), "\n")
cat("gamma   =", round(gamma, 6), "\n")
cat("beta    =", round(beta, 6), "\n")
