
# Ebola Data Clean

# Packages
library(tidyverse) 
library(outbreaks)  
library(janitor)  

# Ebola linelist + contacts

data("ebola_sim", package = "outbreaks")

linelist_raw  <- ebola_sim$linelist %>% as_tibble()
contacts_raw  <- ebola_sim$contacts %>% as_tibble()

# linelist
linelist <- linelist_raw %>%
  clean_names() %>%
  mutate(
    date_of_infection       = as.Date(date_of_infection),
    date_of_onset           = as.Date(date_of_onset),
    date_of_hospitalisation = as.Date(date_of_hospitalisation),
    date_of_outcome         = as.Date(date_of_outcome),
    outcome = factor(outcome, levels = c("Recover", "Death"))
  )

# Derived delay variables
linelist <- linelist %>%
  mutate(
    delay_incubation       = as.numeric(date_of_onset - date_of_infection),
    delay_onset_to_hosp    = as.numeric(date_of_hospitalisation - date_of_onset),
    delay_onset_to_outcome = as.numeric(date_of_outcome - date_of_onset)
  )


# Basic cleaning: remove negative delays
ll_clean <- linelist %>%
  filter(
    is.na(delay_incubation)       | delay_incubation       >= 0,
    is.na(delay_onset_to_hosp)    | delay_onset_to_hosp    >= 0,
    is.na(delay_onset_to_outcome) | delay_onset_to_outcome >= 0
  )

# contacts table
contacts <- contacts_raw %>%
  clean_names() %>%
  rename(infector = infector, infectee = case_id)

# Save cleaned datasets to CSV
write_csv(ll_clean, "linelist_clean.csv")
write_csv(contacts, "contacts_clean.csv")
