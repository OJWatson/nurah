# --------------------------
# Drought Mortality Analysis in Somalia, Ethiopia, Kenya (2021–2022)
# Approach 1: Historical vs Actual CMR
# Approach 2: IPC Phase-based mortality calculations
# --------------------------
library(tidyverse)
# Create helper function to compute total deaths
calculate_deaths <- function(pop_millions, cdr_per_10k_day, days = 365) {
  pop <- pop_millions * 1e6
  cdr <- cdr_per_10k_day / 10000
  deaths <- pop * cdr * days
  return(round(deaths))
}

# Pivot and clean function
piv_clean <- function(res) {

tidy_df <- res %>%
  as.data.frame() %>%
  mutate(id = 1) %>%
  pivot_longer(
    cols = -id,
    names_to = "full_name",
    values_to = "value"
  ) %>%
  separate(
    col = full_name,
    into = c("country", "year", "name"),
    sep = "\\.",
    convert = TRUE
  ) %>%
  select(country, year, name, value) %>%
  mutate(year = gsub("year", "", year) %>% as.numeric)
tidy_df
}

# --------------------------
# Data input
# --------------------------

# ---- SOMALIA ----
# Source: https://reliefweb.int/report/somalia/mortality-patterns-somalia-retrospective-estimates-and-scenario-based-forecasting-report-1-february-2023
# 2021: 1.6M start - 3.5M end in IPC 3+, estimated https://reliefweb.int/report/somalia/somalia-food-security-and-nutrition-analysis-post-gu-2021-technical-release-september-2021
# 2022: 4.1M start - 5.6M end in IPC 3+, estimated from Oct 2022 snapshot
# worst CMR based on 2011 -https://reliefweb.int/report/somalia/mortality-among-populations-southern-and-central-somalia-affected-severe-food
# CMR actual based on SMART analysis in 2021-2022: https://reliefweb.int/report/somalia/mortality-patterns-somalia-retrospective-estimates-and-scenario-based-forecasting-report-1-february-2023
somalia <- list(
  year2021 = list(
    pop = 2.55,
    cmr_actual = 0.325,
    cmr_worst = 2.5,
    cmr_actual_low = 0.275,
    cmr_actual_high = 0.4
  ),
  year2022 = list(
    pop = 4.9,
    cmr_actual = 0.36, # https://reliefweb.int/report/somalia/mortality-patterns-somalia-retrospective-estimates-and-scenario-based-forecasting-report-1-february-2023
    cmr_worst = 2.5,
    cmr_actual_low = 0.3,
    cmr_actual_high = 0.46
  )
)

# ---- ETHIOPIA ----
# Sources: IPC based
# 2021 16.5M & 2022 23.6M -
# IPC Phase populations from 2021-2022 reports
# CMR Worst Based on 2000 in Gode of 3.2 (2.4 - 3.8)
# https://www.cdc.gov/mmwr/preview/mmwrhtml/mm5015a2.htm#:~:text=A%20total%20of%20595%20households,4%20%287%29.%20Of%20the
# Or highs of 4-5 but also IDP Somalia related:
# https://reliefweb.int/report/ethiopia/unicef-ethiopia-horn-africa-crisis-situation-report-11-22-28-september-2011#:~:text=Somali%20Refugee%20Situation,1.1%20deaths%20per%2010%2C000%20daily.
# To be conservative, given the 2000 CMR was in one region and while ago, will go marginally below their lowest and closer to IPC 5
# For actual, no data. But the best study to be found (which U5MR) had 0.38 - 0.47 in drought during 2009-2014, so going marginally as CMR and to align more
# with Somalia, 0.36
# https://www.nature.com/articles/s41598-017-02271-5
ethiopia <- list(
  year2021 = list(
    pop = 16.7,
    cmr_actual = 0.36,
    cmr_worst = 2.4,
    cmr_actual_low = 0.25,
    cmr_actual_high = 0.4
  ),
  year2022 = list(
    pop = 23.6,
    cmr_actual = 0.36,
    cmr_worst = 2.4,
    cmr_actual_low = 0.25,
    cmr_actual_high = 0.4
  )
)

# ---- KENYA ----
# Sources: https://fews.net/east-africa/kenya
# https://humanitarianaction.info/plan/1073
# 2021: 2.9M in IPC 3+
# 2022: 4.4M in IPC 3+
# Very little data directly on CMR - similar mongst IDP from Somalia had up to 1.94
# But this dropped to 0.4 when in camps and managed
# Likely not as high as other countries given broader support and generally managed droughts.
# So just 1.5 chosen to reflect below the worst seen amongst those travellers, and same CMR as Ethiopia assumed

kenya <- list(
  year2021 = list(
    pop = 2.9,
    cmr_actual = 0.32,
    cmr_worst = 1.5,
    cmr_actual_low = 0.25,
    cmr_actual_high = 0.4
  ),
  year2022 = list(
    pop = 4.4,
    cmr_actual = 0.32,
    cmr_worst = 1.5,
    cmr_actual_low = 0.25,
    cmr_actual_high = 0.4
  )
)

# --------------------------
# Approach 1: Total deaths with and without aid
# --------------------------

countries <- list(Somalia = somalia, Ethiopia = ethiopia, Kenya = kenya)
approach1_results <- list()

for (country_name in names(countries)) {
  country <- countries[[country_name]]
  result <- list()

  for (year in names(country)) {
    data <- country[[year]]

    actual_deaths <- calculate_deaths(data$pop, data$cmr_actual)
    no_aid_deaths <- calculate_deaths(data$pop, data$cmr_worst)
    saved_lives <- no_aid_deaths - actual_deaths

    result[[year]] <- list(
      actual_deaths = actual_deaths,
      no_aid_deaths = no_aid_deaths,
      saved_lives = saved_lives,
      actual_low = calculate_deaths(data$pop, data$cmr_actual_low),
      actual_high = calculate_deaths(data$pop, data$cmr_actual_high)
    )
  }

  approach1_results[[country_name]] <- result
}

# Print Approach 1 results
cat("---- Approach 1: Historical vs Actual CMR ----\n")
for (country in names(approach1_results)) {
  cat(paste0("\n", country, "\n"))
  for (year in names(approach1_results[[country]])) {
    vals <- approach1_results[[country]][[year]]
    cat(paste0(year, ": ",
               "Actual = ", vals$actual_deaths,
               ", No Aid = ", vals$no_aid_deaths,
               ", Saved = ", vals$saved_lives,
               " (Low = ", vals$actual_low,
               ", High = ", vals$actual_high, ")\n"))
  }
}

# --------------------------
# Approach 2: IPC phase-based mortality
# --------------------------
# Phase-based CMRs:
# IPC 2: ~0.3/10k/day
# IPC 3: ~0.75/10k/day
# IPC 4: ~1.5/10k/day
# IPC 5: ~2.0/10k/day
# Source: IPC Technical Manual v3.1 https://www.ipcinfo.org/fileadmin/user_upload/ipcinfo/manual/IPC_Technical_Manual_3_Final.pdf

ipc_cmrs <- c("2" = 0.3, "3" = 0.75, "4" = 1.5, "5" = 2.0, "6" = 2.5)

# Population by phase, hardcoded from IPC reports

ipc_data <- list(
  Somalia = list(
    year2021 = list("2" = 3.71, "3" = 2.83, "4" = 0.64, "5" = 0.0),
    year2022 = list("2" = 8.0, "3" = 3.88, "4" = 1.5, "5" = 0.214)
  ),
  Ethiopia = list(
    year2021 = list("2" = 8.0, "3" = 10.0, "4" = 6.0, "5" = 0.5),
    year2022 = list("2" = 12.0, "3" = 15.0, "4" = 8.0, "5" = 0.6)
  ),
  Kenya = list(
    year2021 = list("2" = 4.0, "3" = 1.8, "4" = 0.37, "5" = 0.0),
    year2022 = list("2" = 10.4, "3" = 3.2, "4" = 1.2, "5" = 0.0)
  )
)

# Shift phases up by 1 for counterfactual
shift_phase <- function(phase) {
  if (phase == "2") return("3")
  if (phase == "3") return("4")
  if (phase == "4") return("5")
  return("6")  # IPC 5 goes to 6 which we assign 2.5 CMR to
}

approach2_results <- list()

for (country in names(ipc_data)) {
  years <- ipc_data[[country]]
  result <- list()

  for (year in names(years)) {
    phases <- years[[year]]

    deaths_actual <- 0
    deaths_counterfactual <- 0

    for (phase in names(phases)) {
      pop <- phases[[phase]]
      cdr_actual <- ipc_cmrs[[phase]]
      phase_shifted <- shift_phase(phase)
      cdr_counterfactual <- ipc_cmrs[[phase_shifted]]

      deaths_actual <- deaths_actual + calculate_deaths(pop, cdr_actual)
      deaths_counterfactual <- deaths_counterfactual + calculate_deaths(pop, cdr_counterfactual)
    }

    result[[year]] <- list(
      actual_deaths = deaths_actual,
      no_aid_deaths = deaths_counterfactual,
      saved_lives = deaths_counterfactual - deaths_actual
    )
  }

  approach2_results[[country]] <- result
}

# Print Approach 2 results
cat("\n---- Approach 2: IPC Phase-Based Mortality ----\n")
for (country in names(approach2_results)) {
  cat(paste0("\n", country, "\n"))
  for (year in names(approach2_results[[country]])) {
    vals <- approach2_results[[country]][[year]]
    cat(paste0(year, ": ",
               "Actual = ", vals$actual_deaths,
               ", No Aid = ", vals$no_aid_deaths,
               ", Saved = ", vals$saved_deaths, "\n"))
  }
}


# Collate all and plot

res <- rbind(
  piv_clean(approach1_results) %>% mutate(method = "Post-2000 CMR Extremes"),
  piv_clean(approach2_results) %>% mutate(method = "IPC Phase-Based Mortality")
)

res %>%
  filter(!name %in% c("actual_high", "actual_low")) %>%
  ggplot(aes(year, value, fill = name, alpha = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_alpha_manual(values = c(0.5,1)) +
  facet_wrap(~country)


res %>% filter(name == "saved_lives") %>% group_by(method) %>% summarise(sum(value))
