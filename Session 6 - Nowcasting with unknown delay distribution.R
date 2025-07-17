
  # 2025-07-17
  # Amy Moore
  # SISMID - NFIDD
  # Session 6 - Nowcasting with unknown delay distribution
  
  #Load Packages
library("nfidd")
library("dplyr")
library("tidyr")
library("ggplot2")
library("tidybayes")

  #Set Up
set.seed(123)
options(cmdstanr_print_line_numbers = TRUE)
options(cmdstanr_warn_inits = FALSE)


  #Step 1: Generate Simulated Onset Data
gen_time_pmf <- make_gen_time_pmf()
ip_pmf <- make_ip_pmf()
onset_df <- simulate_onsets(
  make_daily_infections(infection_times), gen_time_pmf, ip_pmf
)
head(onset_df)

  #Add Truncation
cutoff <- 71

  #Simulate Reporting Delays
reporting_delay_pmf <- censored_delay_pmf(
  rlnorm, max = 15, meanlog = 1, sdlog = 0.5
)
plot(reporting_delay_pmf)

  #Now fit the reporting days
reporting_triangle <- onset_df |>
  filter(day < cutoff) |>
  mutate(
    reporting_delay = list(
      tibble(d = 0:15, reporting_delay = reporting_delay_pmf)
    )
  ) |>
  unnest(reporting_delay) |>
  mutate(
    reported_onsets = rpois(n(), onsets * reporting_delay)
  ) |>
  mutate(reported_day = day + d)

tail(reporting_triangle)

  #We can recover the onsets by summing reports by day
noisy_onsets_df <- reporting_triangle |>
  summarise(noisy_onsets = sum(reported_onsets), .by = day)

tail(noisy_onsets_df)

  #filter out events that were reported after the truncated date
filtered_reporting_triangle <- reporting_triangle |>
  filter(reported_day <= max(day))

tail(noisy_onsets_df)

  #Now we recover the onsets reported (with truncation)
available_onsets <- filtered_reporting_triangle |>
  summarise(available_onsets = sum(reported_onsets), .by = day)

tail(available_onsets)



  ##### Fitting the model to jointly estimate delay dist and nowcast #####

  #Model Set Up
joint_mod <- nfidd_cmdstan_model("joint-nowcast")
joint_mod

  #Model Fitting
joint_data <- list(
  n = length(unique(filtered_reporting_triangle$day)), # number of days
  m = nrow(filtered_reporting_triangle),               # number of reports
  p = filtered_reporting_triangle |>
    group_by(day) |>
    filter(d == max(d)) |>
    mutate(d = d + 1) |>
    pull(d),            # number of observations per day
  obs = filtered_reporting_triangle$reported_onsets,   # observed symptom onsets
  d = 16               # number of reporting delays
)
joint_nowcast_fit <- nfidd_sample(joint_mod, data = joint_data)
joint_nowcast_fit

  #Extract the Nowcast estimate
joint_nowcast_onsets <- joint_nowcast_fit |>
  gather_draws(nowcast[day]) |>
  ungroup() |>
  filter(.draw %in% sample(.draw, 100))

  #Plot the Nowcast estimate
ggplot(joint_nowcast_onsets, aes(x = day)) +
  geom_col(
    data = noisy_onsets_df, mapping = aes(y = noisy_onsets), alpha = 0.6
  ) +
  geom_line(mapping = aes(y = .value, group = .draw), alpha = 0.1) +
  geom_point(data = available_onsets, mapping = aes(y = available_onsets))




  ##### Add in estimating Rt too (triple joint model) #####

  #generation time pmf
plot(gen_time_pmf)

  #incubation period pmf
plot(ip_pmf)

  #Set Up for new model with all pieces (estimating delay dist, nowcast, and Rt)
joint_rt_mod <- nfidd_cmdstan_model("joint-nowcast-with-r")
joint_rt_mod


  #Fit the Model
joint_rt_data <- c(joint_data,
                   list(
                     gen_time_max = length(gen_time_pmf),
                     gen_time_pmf = gen_time_pmf,
                     ip_max = length(ip_pmf) - 1,
                     ip_pmf = ip_pmf,
                     h = 0 # this is a small easter egg for the attentive reader
                   )
)
joint_rt_fit <- nfidd_sample(
  joint_rt_mod, data = joint_rt_data,
  adapt_delta = 0.95,
  max_treedepth = 12,
  init = \() list(init_R = 1, rw_sd = 0.01)
)


joint_rt_fit



  #Extract Nowcast Estimate
joint_nowcast_with_r_onsets <- joint_rt_fit |>
  gather_draws(nowcast[day]) |>
  ungroup() |>
  filter(.draw %in% sample(.draw, 100))
ggplot(joint_nowcast_with_r_onsets, aes(x = day)) +
  geom_col(
    data = noisy_onsets_df, mapping = aes(y = noisy_onsets), alpha = 0.6
  ) +
  geom_line(mapping = aes(y = .value, group = .draw), alpha = 0.1) +
  geom_point(data = available_onsets, mapping = aes(y = available_onsets))


#Extract Rt estimates
joint_rt <- joint_rt_fit |>
  gather_draws(R[day]) |>
  ungroup() |>
  filter(.draw %in% sample(.draw, 100))

ggplot(joint_rt, aes(x = day, y = .value, group = .draw)) +
  geom_line(alpha = 0.1)












