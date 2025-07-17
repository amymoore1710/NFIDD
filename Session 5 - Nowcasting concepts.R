
  # 2025-07-17
  # Amy Moore
  # SISMID - NFIDD
  # Session 5 - Nowcasting concepts

  #Load Packages
library("nfidd")
library("dplyr")
library("tidyr")
library("ggplot2")
library("tidybayes")

  #Set Up
set.seed(123)
options(cmdstanr_print_line_numbers = TRUE)




  #Creating Data with a reporting delay
data(infection_times)
df <- infection_times |>
  mutate(
    onset_time = infection_time + rgamma(n(), shape = 5, rate = 1),
    report_time = onset_time + rlnorm(n(), meanlog = 1, sdlog = 0.5)
  )

  #Right truncate data at 70 days
cutoff <- 71
df_co <- df |>
  filter(report_time < cutoff)

## create time series of infections, onsets, and reports
df_co <- df_co |>
  mutate(
    infection_day = floor(infection_time),
    onset_day = floor(onset_time),
    report_day = floor(report_time)
  ) |>
  select(-infection_time, -onset_time, -report_time)

  #Convert to a time series
onset_ts <- df_co |>
  count(day = onset_day, name = "onsets")
reports_ts <- df_co |>
  count(day = report_day, name = "reports")

all_days <- expand_grid(day = seq(0, cutoff - 1)) |>
  full_join(onset_ts, by = "day") |>
  full_join(reports_ts, by = "day") |>
  replace_na(list(onsets = 0, reports = 0))



  #Visualize Data
combined <- all_days |>
  pivot_longer(c(onsets, reports), names_to = "variable")
ggplot(combined, aes(x = day, y = value)) +
  facet_grid(variable ~ .) +
  geom_col()


  #Full outbreak curves (no truncation)
# Use full outbreak dataset
final <- df |>
  mutate(onset_day = floor(onset_time)) |>
  select(-onset_time)
final_onset_ts <- final |>
  count(day = onset_day, name = "onsets")
final_all_days <- expand_grid(day = seq(0, max(final_onset_ts$day))) |>
  full_join(final_onset_ts, by = "day") |>
  replace_na(list(onsets = 0)) |>
  mutate(cutoff = "final")
intermediate <- combined |>
  filter(variable == "onsets") |>
  select(-variable) |>
  rename(onsets = value) |>
  mutate(cutoff = "70 days")
combined_cutoffs <- rbind(
  intermediate,
  final_all_days
)
ggplot(combined_cutoffs, aes(x = day, y = onsets, colour = cutoff)) +
  geom_line() +
  scale_colour_brewer(palette = "Dark2") +
  geom_vline(xintercept = cutoff, linetype = "dashed")





  ##### Simple Model for Nowcasting #####

proportion_reported <- plnorm(1:15, 1, 0.5)
plot(proportion_reported)

gen_time_pmf <- make_gen_time_pmf()
ip_pmf <- make_ip_pmf()
onset_df <- simulate_onsets(
  make_daily_infections(infection_times), gen_time_pmf, ip_pmf
)
reported_onset_df <- onset_df |>
  filter(day < cutoff) |>
  mutate(proportion_reported = c(rep(1, n() - 15), rev(proportion_reported)),
         reported_onsets = rpois(n(), onsets * proportion_reported)
  )
tail(reported_onset_df)

reported_onset_df |>
  ggplot(aes(x = day, y = reported_onsets)) +
  geom_col()


mod <- nfidd_cmdstan_model("simple-nowcast")
mod

data <- list(
  n = nrow(reported_onset_df) - 1,
  obs = reported_onset_df$reported_onsets[-1],
  report_max = length(proportion_reported) - 1,
  report_cdf = proportion_reported 
)
simple_nowcast_fit <- nfidd_sample(mod, data = data)
simple_nowcast_fit

nowcast_onsets <- simple_nowcast_fit |>
  gather_draws(onsets[day]) |>
  ungroup() |>
  filter(.draw %in% sample(.draw, 100)) |>
  mutate(day = day + 1)
ggplot(nowcast_onsets, aes(x = day)) +
  geom_line(mapping = aes(y = .value, group = .draw), alpha = 0.1) +
  geom_col(data = reported_onset_df, mapping = aes(y = onsets), alpha = 0.6) +
  geom_point(data = reported_onset_df, mapping = aes(y = reported_onsets))



  #Adding in a random walk to smooth the now casting estimates
rw_mod <- nfidd_cmdstan_model("simple-nowcast-rw")
rw_mod
rw_nowcast_fit <- nfidd_sample(rw_mod, data = data)
rw_nowcast_fit

rw_nowcast_onsets <- rw_nowcast_fit |>
  gather_draws(onsets[day]) |>
  ungroup() |>
  filter(.draw %in% sample(.draw, 100)) |> ## sample 100 iterations randomly
  mutate(day = day + 1)

ggplot(rw_nowcast_onsets, aes(x = day)) +
  geom_col(data = reported_onset_df, mapping = aes(y = onsets), alpha = 0.6) +
  geom_line(mapping = aes(y = .value, group = .draw), alpha = 0.1) +
  geom_point(data = reported_onset_df, mapping = aes(y = reported_onsets))





    ##### What if we misspecify the delay distribution? #####

  #Now we estimate that the reporting proportion is way different than it actually was simulated
wrong_proportion_reported <- pgamma(1:15, 2, 3)
plot(wrong_proportion_reported)

wrong_delay_data <- data
wrong_delay_data$report_cdf <- wrong_proportion_reported

gamma_nowcast_fit <- nfidd_sample(rw_mod, data = wrong_delay_data)
gamma_nowcast_fit


gamma_nowcast_onsets <- gamma_nowcast_fit |>
  gather_draws(onsets[day]) |>
  ungroup() |>
  filter(.draw %in% sample(.draw, 100)) |>
  mutate(day = day + 1)
ggplot(gamma_nowcast_onsets, aes(x = day)) +
  geom_col(data = reported_onset_df, mapping = aes(y = onsets), alpha = 0.6) +
  geom_line(mapping = aes(y = .value, group = .draw), alpha = 0.1) +
  geom_point(data = reported_onset_df, mapping = aes(y = reported_onsets))


# So now casting works great if we know what the reporting delay was
# but if we don't know the delay the nowcasted estimates can be way off