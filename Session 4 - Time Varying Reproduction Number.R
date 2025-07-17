
  # 2025-07-17
  # Amy Moore
  # SISMID - NFIDD
  # Session 4 - Time Varying Reproduction Number
  
  #Load in Packages
library("nfidd")
library("dplyr")
library("tidyr")
library("ggplot2")
library("tidybayes")
library("purrr")

  #Set Seed
set.seed(123)
options(cmdstanr_print_line_numbers = TRUE)



renewal


inf_ts <- make_daily_infections(infection_times)
head(inf_ts)


  #Creates discrete values for generational time
gen_time_pmf <- censored_delay_pmf(rgamma, max = 14, shape = 4, rate = 1)

gen_time_pmf <- gen_time_pmf[-1] ## remove first element
gen_time_pmf <- gen_time_pmf / sum(gen_time_pmf) ## renormalise

  #Define the model to fit
r_mod <- nfidd_cmdstan_model("estimate-r")
r_mod

  #create dataset for fitting
data <- list(
  n = nrow(inf_ts) - 1,
  obs = inf_ts$infections[-1],
  I0 = inf_ts$infections[1],
  gen_time_max = length(gen_time_pmf),
  gen_time_pmf = gen_time_pmf
)
  #fit the model to the data
r_fit <- nfidd_sample(r_mod, data = data)
r_fit



# Extract posterior draws for estimating Rt at each time point
r_posterior <- r_fit |>
  gather_draws(R[infection_day]) |>
  ungroup() |>
  mutate(infection_day = infection_day - 1) |> 
  filter(.draw %in% sample(.draw, 100))

ggplot(
  data = r_posterior,
  aes(x = infection_day, y = .value, group = .draw))  +
  geom_line(alpha =  0.1) +
  labs(title = "Estimated Rt", 
       subtitle = "Model 1: renewal equation from infections")

  #Posterior draws for estimating # of infections at each time point
inf_posterior <- r_fit |>
  gather_draws(infections[infection_day]) |>
  ungroup() |>
  mutate(infection_day = infection_day - 1) |> 
  mutate(infections = map_dbl(.value, ~ rpois(1, .x))) |>
  filter(.draw %in% sample(.draw, 100))

ggplot(inf_posterior, mapping = aes(x = infection_day)) +
  geom_line(mapping = aes(y = .value, group = .draw), alpha = 0.1) +
  geom_line(
    data = inf_ts, mapping = aes(y = infections), colour = "red"
  ) +
  labs(title = "Infections, estimated (grey) and observed (red)", 
       subtitle = "Model 1: renewal equation from infections")





  ##### Now switching to assuming we don't know infection dates, but we know symptomatic dates #####

  #PMF of delay between infection date to symptom onset
ip_pmf <- censored_delay_pmf(rgamma, max = 14, shape = 5, rate = 1)
  #Convolve infection dates with the pmf above --> gives estimate of when symptoms appear
onsets <- convolve_with_delay(inf_ts$infections, ip_pmf)
  #Convolution gives expected symptom appearance, this adds in individual observations as a layer of uncertainty
obs <- rpois(length(onsets), onsets)


r_inf_mod <- nfidd_cmdstan_model("estimate-inf-and-r")
r_inf_mod

  #Get Bayesian estimates from the model
data <- list(
  n = length(obs) - 1,
  obs = obs[-1],
  I0 = inf_ts$infections[1],
  gen_time_max = length(gen_time_pmf),
  gen_time_pmf = gen_time_pmf,
  ip_max = length(ip_pmf) - 1,
  ip_pmf = ip_pmf
)
  #Fit Bayesian Model -> Estimate Posteriors
r_inf_fit <- nfidd_sample(
  r_inf_mod, data = data, init = \() list(init_R = 1)
)
r_inf_fit

  #Extract estimates of posterior Rt estimates
  # (keeping in mind that the model was fit knowing 
  #  infections on day 1 and symptom dates, not infection dates)
r_inf_posteriors <- r_inf_fit |>
  gather_draws(infections[infection_day], R[infection_day]) |>
  ungroup() |>
  mutate(infection_day = infection_day - 1) |> 
  filter(.draw %in% sample(.draw, 100))

  #Extract estimates of posterior Infection Numbers
inf_posterior <- r_inf_posteriors |>
  filter(.variable == "infections")
ggplot(inf_posterior, mapping = aes(x = infection_day)) +
  geom_line(mapping = aes(y = .value, group = .draw), alpha = 0.1) +
  geom_line(
    data = inf_ts, mapping = aes(y = infections), colour = "red"
  ) +
  labs(title = "Infections, estimated (grey) and observed (red)", 
       subtitle = "Model 2: renewal equation from symptom onsets")

  #Plot Rt estimates
r_inf_posterior <- r_inf_posteriors |>
  filter(.variable == "R")
ggplot(
  r_inf_posterior, mapping = aes(x = infection_day, y = .value, group = .draw)
) +
  geom_line(alpha = 0.1)





  ##### Adding in the assumption that Rt is dependent on R(t-1) - autocorrelation #####

geometric_random_walk

R <- geometric_random_walk(init = 1, noise = rnorm(100), std = 0.1)
data <- tibble(t = seq_along(R), R = exp(R))

# Generate normal(1,1) prior samples for comparison
normal_prior <- rnorm(100, mean = 1, sd = 1)
normal_data <- tibble(t = seq_along(normal_prior), R = normal_prior)

ggplot(data, aes(x = t, y = R)) +
  geom_line() +
  geom_line(data = normal_data, aes(x = t, y = R), colour = "red") +
  labs(title = "Simulated data from a random walk model",
       subtitle = "Random walk (black) vs Normal(1,1) prior (red)",
       x = "Time",
       y = "R")


  #Adding this random walk into the Rt model
rw_mod <- nfidd_cmdstan_model("estimate-inf-and-r-rw")
rw_mod

data <- list(
  n = length(obs) - 1,
  obs = obs[-1],
  I0 = inf_ts$infections[1],
  gen_time_max = length(gen_time_pmf),
  gen_time_pmf = gen_time_pmf,
  ip_max = length(ip_pmf) - 1,
  ip_pmf = ip_pmf
)
r_rw_inf_fit <- nfidd_sample(
  rw_mod, data = data, max_treedepth = 12, 
  init = \() list(init_R = 1, rw_sd = 0.01)
)


rw_posteriors <- r_rw_inf_fit |>
  gather_draws(infections[infection_day], R[infection_day]) |>
  ungroup() |>
  mutate(infection_day = infection_day - 1) |>
  filter(.draw %in% sample(.draw, 100))