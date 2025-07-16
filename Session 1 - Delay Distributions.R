
  # 2025-07-16
  # Amy Moore
  # SISMID - NFIDD
  # Session 1 - Delay Distributions

  #Load in Packages

library("nfidd")
library("ggplot2")
library("dplyr")
library("tidyr")
library("lubridate")
library("posterior")
library("tidybayes")

  #Set Seed

set.seed(1234)
options(cmdstanr_print_line_numbers = TRUE)

  #Load in Data

data(infection_times)
head(infection_times)




  # visualize the infection curve

ggplot(infection_times, aes(x = infection_time)) +
  geom_histogram(binwidth = 1) +
  scale_x_continuous(n.breaks = 10) +
  labs(x = "Infection time (in days)", y = "Number of infections",
       title = "Infections during an outbreak")


  #Incorporate delays into the distribution
df <- add_delays(infection_times)


# convert our data frame to long format
dfl <- df |>
  pivot_longer(
    cols = c(infection_time, onset_time, hosp_time),
    names_to = "type", values_to = "time"
  ) |> 
  mutate(type = ordered(type, 
                        levels = c("infection_time", "onset_time", "hosp_time"), 
                        labels = c("Infections", "Symptom onsets", "Hospitalisations")
  ))

# plot the three event distributions
ggplot(dfl, aes(x = time)) +
  geom_histogram(position = "dodge", binwidth = 1) +
  facet_wrap(~ type, ncol = 1) +
  xlab("Time (in days)") +
  ylab("Count")



  #Creating the Stan Model
mod <- nfidd_cmdstan_model("lognormal")
mod


## Specify the time from onset to hospitalization
df_onset_to_hosp <- df |>
  mutate(onset_to_hosp = hosp_time - onset_time) |> 
  # exclude infections that didn't result in hospitalization
  drop_na(onset_to_hosp)


## Use the data to sample from the model posterior
res <- nfidd_sample(
  mod,
  data = list(
    n = nrow(df_onset_to_hosp),
    y = df_onset_to_hosp$onset_to_hosp
  )
)
res$summary()


  #Summary
res |> summarise_lognormal()

## get shape and rate samples from the posterior
res_df <- res |>
  as_draws_df() |>
  filter(.draw %in% sample(.draw, 100)) # sample 100 draws

## find the value (x) that includes 99% of the cumulative density
max_x <- max(qlnorm(0.99, meanlog = res_df$meanlog, sdlog = res_df$sdlog))

## calculate density on grid of x values
x <- seq(0, max_x, length.out = 100)
res_df <- res_df |>
  crossing(x = x) |> ## add grid to data frame
  mutate(density = dlnorm(x, meanlog, sdlog))

## plot
ggplot(res_df, aes(x = x, y = density, group = .draw)) +
  geom_line(alpha = 0.3)









##### Challenge #1 #####

df_small <- df_onset_to_hosp |> slice_sample(n = 100)

## Use the data to sample from the model posterior
res <- nfidd_sample(
  mod,
  data = list(
    n = nrow(df_small),
    y = df_small$onset_to_hosp
  )
)
res$summary()


#Summary
res |> summarise_lognormal()

## get shape and rate samples from the posterior
res_df <- res |>
  as_draws_df() |>
  filter(.draw %in% sample(.draw, 100)) # sample 100 draws

## find the value (x) that includes 99% of the cumulative density
max_x <- max(qlnorm(0.99, meanlog = res_df$meanlog, sdlog = res_df$sdlog))

## calculate density on grid of x values
x <- seq(0, max_x, length.out = 100)
res_df <- res_df |>
  crossing(x = x) |> ## add grid to data frame
  mutate(density = dlnorm(x, meanlog, sdlog))

## plot
ggplot(res_df, aes(x = x, y = density, group = .draw)) +
  geom_line(alpha = 0.3)


#### When Sample Size decreases, Variability of the estimated model increases








##### Challenge #2 #####

df_long <- add_delays(infection_times, 
                      hosp_params = list(meanlog = 2.5, sdlog = 0.5))

# convert our data frame to long format
dfl <- df_long |>
  pivot_longer(
    cols = c(infection_time, onset_time, hosp_time),
    names_to = "type", values_to = "time"
  ) |> 
  mutate(type = ordered(type, 
                        levels = c("infection_time", "onset_time", "hosp_time"), 
                        labels = c("Infections", "Symptom onsets", "Hospitalisations")
  ))

# plot the three event distributions
ggplot(dfl, aes(x = time)) +
  geom_histogram(position = "dodge", binwidth = 1) +
  facet_wrap(~ type, ncol = 1) +
  xlab("Time (in days)") +
  ylab("Count")



#Creating the Stan Model
mod <- nfidd_cmdstan_model("lognormal")
mod


## Specify the time from onset to hospitalization
df_onset_to_hosp <- df_long |>
  mutate(onset_to_hosp = hosp_time - onset_time) |> 
  # exclude infections that didn't result in hospitalization
  drop_na(onset_to_hosp)


## Use the data to sample from the model posterior
res <- nfidd_sample(
  mod,
  data = list(
    n = nrow(df_onset_to_hosp),
    y = df_onset_to_hosp$onset_to_hosp
  )
)
res$summary()


#Summary
res |> summarise_lognormal()

## get shape and rate samples from the posterior
res_df <- res |>
  as_draws_df() |>
  filter(.draw %in% sample(.draw, 100)) # sample 100 draws

## find the value (x) that includes 99% of the cumulative density
max_x <- max(qlnorm(0.99, meanlog = res_df$meanlog, sdlog = res_df$sdlog))

## calculate density on grid of x values
x <- seq(0, max_x, length.out = 100)
res_df <- res_df |>
  crossing(x = x) |> ## add grid to data frame
  mutate(density = dlnorm(x, meanlog, sdlog))

## plot
ggplot(res_df, aes(x = x, y = density, group = .draw)) +
  geom_line(alpha = 0.3)


#### Looks similar to original --> Data peaks at delay of 10 instead of 5







##### Challenge #3A #####

df_variable <- add_delays(infection_times,
hosp_params = list(meanlog = 1.75, sdlog = 1.0))

# convert our data frame to long format
dfl <- df_variable |>
  pivot_longer(
    cols = c(infection_time, onset_time, hosp_time),
    names_to = "type", values_to = "time"
  ) |> 
  mutate(type = ordered(type, 
                        levels = c("infection_time", "onset_time", "hosp_time"), 
                        labels = c("Infections", "Symptom onsets", "Hospitalisations")
  ))

# plot the three event distributions
ggplot(dfl, aes(x = time)) +
  geom_histogram(position = "dodge", binwidth = 1) +
  facet_wrap(~ type, ncol = 1) +
  xlab("Time (in days)") +
  ylab("Count")



#Creating the Stan Model
mod <- nfidd_cmdstan_model("lognormal")
mod


## Specify the time from onset to hospitalization
df_onset_to_hosp <- df_variable |>
  mutate(onset_to_hosp = hosp_time - onset_time) |> 
  # exclude infections that didn't result in hospitalization
  drop_na(onset_to_hosp)


## Use the data to sample from the model posterior
res <- nfidd_sample(
  mod,
  data = list(
    n = nrow(df_onset_to_hosp),
    y = df_onset_to_hosp$onset_to_hosp
  )
)
res$summary()


#Summary
res |> summarise_lognormal()

## get shape and rate samples from the posterior
res_df <- res |>
  as_draws_df() |>
  filter(.draw %in% sample(.draw, 100)) # sample 100 draws

## find the value (x) that includes 99% of the cumulative density
max_x <- max(qlnorm(0.99, meanlog = res_df$meanlog, sdlog = res_df$sdlog))

## calculate density on grid of x values
x <- seq(0, max_x, length.out = 100)
res_df <- res_df |>
  crossing(x = x) |> ## add grid to data frame
  mutate(density = dlnorm(x, meanlog, sdlog))

## plot
ggplot(res_df, aes(x = x, y = density, group = .draw)) +
  geom_line(alpha = 0.3)



##### Challenge #3B #####

df_precise <- add_delays(infection_times,
                         hosp_params = list(meanlog = 1.75, sdlog = 0.2))

# convert our data frame to long format
dfl <- df_precise |>
  pivot_longer(
    cols = c(infection_time, onset_time, hosp_time),
    names_to = "type", values_to = "time"
  ) |> 
  mutate(type = ordered(type, 
                        levels = c("infection_time", "onset_time", "hosp_time"), 
                        labels = c("Infections", "Symptom onsets", "Hospitalisations")
  ))

# plot the three event distributions
ggplot(dfl, aes(x = time)) +
  geom_histogram(position = "dodge", binwidth = 1) +
  facet_wrap(~ type, ncol = 1) +
  xlab("Time (in days)") +
  ylab("Count")



#Creating the Stan Model
mod <- nfidd_cmdstan_model("lognormal")
mod


## Specify the time from onset to hospitalization
df_onset_to_hosp <- df_precise |>
  mutate(onset_to_hosp = hosp_time - onset_time) |> 
  # exclude infections that didn't result in hospitalization
  drop_na(onset_to_hosp)


## Use the data to sample from the model posterior
res <- nfidd_sample(
  mod,
  data = list(
    n = nrow(df_onset_to_hosp),
    y = df_onset_to_hosp$onset_to_hosp
  )
)
res$summary()


#Summary
res |> summarise_lognormal()

## get shape and rate samples from the posterior
res_df <- res |>
  as_draws_df() |>
  filter(.draw %in% sample(.draw, 100)) # sample 100 draws

## find the value (x) that includes 99% of the cumulative density
max_x <- max(qlnorm(0.99, meanlog = res_df$meanlog, sdlog = res_df$sdlog))

## calculate density on grid of x values
x <- seq(0, max_x, length.out = 100)
res_df <- res_df |>
  crossing(x = x) |> ## add grid to data frame
  mutate(density = dlnorm(x, meanlog, sdlog))

## plot
ggplot(res_df, aes(x = x, y = density, group = .draw)) +
  geom_line(alpha = 0.3)




#### When delays are more variable, spread of posterior curves is wider
#### When delays are less variable, spread of posterior curves in narrower

##### Challenge #4 #####

df_gamma <- add_delays(infection_times,
                       hosp_fun = rgamma,
                       hosp_params = list(shape = 2, rate = 0.3))

# convert our data frame to long format
dfl <- df_precise |>
  pivot_longer(
    cols = c(infection_time, onset_time, hosp_time),
    names_to = "type", values_to = "time"
  ) |> 
  mutate(type = ordered(type, 
                        levels = c("infection_time", "onset_time", "hosp_time"), 
                        labels = c("Infections", "Symptom onsets", "Hospitalisations")
  ))

# plot the three event distributions
ggplot(dfl, aes(x = time)) +
  geom_histogram(position = "dodge", binwidth = 1) +
  facet_wrap(~ type, ncol = 1) +
  xlab("Time (in days)") +
  ylab("Count")



#Creating the Stan Model
mod <- nfidd_cmdstan_model("gamma")
mod


## Specify the time from onset to hospitalization
df_onset_to_hosp <- df_precise |>
  mutate(onset_to_hosp = hosp_time - onset_time) |> 
  # exclude infections that didn't result in hospitalization
  drop_na(onset_to_hosp)


## Use the data to sample from the model posterior
res <- nfidd_sample(
  mod,
  data = list(
    n = nrow(df_onset_to_hosp),
    y = df_onset_to_hosp$onset_to_hosp
  )
)

#hmmmm.... something gets funky over here.
res$summary()


#Summary
res |> summarise_lognormal()

## get shape and rate samples from the posterior
res_df <- res |>
  as_draws_df() |>
  filter(.draw %in% sample(.draw, 100)) # sample 100 draws

## find the value (x) that includes 99% of the cumulative density
max_x <- max(qlnorm(0.99, meanlog = res_df$meanlog, sdlog = res_df$sdlog))

## calculate density on grid of x values
x <- seq(0, max_x, length.out = 100)
res_df <- res_df |>
  crossing(x = x) |> ## add grid to data frame
  mutate(density = dlnorm(x, meanlog, sdlog))

## plot
ggplot(res_df, aes(x = x, y = density, group = .draw)) +
  geom_line(alpha = 0.3)

