  
  # 2025-07-16
  # Amy Moore
  # SISMID - NFIDD
  # Session 2 - Bias in Delay Distributions
  
  #Load in Packages

library("nfidd")
library("ggplot2")
library("dplyr")
library("tidyr")
library("purrr")
library("lubridate")
library("tidybayes")

  #Set Seed

set.seed(123)
options(cmdstanr_print_line_numbers = TRUE)


#Load in Data

data(infection_times)
head(infection_times)


  #Creating Delay Distributions
df <- add_delays(
  infection_times, 
  hosp_params = list(meanlog = 1.0, sdlog = 0.5)
)
head(df)


# Use the floor() function to round down to integers
df_dates <- df |>
  mutate(
    infection_time = floor(infection_time),
    onset_time = floor(onset_time),
    hosp_time = floor(hosp_time)
  )
head(df_dates)


  #Naive Approach --> ignore censoring

df_dates <- df_dates |>
  mutate(
    incubation_period = onset_time - infection_time,
    onset_to_hosp = hosp_time - onset_time
  )

  #Comparing non-censored to censored data summary
summary(df$onset_time - df$infection_time)
summary(df_dates$incubation_period)


  #Fitting the Delay Distribution Model on the Naive Approach Data

#Creating the Stan Model
mod <- nfidd_cmdstan_model("lognormal")
mod


## Specify the time from onset to hospitalization
df_onset_to_hosp <- df_dates |>
  mutate(onset_to_hosp = hosp_time - onset_time) |> 
  # exclude infections that didn't result in hospitalization
  drop_na(onset_to_hosp)


## Use the data to sample from the model posterior
res <- nfidd_sample(
  mod,
  data = list(
    n = nrow(df_onset_to_hosp),
    y = df_onset_to_hosp$onset_to_hosp + 0.01
      #HINT: add the 0.01 to prevent 0 day delays
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




##### Visualizing types of censoring #####

library(ggplot2)
library(dplyr)

# Set parameters for our example
n <- 1e4
meanlog <- 1.0  # Mean ~3 days - try changing this
sdlog <- 0.5
obs_time <- 15

# Generate true delay distribution
true_delays <- rlnorm(n, meanlog = meanlog, sdlog = sdlog)

# Primary censoring only - add uniform uncertainty to primary event timing
primary_censored <- true_delays + runif(n, 0, 1)

# Secondary censoring only (common mistake) - just discretise
secondary_only <- floor(true_delays)

# Double censoring - discretise the primary censored delays
double_censored <- floor(primary_censored)

# Filter to a reasonable range and create PMF for discrete data
keep_range <- function(x) x[x <= obs_time]
primary_filtered <- keep_range(primary_censored)
secondary_filtered <- keep_range(secondary_only)
double_filtered <- keep_range(double_censored)

# Create PMF for discrete data
secondary_pmf <- table(secondary_filtered) / length(secondary_filtered)
double_pmf <- table(double_filtered) / length(double_filtered)

# Create the comparison plot
ggplot() +
  # True distribution (black line)
  geom_function(
    fun = dlnorm,
    args = list(meanlog = meanlog, sdlog = sdlog),
    color = "#252525",
    linewidth = 1.2
  ) +
  # Primary censoring (continuous, blue density)
  geom_density(
    data = data.frame(x = primary_filtered),
    aes(x = x),
    fill = "#4292C6",
    col = "#252525",
    alpha = 0.6
  ) +
  # Secondary censoring only (discrete, coral bars)
  geom_col(
    data = data.frame(
      x = as.numeric(names(secondary_pmf)),
      y = as.numeric(secondary_pmf)
    ),
    aes(x = x, y = y),
    fill = "#E31A1C",
    col = "#252525",
    alpha = 0.6,
    width = 0.9
  ) +
  # Double censoring (discrete, green bars)
  geom_col(
    data = data.frame(
      x = as.numeric(names(double_pmf)),
      y = as.numeric(double_pmf)
    ),
    aes(x = x, y = y),
    fill = "#20b986",
    col = "#252525",
    alpha = 0.4,
    width = 0.9
  ) +
  labs(
    title = "Comparison of Different Censoring Approaches",
    x = "Delay (days)",
    y = "Density / Probability Mass",
    caption = paste0(
      "Black line: True log-normal distribution\n",
      "Blue density: Primary censoring (continuous)\n",
      "Red bars: Secondary censoring only (common mistake)\n",
      "Green bars: Double censoring (both components)"
    )
  ) +
  scale_x_continuous(limits = c(0, 15)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 0)
  )


  ##### Back to the main lecture #####

cmod <- nfidd_cmdstan_model("censored-delay-model")
cmod

cres <- nfidd_sample(cmod,
                     data = list(
                       n = nrow(na.omit(df_dates)),
                       onset_to_hosp = na.omit(df_dates)$onset_to_hosp + 0.01
                     )
)


## Use the data to sample from the model posterior
res <- nfidd_sample(
  mod,
  data = list(
    n = nrow(df_onset_to_hosp),
    y = df_onset_to_hosp$onset_to_hosp + 0.01
    #HINT: add the 0.01 to prevent 0 day delays
  )
)
res$summary()



##### Truncation #####

# Resimulate with longer delays for truncation demonstration
set.seed(123)
df <- add_delays(infection_times, 
                 hosp_params = list(meanlog = 1.75, sdlog = 0.5))

library(ggplot2)
library(dplyr)

# Set parameters for the truncation demo
set.seed(890)
n <- 5000
meanlog <- 1.75
sdlog <- 0.5

# Generate delays
true_delays <- rlnorm(n, meanlog = meanlog, sdlog = sdlog)

# Simulate truncation at different time points
truncation_times <- c(6, 10, 15)
truncated_data <- map_dfr(truncation_times, function(t) {
  truncated_delays <- true_delays[true_delays <= t]
  data.frame(
    delay = truncated_delays,
    truncation = paste("Truncated at", t, "days"),
    mean_delay = mean(truncated_delays)
  )
})

# Create comparison plot
ggplot() +
  # True distribution (black line)
  geom_function(
    fun = dlnorm,
    args = list(meanlog = meanlog, sdlog = sdlog),
    color = "black",
    linewidth = 1.2,
    xlim = c(0, 20)
  ) +
  # Empirical densities for truncated data
  geom_density(
    data = truncated_data,
    aes(x = delay, fill = truncation),
    alpha = 0.6
  ) +
  # Mean lines
  geom_vline(
    data = truncated_data %>% distinct(truncation, mean_delay),
    aes(xintercept = mean_delay, color = truncation),
    linetype = "dashed",
    linewidth = 1
  ) +
  scale_fill_manual(values = c("Truncated at 6 days" = "#E31A1C", 
                               "Truncated at 10 days" = "#FF7F00", 
                               "Truncated at 15 days" = "#1F78B4")) +
  scale_color_manual(values = c("Truncated at 6 days" = "#E31A1C", 
                                "Truncated at 10 days" = "#FF7F00", 
                                "Truncated at 15 days" = "#1F78B4")) +
  labs(
    title = "Effect of Right Truncation on Delay Distributions",
    x = "Delay (days)",
    y = "Density",
    fill = "Truncated data",
    color = "Mean delay",
    caption = "Black line: True log-normal distribution\nColored densities: Truncated data\nDashed lines: Mean delays"
  ) +
  xlim(0, 20) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  ) +
  guides(colour = "none")








df_realtime <- df |>
  mutate(onset_to_hosp = hosp_time - onset_time) |>
  filter(hosp_time <= 70)

# truncated mean delay
mean(df_realtime$onset_to_hosp)

# compare with the mean delay over the full outbreak
mean(df$hosp_time - df$onset_time, na.rm=TRUE)

res <- nfidd_sample(mod,
                    data = list(
                      n = nrow(na.omit(df_realtime)),
                      y = na.omit(df_realtime)$onset_to_hosp
                    )
)
res
res |>
  summarise_lognormal()



tmod <- nfidd_cmdstan_model("truncated-delay-model")
tmod

tres <- nfidd_sample(tmod,
                     data = list(
                       n = nrow(df_realtime),
                       onset_to_hosp = df_realtime$onset_to_hosp, 
                       time_since_onset = 70 - df_realtime$onset_time
                     )
)
tres
tres |>
  summarise_lognormal()




  ##### Varying Time of truncation #####

  #Day 30

df_realtime <- df |>
  mutate(onset_to_hosp = hosp_time - onset_time) |>
  filter(hosp_time <= 30)

# truncated mean delay
mean(df_realtime$onset_to_hosp)

# compare with the mean delay over the full outbreak
mean(df$hosp_time - df$onset_time, na.rm=TRUE)

res <- nfidd_sample(mod,
                    data = list(
                      n = nrow(na.omit(df_realtime)),
                      y = na.omit(df_realtime)$onset_to_hosp
                    )
)
res
res |>
  summarise_lognormal()



tres <- nfidd_sample(tmod,
                     data = list(
                       n = nrow(df_realtime),
                       onset_to_hosp = df_realtime$onset_to_hosp, 
                       time_since_onset = 30 - df_realtime$onset_time
                     )
)
tres
tres |>
  summarise_lognormal()


# --> That blows up the maximum predictions to be large means (yikes!)

#Day 100

df_realtime <- df |>
  mutate(onset_to_hosp = hosp_time - onset_time) |>
  filter(hosp_time <= 100)

# truncated mean delay
mean(df_realtime$onset_to_hosp)

# compare with the mean delay over the full outbreak
mean(df$hosp_time - df$onset_time, na.rm=TRUE)

res <- nfidd_sample(mod,
                    data = list(
                      n = nrow(na.omit(df_realtime)),
                      y = na.omit(df_realtime)$onset_to_hosp
                    )
)
res
res |>
  summarise_lognormal()



tres <- nfidd_sample(tmod,
                     data = list(
                       n = nrow(df_realtime),
                       onset_to_hosp = df_realtime$onset_to_hosp, 
                       time_since_onset = 100 - df_realtime$onset_time
                     )
)
tres
tres |>
  summarise_lognormal()
