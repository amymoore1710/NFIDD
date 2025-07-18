
  # 2025-07-17
  # Amy Moore
  # SISMID - NFIDD
  # Session 10 - Forecasting Ensembles
  
  #Load Packages
library('nfidd')
library("fable")
library("dplyr")
library("ggplot2")
library("epidatr")
library("hubUtils")
library("hubEvals")
library("hubEnsembles")

  #Set Up
theme_set(theme_bw())
set.seed(42) ## for Douglas Adams!

  #Load Data
data(flu_data)

  #New Function for transformation
fourth_root <- function(x) x^0.25
inv_fourth_root <- function(x) x^4
my_fourth_root <- new_transformation(fourth_root, inv_fourth_root)



## make the time-series cross validation splits
flu_data_tscv <- flu_data |> 
  filter(epiweek <= as.Date("2018-07-01")) |> 
  tsibble::stretch_tsibble(
    .init = 732, 
    .step = 4,
    .id = ".split"
  )

## generate the forecasts
cv_models <- flu_data_tscv |> 
  model(
    rw = RW(my_fourth_root(wili)),
    ar2 = ARIMA(my_fourth_root(wili) ~ pdq(2,0,0)),
    fourier = ARIMA(my_fourth_root(wili) ~ pdq(0,0,0) + fourier(period = "year", K=3)),
    fourier_ar2 = ARIMA(my_fourth_root(wili) ~ pdq(2,0,0) + fourier(period = "year", K=3))
  ) 

cv_forecasts <- cv_models |> 
  forecast(h = 8) |> 
  ## the following 3 lines of code ensure that there is a horizon variable in the forecast data
  group_by(.split, .model) |> 
  mutate(h = row_number()) |> 
  ungroup() |> 
  ## this ensures that the output is a fable object
  as_fable(response = "wili", distribution = wili)


  #Plotting the TSCV models

autoplot(
  cv_forecasts,
  flu_data_tscv |> filter(epiweek >= as.Date("2017-07-01")),
  alpha = 0.5
) +
  facet_wrap(~.split) +
  theme(legend.position = "bottom")




  #generate sample forecasts
sampled_forecasts <- generate(
  cv_models, 
  h=8, 
  times=100) |> 
  group_by(.split, region, .model, .rep) |> 
  mutate(h = row_number())
sampled_forecasts

  #reformat to hubverse style
sampled_forecasts_hub <- sampled_forecasts |> 
  rename(value = .sim,
         model_id = .model,
         target_epiweek = epiweek) |> 
  mutate(output_type = "sample",
         output_type_id = stringr::str_c(region, stringr::str_pad(.rep, width = 3, pad = "0")),
         forecast_date = target_epiweek - h*7L) |> 
  ungroup() |> 
  as_tibble() |> 
  select(model_id,
         forecast_date, 
         target_epiweek,
         h,
         output_type,
         output_type_id,
         value)
sampled_forecasts_hub


  #visualize the sampled forecasts
sampled_forecasts_hub |> 
  filter(forecast_date == as.Date("2018-01-14"),
         model_id == "fourier_ar2") |> 
  ggplot() +
  geom_line(aes(x = target_epiweek, y = value, group = output_type_id), 
            color = "red",
            alpha = 0.2) +
  geom_line(data = flu_data |> filter(epiweek >= as.Date("2017-09-01")), 
            aes(x=epiweek, y=wili))






  ##### Creating the ensemble #####

#from samples to quantiles
quantiles_to_save <- c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)
quantile_forecasts <- sampled_forecasts_hub |>
  hubUtils::convert_output_type(
    to = list(quantile = quantiles_to_save)
  )
quantile_forecasts

































