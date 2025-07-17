
  # 2025-07-17
  # Amy Moore
  # SISMID - NFIDD
  # Session 9 - Forecasting Evaluation
  
  #Load Packages
library("nfidd")
library("fable")
library("dplyr")
library("ggplot2")
library("epidatr")

  #Set Up
theme_set(theme_bw())
set.seed(42) # for Jackie Robinson!

  #Load Data
data(flu_data)

## remember to define the transformation
fourth_root <- function(x) x^0.25
inv_fourth_root <- function(x) x^4
my_fourth_root <- new_transformation(fourth_root, inv_fourth_root)

## here is the model
fit_fourier_ar2 <- flu_data |>
  filter(epiweek <= "2017-08-27") |> 
  model(
    arima200_fourier = ARIMA(my_fourth_root(wili) ~ pdq(2,0,0) + fourier(period = "year", K=3))
  )

## make the forecast
forecast_fourier_ar2 <- forecast(fit_fourier_ar2, h=40) 

## plot the forecast against all of the data from the 2017/2018 season
forecast_fourier_ar2 |>
  autoplot(flu_data |> filter(epiweek >= as.Date("2015-09-01"))) +
  facet_wrap(~.model) +
  labs(title = "WILI, US level",
       y="% of visits")

  #Measuring Accuracy on Testing Data
accuracy(forecast_fourier_ar2, 
         data = flu_data, 
         measures = list(mae = MAE, rmse = RMSE, crps = CRPS))

  #Measuring Accuracy on Training Data
accuracy(fit_fourier_ar2, measures = list(mae = MAE, rmse = RMSE))


  #Measuring Accuracy for lots of forecasts
metrics_fourier_ar2 <- forecast_fourier_ar2 |> 
  mutate(h = row_number()) |> 
  accuracy(data = flu_data, 
           measures = list(mae = MAE, rmse = RMSE, crps = CRPS),
           by = c(".model", "h"))
metrics_fourier_ar2

metrics_fourier_ar2 |> 
  tidyr::pivot_longer(
    cols = c("mae", "rmse", "crps"),
    names_to = "metric"
  ) |> 
  ggplot() +
  geom_line(aes(x = h, y = value, color = metric)) +
  facet_wrap(.~metric)




  ##### Time-series CV #####

flu_data_tscv <- flu_data |> 
  filter(epiweek <= as.Date("2018-06-01")) |> 
  tsibble::stretch_tsibble(
    .init = 732, 
    .step = 4,
    .id = ".split"
  )
flu_data_tscv


cv_forecasts <- flu_data_tscv |> 
  model(
    rw = RW(my_fourth_root(wili)),
    ar2 = ARIMA(my_fourth_root(wili) ~ pdq(2,0,0)),
    fourier = ARIMA(my_fourth_root(wili) ~ pdq(0,0,0) + fourier(period = "year", K=3)),
    fourier_ar2 = ARIMA(my_fourth_root(wili) ~ pdq(2,0,0) + fourier(period = "year", K=3))
  ) |> 
  forecast(h = 8) |> 
  ## the following 3 lines of code ensure that there is a horizon variable in the forecast data
  group_by(.split, .model) |> 
  mutate(h = row_number()) |> 
  ungroup() |> 
  ## this ensures that the output is a table object
  as_fable(response = "wili", distribution = wili)

cv_forecasts |>
  filter(.split == 1) |> 
  tsibble::update_tsibble(key=c(.model, region)) |>  
  as_fable(response = "wili", distribution = wili) |> 
  autoplot(flu_data |> filter(epiweek >= as.Date("2016-09-01"))) +
  facet_wrap(~.model) +
  labs(title = "WILI, US level",
       y="% of visits")

cv_forecasts |>
  filter(.split == 4) |> 
  tsibble::update_tsibble(key=c(.model, region)) |>  
  as_fable(response = "wili", distribution = wili) |> 
  autoplot(flu_data |> filter(epiweek >= as.Date("2016-09-01"))) +
  facet_wrap(~.model) +
  labs(title = "WILI, US level",
       y="% of visits")

cv_forecasts |>
  filter(.split == 6) |> 
  tsibble::update_tsibble(key=c(.model, region)) |>  
  as_fable(response = "wili", distribution = wili) |> 
  autoplot(flu_data |> filter(epiweek >= as.Date("2016-09-01"))) +
  facet_wrap(~.model) +
  labs(title = "WILI, US level",
       y="% of visits")

cv_forecasts |>
  filter(.split == 8) |> 
  tsibble::update_tsibble(key=c(.model, region)) |>  
  as_fable(response = "wili", distribution = wili) |> 
  autoplot(flu_data |> filter(epiweek >= as.Date("2016-09-01"))) +
  facet_wrap(~.model) +
  labs(title = "WILI, US level",
       y="% of visits")



  ##### Evaluating Forecasts #####
cv_forecasts |> 
  accuracy(
    flu_data, 
    measures = list(mae = MAE, rmse = RMSE, crps=CRPS)
  ) |> 
  arrange(crps)


cv_forecasts |> 
  accuracy(
    flu_data, 
    by = c("h", ".model"), 
    measures = list(mae = MAE, rmse = RMSE, crps=CRPS)
  ) |> 
  ggplot(aes(x = h, y = crps, color = .model)) +
  geom_point() +
  geom_line()

scores_to_plot <- cv_forecasts |> 
  mutate(forecast_date = epiweek - h*7L) |> 
  accuracy(
    flu_data, 
    by = c("forecast_date", ".model"), 
    measures = list(mae = MAE, rmse = RMSE, crps=CRPS)
  ) 

p1 <- scores_to_plot |> 
  ggplot(aes(x = forecast_date, y = crps, color = .model)) +
  geom_point() +
  geom_line() +
  theme(legend.position = "bottom")

date_range <- range(scores_to_plot$forecast_date)
p2 <- flu_data |> 
  filter(epiweek >= date_range[1],
         epiweek <= date_range[2]) |> 
  autoplot(wili)

gridExtra::grid.arrange(p1, p2)
