  
  # 2025-07-17
  # Amy Moore
  # SISMID - NFIDD
  # Session 8 - Forecasting Improvements
  
  #Load Packages
library("nfidd")
library("fable")
library("dplyr")
library("tidyr")
library("ggplot2")
library("epidatr")

  #Set Up
set.seed(123)


  #load data
data(flu_data)


  #plot data
feasts::gg_tsdisplay(flu_data, y = wili, plot_type='partial', lag_max = 104)

  #
cor(flu_data$wili, dplyr::lag(flu_data$wili, 1), use = "complete.obs")

cor(flu_data$wili, dplyr::lag(flu_data$wili, 12), use = "complete.obs")






##### Fitting some ARIMA models #####

  #ARIMA(1,0,0) -> Just auto-correlation

fit_arima100 <- flu_data |>
  model(ARIMA(wili ~ pdq(1,0,0)))
report(fit_arima100)

forecast(fit_arima100, h=40) |>
  autoplot(flu_data |> filter(epiweek >= as.Date("2015-09-01"))) +
  labs(title = "WILI, US level",
       y="% of visits")

  # ARIMA(1,1,0) --> 1 auto-correlation, 1 differencing

fit_arima110 <- flu_data |>
  model(ARIMA(wili ~ pdq(1,1,0)))
report(fit_arima110)

forecast(fit_arima110, h=40) |>
  autoplot(flu_data |> filter(epiweek >= as.Date("2015-09-01"))) +
  labs(title = "WILI, US level",
       y="% of visits")

  #ARIMA(0,0,1) --> just moving average

fit_arima001 <- flu_data |>
  model(ARIMA(wili ~ pdq(0,0,1)))
report(fit_arima001)

forecast(fit_arima001, h=40) |>
  autoplot(flu_data |> filter(epiweek >= as.Date("2015-09-01"))) +
  labs(title = "WILI, US level",
       y="% of visits")




  #A few ARIMA models that fit this data well (ish)
#original data
feasts::gg_tsdisplay(flu_data, y = wili, plot_type='partial', lag_max = 104)

#1 level of distancing
feasts::gg_tsdisplay(flu_data, y = difference(wili), plot_type='partial', lag_max = 104)


#adding a transformation
fourth_root <- function(x) x^0.25
inv_fourth_root <- function(x) x^4

my_fourth_root <- new_transformation(fourth_root, inv_fourth_root)

autoplot(flu_data, .vars = fourth_root(wili))

  #fit model
fits <- flu_data |>
  model(
    arima200 = ARIMA(wili ~ pdq(2,0,0)),
    arima110 = ARIMA(wili ~ pdq(1,1,0)),
    arima200_trans = ARIMA(my_fourth_root(wili) ~ pdq(2,0,0)),
    arima110_trans = ARIMA(my_fourth_root(wili) ~ pdq(1,1,0))
  )
forecast(fits, h=40) |>
  autoplot(flu_data |> filter(epiweek >= as.Date("2015-09-01"))) +
  facet_wrap(~.model) +
  labs(title = "WILI, US level",
       y="% of visits")









  ##### Doing something about seasonality #####
fits <- flu_data |>
  model(
    arima200_trans = ARIMA(my_fourth_root(wili) ~ pdq(2,0,0)),
    arima110_trans = ARIMA(my_fourth_root(wili) ~ pdq(1,1,0)),
    arima200_fourier = ARIMA(my_fourth_root(wili) ~ pdq(2,0,0) + fourier(period = "year", K=3)))
forecast(fits, h=40) |>
  autoplot(flu_data |> filter(epiweek >= as.Date("2015-09-01"))) +
  facet_wrap(~.model) +
  labs(title = "WILI, US level",
       y="% of visits")

