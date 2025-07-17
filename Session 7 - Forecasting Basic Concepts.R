
  # 2025-07-17
  # Amy Moore
  # SISMID - NFIDD
  # Session 7 - Forecasting Basic Concepts

  #Load Packages
library("nfidd")
library("dplyr")
library("ggplot2")
library("epidatr")
library("feasts")
library("fable")

  #Set Up
theme_set(theme_bw())
set.seed(17)


  #Load in Dataset
data(flu_data)
flu_data <- flu_data |>
  filter(epiweek <= "2017-08-27") # Filter to start of 2017/2018 season

head(flu_data)


ggplot(flu_data) +
  geom_path(aes(x=epiweek, y=wili))

  #Time-series Plot
autoplot(flu_data, .vars = wili) +
  labs(title = "National-level ILI data from the US CDC",
       y = "weighted ILI (% of visits)")

  #Fit ARIMA model
fit <- flu_data |>
  model( arima200 = ARIMA(wili ~ pdq(2,0,0)) )

report(fit)

  #Forecast for next 25 weeks
first_forecast <- forecast(fit, h = 25) 

first_forecast |>
  autoplot(flu_data |> filter(epiweek >= as.Date("2016-09-01"))) +
  labs(title = "WILI, US level",
       y = "% of outpatient visits due to ILI")


