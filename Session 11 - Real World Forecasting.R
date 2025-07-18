
  # 2025-07-18
  # Amy Moore
  # SISMID - NFIDD
  # Session 11 - Real World Forecasting

  #Load Packages
library("nfidd")
library("dplyr")
library("ggplot2")
library("hubData")
library("hubUtils")
library("hubEvals")
library("hubEnsembles")

  #Set Up
theme_set(theme_bw())
set.seed(5050) # for Shohei Ohtani!


  #load Data
data(covid_forecasts)
data(covid_time_series)


unique_per_column <- covid_forecasts  |> 
  select(-value) |> 
  purrr::map(~ sort(unique(.x)))
unique_per_column

  #Visualizing Hub Models
covid_forecasts |> 
  filter(reference_date == "2025-02-15", abbreviation == "GA") |> 
  hubVis::plot_step_ahead_model_output(
    target_data = covid_time_series |> 
      filter(as_of == as.Date("2025-07-09"),
             abbreviation == "GA"),
    use_median_as_point = TRUE,
    x_col_name = "target_end_date",
    x_target_col_name = "date",
    pal_color = "Set3", 
    title = "Weekly hospitalizations due to COVID-19 (data and forecasts)"
  )

  #Summary of info from each model
covid_forecasts |> 
  group_by(model_id) |> 
  summarise(
    n_submissions = n_distinct(reference_date),
    n_rows = n(),
    n_horizons = n_distinct(horizon),
    n_locations = n_distinct(location)
  ) |> 
  arrange(-n_submissions) |> 
  knitr::kable()

  #visualizing the info from each model
covid_forecasts |> 
  group_by(model_id, reference_date) |> 
  summarise(
    n_rows = n(),
    n_locations = n_distinct(location)
  ) |> 
  ungroup() |> 
  mutate(model_id = reorder(model_id, n_rows, FUN = sum)) |> 
  ggplot() +
  geom_tile(aes(x = reference_date, y = model_id, fill = n_locations))



  #building oracle outlook data
covid_oracle_output <- covid_time_series |> 
  filter(as_of == as.Date("2025-07-09")) |> 
  select(target,
         location,
         target_end_date = date,
         oracle_value = observation) 

#visual evaluation
covid_forecasts |> 
  filter(output_type_id == 0.5,
         abbreviation == "GA") |> 
  ggplot() +
  geom_line(aes(x = target_end_date, 
                y = value, 
                group = interaction(reference_date, model_id)),
            alpha = 0.5, 
            color = "red") +
  geom_line(data = filter(covid_time_series, 
                          abbreviation == "GA",
                          as_of == "2025-07-09"), 
            aes(x = date,
                y = observation))+
  ggtitle("Forecasts and data for Georgia") +
  ylab("incident hospital admissions") +
  xlab(NULL)


  #looking just at one state - Virginia
covid_forecasts |> 
  filter(output_type_id == 0.5,
         abbreviation == "VA") |> 
  ggplot() +
  geom_line(aes(x = target_end_date, 
                y = value, 
                group = interaction(reference_date, model_id)),
            alpha = 0.5, 
            color = "red") +
  geom_line(data = filter(covid_time_series, 
                          abbreviation == "VA",
                          as_of == "2025-07-09"), 
            aes(x = date,
                y = observation))+
  ggtitle("Forecasts and data for Virginia") +
  ylab("incident hospital admissions") +
  xlab(NULL)


  #Compute scores
scores <- score_model_out(
  covid_forecasts, 
  covid_oracle_output
)
scores |> 
  arrange(wis) |> 
  knitr::kable(digits = 3)


  #Calculate Relative Scores
score_model_out(
  covid_forecasts, 
  covid_oracle_output,
  metrics = "wis",
  relative_metrics = "wis"
) |> 
  arrange(wis_relative_skill)


  #Remove models that are not 60% complete
## 33 dates, 4 horizons, 53 locations, 60%
threshold_60pct <- 33 * 4 * 53 * 0.6
model_subset <- covid_forecasts |> 
  filter(output_type_id == 0.5) |> 
  group_by(model_id) |> 
  summarize(targets = n()) |> 
  filter(targets > threshold_60pct) |> 
  pull(model_id)

covid_forecasts |> 
  filter(model_id %in% model_subset) |> 
  score_model_out(
    covid_oracle_output,
    metrics = "wis",
    relative_metrics = "wis"
  ) |> 
  arrange(wis_relative_skill)



  #subset to just the state level
covid_forecasts |> 
  filter(model_id %in% model_subset,
         location != "US") |> 
  score_model_out(
    covid_oracle_output,
    metrics = "wis",
    relative_metrics = "wis"
  ) |> 
  arrange(wis_relative_skill)




# Log Calibrated
covid_forecasts |> 
  filter(model_id %in% model_subset) |> 
  mutate(value = log(value+1)) |> 
  score_model_out(
    covid_oracle_output |> 
      mutate(oracle_value = log(oracle_value+1)),
    metrics = "wis",
    relative_metrics = "wis"
  ) |> 
  arrange(wis_relative_skill)



#prediction intervals
covid_forecasts |> 
  filter(model_id %in% model_subset) |> 
  score_model_out(
    covid_oracle_output,
    metrics = c("wis", "interval_coverage_50", "interval_coverage_90")
  ) |> 
  arrange(interval_coverage_90) |> 
  knitr::kable(digits = 3)


#By Horizon
scores_by_horizon <- covid_forecasts |> 
  filter(model_id %in% model_subset,
         location != "US") |> 
  score_model_out(
    covid_oracle_output,
    metrics = "wis",
    relative_metrics = "wis",
    by = c("model_id", "horizon")
  ) 
p <- ggplot(scores_by_horizon)+
  geom_line(aes(x = horizon, y = wis, color = model_id))
plotly::ggplotly(p)


#by reference date
scores_by_reference_date <- covid_forecasts |> 
  filter(model_id %in% model_subset,
         location != "US") |> 
  score_model_out(
    covid_oracle_output,
    metrics = "wis",
    relative_metrics = "wis",
    by = c("model_id", "reference_date")
  ) |> 
  tidyr::pivot_longer(cols = c(wis, wis_relative_skill), 
                      names_to = "metric",
                      values_to = "score")

p <- ggplot(scores_by_reference_date)+
  geom_line(aes(x = reference_date, y = score, color = model_id)) +
  facet_grid(metric~., scales = "free_y")
plotly::ggplotly(p)






  ##### DIY Ensemble Models #####
# Get a list of unique reference dates, removing the last one
reference_dates <- covid_forecasts |> 
  filter(reference_date != as.Date("2025-07-12")) |> 
  pull(reference_date) |> 
  unique()

# Map over each date and apply linear_pool separately
lop_unweighted_all <- purrr::map_dfr(
  reference_dates, 
  function(date) {
    covid_forecasts |> 
      filter(model_id != "CovidHub-ensemble", 
             reference_date == date) |> 
      hubEnsembles::linear_pool(model_id = "lop_unweighted_all")
  })


lop_unweighted_select <- purrr::map_dfr(
  reference_dates, 
  function(date) {
    covid_forecasts |> 
      filter(model_id %in% model_subset,
             model_id != "CovidHub-ensemble", 
             reference_date == date) |> 
      hubEnsembles::linear_pool(model_id = "lop_unweighted_select")
  })
