
  # 2025-07-18
  # Amy Moore
  # SISMID
  # Session 12 - Local Hub Playground

  #Load Packages
library("nfidd")
library("dplyr")
library("ggplot2")
library("epidatr")
library("fable")
library("hubUtils")
library("hubEvals")
library("hubVis")
library("hubData")
library("hubEnsembles")
library("here")

  #Set Up
theme_set(theme_bw())
set.seed(406) # for Ted Williams



  #Loading In the Hub
hub_path <- here::here("sismid-ili-forecasting-sandbox")
hub_con <- connect_hub(hub_path)



  #Seeing Models in the Hub
hub_con |> 
  filter(model_id == "delphi-epicast",
         origin_date == "2015-11-14") |> 
  collect_hub()







  #Setting Location Info
locs <-  c("nat", 
           "hhs1", "hhs2", "hhs3", "hhs4", "hhs5", 
           "hhs6", "hhs7", "hhs8", "hhs9", "hhs10")
location_formal_names <- c("US National", paste("HHS Region", 1:10))
loc_df <- data.frame(region = locs, location = location_formal_names)


  #Load Saved Data
data(flu_data_hhs)



  #Building Simple Model for Multiple locations
fourth_root <- function(x) x^0.25
inv_fourth_root <- function(x) x^4
my_fourth_root <- new_transformation(fourth_root, inv_fourth_root)

fit_arima210 <- flu_data_hhs |>
  filter(origin_date <= "2015-10-24") |> 
  model(ARIMA(my_fourth_root(wili) ~ pdq(2,1,0)))
forecast(fit_arima210, h=4) |>
  autoplot(flu_data_hhs |> 
             filter(
               origin_date <= "2015-10-24",
               origin_date >= "2014-09-01"
             )) +
  facet_wrap(.~location, scales = "free_y") +
  labs(y = "% of visits due to ILI",
       x = NULL)

  #Assessing the model for just 1 location
fit_arima210 |> 
  filter(location == "HHS Region 1") |> 
  report(model)

  #Assessing the model for just 1 location part 2
fit_arima210 |> 
  filter(location == "HHS Region 9") |> 
  report(model)




  ##### Designing Model for validation - Example #####

  #look at all possible prediction dates
origin_dates <- hub_path |> 
  read_config("tasks") |> 
  get_round_ids()
origin_dates

  #define time-series cross-validation data set
flu_data_hhs_tscv_season1 <- flu_data_hhs |> 
  filter(
    origin_date <= "2016-05-07"  ## last 2015/2016 date
  ) |> 
  tsibble::stretch_tsibble(
    .init = 634, ## flu_data_hhs |> filter(location == "HHS Region 1", origin_date <= "2015-10-17") |> nrow()
    .step = 1,
    .id = ".split"
  )
flu_data_hhs_tscv_season1


  #Creating example model and generate predictions
quantile_levels <- c(0.01, 0.025, seq(0.05, 0.95, 0.05), 0.975, 0.99)

cv_forecasts_season1 <- 
  flu_data_hhs_tscv_season1 |> 
  model(
    arima210 = ARIMA(my_fourth_root(wili) ~ pdq(2,1,0))
  ) |> 
  generate(h = 4, times = 100, bootstrap = TRUE) |> 
  ## the following 3 lines of code ensure that there is a horizon variable in the forecast data
  group_by(.split, .rep, location, .model) |> 
  mutate(horizon = row_number()) |> 
  ungroup() |> 
  as_tibble() |> 
  ## make hubverse-friendly names
  rename(
    target_end_date = origin_date,
    value = .sim,
    model_id = .model
  ) |> 
  left_join(loc_df) |> 
  mutate(origin_date = target_end_date - horizon * 7L) |> 
  ## compute the quantiles
  group_by(model_id, location, origin_date, horizon, target_end_date) |> 
  reframe(tibble::enframe(quantile(value, quantile_levels), "quantile", "value")) |> 
  mutate(output_type_id = as.numeric(stringr::str_remove(quantile, "%"))/100, .keep = "unused",
         target = "ili perc",
         output_type = "quantile",
         model_id = "sismid-arima210")


  #plot one set of forecasts
cv_forecasts_season1 |> 
  filter(origin_date == "2015-12-19") |> 
  plot_step_ahead_model_output(
    flu_data_hhs |> 
      filter(origin_date >= "2015-10-01", 
             origin_date <= "2015-12-19") |> 
      rename(observation = wili),
    x_target_col_name = "origin_date",
    x_col_name = "target_end_date",
    use_median_as_point = TRUE,
    facet = "location",
    facet_scales = "free_y",
    facet_nrow=3
  )

  #Continuing with example model --> forecast for season 2

  #create time-series cross-validation data set
first_idx_season2 <- flu_data_hhs |> 
  filter(location == "HHS Region 1", 
         origin_date <= "2016-10-22") |> ## first 2016/2017 date
  nrow()
flu_data_hhs_tscv_season2 <- flu_data_hhs |> 
  filter(
    origin_date <= "2017-05-06"  ## last 2016/2017 date
  ) |> 
  tsibble::stretch_tsibble(
    .init = first_idx_season2,
    .step = 1,
    .id = ".split"
  )


  #generate season 2 predictions
cv_forecasts_season2 <- 
  flu_data_hhs_tscv_season2 |> 
  model(
    arima210 = ARIMA(my_fourth_root(wili) ~ pdq(2,1,0))
  ) |> 
  generate(h = 4, times = 100, bootstrap = TRUE) |> 
  ## the following 3 lines of code ensure that there is a horizon variable in the forecast data
  group_by(.split, .rep, location, .model) |> 
  mutate(horizon = row_number()) |> 
  ungroup() |> 
  as_tibble() |> 
  ## make hubverse-friendly names
  rename(
    target_end_date = origin_date,
    value = .sim,
    model_id = .model
  ) |> 
  left_join(loc_df) |> 
  mutate(origin_date = target_end_date - horizon * 7L) |> 
  ## compute the quantiles
  group_by(model_id, location, origin_date, horizon, target_end_date) |> 
  reframe(tibble::enframe(quantile(value, quantile_levels), "quantile", "value")) |> 
  mutate(output_type_id = as.numeric(stringr::str_remove(quantile, "%"))/100, .keep = "unused",
         target = "ili perc",
         output_type = "quantile",
         model_id = "sismid-arima210")

  #Save the model locally
this_model_id <- "sismid-arima210"

metadata_filepath <- file.path(
  hub_path,
  "model-metadata", 
  paste0(this_model_id, ".yml"))

my_text <- c("team_abbr: \"sismid\"", 
             "model_abbr: \"arima210\"", 
             "designated_model: true")

writeLines(my_text, metadata_filepath)

# Group the forecasts by task id variables
groups <- bind_rows(cv_forecasts_season1, cv_forecasts_season2) |> 
  group_by(model_id, target, location, origin_date, horizon, target_end_date) |> 
  group_split()

# Save each group as a separate CSV
for (i in seq_along(groups)) {
  group_df <- groups[[i]]
  this_model_id <- group_df$model_id[1]
  this_origin_date <- group_df$origin_date[1]
  
  ## remove model_id from saved data, as it is implied from filepath
  group_df <- select(group_df, -model_id)
  
  ## path to the file from the working directory of the instructional repo
  model_folder <- file.path(
    hub_path,
    "model-output", 
    this_model_id)
  
  ## just the filename, no path
  filename <- paste0(this_origin_date, "-", this_model_id, ".csv")
  
  ## path to the file
  results_path <- file.path(
    model_folder, 
    filename)
  
  ## if this model's model-out directory doesn't exist yet, make it
  if (!file.exists(model_folder)) {
    dir.create(model_folder, recursive = TRUE)
  }
  
  write.csv(group_df, file = results_path, row.names = FALSE)
  
  ### if you run into errors, the code below could help trouble-shoot validations
  # hubValidations::validate_submission(
  #   hub_path = hub_path,
  #   file_path = file.path(this_model_id, filename)
  # )
}


  #local hub validation

  #collect the models from the hub
validation_origin_dates <- origin_dates[which(as.Date(origin_dates) <= as.Date("2017-05-06"))]

new_hub_con <- connect_hub(hub_path)

validation_forecasts <- new_hub_con |> 
  filter(origin_date %in% validation_origin_dates) |> 
  collect_hub()

  #compare to oracle
oracle_output <- connect_target_oracle_output(hub_path) |> 
  collect()

hubEvals::score_model_out(
  validation_forecasts,
  oracle_output
) |> 
  arrange(wis) |> 
  knitr::kable(digits = 2)























  ##### My model building - Fourier AR2 #####


cv_forecasts_season1 <- 
  flu_data_hhs_tscv_season1 |> 
  model(
    fourier_ar2 = ARIMA(my_fourth_root(wili) ~ pdq(2,0,0) + fourier(period = "year", K=3))
  ) |> 
  generate(h = 4, times = 100, bootstrap = TRUE) |> 
  ## the following 3 lines of code ensure that there is a horizon variable in the forecast data
  group_by(.split, .rep, location, .model) |> 
  mutate(horizon = row_number()) |> 
  ungroup() |> 
  as_tibble() |> 
  ## make hubverse-friendly names
  rename(
    target_end_date = origin_date,
    value = .sim,
    model_id = .model
  ) |> 
  left_join(loc_df) |> 
  mutate(origin_date = target_end_date - horizon * 7L) |> 
  ## compute the quantiles
  group_by(model_id, location, origin_date, horizon, target_end_date) |> 
  reframe(tibble::enframe(quantile(value, quantile_levels), "quantile", "value")) |> 
  mutate(output_type_id = as.numeric(stringr::str_remove(quantile, "%"))/100, .keep = "unused",
         target = "ili perc",
         output_type = "quantile",
         model_id = "AEM-fourier_ar2")


#plot one set of forecasts
cv_forecasts_season1 |> 
  filter(origin_date == "2015-12-19") |> 
  plot_step_ahead_model_output(
    flu_data_hhs |> 
      filter(origin_date >= "2015-10-01", 
             origin_date <= "2015-12-19") |> 
      rename(observation = wili),
    x_target_col_name = "origin_date",
    x_col_name = "target_end_date",
    use_median_as_point = TRUE,
    facet = "location",
    facet_scales = "free_y",
    facet_nrow=3
  )

#generate season 2 predictions
cv_forecasts_season2 <- 
  flu_data_hhs_tscv_season2 |> 
  model(
    fourier_ar2 = ARIMA(my_fourth_root(wili) ~ pdq(2,0,0) + fourier(period = "year", K=3))
  ) |> 
  generate(h = 4, times = 100, bootstrap = TRUE) |> 
  ## the following 3 lines of code ensure that there is a horizon variable in the forecast data
  group_by(.split, .rep, location, .model) |> 
  mutate(horizon = row_number()) |> 
  ungroup() |> 
  as_tibble() |> 
  ## make hubverse-friendly names
  rename(
    target_end_date = origin_date,
    value = .sim,
    model_id = .model
  ) |> 
  left_join(loc_df) |> 
  mutate(origin_date = target_end_date - horizon * 7L) |> 
  ## compute the quantiles
  group_by(model_id, location, origin_date, horizon, target_end_date) |> 
  reframe(tibble::enframe(quantile(value, quantile_levels), "quantile", "value")) |> 
  mutate(output_type_id = as.numeric(stringr::str_remove(quantile, "%"))/100, .keep = "unused",
         target = "ili perc",
         output_type = "quantile",
         model_id = "AEM-fourier_ar2")


#Save the model locally
this_model_id <- "AEM-fourier_ar2"

metadata_filepath <- file.path(
  hub_path,
  "model-metadata", 
  paste0(this_model_id, ".yml"))

my_text <- c("team_abbr: \"AEM\"", 
             "model_abbr: \"fourier_ar2\"", 
             "designated_model: true")

writeLines(my_text, metadata_filepath)

# Group the forecasts by task id variables
groups <- bind_rows(cv_forecasts_season1, cv_forecasts_season2) |> 
  group_by(model_id, target, location, origin_date, horizon, target_end_date) |> 
  group_split()

# Save each group as a separate CSV
for (i in seq_along(groups)) {
  group_df <- groups[[i]]
  this_model_id <- group_df$model_id[1]
  this_origin_date <- group_df$origin_date[1]
  
  ## remove model_id from saved data, as it is implied from filepath
  group_df <- select(group_df, -model_id)
  
  ## path to the file from the working directory of the instructional repo
  model_folder <- file.path(
    hub_path,
    "model-output", 
    this_model_id)
  
  ## just the filename, no path
  filename <- paste0(this_origin_date, "-", this_model_id, ".csv")
  
  ## path to the file
  results_path <- file.path(
    model_folder, 
    filename)
  
  ## if this model's model-out directory doesn't exist yet, make it
  if (!file.exists(model_folder)) {
    dir.create(model_folder, recursive = TRUE)
  }
  
  write.csv(group_df, file = results_path, row.names = FALSE)
  
  ### if you run into errors, the code below could help trouble-shoot validations
  # hubValidations::validate_submission(
  #   hub_path = hub_path,
  #   file_path = file.path(this_model_id, filename)
  # )
}


#local hub validation

#collect the models from the hub
validation_origin_dates <- origin_dates[which(as.Date(origin_dates) <= as.Date("2017-05-06"))]

new_hub_con <- connect_hub(hub_path)

validation_forecasts <- new_hub_con |> 
  filter(origin_date %in% validation_origin_dates) |> 
  collect_hub()

#compare to oracle
oracle_output <- connect_target_oracle_output(hub_path) |> 
  collect()

hubEvals::score_model_out(
  validation_forecasts,
  oracle_output
) |> 
  arrange(wis) |> 
  knitr::kable(digits = 2)
































##### My model building - Model 2 #####

#fourier ar2 with 1 difference

# log_transform <- function(x) log(x)
# inv_log_transform <- function(x) exp(x)
# my_log_transform <- new_transformation(log_transform, inv_log_transform)
# 
# flu_data_hhs_tscv_season1 %>%
#   ggplot( aes(x=origin_date, y=wili)) +
#   geom_line()
# 
# flu_data_hhs_tscv_season1 %>%
#   ggplot( aes(x=origin_date, y=log_transform(wili))) +
#   geom_line()


cv_forecasts_season1 <- 
  flu_data_hhs_tscv_season1 |> 
  model(
    fourier_ar210 = ARIMA(my_fourth_root(wili) ~ pdq(2,1,0) + fourier(period = "year", K=3))
  ) |> 
  generate(h = 4, times = 100, bootstrap = TRUE) |> 
  ## the following 3 lines of code ensure that there is a horizon variable in the forecast data
  group_by(.split, .rep, location, .model) |> 
  mutate(horizon = row_number()) |> 
  ungroup() |> 
  as_tibble() |> 
  ## make hubverse-friendly names
  rename(
    target_end_date = origin_date,
    value = .sim,
    model_id = .model
  ) |> 
  left_join(loc_df) |> 
  mutate(origin_date = target_end_date - horizon * 7L) |> 
  ## compute the quantiles
  group_by(model_id, location, origin_date, horizon, target_end_date) |> 
  reframe(tibble::enframe(quantile(value, quantile_levels), "quantile", "value")) |> 
  mutate(output_type_id = as.numeric(stringr::str_remove(quantile, "%"))/100, .keep = "unused",
         target = "ili perc",
         output_type = "quantile",
         model_id = "AEM_fourier_ar210")


#plot one set of forecasts
cv_forecasts_season1 |> 
  filter(origin_date == "2015-12-19") |> 
  plot_step_ahead_model_output(
    flu_data_hhs |> 
      filter(origin_date >= "2015-10-01", 
             origin_date <= "2015-12-19") |> 
      rename(observation = wili),
    x_target_col_name = "origin_date",
    x_col_name = "target_end_date",
    use_median_as_point = TRUE,
    facet = "location",
    facet_scales = "free_y",
    facet_nrow=3
  )

#generate season 2 predictions
cv_forecasts_season2 <- 
  flu_data_hhs_tscv_season2 |> 
  model(
    fourier_ar210 = ARIMA(my_fourth_root(wili) ~ pdq(2,1,0) + fourier(period = "year", K=3))
  ) |> 
  generate(h = 4, times = 100, bootstrap = TRUE) |> 
  ## the following 3 lines of code ensure that there is a horizon variable in the forecast data
  group_by(.split, .rep, location, .model) |> 
  mutate(horizon = row_number()) |> 
  ungroup() |> 
  as_tibble() |> 
  ## make hubverse-friendly names
  rename(
    target_end_date = origin_date,
    value = .sim,
    model_id = .model
  ) |> 
  left_join(loc_df) |> 
  mutate(origin_date = target_end_date - horizon * 7L) |> 
  ## compute the quantiles
  group_by(model_id, location, origin_date, horizon, target_end_date) |> 
  reframe(tibble::enframe(quantile(value, quantile_levels), "quantile", "value")) |> 
  mutate(output_type_id = as.numeric(stringr::str_remove(quantile, "%"))/100, .keep = "unused",
         target = "ili perc",
         output_type = "quantile",
         model_id = "AEM_fourier_ar210")


#Save the model locally
this_model_id <- "AEM_fourier_ar210"

metadata_filepath <- file.path(
  hub_path,
  "model-metadata", 
  paste0(this_model_id, ".yml"))

my_text <- c("team_abbr: \"AEM\"", 
             "model_abbr: \"fourier_ar210\"", 
             "designated_model: true")

writeLines(my_text, metadata_filepath)

# Group the forecasts by task id variables
groups <- bind_rows(cv_forecasts_season1, cv_forecasts_season2) |> 
  group_by(model_id, target, location, origin_date, horizon, target_end_date) |> 
  group_split()

# Save each group as a separate CSV
for (i in seq_along(groups)) {
  group_df <- groups[[i]]
  this_model_id <- group_df$model_id[1]
  this_origin_date <- group_df$origin_date[1]
  
  ## remove model_id from saved data, as it is implied from filepath
  group_df <- select(group_df, -model_id)
  
  ## path to the file from the working directory of the instructional repo
  model_folder <- file.path(
    hub_path,
    "model-output", 
    this_model_id)
  
  ## just the filename, no path
  filename <- paste0(this_origin_date, "-", this_model_id, ".csv")
  
  ## path to the file
  results_path <- file.path(
    model_folder, 
    filename)
  
  ## if this model's model-out directory doesn't exist yet, make it
  if (!file.exists(model_folder)) {
    dir.create(model_folder, recursive = TRUE)
  }
  
  write.csv(group_df, file = results_path, row.names = FALSE)
  
  ### if you run into errors, the code below could help trouble-shoot validations
  # hubValidations::validate_submission(
  #   hub_path = hub_path,
  #   file_path = file.path(this_model_id, filename)
  # )
}


#local hub validation

#collect the models from the hub
validation_origin_dates <- origin_dates[which(as.Date(origin_dates) <= as.Date("2017-05-06"))]

new_hub_con <- connect_hub(hub_path)

validation_forecasts <- new_hub_con |> 
  filter(origin_date %in% validation_origin_dates) |> 
  collect_hub()

#compare to oracle
oracle_output <- connect_target_oracle_output(hub_path) |> 
  collect()

hubEvals::score_model_out(
  validation_forecasts,
  oracle_output
) |> 
  arrange(wis) |> 
  knitr::kable(digits = 2)


