#### This script performs the analysis on the simulation output ####
library(tidyverse)
#### analysis helper functions ####
get_performance <- function(
  results,
  aspect_var = NULL
) {
  # custom regex for this simulation study (indicates structure of the variable
  # names in the results data frame)
  custom_regex <- paste0(
    "^",
    "(phi|zeta)",
    "\\d{1,2}_",
    "k\\d{1}",
    "_",
    "(est|pop|se)"
  )

  ## reorganize results data frame
  long_results <- results |>
    # bring results in long format: 1 row per parameter (e.g., phi11), group,
    # and stat (est(imate), pop(ulation) value, or se(standard error))
    pivot_longer(
      cols = matches(custom_regex),
      names_to = c("parameter", "cluster", "stat"),
      names_sep = "_"
    ) |>
    # and then back into wide format with est, pop, and se in different columns
    pivot_wider(names_from = stat, values_from = value) |>
    mutate(param_id = paste(parameter, cluster, sep = "_"))

  # compute bias per parameter (grouped by aspect_var)
  bias <- long_results |>
    filter(!is.na(est), !is.na(pop)) |>
    mutate(difference = est - pop) |>
    group_by(across(all_of(aspect_var)), param_id) |>
    summarize(
      bias = mean(difference),
      .groups = "drop"
    ) |>
    pivot_wider(
      names_from = param_id,
      values_from = bias,
      names_prefix = "bias_"
    )

  # compute SE recovery per parameter (grouped by aspect_var)
  se_recovery <- long_results |>
    filter(!is.na(est), !is.na(se)) |>
    group_by(across(all_of(aspect_var)), param_id) |>
    summarize(
      SE_recovery = mean(se) / sd(est),
      .groups = "drop"
    ) |>
    pivot_wider(
      names_from = param_id,
      values_from = SE_recovery,
      names_prefix = "SErec_"
    )

  ## compute ARI, convergence rate, and computation time
  outcomes <- results |>
    mutate(
      convergence_rate = (15 - nonconvergences) / 15
    ) |>
    group_by(across(all_of(aspect_var))) |>
    summarize(
      mean_ARI = mean(ARI),
      mean_convrate = mean(convergence_rate),
      mean_comptime = mean(duration / 3600)
    )

  ## combine
  if (is.null(aspect_var)) {
    outcomes <- bind_cols(bias, se_recovery, outcomes)
  } else {
    outcomes <- bias |>
      full_join(se_recovery, by = aspect_var) |>
      full_join(outcomes, by = aspect_var)
  }

  # add "aspect" and "level" columns
  if (is.null(aspect_var)) {
    outcomes <- outcomes |>
      mutate(aspect = "overall", level = NA_character_, .before = everything())
  } else {
    outcomes <- outcomes |>
      rename(level = !!aspect_var) |>
      mutate(aspect = aspect_var, .before = level) |>
      mutate(level = as.character(level))
  }

  return(outcomes)
}

# read data and sort by iteration
results <- read_csv("output_sim.csv") |>
  arrange(iteration)

# names of the condition columns
cond_cols <- c("n_obs", "n_clusters", "n_factors")

#### errors, warnings, non converged starts ####
# display number of warnings and errors in each step
results |>
  summarize(across(step1_warning:step3_error, ~ sum(.x)))
# 60 warnings/errors in step3
unique(results$step3_warning_text)
unique(results$step3_error_text)
# --> warnings and errors refer to non-zero status codes

# display number of non-converged starts
results$nonconvergences |> table()

# check codes:
table(results$status_k1)
table(results$status_k2)
table(results$status_k3)
table(results$status_k4)

# remove all rows where at least one status code is not 0:
results_clean <- results %>%
  filter(if_all(
    c(status_k1, status_k2, status_k3, status_k4),
    ~ .x == 0 | is.na(.x)
  ))

# remove all rows where one phi or zeta estimate or se is NA when it should exist:
results_clean <- results_clean |>
  filter(
    # factors 1–2, clusters 1–2: always required
    if_all(
      matches("^(phi|zeta)[12][12]_k[12]_(est|se)$"),
      ~ !is.na(.x)
    ),

    # factors 1–2, clusters 3–4: required only with 4 clusters
    n_clusters == 2 |
      if_all(
        matches("^(phi|zeta)[12][12]_k[34]_(est|se)$"),
        ~ !is.na(.x)
      ),

    # parameters involving factor 3 or 4, clusters 1–2
    n_factors == 2 |
      if_all(
        matches("^(phi|zeta)([34][1-4]|[1-4][34])_k[12]_(est|se)$"),
        ~ !is.na(.x)
      ),

    # parameters involving factor 3 or 4, clusters 3–4
    n_factors == 2 |
      n_clusters == 2 |
      if_all(
        matches("^(phi|zeta)([34][1-4]|[1-4][34])_k[34]_(est|se)$"),
        ~ !is.na(.x)
      )
  )


# remove rows with SEs greater than 1
results_clean <- results_clean |>
  filter(if_all(matches("^phi\\d{1,2}_k\\d{1}_se$"), ~ is.na(.x) | .x <= 1))

#### duration in hours (detailed in unique conditions) ####
mean(results_clean$duration / 3600)
results_clean |>
  group_by(across(all_of(cond_cols))) |>
  summarize(
    duration_avg = mean(duration / 3600, na.rm = TRUE),
    min = min(duration),
    max = max(duration)
  )


#### outcomes (bias, ARI, convergence rate, local maxima, duration) ####
performance <- map(
  c(list(NULL), as.list(cond_cols)),
  ~ get_performance(results_clean, aspect_var = .x)
) |>
  list_rbind()
performance |>
  mutate(across(where(is.numeric), ~ round(.x, 2))) |>
  print(width = Inf)
