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
      bias = mean(difference, na.rm = TRUE),
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
      SE_recovery = mean(se, na.rm = TRUE) / sd(est, na.rm = TRUE),
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
      mean_ARI = mean(ARI, na.rm = TRUE),
      mean_convrate = mean(convergence_rate, na.rm = TRUE),
      mean_comptime = mean(duration / 3600, na.rm = TRUE)
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

results <- results |>
  filter(if_all(matches("^phi\\d{1,2}_k\\d{1}_se$"), ~ is.na(.x) | .x <= 1))

# rename estimate columns (add _est suffix):
est_cols <- names(results)[
  str_detect(names(results), "^(phi|zeta)\\d{1,2}_k\\d{1}$")
]
names(est_cols) <- paste0(est_cols, "_est")
results <- results |>
  rename(all_of(est_cols))


# names of the condition columns
cond_cols <- c("n_obs", "n_clusters", "n_factors")

#### errors, warnings, non converged starts ####
results |>
  summarize(across(step1_warning:step3_error, ~ sum()))
# no warnings or errors

results$nonconvergences |> table()

# maximum of 6 non-convergences

#### duration in hours (detailed in unique conditions) ####
mean(results$duration / 3600)
results |>
  group_by(across(all_of(cond_cols))) |>
  summarize(
    duration_avg = mean(duration / 3600, na.rm = TRUE),
    min = min(duration),
    max = max(duration)
  )


#### outcomes (bias, ARI, convergence rate, local maxima, duration) ####
performance <- map(
  c(list(NULL), as.list(cond_cols)),
  ~ get_performance(results, aspect_var = .x)
) |>
  list_rbind()
print(performance, width = Inf)
