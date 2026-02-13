#### This script performs the analysis on the simulation output ####
library(tidyverse)
#### analysis helper function ####
get_performance <- function(results,
                            aspect_var = NULL,
                            param_regex = "(phi|zeta)",
                            grouping_regex = "k\\d{1}",
                            stat_regex = "(est|pop)"
) {

  custom_regex <- paste0("^", param_regex, "\\d{1,2}_", grouping_regex, "_", stat_regex)

  ## compute bias
  bias <- results |>
    # bring results in wide format: 1 row per parameter (e.g., phi11), group, and stat (est(imate) or pop(ulation) value)
    pivot_longer(
      cols = matches(custom_regex),
      names_to = c("parameter", "cluster", "stat"),
      names_sep = "_"
    ) |>
    # and then back into wide format with est and pop in different columns
    pivot_wider(names_from = stat, values_from = value) |>
    mutate(
      # compute bias
      difference = est - pop,
      # classify as either AR or CR:
      ij     = str_extract(parameter, "(?<=phi)\\d{2}"),
      i      = as.integer(str_sub(ij, 1, 1)),
      j      = as.integer(str_sub(ij, 2, 2)),
      phi_type = if_else(i == j, "AR", "CR")) |>
    # only keep phi (not interested in zeta)
    filter(phi_type %in% c("AR", "CR")) |>
    # compute bias per parameter (grouped by aspect_var)
    group_by(across(all_of(aspect_var)), parameter, cluster, phi_type) |>
    summarize(AB = mean(difference, na.rm = TRUE),
              .groups = "drop") |>
    group_by(across(all_of(aspect_var)), phi_type) |>
    summarize(mean_bias = mean(AB, na.rm = TRUE),
              .groups = "drop") |>
    pivot_wider(names_from = phi_type,
                names_prefix = "meanbias_",
                values_from = mean_bias)

  ## compute ARI, local max percentage, computation time
  outcomes <- results |>
    mutate(convergence_rate = (15-nonconvergences)/15) |>
    group_by(across(all_of(aspect_var))) |>
    summarize(mean_ARI = mean(ARI, na.rm = TRUE),
              mean_convrate = mean(convergence_rate, na.rm = TRUE),
              pct_localmax = mean(local_max, na.rm = TRUE),
              mean_comptime = mean(duration/3600, na.rm = TRUE))

  ## combine
  if(is.null(aspect_var)) {
    outcomes <- bind_cols(bias, outcomes)
  } else {
    outcomes <- bias |>
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
mean(results$duration/3600)
results |>
  group_by(across(all_of(cond_cols))) |>
  summarize(duration = mean(duration/3600, na.rm = TRUE))


#### outcomes (bias, ARI, convergence rate, local maxima, duration) ####
performance <- map(c(list(NULL), as.list(cond_cols)),
                   ~ get_performance(results, aspect_var = .x)) |>
  list_rbind()
