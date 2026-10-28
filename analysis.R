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
  # Compute bias, its Monte Carlo standard error, and a 95% CI
  bias <- long_results |>
    filter(!is.na(est), !is.na(pop)) |>
    mutate(difference = est - pop) |>
    group_by(across(all_of(aspect_var)), param_id) |>
    summarize(
      n_iterations = n(),
      bias = mean(difference),
      sd_difference = sd(difference),
      .groups = "drop"
    ) |>
    mutate(
      MCSE_bias = sd_difference / sqrt(n_iterations),
      ci_lower_bias = bias -
        qt(0.975, df = pmax(n_iterations - 1L, 1L)) * MCSE_bias,
      ci_upper_bias = bias +
        qt(0.975, df = pmax(n_iterations - 1L, 1L)) * MCSE_bias
    ) |>
    select(-sd_difference, -n_iterations) |>
    pivot_wider(
      names_from = param_id,
      values_from = c(bias, MCSE_bias, ci_lower_bias, ci_upper_bias),
      names_glue = "{.value}_{param_id}"
    )

  # Compute Relative % error in ModSE per parameter (grouped by aspect_var)
  # following Morris et al. (2019), Table 6
  se_recovery <- long_results |>
    filter(!is.na(est), !is.na(se)) |>
    group_by(across(all_of(aspect_var)), param_id) |>
    summarize(
      n_iterations = n(),
      # EmpSE: empirical SD of the estimates
      EmpSE = sd(est),
      # ModSE: root-mean of the squared model SEs  [i.e. sqrt(mean(se^2))]
      ModSE = sqrt(mean(se^2)),
      # Var[Var_hat(theta_hat)]: variance of the per-iteration squared SEs
      # (se_i^2 is the per-iteration variance estimate)
      var_of_var = var(se^2),
      .groups = "drop"
    ) |>
    mutate(
      # Relative % error in ModSE  (Morris 2019, Table 6)
      relpcte = 100 * (ModSE / EmpSE - 1),
      # Monte Carlo SE for Relative % error in ModSE (Morris 2019, Table 6)
      MCSE_relpcte = 100 *
        (ModSE / EmpSE) *
        sqrt(
          var_of_var /
            (4 * n_iterations * ModSE^4) +
            1 / (2 * (n_iterations - 1))
        ),
      # 95% CI
      ci_lower_relpcte = relpcte -
        qt(0.975, df = pmax(n_iterations - 1L, 1L)) * MCSE_relpcte,
      ci_upper_relpcte = relpcte +
        qt(0.975, df = pmax(n_iterations - 1L, 1L)) * MCSE_relpcte
    ) |>
    select(-n_iterations, -EmpSE, -ModSE, -var_of_var) |>
    pivot_wider(
      names_from = param_id,
      values_from = c(
        relpcte,
        MCSE_relpcte,
        ci_lower_relpcte,
        ci_upper_relpcte
      ),
      names_glue = "{.value}_{param_id}"
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
      mean_comptime = mean(duration / 3600),
      .groups = "drop"
    )

  #
  # combine
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
      mutate(
        # paste the variable names together, e.g. "n_factors:n_obs"
        aspect = paste(aspect_var, collapse = ":"),
        # paste the corresponding values together, e.g. "100:20"
        level = paste(!!!syms(aspect_var), sep = ":"),
        .before = everything()
      ) |>
      # drop the original individual aspect columns (now redundant)
      select(-all_of(aspect_var)) |>
      mutate(level = as.character(level))
  }

  # reorder columns: group all metrics per parameter together
  # get all param_ids from the column names (using bias_ prefix as reference)
  param_ids <- unique(long_results$param_id)

  param_col_order <- purrr::map(param_ids, function(pid) {
    c(
      paste0("bias_", pid),
      paste0("MCSE_bias_", pid),
      paste0("ci_lower_bias_", pid),
      paste0("ci_upper_bias_", pid),
      paste0("relpcte_", pid),
      paste0("MCSE_relpcte_", pid),
      paste0("ci_lower_relpcte_", pid),
      paste0("ci_upper_relpcte_", pid)
    )
  }) |>
    purrr::list_c()

  # reorder: bookkeeping columns first, then parameter columns, then other outcomes
  outcomes <- outcomes |>
    select(
      any_of(c("aspect", "level")),
      # bookkeeping
      any_of(param_col_order),
      # parameter metrics, grouped by param_id
      everything()
      # remaining cols (mean_ARI, mean_convrate, etc.)
    )
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


#### split into two data frames: one with two clusters and one with four clusters ####
results_2k <- results_clean |>
  filter(n_clusters == 2)

results_4k <- results_clean |>
  filter(n_clusters == 4)

#### outcomes (bias, ARI, convergence rate, local maxima, duration) ####
performance_2k <- map(
  c(list(NULL), as.list(cond_cols[-3])),
  ~ get_performance(results_2k, aspect_var = .x)
) |>
  list_rbind()
performance_2k |>
  mutate(across(where(is.numeric), ~ round(.x, 2))) |>
  print(width = Inf)
