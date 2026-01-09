library(tidyverse)

# read data and sort by iteration
results <- read_csv("output_sim.csv") |>
  arrange(iteration)

# names of the condition columns
cond_cols <- c("n_obs", "n_clusters", "n_factors")

#### errors, warnings, non converged starts ####
results |>
  summarize(across(step1_warning:step3_error, ~ sum()))
# no warnings or errors

results$nonconvergences |> table()
# 1 non-converged start in 6 iterations

#### duration and sensitivity to local maxima ####
duration <- results |>
  group_by(across(all_of(cond_cols))) |>
  summarize(duration_hours = mean(duration)/3600, .groups = "drop")

localmax <- results |>
  group_by(across(all_of(cond_cols))) |>
  summarize(pct_locmax = mean(local_max), .groups = "drop")

#### bias and ARI ####
# transform into long format (1 row per parameter)
phi_long <- results |>
  pivot_longer(
    # get both estimates (e.g., phi_11_k1), and population parameters (phi_11_k1_pop)
    cols = matches("^phi\\d{2}_k\\d+(?:_pop)?$"),
    names_to  = "name",
    values_to = "value"
  ) |>
  # remove NAs
  drop_na(value)  |>
  mutate(
    # get source (estimate or population)
    source = if_else(str_ends(name, "_pop"), "pop", "est"),
    # which parameter?
    param  = str_remove(name, "_pop$"),
    # classify as either AR or CR:
    ij     = str_extract(param, "(?<=phi)\\d{2}"),
    i      = as.integer(str_sub(ij, 1, 1)),
    j      = as.integer(str_sub(ij, 2, 2)),
    phi_type = if_else(i == j, "AR", "CR")
  )

# extract population values per parameter (equal across all replications within a condition)
pop <- phi_long |>
  filter(source == "pop") |>
  group_by(across(all_of(cond_cols)), param, phi_type) |>
  summarise(pop_value = first(value), .groups = "drop")

# mean estimates for each parameter across replications within a condition
est <- phi_long |>
  filter(source == "est") |>
  group_by(across(all_of(cond_cols)), param, phi_type) |>
  summarise(mean_estimate = mean(value), .groups = "drop")

# combine estimates and population values
combined <- est |>
  full_join(pop, by = c(cond_cols, "param", "phi_type"))

# bias per condition & phi_type (average across parameters of same type (AR/CR))
bias <- combined |>
  # compute bias by condition:
  mutate(bias = mean_estimate - pop_value) |>
  # average across parameters of the same type:
  group_by(across(all_of(cond_cols)), phi_type) |>
  summarise(bias = mean(bias), .groups = "drop") |>
  # transform into wide format, with one column per type (each row corresponds to one condition)
  pivot_wider(names_from = phi_type, values_from = bias, names_prefix = "bias_")  # -> bias_AR, bias_CR

# ARI per condition (single ARI column)
ari <- results |>
  group_by(across(all_of(cond_cols))) |>
  summarise(ARI = mean(ARI), .groups = "drop")   # uses 'ARI' column

#### combine all results ####
cond_summary <- ari |>
  full_join(bias, by = cond_cols) |>
  full_join(duration, by = cond_cols) |>
  full_join(localmax, by = cond_cols)


# 7) Aspect-level summary: for each aspect, average over the other aspects
overall <- cond_summary |>
  summarise(
    ARI     = mean(ARI),
    bias_AR = mean(bias_AR),
    bias_CR = mean(bias_CR),
    duration_hours = mean(duration_hours),
    pct_locmax = mean(pct_locmax)
  ) |>
  mutate(aspect = "overall", level = "overall") |>
  select(aspect, level, bias_AR, bias_CR, ARI, pct_locmax, duration_hours)

aspects <- map_dfr(cond_cols, function(aspect) {
  cond_summary |>
    group_by(.data[[aspect]]) |>
    summarise(
      ARI     = mean(ARI, na.rm = TRUE),
      bias_AR = mean(bias_AR, na.rm = TRUE),
      bias_CR = mean(bias_CR, na.rm = TRUE),
      duration_hours = mean(duration_hours),
      pct_locmax = mean(pct_locmax),
      .groups = "drop"
    ) |>
    rename(level = !!aspect) |>
    mutate(aspect = aspect,
           level = as.character(level)) |>
    select(aspect, level, bias_AR, bias_CR, ARI, pct_locmax, duration_hours)
})

final_table <- bind_rows(overall, aspects)