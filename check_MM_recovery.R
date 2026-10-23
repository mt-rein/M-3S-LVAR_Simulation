mm_est <- grep(
  "^(lambda|theta|tau)_f[1-4]_v[1-4]_est$",
  names(results_clean),
  value = TRUE
)

mm_bias <- do.call(
  rbind,
  lapply(mm_est, function(est_col) {
    pop_col <- sub("_est$", "_pop", est_col)
    factor_number <- as.integer(
      sub(".*_f([1-4])_v.*", "\\1", est_col)
    )

    # f3 and f4 are not expected when n_factors == 2
    eligible <- !is.na(results_clean$n_factors) &
      results_clean$n_factors >= factor_number

    error <- results_clean[[est_col]][eligible] -
      results_clean[[pop_col]][eligible]
    error <- error[is.finite(error)]

    n <- length(error)
    bias <- if (n > 0) mean(error) else NA_real_
    mcse <- if (n > 1) sd(error) / sqrt(n) else NA_real_

    data.frame(
      parameter = sub("_est$", "", est_col),
      n_eligible = sum(eligible),
      n_used = n,
      bias = bias,
      CI_lower = if (n > 1) bias - qt(0.975, n - 1) * mcse else NA_real_,
      CI_upper = if (n > 1) bias + qt(0.975, n - 1) * mcse else NA_real_
    )
  })
)

print(mm_bias, row.names = FALSE)
mm_bias |>
  mutate(across(c(bias, CI_lower, CI_upper), ~ sprintf("%.2f", .x))) |>
  print(n = Inf, width = Inf)

mm_bias |>
  mutate(across(c(bias, CI_lower, CI_upper), ~ sprintf("%.2f", .x))) |>
  print(row.names = FALSE)
