#### This script defines auxiliary functions for the simulation ####

#### sim_VAR() ####
# this function generates data for a single individual
# according to a vector autoregressive model
sim_VAR <- function(factors, obs, phi, zeta, mu, burn_in = 0) {
  # factors = number of factors
  # obs = number of observations
  # phi = auto-regressive effect (a matrix in case of multiple constructs)
  # zeta = innovation variance (a matrix in case of multiple constructs)
  # mu = latent means (a vector in case of multiple constructs)
  # burn_in = length of burn in (remove influence of initial random draw)

  # create empty dataframe of length obs + burn_in
  data <- as.data.frame(matrix(NA, nrow = burn_in + obs, ncol = factors))
  names(data) <- paste0("eta", 1:factors)

  for (i in seq_len(nrow(data))) {
    innovation <- MASS::mvrnorm(
      1,
      mu = rep(0, factors),
      Sigma = zeta,
      empirical = FALSE
    )
    # simulate the first deviation (delta) only from the innovation
    if (i == 1) {
      delta <- innovation
    }

    # loop through all the rows: predict the current temporal deviation (delta)
    # from the previous deviation then add random innovation
    if (i > 1) {
      delta <- phi %*% delta + innovation
    }
    # sum stable base-line (mu) and temporal deviation (delta):
    data[i, ] <- mu + delta
  }

  # remove the first rows, depending on length of burn in
  if (burn_in > 0) {
    data <- dplyr::slice(data, burn_in + seq_len(dplyr::n()))
  }

  data$obs <- seq_len(nrow(data))

  return(data)
}

#### adjust_labels() ####
adjust_labels <- function(estimated_clusters, true_clusters, n_clusters) {
  estimated_clusters <- factor(
    estimated_clusters,
    levels = levels(true_clusters)
  )
  # create matrix with all permutations of labels
  combinations <- levels(true_clusters) |>
    RcppAlgos::permuteGeneral() |>
    as.data.frame()
  combinations$diagsum <- 0

  for (i in seq_len(nrow(combinations))) {
    temp <- estimated_clusters |>
      factor(labels = combinations[i, 1:n_clusters]) |>
      # turn the factor into a character vector and then back into a factor with
      # the same levels as true_clusters:
      as.character() |>
      factor(levels = levels(true_clusters))

    # creates a cross table of estimated and true cluster assignments:
    crosstable <- table(temp, true_clusters)
    # compute the sum of the diagonal of the cross table and save it
    combinations$diagsum[i] <- sum(diag(crosstable))
  }
  # choose the permutation with the largest sum of the diagonal:
  max_diagsum <- which.max(combinations$diagsum)
  new_labels <- combinations[max_diagsum, 1:n_clusters] |>
    as.character()

  return(new_labels)
}

#### safely/quietly functions ####
run_step1 <- quietly(safely(ezLVAR::step1))
run_step2 <- quietly(safely(ezLVAR::step2))
run_step3 <- quietly(safely(ezLVAR::step3))
