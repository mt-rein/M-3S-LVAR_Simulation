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
adjust_labels <- function(estimated_clusters, true_clusters) {
  # turn estimated_clusters into factor with levels equal to true_clusters
  estimated_clusters <- factor(
    estimated_clusters,
    levels = levels(true_clusters)
  )
  # create contingency table
  crosstable <- table(estimated_clusters, true_clusters)

  # solve the linear sum assignment problem to find the best matching between
  # estimated and true clusters
  assignment <- clue::solve_LSAP(crosstable, maximum = TRUE)

  # create a mapping of old estimated cluster labels to new labels based on the
  # assignment
  mapping <- setNames(
    colnames(crosstable)[as.integer(assignment)],
    rownames(crosstable)
  )

  return(mapping)
}

#### safely/quietly functions ####
run_step1 <- quietly(safely(ezLVAR::step1))
run_step2 <- quietly(safely(ezLVAR::step2))
run_step3 <- quietly(safely(ezLVAR::step3))
