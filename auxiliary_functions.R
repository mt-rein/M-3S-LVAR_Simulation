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
    innovation <-  MASS::mvrnorm(1, mu = rep(0, factors),
                                 Sigma = zeta,
                                 empirical = FALSE)
    # simulate the first deviation (delta) only from the innovation
    if (i == 1) {
      delta <- innovation
    }

    # loop through all the rows:
    # predict the current temporal deviation (delta) from the previous deviation
    # then add random innovation
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

#### EStep() ####
# perform the E-Step of the EM algorithm
EStep <- function(pi_ks, ngroup, nclus, loglik) {

  max_g <- rep(0, ngroup)
  z_gks <- matrix(NA, nrow = ngroup, ncol = nclus)

  for (g in 1:ngroup) {
    for (k in 1:nclus) {
      z_gks[g, k] <- log(pi_ks[k]) + loglik[g, k]
    }
    max_g[g] <- max(z_gks[g, ]) # prevent arithmetic underflow
    z_gks[g, ] <- exp(z_gks[g, ] - rep(max_g[g], nclus))
  }

  # divide by the rowwise sum of the above calculated part
  z_gks <- diag(1 / apply(z_gks, 1, sum)) %*% z_gks

  return(z_gks)
}
# taken from https://github.com/AndresFPA/mmgsem/blob/main/R/E_Step.R

#### MStep() ####
MStep <- function(n_clusters,
                  weights,
                  objectives,
                  model_list,
                  startvalues = startvalues,
                  verbose = FALSE) {
  # create one model per cluster
  # each model is a multi-group model
  # where the person-models (each person is a "group") are weighted
  # by their posterior probabilities
  clustermodels <- vector(mode = "list", length = n_clusters)
  for (k in 1:n_clusters) {
    clustername <- paste0("model_k", k)
    weighted_objectives <- paste(weights[, k], "*", objectives,
                                 collapse = " + ")
    model <- OpenMx::mxModel(clustername, model_list,
                             OpenMx::mxAlgebraFromString(weighted_objectives,
                                                         name = "weightedfit"),
                             OpenMx::mxFitFunctionAlgebra("weightedfit"))
    # add start values:
    model <- OpenMx::omxSetParameters(model,
                                      labels = names(startvalues[[k]]),
                                      values = startvalues[[k]])

    clustermodels[[k]] <- model
  }
  names(clustermodels) <- paste0("model_k", 1:n_clusters)

  # run the models
  clustermodels_run <- clustermodels |>
    purrr::map(OpenMx::mxRun,
               silent = !verbose,
               suppressWarnings = TRUE)

  return(clustermodels_run)
}

#### generate_startvalues() ####
# generate random starting values for OpenMx
generate_startvalues <- function(n_clusters, n_factors,
                                 labels_phi, labels_zeta) {
  startvalues <- vector(mode = "list", length = n_clusters)
  for (k in 1:n_clusters) {
    # generate a stationary matrix of regression coefficients
    phistart <- matrix(runif(n_factors * n_factors, 0.05, .6), nrow = n_factors)
    ev <- eigen(phistart)$values
    phistart_scaled <- phistart * (.9 / max(Mod(ev)))

    # generate a positive definitive matrix of innovation (co)variances
    zetastart <- matrix(runif(n_factors * n_factors, 0.3, 1.5),
                        nrow = n_factors)
    zetastartPD <- Matrix::nearPD(zetastart)$mat |> as.matrix()

    startvalues[[k]] <- c(
      c(phistart_scaled),
      zetastartPD[lower.tri(zetastartPD, diag = TRUE)]
    )
    phi_names <- labels_phi[!is.na(labels_phi)]
    zeta_names <- labels_zeta[!is.na(labels_zeta)]
    names(startvalues[[k]]) <- c(phi_names, unique(zeta_names))
  }
  return(startvalues)
}

#### compute_observed_data_LL() ####
# compute observed LL from personwise LL and class proportions
compute_observed_data_LL <- function(personLL, class_proportions) {
  # sum with the log of class proportions:
  personLL_weighted <- sweep(personLL, 2, log(class_proportions), "+")
  # get the max value per row:
  max_i <- apply(personLL_weighted, 1, max)
  # substract max value from each row (prevent arithmetic underflow):
  minus_max <- sweep(personLL_weighted, 1, max_i, "-")
  # exp to get out of log space:
  personLL_exp <- exp(minus_max)
  # sum per row and then take the log again:
  personLL_summed <- log(rowSums(personLL_exp))
  # re-add the max value and then sum to obtain observed data LL
  observed_data_LL <- sum(personLL_summed + max_i)

  return(observed_data_LL)
}

#### adjust_labels() ####
adjust_labels <- function(modal_matrix, true_clusters) {
  # create matrix with all permutations of labels
  n_clusters <- ncol(modal_matrix)
  combinations <- paste0("cluster", 1:n_clusters) |>
    RcppAlgos::permuteGeneral() |>
    as.data.frame()
  combinations$diagsum <- 0

  for (i in seq_len(nrow(combinations))) {
    temp <- modal_matrix
    # relabel the clusters according to this permutation:
    colnames(temp) <- combinations[i, 1:n_clusters]
    # turn the clusterassignment matrix (0s and 1s) into a vector (factor):
    clusterassignment_estimated <- colnames(temp)[
      max.col(temp)
    ] |> as.factor()
    # creates a cross table of estimated and true cluster assignments:
    crosstable <- table(clusterassignment_estimated, true_clusters)
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
run_step1 <- quietly(safely(step1))
run_step2 <- quietly(safely(step2))
run_step3 <- quietly(safely(step3))