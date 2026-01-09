#### This script defines the simulation function do_sim() ####

do_sim <- function(pos, cond, outputfile, verbose = FALSE) {
  # pos = position in the condition grid
  # cond = the condition grid
  # outputfile = file name for the output CSV file
  # verbose = if TRUE, prints a message after the iteration is finished


  output_list <- list("iteration" = cond$iteration[pos],
                      "replication" = cond$replication[pos])

  # get condition levels and set seed:
  n_obs <- output_list[["n_obs"]] <- cond$n_obs[pos]
  n_clusters <- output_list[["n_clusters"]] <- cond$n_clusters[pos]
  n_factors <- output_list[["n_factors"]] <- cond$n_factors[pos]
  seed <- output_list[["seed"]] <- cond$seed[pos]
  set.seed(seed)

  # change some global options in OpenMx:
  OpenMx::mxOption(key = "Calculate Hessian", value = "No")
  OpenMx::mxOption(key = "Standard Errors", value = "No")

  #### set data generation parameters ####
  # number of individuals:
  n_persons <- 60
  ## regression parameters:
  # if two factors:
  if (n_factors == 2) {
    # Cluster 1: stable unconnected
    phimat_k1_pop <- matrix(c(
      0.7,  0.0,
      0.0,  0.7
    ), ncol = 2, byrow = TRUE)

    # Cluster 2: interconnected
    phimat_k2_pop <- matrix(c(
      0.3,  0.3,
      0.3,  0.3
    ), ncol = 2, byrow = TRUE)

    if (n_clusters == 4) {
      # Cluster 3: mixed
      phimat_k3_pop <- matrix(c(
        0.6, -0.3,
        0.2,  0.1
      ), ncol = 2, byrow = TRUE)

      # Cluster 4: antagonistic
      phimat_k4_pop <- matrix(c(
        0.3, -0.2,
        -0.2,  0.3
      ), ncol = 2, byrow = TRUE)
    } else {
      # create matrices with NA if there's only two clusters:
      phimat_k3_pop <- phimat_k4_pop <- matrix(NA, nrow = 2, ncol = 2)
    }
  }

  # if four factors:
  if (n_factors == 4) {
    # Cluster 1: stable unconnected
    phimat_k1_pop <- matrix(c(
      0.68,  0.00,  0.00,  0.00,
      0.00,  0.68,  0.00,  0.00,
      0.00,  0.00,  0.68,  0.00,
      0.00,  0.00,  0.00,  0.68
    ), ncol = 4, byrow = TRUE)

    # Cluster 2: interconnected
    phimat_k2_pop <- matrix(c(
      0.28,  0.25,  0.20,  0.18,
      0.25,  0.28,  0.18,  0.20,
      0.20,  0.18,  0.28,  0.25,
      0.18,  0.20,  0.25,  0.28
    ), ncol = 4, byrow = TRUE)

    if (n_clusters == 4) {
      # Cluster 3: mixed
      phimat_k3_pop <- matrix(c(
        0.60, -0.25,  0.18,  0.05,
        0.18,  0.12, -0.15,  0.22,
        -0.10,  0.22,  0.48, -0.18,
        0.12, -0.05,  0.08,  0.22
      ), ncol = 4, byrow = TRUE)

      # Cluster 4: antagonistic
      phimat_k4_pop <- matrix(c(
        0.30, -0.12, -0.08,  0.00,
        -0.12,  0.30,  0.00, -0.08,
        -0.08,  0.00,  0.30, -0.12,
        0.00, -0.08, -0.12,  0.30
      ), ncol = 4, byrow = TRUE)
    } else {
      # create matrices with NA if there's only two clusters:
      phimat_k3_pop <- phimat_k4_pop <- matrix(NA, nrow = 4, ncol = 4)
    }
  }

  ## innovation variances
  # fixed effect innovation covariance matrix is invariant across clusters
  zetamat_pop <- matrix(.5, nrow = n_factors, ncol = n_factors)
  diag(zetamat_pop) <- 1.5

  ## grand means
  if (n_factors == 2) {
    grandmeans <- c(5, 3)
  } else {
    grandmeans <- c(5, 3, 6, 2)
  }
  mu_variance <- matrix(.3, nrow = n_factors, ncol = n_factors)
  diag(mu_variance) <- 1

  ## measurement model parameters
  # loadings:
  lambda <- replicate(n_factors, rep(.8, 4), simplify = FALSE) |>
    lavaan::lav_matrix_bdiag()

  # residual variances:
  theta <- rep(.36, n_factors * 4) |>
    diag()

  # intercepts:
  tau <- rep(1, n_factors * 4)

  #### generate data ####
  ## create cluster assignment vector:
  clusterassignment_true <- rep(paste0("cluster", 1:n_clusters),
                                length.out = n_persons) |>
    # shuffle:
    sample(n_persons) |>
    factor()

  # create empty data frame for all (observed) items:
  eta_vars <- paste0("eta", 1:n_factors)
  eta <- data.frame(id = integer(),
                    obs = integer(),
                    k_true = factor(levels = levels(clusterassignment_true)))
  # add eta variables (depending on number of factors)
  for (var in eta_vars) {
    eta[[var]] <- numeric()
  }

  for (i in 1:n_persons) {
    # get cluster membership:
    k_i <- clusterassignment_true[i]

    ## create (partly person-specific) data-generating parameter values
    ## from population values
    # get correct phi matrix:
    if (k_i == "cluster1") {
      phimat <- phimat_k1_pop
    }
    if (k_i == "cluster2") {
      phimat <- phimat_k2_pop
    }
    if (k_i == "cluster3") {
      phimat <- phimat_k3_pop
    }
    if (k_i == "cluster4") {
      phimat <- phimat_k4_pop
    }

    # generate person-specific latent means:
    mu_i <- MASS::mvrnorm(1, mu = grandmeans,
                          Sigma = mu_variance)

    eta_i <- sim_VAR(factors = n_factors, obs = n_obs,
                     phi = phimat, zeta = zetamat_pop,
                     mu = mu_i,
                     burn_in = 10)

    eta_i$id <- i
    eta_i$k_true <- k_i

    # add person-data to full data frame:
    eta <- dplyr::bind_rows(eta, eta_i)
  }

  # generate measurement error from residual variances:
  epsilon <- MASS::mvrnorm(nrow(eta), mu = rep(0, n_factors * 4),
                           Sigma = theta, empirical = FALSE)

  # transform factor scores into observed scores:
  data <- t(tau + lambda %*% t(eta[, eta_vars])) + epsilon |>
    as.data.frame()
  colnames(data) <- paste0("v", 1:(n_factors * 4))
  # add id, obs and true cluster variable:
  data$id <- eta$id
  data$obs <- eta$obs
  data$k_true <- eta$k_true

  #### add population parameters to output ####
  # phi:
  for (k in 1:4) {
    # get phimat of cluster k:
    phimat <- get(paste0("phimat_k", k, "_pop"))
    # loop over all entries, check if they exist, and extract value/set to NA:
    for (i in 1:4) {
      for (j in 1:4) {
        out_name <- paste0("phi", i, j, "_k", k, "_pop")
        value <- ifelse(i <= n_factors && j <= n_factors, phimat[i, j], NA)
        output_list[[out_name]] <- value
      }
    }
  }

  # zeta
  for (i in 1:4) {
    for (j in i:4) {
      # (j in i:4) ensures j >= i (diagonal + lower triangle)
      out_name <- paste0("zeta", i, j, "_pop")
      value <- ifelse(i <= n_factors && j <= n_factors, zetamat_pop[i, j], NA)
      output_list[[out_name]] <- value
    }
  }


  #### Step 1 ####
  start <- Sys.time()
  # MM in lavaan syntax:
  if (n_factors == 2) {
    model_step1 <- list(
      "
      f1 =~ v1 + v2 + v3 + v4
      v1 ~ 1*1
      f1 ~ NA*1
      ",
      "
      f2 =~ v5 + v6 + v7 + v8
      v5 ~ 1*1
      f2 ~ NA*1
      "
    )
  }
  if (n_factors == 4) {
    model_step1 <- list("
      f1 =~ v1 + v2 + v3 + v4
      v1 ~ 1*1
      f1 ~ NA*1
      ",
                        "
      f2 =~ v5 + v6 + v7 + v8
      v5 ~ 1*1
      f2 ~ NA*1
      ",
                        "
      f3 =~ v9 + v10 + v11 + v12
      v9 ~ 1*1
      f3 ~ NA*1
      ",
                        "
      f4 =~ v13 + v14 + v15 + v16
      v13 ~ 1*1
      f4 ~ NA*1
      "
    )
  }

  # run step 1:
  output_step1 <- run_step1(data = data,
                            measurementmodel = model_step1,
                            id = "id")
  # extract error/warning messages (if applicable):
  step1_warning <- ifelse(rlang::is_empty(output_step1$warnings),
                          FALSE, TRUE)
  step1_warning_text <- ifelse(step1_warning,
                               paste(c(output_step1$warnings),
                                     collapse = "; "),
                               "")
  step1_error <- ifelse(rlang::is_empty(output_step1$result$error),
                        FALSE, TRUE)
  step1_error_text <- ifelse(step1_error,
                             paste(c(output_step1$result$error),
                                   collapse = "; "),
                             "")

  #### Step 2 ####
  # only proceed if there is no error in step 1:
  if (!step1_error) {
    output_step2 <- run_step2(step1output = output_step1$result$result,
                              id = "id")
    # extract error/warning messages (if applicable):
    step2_warning <- ifelse(rlang::is_empty(output_step2$warnings),
                            FALSE, TRUE)
    step2_warning_text <- ifelse(step2_warning,
                                 paste(c(output_step2$warnings),
                                       collapse = "; "),
                                 "")
    step2_error <- ifelse(rlang::is_empty(output_step2$result$error),
                          FALSE, TRUE)
    step2_error_text <- ifelse(step2_error,
                               paste(c(output_step2$result$error),
                                     collapse = "; "),
                               "")
  } else {
    step2_warning <- FALSE
    step2_warning_text <- "step1 not successful"
    step2_error <- FALSE
    step2_error_text <- "step1 not successful"
  }

  #### Step 3 ####
  # only proceed if there is no error in step 1 as well as step 2
  if (!step1_error && !step2_error) {
    output_step3 <- run_step3(step2output = output_step2$result$result,
                              n_clusters = n_clusters,
                              n_starts = 25, n_best_starts = 15,
                              maxit = 100,
                              true_clusters = clusterassignment_true,
                              verbose = FALSE)
    duration <- difftime(Sys.time(), start, unit = "s")
    # extract error/warning messages (if applicable):
    step3_warning <- ifelse(rlang::is_empty(output_step3$warnings),
                            FALSE, TRUE)
    step3_warning_text <- ifelse(step3_warning,
                                 paste(c(output_step3$warnings),
                                       collapse = "; "),
                                 "")
    step3_error <- ifelse(rlang::is_empty(output_step3$result$error),
                          FALSE, TRUE)
    step3_error_text <- ifelse(step3_error,
                               paste(c(output_step3$result$error),
                                     collapse = "; "),
                               "")
  } else {
    step3_warning <- FALSE
    step3_warning_text <- "step1 or step2 not successful"
    step3_error <- FALSE
    step3_error_text <- "step1 or step2 not successful"
  }

  #### extract results of three-step estimation ####
  if (!step1_error && !step2_error && !step3_error) {
    ## if step 3 was successful:
    results <- output_step3$result$result

    ## adjust potential label switching:
    new_labels <- results$clustering$modal_assignment |>
      adjust_labels(true_clusters = clusterassignment_true)
    # swap the labels accordingly in the output of step 3:
    names(results$estimates) <-
      colnames(results$clustering$posterior_prob) <-
      names(results$clustering$class_proportions) <-
      colnames(results$clustering$modal_assignment) <-
      new_labels

    ## duration and number of non-convergences:
    output_list[["duration"]] <- duration |>
      as.numeric()
    output_list[["nonconvergences"]] <- results$other$nonconvergences

    ## ARI and local maximum:
    # ARI:
    # turn the modal clustermembership matrix (0s and 1s) into a vector (factor):
    modal_mat <- results$clustering$modal_assignment
    clusterassignment_estimated <- colnames(modal_mat)[
      max.col(modal_mat)
    ] |> as.factor()
    output_list[["ARI"]] <- mcclust::arandi(clusterassignment_true,
                                                      clusterassignment_estimated,
                                                      adjust = TRUE)

    # local maximum:
    proxy_max <- results$other$proxy_maximum
    best_loglik <- results$other$loglik
    output_list[["local_max"]] <- abs(proxy_max - best_loglik) > .001

    ## parameter estimates:
    estimates <- results$estimates
    # phi:
    for (k in 1:4) {
      for (i in 1:4) {
        for (j in 1:4) {
          out_name <- paste0("phi", i, j, "_k", k)
          # check if cluster AND factors exist:
          if (k <= n_clusters && i <= n_factors && j <= n_factors) {
            cluster_name <- paste0("cluster", k)
            param_name <- paste0("phi_f", i, "_f", j)
            output_list[[out_name]] <- estimates[[cluster_name]][[param_name]]
          } else {
            output_list[[out_name]] <- NA
          }
        }
      }
    }

    # zeta:
    for (k in 1:4) {
      for (i in 1:4) {
        for (j in i:4) {
          # (j in i:4) ensures j >= i (diagonal + lower triangle)
          out_name <- paste0("zeta", i, j, "_k", k)
          # check if cluster AND factors exist:
          if (k <= n_clusters && i <= n_factors && j <= n_factors) {
            cluster_name <- paste0("cluster", k)
            param_name <- paste0("zeta_f", i, "_f", j)
            output_list[[out_name]] <- estimates[[cluster_name]][[param_name]]
          } else {
            output_list[[out_name]] <- NA
          }
        }
      }
    }
  } else {
    # if step 3 was not successful, set all values to NA
    output_list[["duration"]] <- NA
    output_list[["nonconvergences"]] <- NA
    output_list[["ARI"]] <- NA
    output_list[["local_max"]] <- NA
    # phis:
    for (k in 1:4) {
      for (i in 1:4) {
        for (j in 1:4) {
          out_name <- paste0("phi", i, j, "_k", k)
          output_list[[out_name]] <- NA
        }
      }
    }

    # zetas:
    for (k in 1:4) {
      for (i in 1:4) {
        for (j in i:4) {  # ensures j >= i (diagonal + lower triangle)
          out_name <- paste0("zeta", i, j, "_k", k)
          output_list[[out_name]] <- NA
        }
      }
    }
  }

  #### add warnings and errors to output ####
  warnings_errors <- c("step1_warning", "step2_warning", "step3_warning",
                       "step1_error", "step2_error", "step3_error",
                       "step1_warning_text", "step2_warning_text", "step3_warning_text",
                       "step1_error_text", "step2_error_text", "step3_error_text")
  for (x in warnings_errors){
    output_list[[x]] <- get(x)
  }

  # reset the OpenMx parameters:
  mxOption(key = "Major Iterations", reset = TRUE)
  mxOption(key = "Calculate Hessian", reset = TRUE)
  mxOption(key = "Standard Errors", reset = TRUE)

  # remove all whitespace, linebreaks, and commata from error and warning strings
  text_elements <- grep("_text$", names(output_list), value = TRUE)
  output_list[text_elements] <- output_list[text_elements] |>
    map(~ {
      .x |>
        str_squish() |>
        str_replace_all(",", "")
    })

  #### write output file ####
  output_vector <- unlist(output_list, use.names = TRUE)

  # check if file exists
  if (!file.exists(outputfile)) {
    # if file does not yet exist
    write.table(t(output_vector), file = outputfile, append = FALSE, quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)
  } else {
    # lock the file to prevent multiple processes accessing it simultaneously
    lock <- flock::lock(outputfile)
    write.table(t(output_vector), file = outputfile, append = TRUE, quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)
    # unlock the file
    flock::unlock(lock)
  }

  if (verbose == TRUE) {
    print(paste("Simulation", pos, "completed at", Sys.time()))                 # prints a message when a replication is done (as a sign that R did not crash)
  }

  return(output_vector)
}