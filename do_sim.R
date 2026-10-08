#### This script defines the simulation function do_sim() ####

do_sim <- function(pos, cond, outputfile, verbose = FALSE) {
  # pos = position in the condition grid
  # cond = the condition grid
  # outputfile = file name for the output CSV file
  # verbose = if TRUE, prints a message after the iteration is finished

  output_list <- list(
    "iteration" = cond$iteration[pos],
    "replication" = cond$replication[pos]
  )

  # get condition levels and set seed:
  n_obs <- output_list[["n_obs"]] <- cond$n_obs[pos]
  n_clusters <- output_list[["n_clusters"]] <- cond$n_clusters[pos]
  n_factors <- output_list[["n_factors"]] <- cond$n_factors[pos]
  seed <- output_list[["seed"]] <- cond$seed[pos]
  set.seed(seed)

  #### set data generation parameters ####
  # number of individuals:
  n_persons <- 60
  ## regression parameters:
  # if two factors:
  if (n_factors == 2) {
    # Cluster 1: interconnected
    phimat_k1_pop <- matrix(c(0.3, 0.2, 0.2, 0.3), ncol = 2, byrow = TRUE)

    # Cluster 2: antagonistic
    phimat_k2_pop <- matrix(c(0.1, -0.2, -0.1, -0.1), ncol = 2, byrow = TRUE)

    if (n_clusters == 4) {
      # Cluster 3: mixed
      phimat_k3_pop <- matrix(c(0.3, -0.3, 0.1, 0.9), ncol = 2, byrow = TRUE)

      # Cluster 4: unconnected
      phimat_k4_pop <- matrix(c(0.6, 0.1, 0.1, 0.6), ncol = 2, byrow = TRUE)
    } else {
      # create matrices with NA if there's only two clusters:
      phimat_k3_pop <- phimat_k4_pop <- matrix(NA, nrow = 2, ncol = 2)
    }
  }

  # if four factors:
  if (n_factors == 4) {
    # Cluster 1: interconnected
    phimat_k1_pop <- matrix(
      c(
        0.3,
        0.2,
        0.2,
        0.2,
        0.2,
        0.4,
        0.2,
        0.2,
        0.2,
        0.2,
        0.3,
        0.2,
        0.2,
        0.2,
        0.2,
        0.4
      ),
      ncol = 4,
      byrow = TRUE
    )

    # Cluster 2: antagonistic
    phimat_k2_pop <- matrix(
      c(
        0.1,
        -0.2,
        -0.1,
        0.2,
        -0.1,
        -0.1,
        -0.1,
        -0.1,
        -0.2,
        -0.1,
        0.2,
        -0.1,
        -0.1,
        -0.2,
        -0.1,
        0.1
      ),
      ncol = 4,
      byrow = TRUE
    )

    if (n_clusters == 4) {
      # Cluster 3: mixed
      phimat_k3_pop <- matrix(
        c(
          0.3,
          -0.3,
          0.2,
          -0.2,
          0.1,
          0.9,
          -0.3,
          0,
          0.2,
          -0.2,
          0.4,
          0.2,
          -0.3,
          0.1,
          0.2,
          0.7
        ),
        ncol = 4,
        byrow = TRUE
      )

      # Cluster 4: unconnected
      phimat_k4_pop <- matrix(
        c(
          0.6,
          0.1,
          0,
          0,
          0.1,
          0.6,
          0.1,
          0.1,
          0,
          0.1,
          0.6,
          0.1,
          0.1,
          0,
          0,
          0.6
        ),
        ncol = 4,
        byrow = TRUE
      )
    } else {
      # create matrices with NA if there's only two clusters:
      phimat_k3_pop <- phimat_k4_pop <- matrix(NA, nrow = 4, ncol = 4)
    }
  }

  ## innovation variances
  # if two factors:
  if (n_factors == 2) {
    zetamat_k1_pop <- zetamat_k2_pop <- matrix(
      c(292, -104, -104, 46),
      ncol = 2,
      byrow = TRUE
    )

    if (n_clusters == 4) {
      zetamat_k3_pop <- zetamat_k4_pop <- zetamat_k1_pop
    } else {
      # create matrices with NA if there's only two clusters:
      zetamat_k3_pop <- zetamat_k4_pop <- matrix(NA, nrow = 2, ncol = 2)
    }
  }

  # if four factors:
  if (n_factors == 4) {
    zetamat_k1_pop <- zetamat_k2_pop <- matrix(
      c(
        293,
        -103,
        47,
        -22,
        -103,
        157,
        -27,
        19,
        47,
        -27,
        163,
        -7,
        -22,
        19,
        -7,
        47
      ),
      ncol = 4,
      byrow = TRUE
    )

    if (n_clusters == 4) {
      zetamat_k3_pop <- zetamat_k3_pop <- zetamat_k1_pop
    } else {
      # create matrices with NA if there's only two clusters:
      zetamat_k3_pop <- zetamat_k4_pop <- matrix(NA, nrow = 4, ncol = 4)
    }
  }

  ## grand means
  if (n_factors == 2) {
    grandmeans <- rep(0, 2)
  } else {
    grandmeans <- rep(0, 4)
  }
  mu_variance <- matrix(-15, nrow = n_factors, ncol = n_factors)
  if (n_factors == 2) {
    diag(mu_variance) <- c(60, 60)
  }
  if (n_factors == 4) {
    diag(mu_variance) <- c(60, 60, 60, 60)
  }

  ## measurement model parameters
  # loadings:
  lambda_f1_pop <- c(1, 0.9, 0.6, 0.7)
  lambda_f2_pop <- c(1, 1.5, 1.4, 1.6)
  if (n_factors == 4) {
    lambda_f3_pop <- c(1, 0.8, 0.5, 0.9)
    lambda_f4_pop <- c(1, 1.2, 1.3, 1.5)
  } else {
    lambda_f3_pop <- lambda_f4_pop <- rep(NA, 4)
  }

  if (n_factors == 2) {
    lambda <- list(
      f1 = lambda_f1_pop,
      f2 = lambda_f2_pop
    ) |>
      lavaan::lav_matrix_bdiag()
  }
  if (n_factors == 4) {
    lambda <- list(
      f1 = lambda_f1_pop,
      f2 = lambda_f2_pop,
      f3 = lambda_f3_pop,
      f4 = lambda_f4_pop
    ) |>
      lavaan::lav_matrix_bdiag()
  }

  # residual variances:
  theta_f1_pop <- c(67, 195, 428, 351)
  theta_f2_pop <- c(449, 204, 304, 369)
  if (n_factors == 4) {
    theta_f3_pop <- c(85, 187, 303, 330)
    theta_f4_pop <- c(361, 201, 304, 351)
  } else {
    theta_f3_pop <- theta_f4_pop <- rep(NA, 4)
  }

  if (n_factors == 2) {
    theta <- diag(c(theta_f1_pop, theta_f2_pop))
  }
  if (n_factors == 4) {
    theta <- diag(c(
      theta_f1_pop,
      theta_f2_pop,
      theta_f3_pop,
      theta_f4_pop
    ))
  }

  # intercepts:
  tau_f1_pop <- c(55, 53, 57, 59)
  tau_f2_pop <- c(25, 19, 22, 31)
  if (n_factors == 4) {
    tau_f3_pop <- c(43, 41, 45, 47)
    tau_f4_pop <- c(15, 9, 12, 21)
  } else {
    tau_f3_pop <- tau_f4_pop <- rep(NA, 4)
  }

  if (n_factors == 2) {
    tau <- c(tau_f1_pop, tau_f2_pop)
  }
  if (n_factors == 4) {
    tau <- c(tau_f1_pop, tau_f2_pop, tau_f3_pop, tau_f4_pop)
  }

  #### generate data ####
  ## create cluster assignment vector:
  clusterassignment_true <- rep(
    paste0("cluster", 1:n_clusters),
    length.out = n_persons
  ) |>
    # shuffle:
    sample(n_persons) |>
    factor()

  # create empty data frame for all (observed) items:
  eta_vars <- paste0("eta", 1:n_factors)
  eta <- data.frame(
    id = integer(),
    obs = integer(),
    k_true = factor(levels = levels(clusterassignment_true))
  )
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
      zetamat <- zetamat_k1_pop
    }
    if (k_i == "cluster2") {
      phimat <- phimat_k2_pop
      zetamat <- zetamat_k2_pop
    }
    if (k_i == "cluster3") {
      phimat <- phimat_k3_pop
      zetamat <- zetamat_k3_pop
    }
    if (k_i == "cluster4") {
      phimat <- phimat_k4_pop
      zetamat <- zetamat_k4_pop
    }

    # generate person-specific latent means:
    mu_i <- MASS::mvrnorm(1, mu = grandmeans, Sigma = mu_variance)

    eta_i <- sim_VAR(
      factors = n_factors,
      obs = n_obs,
      phi = phimat,
      zeta = zetamat,
      mu = mu_i,
      burn_in = 10
    )

    eta_i$id <- i
    eta_i$k_true <- k_i

    # add person-data to full data frame:
    eta <- dplyr::bind_rows(eta, eta_i)
  }

  # generate measurement error from residual variances:
  epsilon <- MASS::mvrnorm(
    nrow(eta),
    mu = rep(0, n_factors * 4),
    Sigma = theta,
    empirical = FALSE
  )

  # transform factor scores into observed scores:
  data <- t(tau + lambda %*% t(eta[, eta_vars])) +
    epsilon |>
      as.data.frame()
  colnames(data) <- paste0(
    "f",
    rep(1:n_factors, each = 4),
    "_v",
    rep(1:4, times = n_factors)
  )
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

  # zeta:
  for (k in 1:4) {
    # get zetamat of cluster k:
    zetamat <- get(paste0("zetamat_k", k, "_pop"))
    # loop over all entries, check if they exist, and extract value/set to NA:
    for (i in 1:4) {
      for (j in i:4) {
        # (j in i:4) ensures j >= i (diagonal + upper triangle)
        out_name <- paste0("zeta", i, j, "_k", k, "_pop")
        value <- ifelse(i <= n_factors && j <= n_factors, zetamat[i, j], NA)
        output_list[[out_name]] <- value
      }
    }
  }

  # loadings:
  for (f in 1:4) {
    lambda_f <- get(paste0("lambda_f", f, "_pop"))
    for (v in 1:4) {
      out_name <- paste0("lambda_f", f, "_v", v, "_pop")
      value <- ifelse(f <= n_factors, lambda_f[v], NA)
      output_list[[out_name]] <- value
    }
  }

  # residual variances:
  for (f in 1:4) {
    theta_f <- get(paste0("theta_f", f, "_pop"))
    for (v in 1:4) {
      out_name <- paste0("theta_f", f, "_v", v, "_pop")
      value <- ifelse(f <= n_factors, theta_f[v], NA)
      output_list[[out_name]] <- value
    }
  }

  # intercepts:
  for (f in 1:4) {
    tau_f <- get(paste0("tau_f", f, "_pop"))
    for (v in 1:4) {
      out_name <- paste0("tau_f", f, "_v", v, "_pop")
      value <- ifelse(f <= n_factors, tau_f[v], NA)
      output_list[[out_name]] <- value
    }
  }

  #### Step 1 ####
  start <- Sys.time()
  # MM in lavaan syntax:
  if (n_factors == 2) {
    model_step1 <- list(
      "
      f1 =~ f1_v1 + f1_v2 + f1_v3 + f1_v4
      ",
      "
      f2 =~ f2_v1 + f2_v2 + f2_v3 + f2_v4
      "
    )
  }
  if (n_factors == 4) {
    model_step1 <- list(
      "
      f1 =~ f1_v1 + f1_v2 + f1_v3 + f1_v4
      ",
      "
      f2 =~ f2_v1 + f2_v2 + f2_v3 + f2_v4
      ",
      "
      f3 =~ f3_v1 + f3_v2 + f3_v3 + f3_v4
      ",
      "
      f4 =~ f4_v1 + f4_v2 + f4_v3 + f4_v4
      "
    )
  }

  # run step 1:
  output_step1 <- run_step1(
    data = data,
    measurementmodel = model_step1,
    id = "id"
  )
  # extract error/warning messages (if applicable):
  step1_warning <- ifelse(rlang::is_empty(output_step1$warnings), FALSE, TRUE)
  step1_warning_text <- ifelse(
    step1_warning,
    paste(c(output_step1$warnings), collapse = "; "),
    ""
  )
  step1_error <- ifelse(rlang::is_empty(output_step1$result$error), FALSE, TRUE)
  step1_error_text <- ifelse(
    step1_error,
    paste(c(output_step1$result$error), collapse = "; "),
    ""
  )

  #### Step 2 ####
  # only proceed if there is no error in step 1:
  if (!step1_error) {
    ## extract MM estimates:
    results_MM <- output_step1$result$result$mm_output
    # loadings:
    for (f in 1:4) {
      if (f <= n_factors) {
        PE_f <- lavaan::parameterEstimates(results_MM[[f]])
        for (v in 1:4) {
          item_name <- paste0("f", f, "_v", v)
          out_name <- paste0("lambda_", item_name, "_est")
          output_list[[out_name]] <- PE_f$est[
            PE_f$lhs == paste0("f", f) &
              PE_f$op == "=~" &
              PE_f$rhs == item_name
          ]
        }
      } else {
        for (v in 1:4) {
          out_name <- paste0("lambda_f", f, "_v", v, "_est")
          output_list[[out_name]] <- NA
        }
      }
    }

    # residual variances:
    for (f in 1:4) {
      if (f <= n_factors) {
        PE_f <- lavaan::parameterEstimates(results_MM[[f]])
        for (v in 1:4) {
          item_name <- paste0("f", f, "_v", v)
          out_name <- paste0("theta_f", f, "_v", v, "_est")
          output_list[[out_name]] <- PE_f$est[
            PE_f$lhs == item_name &
              PE_f$op == "~~" &
              PE_f$rhs == item_name
          ]
        }
      } else {
        for (v in 1:4) {
          out_name <- paste0("theta_f", f, "_v", v, "_est")
          output_list[[out_name]] <- NA
        }
      }
    }

    # intercepts:
    for (f in 1:4) {
      if (f <= n_factors) {
        PE_f <- lavaan::parameterEstimates(results_MM[[f]])
        for (v in 1:4) {
          item_name <- paste0("f", f, "_v", v)
          out_name <- paste0("tau_f", f, "_v", v, "_est")
          output_list[[out_name]] <- PE_f$est[
            PE_f$lhs == item_name &
              PE_f$op == "~1"
          ]
        }
      } else {
        for (v in 1:4) {
          out_name <- paste0("tau_f", f, "_v", v, "_est")
          output_list[[out_name]] <- NA
        }
      }
    }

    output_step2 <- run_step2(
      step1output = output_step1$result$result
    )
    # extract error/warning messages (if applicable):
    step2_warning <- ifelse(rlang::is_empty(output_step2$warnings), FALSE, TRUE)
    step2_warning_text <- ifelse(
      step2_warning,
      paste(c(output_step2$warnings), collapse = "; "),
      ""
    )
    step2_error <- ifelse(
      rlang::is_empty(output_step2$result$error),
      FALSE,
      TRUE
    )
    step2_error_text <- ifelse(
      step2_error,
      paste(c(output_step2$result$error), collapse = "; "),
      ""
    )
  } else {
    step2_warning <- FALSE
    step2_warning_text <- "step1 not successful"
    step2_error <- FALSE
    step2_error_text <- "step1 not successful"
  }

  #### Step 3 ####
  # only proceed if there is no error in step 1 as well as step 2
  if (!step1_error && !step2_error) {
    A_matrix <- create_A(
      step2output = output_step2$result$result,
      random_intercept = TRUE
    )
    Q_matrix <- create_Q(
      step2output = output_step2$result$result,
      random_intercept = TRUE
    )

    output_step3 <- run_step3(
      step2output = output_step2$result$result,
      A = A_matrix,
      Q = Q_matrix,
      mixture = TRUE,
      n_clusters = n_clusters,
      n_starts = 25,
      n_best_starts = 15,
      verbose = FALSE
    )
    duration <- difftime(Sys.time(), start, unit = "s")
    # extract error/warning messages (if applicable):
    step3_warning <- ifelse(rlang::is_empty(output_step3$warnings), FALSE, TRUE)
    step3_warning_text <- ifelse(
      step3_warning,
      paste(c(output_step3$warnings), collapse = "; "),
      ""
    )
    step3_error <- ifelse(
      rlang::is_empty(output_step3$result$error),
      FALSE,
      TRUE
    )
    step3_error_text <- ifelse(
      step3_error,
      paste(c(output_step3$result$error), collapse = "; "),
      ""
    )
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
    new_labels <- adjust_labels(
      estimated_clusters = results$posterior_probabilities$modal,
      true_clusters = clusterassignment_true,
      n_clusters = n_clusters
    )
    # swap the labels accordingly in the output of step 3:
    rownames(results$estimates) <-
      rownames(results$standarderrors) <-
        names(results$model) <-
          new_labels
    clusterassignment_estimated <- results$posterior_probabilities$modal |>
      factor(labels = new_labels) |>
      as.character() |>
      factor(levels = levels(clusterassignment_true))

    ## duration, number of non-convergences, status check:
    output_list[["duration"]] <- duration |>
      as.numeric()
    output_list[["nonconvergences"]] <- results$n_nonconverged
    output_list[["status_k1"]] <- results$model$cluster1$output$status$code
    output_list[["status_k2"]] <- results$model$cluster2$output$status$code
    if (n_clusters == 4) {
      output_list[["status_k3"]] <- results$model$cluster3$output$status$code
      output_list[["status_k4"]] <- results$model$cluster4$output$status$code
    } else {
      output_list[["status_k3"]] <- NA
      output_list[["status_k4"]] <- NA
    }

    ## ARI:
    output_list[["ARI"]] <- mcclust::arandi(
      clusterassignment_true,
      clusterassignment_estimated,
      adjust = TRUE
    )

    ## parameter estimates:
    estimates <- results$estimates
    # phi:
    for (k in 1:4) {
      for (i in 1:4) {
        for (j in 1:4) {
          out_name <- paste0("phi", i, j, "_k", k, "_est")
          # check if cluster AND factors exist:
          if (k <= n_clusters && i <= n_factors && j <= n_factors) {
            cluster_name <- paste0("cluster", k)
            param_name <- paste0("phi_f", i, "_f", j)
            output_list[[out_name]] <- estimates[cluster_name, param_name]
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
          # (j in i:4) ensures j >= i (diagonal + upper triangle)
          out_name <- paste0("zeta", i, j, "_k", k, "_est")
          # check if cluster AND factors exist:
          if (k <= n_clusters && i <= n_factors && j <= n_factors) {
            cluster_name <- paste0("cluster", k)
            param_name <- paste0("zeta_f", i, "_f", j)
            output_list[[out_name]] <- estimates[cluster_name, param_name]
          } else {
            output_list[[out_name]] <- NA
          }
        }
      }
    }

    ## standard errors:
    standarderrors <- results$standarderrors
    # phi:
    for (k in 1:4) {
      for (i in 1:4) {
        for (j in 1:4) {
          out_name <- paste0("phi", i, j, "_k", k, "_se")
          # check if cluster AND factors exist:
          if (k <= n_clusters && i <= n_factors && j <= n_factors) {
            cluster_name <- paste0("cluster", k)
            param_name <- paste0("phi_f", i, "_f", j)
            output_list[[out_name]] <- standarderrors[cluster_name, param_name]
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
          # (j in i:4) ensures j >= i (diagonal + upper triangle)
          out_name <- paste0("zeta", i, j, "_k", k, "_se")
          # check if cluster AND factors exist:
          if (k <= n_clusters && i <= n_factors && j <= n_factors) {
            cluster_name <- paste0("cluster", k)
            param_name <- paste0("zeta_f", i, "_f", j)
            output_list[[out_name]] <- standarderrors[cluster_name, param_name]
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
    output_list[["status_k1"]] <- NA
    output_list[["status_k2"]] <- NA
    output_list[["status_k3"]] <- NA
    output_list[["status_k4"]] <- NA
    output_list[["ARI"]] <- NA
    output_list[["solution_loglik"]] <- NA
    output_list[["proxy_loglik"]] <- NA

    ## estimates:
    # phis:
    for (k in 1:4) {
      for (i in 1:4) {
        for (j in 1:4) {
          out_name <- paste0("phi", i, j, "_k", k, "_est")
          output_list[[out_name]] <- NA
        }
      }
    }

    # zetas:
    for (k in 1:4) {
      for (i in 1:4) {
        for (j in i:4) {
          # ensures j >= i (diagonal + upper triangle)
          out_name <- paste0("zeta", i, j, "_k", k, "_est")
          output_list[[out_name]] <- NA
        }
      }
    }

    ## standard errors:
    # phis:
    for (k in 1:4) {
      for (i in 1:4) {
        for (j in 1:4) {
          out_name <- paste0("phi", i, j, "_k", k, "_se")
          output_list[[out_name]] <- NA
        }
      }
    }

    # zetas:
    for (k in 1:4) {
      for (i in 1:4) {
        for (j in i:4) {
          # ensures j >= i (diagonal + upper triangle)
          out_name <- paste0("zeta", i, j, "_k", k, "_se")
          output_list[[out_name]] <- NA
        }
      }
    }
  }

  #### add warnings and errors to output ####
  warnings_errors <- c(
    "step1_warning",
    "step2_warning",
    "step3_warning",
    "step1_error",
    "step2_error",
    "step3_error",
    "step1_warning_text",
    "step2_warning_text",
    "step3_warning_text",
    "step1_error_text",
    "step2_error_text",
    "step3_error_text"
  )
  for (x in warnings_errors) {
    output_list[[x]] <- get(x)
  }

  # remove all whitespace, linebreaks, and commata from error and warning strings
  text_elements <- grep("_text$", names(output_list), value = TRUE)
  output_list[text_elements] <- output_list[text_elements] |>
    map(
      ~ {
        .x |>
          str_squish() |>
          str_replace_all(",", "")
      }
    )

  #### write output file ####
  output_vector <- unlist(output_list, use.names = TRUE)

  # check if file exists
  if (!file.exists(outputfile)) {
    # if file does not yet exist
    write.table(
      t(output_vector),
      file = outputfile,
      append = FALSE,
      quote = FALSE,
      sep = ",",
      row.names = FALSE,
      col.names = TRUE
    )
  } else {
    # lock the file to prevent multiple processes accessing it simultaneously
    lock <- flock::lock(outputfile)
    write.table(
      t(output_vector),
      file = outputfile,
      append = TRUE,
      quote = FALSE,
      sep = ",",
      row.names = FALSE,
      col.names = FALSE
    )
    # unlock the file
    flock::unlock(lock)
  }

  if (verbose == TRUE) {
    print(paste("Simulation", pos, "completed at", Sys.time())) # prints a message when a replication is done (as a sign that R did not crash)
  }

  return(output_vector)
}
