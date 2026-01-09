step3 <- function(step2output, n_clusters, n_starts = 25, n_best_starts = 5,
                  maxit = 100,
                  convergence_criterion = 1e-6,
                  true_clusters, verbose = FALSE) {

  #### 1) Preparations ####
  data <- step2output$data
  lambda_star <- step2output$lambda_star
  theta_star <- step2output$theta_star
  # names and number of latent factors:
  factors <- step2output$other$factors
  n_factors <- length(factors)
  # create names of factor score variables (single indicators)
  factors_ind <- paste0(factors, "_ind")
  # vector and number of unique ids
  id <- step2output$other$id
  unique_ids <- unique(data[, id])
  n_persons <- length(unique_ids)

  #### 2) data manipulation ####
  ## rename the factor score variables in the data
  # to use them as indicators of the latent variables
  data <- data |>
    dplyr::rename_with(~ factors_ind, tidyselect::all_of(factors))

  #### 3) create OpenMx matrices ####
  # number of latent variables in model is number of factors times 2
  # due to the random intercepts:
  xdim <- n_factors * 2
  # exogenous covariates (ignored)
  udim <- 1
  # number of indicators (i.e., factor score variables)
  ydim <- n_factors

  ## A matrices (= dynamics)
  # create matrix that indicates free parameters:
  free_phi <- matrix(TRUE, nrow = n_factors, ncol = n_factors)
  # expand free matrix with random intercept specification
  others <- matrix(FALSE, nrow = n_factors, ncol = n_factors)
  free_phi <- rbind(cbind(free_phi, others),
                    cbind(others, others))

  # create matrix of labels for regression coefficients:
  labels_phi <- outer(factors, factors,
                      FUN = function(i, j) paste0("phi_", i, "_", j))
  # expand label matrix with random intercept specification
  others <- matrix(NA, nrow = n_factors, ncol = n_factors)
  labels_phi <- rbind(cbind(labels_phi, others),
                      cbind(others, others))

  # create matrix of starting values for the free parameters
  values_phi <- matrix(.1, nrow = n_factors, ncol = n_factors)
  # expand values matrix with random intercept specification:
  values_phi <- rbind(cbind(values_phi,
                            matrix(0, nrow = n_factors, ncol = n_factors)),
                      cbind(matrix(0, nrow = n_factors, ncol = n_factors),
                            diag(n_factors)))

  amat <- OpenMx::mxMatrix("Full", name = "A",
                           nrow = xdim, ncol = xdim,
                           free = free_phi,
                           values = values_phi,
                           labels = labels_phi,
                           lbound = NA,
                           ubound = NA)

  ## B matrix (= exogenous covariates on latent constructs)
  bmat <- OpenMx::mxMatrix("Zero", name = "B",
                           nrow = xdim, ncol = udim)

  ## D matrix (= exogenous covariates on observed variables)
  dmat <- OpenMx::mxMatrix("Zero", name = "D",
                           nrow = ydim, ncol = udim)

  ## Q matrix (= innovation (co)variances)
  # create matrix that indicates free parameters:
  free_zeta <- matrix(TRUE, nrow = n_factors, ncol = n_factors)
  # expand free matrix with random intercept specification
  others <- matrix(FALSE, nrow = n_factors, ncol = n_factors)
  free_zeta <- rbind(cbind(free_zeta, others),
                     cbind(others, others))

  # create matrix of labels for innovation covariances:
  labels_zeta <- matrix(NA, n_factors, n_factors)
  # fill the matrix with symmetric labels
  for (i in 1:n_factors) {
    for (j in 1:n_factors) {
      labels_zeta[i, j] <- paste0("zeta_",
                                  factors[min(i, j)],
                                  "_",
                                  factors[max(i, j)])
    }
  }
  # expand label matrix with random intercept specification
  others <- matrix(NA, nrow = n_factors, ncol = n_factors)
  labels_zeta <- rbind(cbind(labels_zeta, others),
                       cbind(others, others))

  # create matrix of starting values for the free parameters
  values_zeta <- matrix(.3, nrow = n_factors, ncol = n_factors)
  diag(values_zeta) <- 1
  # expand startvalue matrix with random intercept specification
  others <- matrix(0, nrow = n_factors, ncol = n_factors)
  values_zeta <- rbind(cbind(values_zeta, others),
                       cbind(others, others))

  qmat <- OpenMx::mxMatrix("Full", name = "Q",
                           nrow = xdim, ncol = xdim,
                           free = free_zeta,
                           values = values_zeta,
                           labels = labels_zeta,
                           lbound = NA,
                           ubound = NA)

  ## x0 and P0 (= initial means and (co)variances of the latent states)
  # x0:
  xmat <- OpenMx::mxMatrix("Full", name = "x0",
                           nrow = xdim, ncol = 1,
                           free = FALSE,
                           values = 0)

  #P0:
  # create matrix that indicates free parameters
  free_P0 <- matrix(TRUE, nrow = n_factors, ncol = n_factors)
  # expand with random intercept specification:
  others <- matrix(FALSE, nrow = n_factors, ncol = n_factors)
  free_P0 <- rbind(cbind(free_P0, others),
                   cbind(others, free_P0))

  # create matrix of labels for P0:
  labels_P0 <- matrix(NA, n_factors, n_factors)
  # fill the matrix with symmetric labels
  for (i in 1:n_factors) {
    for (j in 1:n_factors) {
      labels_P0[i, j] <- paste0("P0_",
                                factors[min(i, j)],
                                "_",
                                factors[max(i, j)])
    }
  }
  # expand with random intercept specification:
  others <- matrix(NA, nrow = n_factors, ncol = n_factors)
  intercepts <- matrix(NA, n_factors, n_factors)
  for (i in 1:n_factors) {
    for (j in 1:n_factors) {
      intercepts[i, j] <- paste0("P0_icp_",
                                 factors[min(i, j)],
                                 "_",
                                 factors[max(i, j)])
    }
  }
  labels_P0 <- rbind(cbind(labels_P0, others),
                     cbind(others, intercepts))

  # create matrix of starting values:
  values_P0 <- matrix(10, nrow = n_factors, ncol = n_factors)
  diag(values_P0) <- 100
  others <- matrix(0, nrow = n_factors, ncol = n_factors)
  values_P0 <- rbind(cbind(values_P0, others),
                     cbind(others, values_P0))

  pmat <- OpenMx::mxMatrix("Full", name = "P0",
                           nrow = xdim, ncol = xdim,
                           free = free_P0,
                           values = values_P0,
                           labels = labels_P0)

  # u (= covariates)
  umat <- OpenMx::mxMatrix("Zero", name = "u",
                           nrow = udim, ncol = 1)

  #### 4) create OpenMx models ####
  # create a list of models (one for each individual):
  personmodelnames <- paste0("id_", unique_ids)
  # vector of objective strings, needed later to create custom fit functions:
  objectives <- paste0(personmodelnames, ".objective")

  personmodel_list <- vector(mode = "list", length = n_persons)
  for (i in 1:n_persons) {
    # C matrix (= factor loadings, here fixed to lambda_star)
    values <- matrix(0, nrow = n_factors, ncol = n_factors)
    diag(values) <- lambda_star[i, ]
    cmat <- OpenMx::mxMatrix("Full", name = "C",
                             nrow = ydim, ncol = xdim,
                             free = FALSE,
                             values = values,
                             dimnames = list(factors_ind,
                                             c(paste0(factors),
                                               paste0("intercept_", factors))
                             )
    )
    # R matrix (= measurement noise, here fixed to theta_star)
    rmat <- OpenMx::mxMatrix("Diag", name = "R",
                             nrow = ydim, ncol = ydim,
                             free = FALSE,
                             values = theta_star[i, ]
    )

    # crate model and store in list
    personmodel_list[[i]] <- OpenMx::mxModel(name = personmodelnames[i],
                                             amat, bmat, cmat, dmat,
                                             qmat, rmat, xmat, pmat,
                                             umat,
                                             OpenMx::mxExpectationStateSpace(
                                               "A", "B", "C", "D",
                                               "Q", "R", "x0", "P0",
                                               "u"),
                                             OpenMx::mxFitFunctionML(),
                                             OpenMx::mxData(
                                               data[data[, id] == i,
                                                    factors_ind],
                                               "raw")
    )
  }
  names(personmodel_list) <- personmodelnames

  #### 5) mixture modeling ####
  # provide seeds for the multiple starts (for replicability)
  all_seeds <- sample(1:100000000, n_starts)
  estimation_start <- Sys.time()
  # The estimation happens in two parts:
  # First, the algorithm aims for a stable clustering with an imprecise M-Step
  # after a stable clustering is achieved, a final precise M-Step is performed
  # Then, the n_best_starts best starts are completed

  ## Part 1: all random starts
  all_starts <- vector(mode = "list", length = n_starts)
  default_iterations <- OpenMx::mxOption(key = "Major iterations")
  for (random_start in 1:n_starts) {
    # set seed:
    set.seed(all_seeds[random_start])
    # create random hard cluster assignment:
    post <- matrix(0, nrow = n_persons, ncol = n_clusters)
    for (person in 1:n_persons) {
      if (person <= n_clusters) {
        # assign the first persons to certain clusters to avoid empty clusters
        # (first person in first cluster, second person in second cluster, ...)
        # until every cluster has one individual
        post[person, person] <- 1
      } else {
        post[person, sample(1:n_clusters, 1)] <- 1
      }
    }

    # reset the observed data LL:
    observed_data_LL0 <- -Inf
    # reset the check for stable clustering:
    stable_clustering <- FALSE
    # reduce the number of iterations in OpenMx (to speed it up):
    OpenMx::mxOption(key = "Major iterations", value = 3)

    # loop over iterations until convergence or max iterations are reached
    for (it in 1:maxit) {
      ## E-step: update class membership and class proportions
      if (it > 1) {
        post0 <- post
        # update posteriors:
        post <- EStep(pi_ks = class_proportions, ngroup = n_persons,
                      nclus = n_clusters, loglik = personLL)
        # compute (absolute) differences in posteriors
        delta <- abs(post - post0)

        if (verbose) {
          print(paste0("Highest person-delta in posterior probabilities: ",
                       max(rowSums(delta)) |> round(3),
                       ". Average person-delta: ",
                       mean(rowSums(delta)) |> round(3),
                       "."))
        }

        # check if stable clustering has been achieved:
        if (mean(rowSums(delta)) < .1) {
          stable_clustering <- TRUE
          OpenMx::mxOption(key = "Major iterations", value = default_iterations)
          if (verbose) {
            print("Switched to full M-Step.")
          }
        }
      }

      # compute class proportions:
      class_proportions <- colMeans(post)

      ## M-step: fitting SSM model and update parameter estimates
      # get start values from the coefficients from previous iteration:
      if (it > 1) {
        startvalues <- clustermodels_run |>
          purrr::map(coef)
      } else {
        # in the first iteration, generate random start values in each cluster
        startvalues <- generate_startvalues(n_clusters = n_clusters,
                                            n_factors = n_factors,
                                            labels_phi = labels_phi,
                                            labels_zeta = labels_zeta)
      }

      clustermodels_run <- MStep(n_clusters = n_clusters,
                                 weights = post,
                                 objectives = objectives,
                                 model_list = personmodel_list,
                                 startvalues = startvalues,
                                 verbose = verbose)

      # obtain person-wise LL in a n_persons x n_clusters matrix:
      personLL <- clustermodels_run |>
        purrr::map(~ purrr::map_dbl(.x$submodels, ~ .x$fitfunction$result)) |>
        simplify2array()
      personLL <- personLL / (-2)

      # compute observed-data log likelihood from person-wise LL:
      observed_data_LL <- compute_observed_data_LL(personLL = personLL,
                                                   class_proportions = class_proportions)

      if (verbose) {
        if (it != 1) {
          print(paste0("Iteration: ", it,
                       ". Log Likelihood: ", round(observed_data_LL, 4),
                       ". Change: ", round(observed_data_LL - observed_data_LL0, 6),
                       "."))
        } else {
          print(paste0("Iteration: ", it,
                       ". Log Likelihood: ", round(observed_data_LL, 4), "."))
        }
      }

      # check stable clustering and break loop if applicable:
      if (stable_clustering) {
        if (verbose) {
          print(paste("Start", random_start, " has a stable clustering."))
          if (random_start < n_starts) {
            print("Proceeding to next start.")
          }
        }
        break
      }

      observed_data_LL0 <- observed_data_LL

      # check if maxit has been reached:
      if (verbose && it == maxit) {
        print(paste("Start", random_start,
                    " did not arrive at a stable clustering."))
        if (random_start < n_starts) {
          print("Proceeding to next start.")
        }
      }

    }

    # save relevant objects in list:
    all_starts[[random_start]] <- list("observed_data_LL" = observed_data_LL,
                                       "class_proportions" = class_proportions,
                                       "personLL" = personLL,
                                       "clustermodels_run" = clustermodels_run,
                                       "seed" = all_seeds[random_start])
  }

  if (verbose) {
    print(paste("Setup finished. Finishing", n_best_starts, "best starts now."))
  }

  ## Part 2: complete the best starts
  # extract the best starts:
  best_starts <- all_starts |>
    purrr::map_dbl(~ .x$observed_data_LL) |>
    order(decreasing = TRUE) |>
    head(n_best_starts)

  best_loglik <- -Inf
  nonconvergences <- 0
  # reset the maximum number of iterations in OpenMx (for greater precision):
  OpenMx::mxOption(key = "Major iterations", value = default_iterations)
  for (random_start in 1:n_best_starts) {
    start_number <- best_starts[random_start]

    # extract outputs from corresponding start in Part 1:
    observed_data_LL0 <- all_starts[[start_number]]$observed_data_LL
    class_proportions <- all_starts[[start_number]]$class_proportions
    personLL <- all_starts[[start_number]]$personLL
    clustermodels_run <- all_starts[[start_number]]$clustermodels_run


    # loop over iterations until convergence or max iterations are reached
    for (it in 1:maxit) {
      post <- EStep(pi_ks = class_proportions, ngroup = n_persons,
                    nclus = n_clusters, loglik = personLL)

      # compute class proportions:
      class_proportions <- colMeans(post)

      ## M-step: fitting SSM model and update parameter estimates
      # get start values from the coefficients from previous iteration:
      startvalues <- clustermodels_run |>
        purrr::map(coef)

      clustermodels_run <- MStep(n_clusters = n_clusters,
                                 weights = post,
                                 objectives = objectives,
                                 model_list = personmodel_list,
                                 startvalues = startvalues,
                                 verbose = verbose)

      # obtain person-wise LL in a n_persons x n_clusters matrix:
      personLL <- clustermodels_run |>
        purrr::map(~ purrr::map_dbl(.x$submodels, ~ .x$fitfunction$result)) |>
        simplify2array()
      personLL <- personLL / (-2)

      # compute observed-data log likelihood from person-wise LL:
      observed_data_LL <- compute_observed_data_LL(personLL = personLL,
                                                   class_proportions = class_proportions)
      if (verbose) {
        print(paste0("Iteration: ", it,
                     ". Log Likelihood: ", round(observed_data_LL, 4),
                     ". Change: ", round(observed_data_LL - observed_data_LL0, 6),
                     "."))
      }

      # check convergence and break loop if applicable:
      if ((observed_data_LL - observed_data_LL0) < convergence_criterion) {
        if (verbose) {
          print("Convergence achieved.")
        }
        break
      }

      observed_data_LL0 <- observed_data_LL

      # check if maximum number of iterations has been reached:
      if (it == maxit) {
        if (verbose) {
          print(paste("Max iterations reached without convergence. Start:", random_start))
          if (random_start < n_best_starts) {
            print("Proceeding to next start.")
          }
        }
        nonconvergences <- nonconvergences + 1
      }
    }

    # check if the new fit is better than the previous best one
    if (observed_data_LL > best_loglik) {
      best_startnumber <- start_number
      best_seed <- all_seeds[start_number]
      best_loglik <- observed_data_LL
      best_post <- post
      best_models <- clustermodels_run
    }

    if (verbose) {
      print(paste("Start", random_start, "out of", n_best_starts, "best starts completed."))
      if (random_start < n_best_starts) {
        print("Proceeding to next start.")
      }
    }

  }

  #### 6) find proxy maximum ####
  # --> fit the model with true cluster memberships as starting "post"

  # create a matrix of 0s
  # then replace the value in the column of the cluster to which an individual belongs with 1
  post <- matrix(0, nrow = n_persons, ncol = n_clusters)
  post[cbind(1:n_persons, as.integer(true_clusters))] <- 1


  observed_data_LL0 <- -Inf
  for (it in 1:maxit) {
    ## E-step: update class membership and class proportions
    if (it > 1) {
      # update posteriors:
      post <- EStep(pi_ks = class_proportions, ngroup = n_persons,
                    nclus = n_clusters, loglik = personLL)
    }

    # update class proportions:
    class_proportions <- colMeans(post)

    ## M-step: fitting SSM model and update parameter estimates
    # get start values from the coefficients from previous iteration:
    if (it > 1) {
      startvalues <- clustermodels_run |>
        purrr::map(coef)
    } else {
      # in the first iteration, generate random start values in each cluster
      startvalues <- generate_startvalues(n_clusters = n_clusters,
                                          n_factors = n_factors,
                                          labels_phi = labels_phi,
                                          labels_zeta = labels_zeta)
    }

    clustermodels_run <- MStep(n_clusters = n_clusters,
                               weights = post,
                               objectives = objectives,
                               model_list = personmodel_list,
                               startvalues = startvalues)

    # obtain person-wise LL in a n_persons x n_clusters matrix:
    personLL <- clustermodels_run |>
      purrr::map(~ purrr::map_dbl(.x$submodels, ~ .x$fitfunction$result)) |>
      simplify2array()
    personLL <- personLL / (-2)

    # compute observed-data log likelihood from person-wise LL:
    observed_data_LL <- compute_observed_data_LL(personLL = personLL,
                                                 class_proportions = class_proportions)

    if (verbose) {
      if (it != 1) {
        print(paste0("Iteration: ", it,
                     ". Log Likelihood: ", round(observed_data_LL, 4),
                     ". Change: ", round(observed_data_LL - observed_data_LL0, 6),
                     "."))
      } else {
        print(paste0("Iteration: ", it,
                     ". Log Likelihood: ", round(observed_data_LL, 4),
                     "."))
      }

    }

    # check convergence and break loop if applicable:
    if (it > 1 && (observed_data_LL - observed_data_LL0) < 1.0e-6) {
      proxy_maximum <- observed_data_LL
      if (verbose) {
        print("Convergence achieved.")
      }
      break
    }

    observed_data_LL0 <- observed_data_LL
    # check if maximum number of iterations has been reached:
    if (it == maxit) {
      proxy_maximum <- NA
      if (verbose) {
        print(paste("Max iterations reached without convergence."))
      }
    }
  }


  #### 7) build output ####
  estimates <- purrr::map(best_models, coef)
  loglik <- best_loglik
  post <- as.data.frame(best_post)
  rownames(post) <- unique_ids
  colnames(post) <- paste0("cluster", 1:n_clusters)
  modal_assignment <- t(apply(post, MARGIN = 1, function(x) ifelse(x == max(x), 1, 0)))
  class_proportions <- colMeans(post)

  clustering <- list("class_proportions" = class_proportions,
                     "posterior_prob" = post,
                     "modal_assignment" = modal_assignment)

  other <- list("loglik" = loglik,
                "nonconvergences" = nonconvergences,
                "proxy_maximum" = proxy_maximum,
                "best_startnumber" = best_startnumber,
                "best_seed" = best_seed,
                "all_seeds" = all_seeds)

  output <- list("data" = data,
                 "estimates" = estimates,
                 "clustering" = clustering,
                 "other" = other)


  return(output)
}