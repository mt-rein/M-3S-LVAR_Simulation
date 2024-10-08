step3 <- function(step2output, n_clusters, n_starts = 25, n_best_starts = 5,
                  maxit = 100,
                  true_clusters, verbose = FALSE){
  # step2output:
  #   the object that was generated using the step2() function
  # ...
  
  
  #### FOR TESTING ####
  # step2output = output_step2$result$result
  # n_clusters = n_k
  # n_starts = 20
  # n_best_starts = 5
  # maxit = 100
  # true_clusters = clusterassignment_true
  # verbose = TRUE
  
  #### 1) Preparations ####
  ## extract objects from step 1 output:
  data <- step2output$data
  lambda_star <- step2output$lambda_star
  theta_star <- step2output$theta_star
  
  ## generate objects for later use
  factors <- step2output$other$factors                                          # names of latent factors
  factors_ind <- paste0(factors, "_ind")                                        # names of factor score variables (single indicators)
  id <- step2output$other$id                                                    # name of the variable that served as cluster variable in step1
  unique_ids <- unique(data[, id])                                              # vector of unique ids
  n_persons <- length(unique_ids)                                               # number of individuals
  
  #### 2) data manipulation ####
  ## rename the factor score variables in the data
  # to use them as indicators of the latent variables
  data <- data |> 
    dplyr::rename_with(~ factors_ind, all_of(factors))
  
  #### 3) create OpenMx matrices ####
  xdim <- length(factors)*2 # number of latent constructs in model is number of factors times 2 due to the random intercepts
  udim <- 1 # exogenous covariates (ignored so far)
  ydim <- length(factors) # number of indicators (i.e., factor score variables)
  
  ## A matrices (= dynamics)
  amat <- mxMatrix(type = "Full", nrow = xdim, ncol = xdim,
                   name = "A",
                   free = c(TRUE, TRUE, FALSE, FALSE,
                            TRUE, TRUE, FALSE, FALSE,
                            FALSE, FALSE, FALSE, FALSE,
                            FALSE, FALSE, FALSE, FALSE),
                   values = c(.1, .1, 1, 0,
                              .1, .1, 0, 1,
                              0, 0, 1, 0,
                              0, 0, 0, 1),
                   labels = c("phi11", "phi12", NA, NA,
                              "phi21", "phi22", NA, NA,
                              NA, NA, NA, NA,
                              NA, NA, NA, NA),
                   lbound= c(-.9, -.9, NA, NA,
                             -.9, -.9, NA, NA,
                             NA, NA, NA, NA,
                             NA, NA, NA, NA),
                   ubound= c(.9, .9, NA, NA,
                             .9, .9, NA, NA,
                             NA, NA, NA, NA,
                             NA, NA, NA, NA),
                   byrow = TRUE)
  
  # B matrix (= exogenous covariates on latent constructs)
  bmat <- mxMatrix('Zero', nrow = xdim, ncol = udim,
                   name='B')
  
  
  
  # D matrix (= exogenous covariates on observed variables)
  dmat <- mxMatrix('Zero', nrow = ydim, ncol = udim,
                   name='D')
  
  # Q matrix (= innovation (co)variances)
  qmat <- mxMatrix("Symm", xdim, xdim,
                   name = "Q",
                   free = c(TRUE,
                            TRUE, TRUE,
                            FALSE, FALSE, FALSE,
                            FALSE, FALSE, FALSE, FALSE),
                   values = c(1,
                              .3, 1,
                              0, 0, 0,
                              0, 0, 0, 0),
                   labels = c("zeta1",
                              "zeta12", "zeta2",
                              NA, NA, NA,
                              NA, NA, NA, NA),
                   lbound= c(1e-6,
                             NA, 1e-6,
                             NA, NA, NA,
                             NA, NA, NA, NA),
                   byrow = TRUE)
  
  
  
  
  # x0 and P0 (= initial values and (co)variances of the latent constructs)
  xmat <- mxMatrix('Full', nrow = xdim, ncol = 1,
                   name='x0',
                   free = FALSE,
                   values = c(0, 0, 0, 0),
                   labels = c(paste0("ini_", factors),
                              paste0("m_intercept_", factors))
  )
  
  pmat <- mxMatrix('Symm', nrow = xdim, ncol = xdim,
                   name='P0',
                   free = FALSE,
                   values = c(1,
                              .5, 1,
                              0, 0, 1,
                              0, 0, .5, 1),
                   # lbound = c(1e-2,
                   #            NA, 1e-2,
                   #            NA, NA,  1e-2,
                   #            NA, NA, NA,  1e-2),
                   labels = c("P0_f1",
                              "P0_f1f2", "P0_f2",
                              NA, NA, "P0_icp1",
                              NA, NA, "P0_icp12", "P0_icp2"),
                   byrow = TRUE)
  
  # u (= covariates)
  umat <- mxMatrix('Zero', nrow = udim, ncol = 1, name='u')
  
  #### 4) create OpenMx models ####
  # create a list of models (one for each individual):
  personmodelnames <- paste0("id_", unique_ids)
  objectives <- paste0(personmodelnames, ".objective")                                # vector of objective strings, needed later to create custom fit functions
  
  personmodel_list <- vector(mode = "list", length = n_persons)
  for(i in 1:n_persons){
    # C matrix (= factor loadings, here fixed to lambda_star)
    cmat <- mxMatrix('Full', nrow = ydim, ncol = xdim,
                     name='C',
                     free = FALSE,
                     values = c(lambda_star[i, 1], 0, 0, 0,
                                0, lambda_star[i, 2], 0, 0),
                     labels = NA,
                     byrow = TRUE,
                     dimnames = list(factors_ind, c(paste0(factors),
                                                    paste0("intercept_", factors))
                     )
    )
    # R matrix (= measurement noise, here fixed to theta_star)
    rmat <- mxMatrix('Diag', nrow = ydim, ncol = ydim,
                     name = 'R',
                     free = FALSE,
                     values = theta_star[i,]
    )
    
    personmodel_list[[i]] <- mxModel(name = personmodelnames[i],
                                     amat, bmat, cmat, dmat,
                                     qmat, rmat, xmat, pmat, 
                                     umat,
                                     mxExpectationStateSpace('A', 'B', 'C', 'D', 
                                                             'Q', 'R', 'x0', 'P0', 
                                                             'u'),
                                     mxFitFunctionML(),
                                     mxData(data[data[, id] == i, factors_ind], 
                                            'raw'))
  }
  names(personmodel_list) <- personmodelnames
  
  #### 5) mixture modeling ####
  #!!!! work on the whole replicability thing!!!!
  # provide seeds for the multiple starts (for replicability)
  seeds <- sample(1:100000000, n_starts)
  estimation_start <- Sys.time()
  # The estimation happens in two parts: First, a single iteration of the EM algorithm is run for each random start
  # the n_best_starts best starts are then completed
  
  ## Part 1: all random starts
  all_starts <- vector(mode = "list", length = n_starts)                        # create empty vector for the output
  mxOption(key = "Major iterations", value = 2)                                 # reduce the number of iterations when estimating the parameters in OpenMx (to speed it up)
  for(random_start in 1:n_starts){
    # set seed:
    set.seed(seeds[random_start])
    # create random cluster assignment:
    post <- matrix(0, nrow = n_persons, ncol = n_clusters)
    for(person in 1:n_persons){
      if(person <= n_clusters){
        # assign the first persons to certain clusters to avoid empty clusters
        # (i.e., first person in first cluster, second person in second cluster, etc)
        # until every cluster has one individual
        post[person, person] <- 1
      } else {
        post[person, sample(1:n_clusters, 1)] <- 1
      }
    }
    observed_data_LL0 <- -Inf                                                   # reset the observed data LL
    # loop over iterations until convergence or max iterations are reached
    for(it in 1:maxit){
      ## E-step: update class membership (skip in first iteration) and class proportions
      if(it >1){
        # update posteriors:
        post <- EStep(pi_ks = class_proportions, ngroup = n_persons, 
                      nclus = n_clusters, loglik = casewiseLL)
      }
      
      # compute class proportions:
      class_proportions <- colMeans(post)
      
      ## M-step: fitting SSM model and update parameter estimates
      clustermodels <- vector(mode = "list", length = n_clusters)
      for (i in 1:n_clusters){
        clustername <- paste0("model_k", i)
        model <- create_model(clustername,
                              weights = post[, i],
                              objectives = objectives,
                              model_list = personmodel_list)
        if(it == 1){
          # in first iteration, generate random starting values
          model <- generate_startval(model)
        } else {
          # if it's not the first iteration, set starting values to estimates 
          # from previous iteration:
          model <- omxSetParameters(model,
                                    labels = names(coef(clustermodels_run[[i]])),
                                    values = coef(clustermodels_run[[i]]))
        }
        
        clustermodels[[i]] <- model
      }
      names(clustermodels) <- paste0("model_k", 1:n_clusters)
      
      clustermodels_run <- purrr::map(clustermodels, 
                                      mxRun, 
                                      silent = !verbose, suppressWarnings = TRUE)
      
      casewiseLL <- get_casewiseLL(clustermodels_run, n_clusters = n_clusters, n = n_persons)
      
      # compute complete-data log likelihood
      # complete_data_loglik <-  sum(post*log(class_proportions)-post*casewise.loglik)
      
      # compute observed-data log likelihood
      observed_data_LL <- compute_observed_data_LL(casewiseLL = casewiseLL, 
                                                   class_proportions = class_proportions)
      
      if(verbose){
        if(it != 1){
          print(paste0("Iteration: ", it, ". Log Likelihood: ", round(observed_data_LL, 4), ". Change: ", round(observed_data_LL - observed_data_LL0, 6), "."))
        } else {
          print(paste0("Iteration: ", it, ". Log Likelihood: ", round(observed_data_LL, 4), "."))
        }
      }
      
      ##!!! check for negative change##
      if(it > 2 & (observed_data_LL - observed_data_LL0) < 0){
        warning(paste0("Change in LL is negative. Part 1, random start: ", random_start))
      }
      
      # check close convergence and break loop if applicable:
      if(it > 2 & (observed_data_LL - observed_data_LL0) < 5){
        if(verbose){
          print(paste("Start", random_start, "close to convergence. Proceeding to next start."))
        }
        break
      }
      
      observed_data_LL0 = observed_data_LL
      
      if(it == maxit & verbose){
        print(paste("Start", random_start, " did not come close to convergence. Proceeding to next start."))
      }
      
    }
    
    # save observed_data_loglik, class_proportions, casewise.loglik, and parameter estimates for later:
    all_starts[[random_start]] <- list("observed_data_LL" = observed_data_LL,
                                       "class_proportions" = class_proportions,
                                       "casewiseLL" = casewiseLL,
                                       "clustermodels_run" = clustermodels_run)
  }
  
  if(verbose){
    print(paste("Setup finished. Finishing", n_best_starts, "best starts now."))
  }
  
  ## Part 2: complete the best starts
  # extract the best starts:
  best_starts <- purrr::map_dbl(all_starts, ~ .x$observed_data_LL) |> 
    order(decreasing = TRUE) |> 
    head(n_best_starts)
  
  # free the x0 and P0 parameters (for greater precision):
  personmodel_list <- purrr::map(personmodel_list,
                                 omxSetParameters,
                                 labels = c(paste0("ini_", factors),
                                            paste0("m_intercept_", factors),
                                            "P0_f1", "P0_f1f2", "P0_f2",
                                            "P0_icp1", "P0_icp12", "P0_icp2"),
                                 free = TRUE)
  
  best_loglik <- -Inf
  nonconvergences <- 0
  mxOption(key = "Major iterations", value = 1000)                              # reset the maximum number of iterations in OpenMx to 1000 (for greater precision)
  for(random_start in 1:n_best_starts){
    start_number <- best_starts[random_start]
    
    # extract outputs from corresponding start in Part 1:
    observed_data_LL0 <- all_starts[[start_number]]$observed_data_LL
    class_proportions <- all_starts[[start_number]]$class_proportions
    casewiseLL = all_starts[[start_number]]$casewiseLL
    clustermodels_run <- all_starts[[start_number]]$clustermodels_run
    
    
    # loop over iterations until convergence or max iterations are reached
    for(it in 1:maxit){
      post <- EStep(pi_ks = class_proportions, ngroup = n_persons, 
                    nclus = n_clusters, loglik = casewiseLL)
      
      # compute class proportions:
      class_proportions <- colMeans(post)
      
      ## M-step: fitting SSM model and update parameter estimates
      clustermodels <- vector(mode = "list", length = n_clusters)
      for (i in 1:n_clusters){
        clustername <- paste0("model_k", i)
        model <- create_model(clustername,
                              weights = post[, i],
                              objectives = objectives,
                              model_list = personmodel_list)
        model <- omxSetParameters(model,
                                  labels = names(coef(clustermodels_run[[i]])),
                                  values = coef(clustermodels_run[[i]]))
        clustermodels[[i]] <- model
      }
      names(clustermodels) <- paste0("model_k", 1:n_clusters)
      
      clustermodels_run <- purrr::map(clustermodels, 
                                      mxRun, 
                                      silent = !verbose, suppressWarnings = TRUE)      
      casewiseLL <- get_casewiseLL(clustermodels_run, n = n_persons, n_clusters = n_clusters)
      
      # compute complete-data log likelihood
      # complete_data_loglik <-  sum(post*log(class_proportions)-post*casewise.loglik)
      
      # compute observed-data log likelihood
      observed_data_LL <- compute_observed_data_LL(casewiseLL = casewiseLL, 
                                                   class_proportions = class_proportions)
      
      if(verbose){
        print(paste0("Iteration: ", it, ". Log Likelihood: ", round(observed_data_LL, 4), ". Change: ", round(observed_data_LL - observed_data_LL0, 6), "."))
      }
      
      ##!!! check for negative change##
      if((observed_data_LL - observed_data_LL0) < 0){
        warning(paste0("Change in LL is negative. Part 2, random start: ", random_start))
      }
      
      # check convergence and break loop if applicable:
      if((observed_data_LL - observed_data_LL0) < 1.0e-6){
        if(verbose){
          print("Convergence achieved.")
        }
        break
      }
      
      observed_data_LL0 = observed_data_LL
      
      # check if maximum number of iterations has been reached:
      if(it == maxit){
        if(verbose){
          print(paste("Max iterations reached without convergence. Start:", random_start))
        }
        nonconvergences <- nonconvergences + 1
      }
    }
    
    # check if the new fit is better than the best one from previous starts
    if(observed_data_LL > best_loglik){
      best_loglik <- observed_data_LL
      best_post <- post
      best_models <- clustermodels_run
    }
    
    if(verbose){
      print(paste("Start", random_start, "out of", n_best_starts, "best starts completed."))
    }
    
  }
  duration <- difftime(Sys.time(), estimation_start, unit = "s")
  
  #### find proxy maximum ####
  # create matrix with true cluster assignments
  post <- matrix(0, nrow = n_persons, ncol = n_clusters)
  for (i in 1:n_clusters) {
    post[true_clusters == i, i] <- 1
  }
  observed_data_LL0 <- -Inf
  mxOption(key = "Major iterations", value = 1000)                              # reset the maximum number of iterations in OpenMx to 1000
  
  for(it in 1:maxit){
    ## E-step: update class membership (skip in first iteration) and class proportions
    if(it >1){
      # update posteriors:
      post <- EStep(pi_ks = class_proportions, ngroup = n_persons,
                    nclus = n_clusters, loglik = casewiseLL)
    }
    
    # update class proportions:
    class_proportions <- colMeans(post)

    ## M-step: fitting SSM model and update parameter estimates
    clustermodels <- vector(mode = "list", length = n_clusters)
    for (i in 1:n_clusters){
      clustername <- paste0("model_k", i)
      model <- create_model(clustername,
                            weights = post[, i],
                            objectives = objectives,
                            model_list = personmodel_list)
      if(it == 1){
        model <- generate_startval(model)
      } else {
        model <- omxSetParameters(model,
                                  labels = names(coef(clustermodels_run[[i]])),
                                  values = coef(clustermodels_run[[i]]))
      }
      
      clustermodels[[i]] <- model
    }
    names(clustermodels) <- paste0("model_k", 1:n_clusters)
    
    clustermodels_run <- purrr::map(clustermodels, 
                                    mxRun, 
                                    silent = !verbose, suppressWarnings = TRUE)    
    casewiseLL <- get_casewiseLL(clustermodels_run, n_clusters = n_clusters, n = n_persons)
    
    # compute complete-data log likelihood
    # complete_data_loglik <-  sum(post*log(class_proportions)-post*casewise.loglik)
    
    # compute observed-data log likelihood
    observed_data_LL <- compute_observed_data_LL(casewiseLL = casewiseLL, 
                                                 class_proportions = class_proportions)
    
    if(verbose){
      if(it != 1){
        print(paste0("Iteration: ", it, ". Log Likelihood: ", round(observed_data_LL, 4), ". Change: ", round(observed_data_LL - observed_data_LL0, 6), "."))
      } else {
        print(paste0("Iteration: ", it, ". Log Likelihood: ", round(observed_data_LL, 4), "."))
      }
      
    }
    
    ##!!! check for negative change##
    if(it > 1 & (observed_data_LL - observed_data_LL0) < 0){
      warning("Change in LL is negative. Proxy maximum finding.")
    }
    
    # check convergence and break loop if applicable:
    if(it > 1 & (observed_data_LL - observed_data_LL0) < 1.0e-6){
      if(verbose){
        print("Convergence achieved.")
      }
      break
    }
    
    observed_data_LL0 = observed_data_LL
  }
  
  proxy_maximum <- observed_data_LL
  
  
  #### 6) extract estimates ####
  estimates <- lapply(best_models, coef)
  loglik <- best_loglik
  post <- as.data.frame(best_post)
  names(estimates) <- colnames(post) <- paste0("cluster", 1:n_clusters)
  modal_assignment <- t(apply(post, MARGIN = 1, function(x) ifelse(x == max(x), 1, 0)))
  class_proportions <- colMeans(post)
  
  
  
  clustering <- list("class_proportions" = class_proportions,
                     "posterior_prob" = post,
                     "modal_assignment" = modal_assignment)
  
  other <- list("loglik" = loglik,
                "nonconvergences" = nonconvergences,
                "proxy_maximum" = proxy_maximum,
                "duration" = duration)
  
  #### 7) build the output ####
  output <- list("data" = data,
                 "estimates" = estimates,
                 "clustering" = clustering,
                 "other" = other)
  
  
  return(output)
}