step2 <- function(step1output, id) {
  #### Preparations ####
  fit_step1 <- step1output$MMoutput
  data <- step1output$data
  measurementmodel <- step1output$measurementmodel
  n_mb <- length(measurementmodel)
  indicators <- lavaan::lavNames(lavaan::lavaanify(measurementmodel), "ov")
  factors <- lavaan::lavNames(lavaan::lavaanify(measurementmodel), "lv")
  n_factors <- length(factors)
  unique_ids <- unique(data[, id])
  n_persons <- length(unique_ids)

  #### compute factor scores ####
  for (m in 1:n_mb) {
    # if there are measurement blocks, compute factor scores in each block
    temp <- lavaan::lavPredict(fit_step1[[m]],
                               assemble = TRUE,
                               append.data = TRUE) |>
      as.data.frame()
    # and append to original data
    data <- dplyr::full_join(data, temp,
                             by = dplyr::intersect(colnames(data),
                                                   colnames(temp)))
  }

  #### compute lambda_star and theta_star ####


  # create lists to store MM parameter values per block
  # (lists with n_mb elements)
  psi_block <-
    lambda_block <-
    theta_block <-
    vector(mode = "list", length = n_mb)
  for (m in 1:n_mb) {
    EST_block <- lavaan::lavInspect(fit_step1[[m]], "est")
    psi_block[[m]] <- EST_block[["psi"]]
    lambda_block[[m]] <- EST_block[["lambda"]]
    theta_block[[m]] <- EST_block[["theta"]]
  }

  # combine the matrices of different measurement blocks into a single matrix
  psi <- lavaan::lav_matrix_bdiag(psi_block)
  lambda <- lavaan::lav_matrix_bdiag(lambda_block)
  theta  <- lavaan::lav_matrix_bdiag(theta_block)

  # name the matrices' rows and columns
  rownames(psi) <- colnames(psi) <- factors
  rownames(lambda) <- indicators
  colnames(lambda) <- factors
  rownames(theta) <- colnames(theta) <- indicators

  # compute lambda_star and theta_star
  sigma <- lambda %*% psi %*% t(lambda) + theta
  A <- psi %*% t(lambda) %*% solve(sigma)
  lambda_star <- matrix(diag(A %*% lambda),
                        nrow = n_persons,
                        ncol = n_factors,
                        byrow = TRUE)
  theta_star <- matrix(diag(A %*% theta %*% t(A)),
                       nrow = n_persons,
                       ncol = n_factors,
                       byrow = TRUE)
  colnames(lambda_star) <- colnames(theta_star) <- factors
  rownames(lambda_star) <- rownames(theta_star) <- unique_ids

  # assemble output
  output <- list("data" = data,
                 "lambda_star" = lambda_star,
                 "theta_star" = theta_star,
                 "other" = list("factors" = factors,
                                "indicators" =  indicators,
                                "id" = id))

  return(output)
}