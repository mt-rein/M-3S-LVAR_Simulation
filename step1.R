step1 <- function(data, measurementmodel, id) {
  # How many measurement blocks do we have?
  n_mb <- length(measurementmodel)

  # create list, with elements equal to number of measurement blocks:
  MMoutput <- vector("list", length = n_mb)

  # estimate MM in each block:
  for (m in 1:n_mb) {
    MMoutput[[m]] <- lavaan::cfa(measurementmodel[[m]],
                                 data = data)
  }

  output <- list("MMoutput" = MMoutput,
                 "data" = data,
                 "measurementmodel" = measurementmodel)

  return(output)
}