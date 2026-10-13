#### This script performs the simulation and saves the output ####

#### load packages and functions ####
library(conflicted)
library(dplyr)
library(ezLVAR)
library(flock)
library(lavaan)
library(MASS)
library(mcclust)
library(OpenMx)
library(purrr)
library(tidyr)
library(truncnorm)
library(RcppAlgos)
library(stringr)

# packages for parallel processing:
library(parabar)
library(parallel)

# load functions
source("do_sim.R")
source("auxiliary_functions.R")

#### define condition grid
cond <- expand.grid(
  replication = 1:50,
  n_obs = c(56, 112),
  n_clusters = c(2, 4),
  n_factors = c(2, 4)
)

# add seeds:
set.seed(123)
cond$seed <- sample(1:(nrow(cond) * 5), size = nrow(cond), replace = FALSE)
# add iteration number:
cond$iteration <- 1:nrow(cond) # (unique) number of each iteration
cond <- cond %>% relocate(iteration) # move iteration column to beginning

#### set up parallel computing ####
## open cluster
numCores <- parallel::detectCores() - 1
backend <- start_backend(
  numCores,
  cluster_type = "psock",
  backend_type = "async"
)


## load libraries in cluster
parabar::evaluate(backend, {
  library(conflicted)
  library(dplyr)
  library(ezLVAR)
  library(flock)
  library(lavaan)
  library(MASS)
  library(mcclust)
  library(OpenMx)
  library(purrr)
  library(tidyr)
  library(truncnorm)
  library(RcppAlgos)
  library(stringr)

  source("do_sim.R")
  source("auxiliary_functions.R")
})

## load objects in cluster
export(backend, "cond")

#### prepare error file ####
error_template <- cond[0, ]
error_template$error_message <- character()
write.table(
  error_template,
  file = "error_sim.csv",
  sep = ",",
  row.names = FALSE,
  col.names = TRUE
)

#### perform simulation ####
start <- Sys.time()
output <- par_lapply(
  backend,
  1:nrow(cond),
  run_one_sim,
  cond = cond,
  outputfile = "output_sim.csv",
  errorfile = "error_sim.csv",
  verbose = FALSE
)
end <- Sys.time()

# close cluster:
stop_backend(backend)
