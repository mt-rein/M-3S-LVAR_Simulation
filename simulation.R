#### This script performs the simulation and saves the output ####

#### load packages and functions ####
library(conflicted)
library(dplyr)
library(flock)
library(lavaan)
library(MASS)
library(mcclust)
library(OpenMx)
library(EasyMx)
library(purrr)
library(tidyr)
library(truncnorm)
library(RcppAlgos)
library(stringr)

# packages for parallel processing:
library(parabar)
library(parallel)

# load functions
source("step1.R")
source("step2.R")
source("step3.R")
source("do_sim.R")
source("auxiliary_functions.R")

#### define condition grid
cond <- expand.grid(replication = 1:50,
                    n_obs = c(56, 112),
                    n_clusters = c(2, 4),
                    n_factors = c(2, 4))

# add seeds:
set.seed(123)
cond$seed <- sample(1:(nrow(cond)*5), size = nrow(cond), replace = FALSE)
# add iteration number:
cond$iteration <- 1:nrow(cond)                                                  # (unique) number of each iteration

# # split off rows where n_clusters and n_factors are 4
# cond_subset <- cond %>%
#   dplyr::filter(n_clusters == 4, n_factors == 4)
# # remove those rows from original
# cond <- cond %>%
#   dplyr::filter(!(n_clusters == 4 & n_factors == 4))

#### set up parallel computing ####
## open cluster
numCores <- parallel::detectCores() - 1
backend <- start_backend(numCores, cluster_type = "psock", backend_type = "async")


## load libraries in cluster
parabar::evaluate(backend, {
  library(conflicted)
  library(dplyr)
  library(flock)
  library(lavaan)
  library(MASS)
  library(mcclust)
  library(OpenMx)
  library(EasyMx)
  library(purrr)
  library(tidyr)
  library(truncnorm)
  library(RcppAlgos)
  library(stringr)

  source("step1.R")
  source("step2.R")
  source("step3.R")
  source("do_sim.R")
  source("auxiliary_functions.R")
})

## load objects in cluster
export(backend, "cond")

#### perform simulation ####
start  <- Sys.time()
output <- par_lapply(backend, 1:nrow(cond), do_sim, cond = cond, outputfile = "output_sim.csv", verbose = FALSE)
end <- Sys.time()

# close cluster:
stop_backend(backend)