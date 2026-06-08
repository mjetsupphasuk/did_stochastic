
# set directories
setwd("/work/users/m/j/mjets/dissertation/paper2")
.libPaths("/nas/longleaf/home/mjets/RLibs")

# load libraries
library(BART)
library(SuperLearner)
library(hal9001)
library(caret)
library(here)
library(tidyverse)

# set ggplot theme
theme_set(theme_bw())

# import functions
source(here('code', 'utils.R'))

# data generating scenario
scenario = 2
source(here('code', 'sims_scenarios.R'))

# collect arguments passed in from job array
args = commandArgs(trailingOnly=TRUE)
sim = as.numeric(args[1])

# read data
data_file = here('data', 'sims', paste0('sim_data_sc_', scenario, '.rds'))
data = readRDS(data_file)[[sim]]

# delta values want
delta_vals = seq(-100, 100, length.out = 201)

# cross-fitting
cf_folds = 5

# # nuisance function estimator
# nuisance_estimator = 'bart'

# place to save results
save_string = here('results', 'sims', 'interim', 
                   paste0('sim_res_sc', scenario, 
                          '_cfK', cf_folds, 
                          '_sim', sim,
                          '.rds'))


# Perform estimation ------------------------------------------------------

# estimates
estimates.glm = estimate_psi(data$deltaY, data$treat,
                             data %>% select(starts_with('X')),
                             stochastic_policy = 'exp_tilt',
                             policy_params = list('delta'=delta_vals), 
                             K = cf_folds,
                             nuisance_estimator = 'glm')

estimates.bart = estimate_psi(data$deltaY, data$treat,
                              data %>% select(starts_with('X')),
                              stochastic_policy = 'exp_tilt',
                              policy_params = list('delta'=delta_vals), 
                              K = cf_folds,
                              nuisance_estimator = 'bart')

estimates.oracle = estimate_psi(data$deltaY, data$treat,
                                data %>% select(starts_with('X')),
                                stochastic_policy = 'exp_tilt',
                                policy_params = list('delta'=delta_vals), 
                                K = cf_folds,
                                nuisance_estimator = 'oracle')

# remove IF matrix (to save memory)
estimates.glm$IF_matrix = NULL
estimates.bart$IF_matrix = NULL
estimates.oracle$IF_matrix = NULL

# save results
saveRDS(list('glm' = estimates.glm,
             'bart' = estimates.bart,
             'oracle' = estimates.oracle),
        file = save_string)

