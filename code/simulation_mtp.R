
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
source(here('code', 'sims_scenarios_mtp.R'))

# collect arguments passed in from job array
args = commandArgs(trailingOnly=TRUE)
sim = as.numeric(args[1])

# read data
data_file = here('data', 'sims', paste0('sim_data_mtp_sc_', scenario, '.rds'))
panel_data = readRDS(data_file)[[sim]]

# delta values want
delta_vals = seq(-100, 100, length.out = 201)

# cross-fitting
cf_folds = 5

# place to save results
save_string_disagg = here('results', 'sims', 'interim', 
                          paste0('sim_res_mtp_sc', scenario, 
                                 '_cfK', cf_folds, 
                                 '_sim', sim,
                                 '.rds'))

save_string_agg = here('results', 'sims', 'interim', 
                       paste0('sim_res_mtp-agg_sc', scenario, 
                              '_cfK', cf_folds, 
                              '_sim', sim,
                              '.rds'))



# Perform estimation ------------------------------------------------------

# initiate lists/counters to store results
results.all <- list()
IF_list <- list()
counter <- 1

# cohorts and calendar times
cohort_use <- 1:T_periods
post_times <- 1:T_periods

for (g in cohort_use) {
  for (t in post_times) {
    if (t < g) next
    
    cat("Processing sim=", sim, " g=", g, " t=", t, "\n")
    
    # create a "wide" cross-sectional dataset out of the panel data
    cs <- make_gt_cross_section(panel_data, g, t)
    
    # main estimation
    
    # glm nuisance functions
    estimates.glm <- tryCatch(
      estimate_psi_mtp(
        Y_t = cs$Y_t,
        Y_g1 = cs$Y_g1,
        treat = cs$treat,
        G = cs$G,
        x = cs$x,
        stochastic_policy = 'exp_tilt',
        policy_params = list('delta' = delta_vals),
        nuisance_estimator = 'glm',
        g = g,
        t = t,
        T_periods = T_periods,
        K = cf_folds,
        design_d = seq(0.001, 0.999, length.out = 500)
      ),
      error = function(e) {
        cat("  -> ERROR at g=", g, " t=", t, ": ", e$message, "\n")
        return(NULL)
      }
    )
    
    # bart nuisance functions
    estimates.bart <- tryCatch(
      estimate_psi_mtp(
        Y_t = cs$Y_t,
        Y_g1 = cs$Y_g1,
        treat = cs$treat,
        G = cs$G,
        x = cs$x,
        stochastic_policy = 'exp_tilt',
        policy_params = list('delta' = delta_vals),
        nuisance_estimator = 'bart',
        g = g,
        t = t,
        T_periods = T_periods,
        K = cf_folds,
        design_d = seq(0.001, 0.999, length.out = 500)
      ),
      error = function(e) {
        cat("  -> ERROR at g=", g, " t=", t, ": ", e$message, "\n")
        return(NULL)
      }
    )
    
    # oracle nuisance functions
    estimates.oracle <- tryCatch(
      estimate_psi_mtp(
        Y_t = cs$Y_t,
        Y_g1 = cs$Y_g1,
        treat = cs$treat,
        G = cs$G,
        x = cs$x,
        stochastic_policy = 'exp_tilt',
        policy_params = list('delta' = delta_vals),
        nuisance_estimator = 'oracle',
        g = g,
        t = t,
        T_periods = T_periods,
        K = cf_folds,
        design_d = seq(0.001, 0.999, length.out = 500)
      ),
      error = function(e) {
        cat("  -> ERROR at g=", g, " t=", t, ": ", e$message, "\n")
        return(NULL)
      }
    )
    
    # save the results
    results.all[[paste(g, t, sep="_")]] <- list('glm' = estimates.glm$main,
                                                'bart' = estimates.bart$main,
                                                'oracle' = estimates.oracle$main)
    
    # save the influence matrix separately
    IF_list[[paste(g, t, sep="_")]] <- list('glm' = estimates.glm$IF_matrix,
                                            'bart' = estimates.bart$IF_matrix,
                                            'oracle' = estimates.oracle$IF_matrix)
    
    counter <- counter + 1
  }
}

# save results for disaggregated parameters
saveRDS(results.all, save_string_disagg)



# estimate aggregated parameters
results.df = format_results_mtp(results.all, sim)
results.all.agg = list()
for (nuisance in c('glm', 'bart', 'oracle')) {
  
  results.all.agg[[nuisance]] <- aggregate_continuous_did(
    results_df = results.df, 
    IF_list = IF_list, 
    G = panel_data %>% filter(periods==0) %>% pull(G), 
    T_periods = T_periods, 
    type = "all",
    nuisance_estimator = nuisance
  )
  
}

# save results for aggregated parameters
saveRDS(results.all.agg, save_string_agg)
