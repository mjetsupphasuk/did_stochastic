
# set directories
setwd("/work/users/m/j/mjets/dissertation/paper2")
.libPaths("/nas/longleaf/home/mjets/RLibs")

# load libraries
library(here)
library(tidyverse)

# set ggplot theme
theme_set(theme_bw(base_size = 14))

# import functions
source(here('code', 'utils.R'))

# data generating scenario
scenario = 1
source(here('code', 'sims_scenarios_mtp.R'))

# create very large dataset
set.seed(993811 * scenario)
n_big = 100000
panel_data <- gen_staggered_data(
  n = n_big, 
  shape1_alpha_f = shape1_alpha_f,
  shape2_beta_f  = shape2_beta_f,
  T_periods = T_periods,
  c_fun = c_fun,
  theta_t = theta_t,
  cohort_predictor = cohort_predictor
)

# delta values want
delta = seq(-100, 100, length.out = 201)


# Compute disaggregated asdt (for each g, t pair) -------------------------

cohort_use <- 1:T_periods
post_times <- 1:T_periods

asdt_true_list <- list()
counter <- 1

for (g in cohort_use) {
  for (t in post_times) {
    if (t < g) next # Only calculate effects for post-adoption periods
    
    # Calculate the true Average Shifted Dose Effect (ASDT) utilizing the 
    # Oracle true outcome and propensity functions
    asdt_vec <- compute_true_stag_asdt(
      panel_data = panel_data,
      delta      = delta,
      g_spec     = g,
      t_spec     = t,
      shape1_alpha_f = shape1_alpha_f,
      shape2_beta_f  = shape2_beta_f,
      c_fun      = c_fun,
      mu_0_gt_f  = mu_0_gt_f,
      mu_d_gt_f  = mu_d_gt_f,
      gamma_t    = gamma_t,
      theta_t    = theta_t,
      T_periods  = T_periods
    )
    
    asdt_true_list[[counter]] <- data.frame(
      g     = g,
      t     = t,
      e     = t - g,
      delta = delta,
      asdt  = as.numeric(asdt_vec)
    )
    counter <- counter + 1
    
    print(sprintf("Completed t %i; Time: %s", t, Sys.time()))
  }
  
  print(sprintf("Completed g %i; Time: %s", g, Sys.time()))
}

asdt_true <- bind_rows(asdt_true_list)

# save true asdt
truth_file_path <- here("results", "sims", paste0("asdt_mtp_true_sc", scenario, ".rds"))
saveRDS(asdt_true, truth_file_path)



# Compute aggregated ASDT parameters --------------------------------------

# Calculate the baseline population probabilities for each cohort: P(G=g | G <= T_periods)
# These act as the fundamental weights for all subsequent aggregations
pop_units <- panel_data %>% distinct(id, .keep_all = TRUE)
prob_g <- pop_units %>%
  filter(G <= T_periods) %>% 
  group_by(G) %>%
  summarise(p_g = n(), .groups = "drop") %>%
  mutate(p_g = p_g / sum(p_g)) %>%  
  rename(g = G)

# Merge cohort probabilities onto the disaggregated results
asdt_true_agg_base <- asdt_true %>%
  left_join(prob_g, by = "g") %>%
  filter(g <= t) 

#Aggregation Type 1: Simple (Overall)
# A straightforward weighted average of all post-treatment effects across all g and t
true_simple <- asdt_true_agg_base %>%
  group_by(delta) %>%
  summarise(
    true_estimate = sum(asdt * p_g) / sum(p_g), 
    .groups = "drop"
  ) %>%
  mutate(type = "simple", target = "Overall")

# --- Aggregation Type 2: Group ---
# Averages across time for each cohort first, then computes an overall average weighted by cohort size
true_cohort <- asdt_true_agg_base %>%
  group_by(delta, g) %>%
  summarise(true_estimate = mean(asdt), .groups = "drop") %>% 
  mutate(type = "group", target = paste0("g=", g))

true_group_overall <- true_cohort %>%
  left_join(prob_g, by = "g") %>%
  group_by(delta) %>%
  summarise(true_estimate = sum(true_estimate * p_g) / sum(p_g), .groups = "drop") %>% 
  mutate(type = "group", target = "Overall")

#Aggregation Type 3: Dynamic (Unbalanced Event Study)
# Averages effects across cohorts for each specific relative time 'e' (event time)
true_dynamic <- asdt_true_agg_base %>%
  group_by(delta, e) %>%
  summarise(
    true_estimate = sum(asdt * p_g) / sum(p_g),
    .groups = "drop"
  ) %>%
  mutate(type = "dynamic", target = paste0("e=", e))

true_dynamic_overall <- true_dynamic %>%
  filter(e >= 0) %>%
  group_by(delta) %>%
  summarise(true_estimate = mean(true_estimate), .groups = "drop") %>% 
  mutate(type = "dynamic", target = "Overall")

# Aggregation Type 4: Balanced Dynamic (Balanced Event Study)
# Restricts the dynamic aggregation only to cohorts observed for at least 'balance_e' periods
# This prevents compositional changes in the cohorts driving the dynamic effects over time
balance_e <- 1

balanced_groups <- prob_g %>%
  filter(g + balance_e <= T_periods)

# Re-weight probabilities among the restricted subset of balanced cohorts
balanced_prob_g <- balanced_groups %>%
  mutate(p_g_balanced = p_g / sum(p_g)) %>% 
  select(g, p_g_balanced)

true_dynamic_balanced <- asdt_true_agg_base %>%
  filter(g %in% balanced_groups$g) %>% 
  filter(e <= balance_e) %>%           
  left_join(balanced_prob_g, by = "g") %>%
  group_by(delta, e) %>%
  summarise(
    true_estimate = sum(asdt * p_g_balanced), 
    .groups = "drop"
  ) %>%
  mutate(type = "dynamic_balanced", target = paste0("e=", e))

true_dynamic_balanced_overall <- true_dynamic_balanced %>%
  filter(e >= 0) %>%
  group_by(delta) %>%
  summarise(true_estimate = mean(true_estimate), .groups = "drop") %>%
  mutate(type = "dynamic_balanced", target = "Overall")

# Assembly and Save
truth_aggregated <- bind_rows(
  true_simple,
  true_cohort, true_group_overall,
  true_dynamic, true_dynamic_overall,
  true_dynamic_balanced, true_dynamic_balanced_overall
)

agg_truth_path <- here("results", "sims", paste0("asdt_mtp-agg_true_sc", scenario, ".rds"))
saveRDS(truth_aggregated, agg_truth_path)

