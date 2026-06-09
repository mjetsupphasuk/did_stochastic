
# load libraries
library(here)
library(tidyverse)

# set ggplot theme
theme_set(theme_bw(base_size = 14))

# import functions
source(here('code', 'utils.R'))

# data generating scenario
scenario = 2
source(here('code', 'sims_scenarios.R'))


# Compute the true quantities ---------------------------------------------

# create very large dataset
set.seed(417816 * scenario)
n_big = 100000
data = gen_data(n_big, pi_bin, mu.0_f, shape1_alpha_f, shape2_beta_f, mu.d_f)

# delta values want
delta = seq(-100, 100, length.out = 201)

# compute asdt
asdt.true = data.frame(delta = delta,
                       asdt = compute_true_asdt(data, delta, shape1_alpha_f, shape2_beta_f))

# plot it
ggplot(asdt.true, aes(x=delta, y=asdt)) +
  geom_point() +
  geom_line()

# save true asdt
save_path = here('results', 'sims', paste0('asdt_true_sc', scenario, '.rds'))
saveRDS(asdt.true, save_path)