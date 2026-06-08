
# set directories
setwd("/work/users/m/j/mjets/dissertation/paper2")
.libPaths("/nas/longleaf/home/mjets/RLibs")

# load libraries
library(here)
library(tidyverse)

# set ggplot theme
theme_set(theme_bw(base_size = 18))

# import functions
source(here('code', 'utils.R'))

# data generating scenario
scenario = 2
source(here('code', 'sims_scenarios.R'))

# delta values want
delta_vals = seq(-100, 100, length.out = 201)

# cross-fitting
cf_folds = 5

# # nuisance function estimator
# nuisance_estimator = 'bart'


# Evaluate sims -----------------------------------------------------------

num_sims = 1000
results_path = here('results', 'sims', 'interim')

# read in and aggregate results
results.all = list()
for (sim in 1:num_sims) {
  read_string = file.path(results_path, paste0('sim_res_sc', scenario,
                                               '_cfK', cf_folds,
                                               '_sim', sim, '.rds'))

  results.tmp = readRDS(read_string)
  results.all[[sim]] = format_results(results.tmp, sim)
  
}
results.all = bind_rows(results.all)


# get truth
asdt.true = readRDS(here('results', 'sims', paste0('asdt_true_sc', scenario, '.rds')))
results.all = left_join(results.all, asdt.true, by='delta')

# average psi hat
results.all = results.all %>%
  group_by(nuisance, delta) %>%
  mutate(psihat_avg = mean(psihat)) %>%
  ungroup()

# nicer nuisance names
results.all = results.all %>%
  mutate(nuisance2 = case_when(nuisance=='bart' ~ 'BART',
                               nuisance=='glm' ~ 'GLMs',
                               nuisance=='oracle' ~ 'Oracle'))

# plot
ggplot(results.all %>% filter(abs(delta) <= 10)) +
  geom_line(mapping = aes(x=delta, y=psihat, group=sim), color='grey', alpha=0.1) +
  geom_line(mapping = aes(x=delta, y=asdt), inherit.aes=FALSE, color='red') +
  geom_line(mapping = aes(x=delta, y=psihat_avg), inherit.aes=FALSE, color='blue') +
  facet_grid(cols = vars(nuisance2)) +
  labs(x = 'δ', y = 'ASDT') +
  theme(panel.grid = element_blank())
ggsave(here('results', 'sims', paste0('result_plot_sc', scenario, '_cfK', cf_folds, '.png')),
       width=8, height=4.5, units='in')


# compute coverage
results.all = results.all %>%
  mutate(cover = (asdt >= ci_lower) & (asdt <= ci_upper))
coverage = results.all %>%
  group_by(nuisance2, delta) %>%
  summarize(coverage = mean(cover))

# plot coverage
ggplot(coverage %>% filter(abs(delta) <= 10), 
       aes(x=delta, y=coverage)) +
  geom_line() +
  geom_abline(intercept=0.95, slope=0, color='red', linetype='dashed') +
  facet_grid(cols = vars(nuisance2)) +
  ylim(c(0,1)) +
  labs(x = 'δ', y = 'Coverage') +
  theme(panel.grid = element_blank())
ggsave(here('results', 'sims', paste0('result_coverage_sc', scenario, '_cfK', cf_folds, '.png')),
       width=8, height=4.5, units='in')

# summary stats
sim_stats = results.all %>%
  group_by(nuisance, delta) %>%
  summarize(bias = mean(psihat - asdt),
            mse = mean((psihat - asdt)^2),
            emp_se = sd(psihat),
            avg_se = mean(sqrt(psihat_var)),
            coverage = mean(cover)) %>%
  mutate_all(round, 3)
sim_stats %>% filter(delta %in% c(-10,-5,0,5,10))
write_csv(sim_stats %>% filter(delta %in% c(-10,-5,0,5,10)),
          here('results', 'sims', paste0('result_stats_sc', scenario, '_cfK', cf_folds, '.csv')))

# check empirical vs average variance for psi1 and psi2 separately
results.all %>% 
  filter(abs(delta) <= 10) %>%
  group_by(nuisance, delta) %>% 
  summarize(psi1_emp_se = sd(psi1_hat),
            psi1_avg_se = mean(sqrt(psi1_hat_var)),
            psi2_emp_se = sd(psi2_hat),
            psi2_avg_se = mean(sqrt(psi2_hat_var)),
            psi_emp_se = sd(psihat),
            psi_avg_se = mean(sqrt(psihat_var))) %>%
  mutate_all(round, 3) %>%
  as.data.frame() 

# check psi2
results.all %>%
  filter(delta==0) %>%
  group_by(nuisance) %>%
  summarize(mean(psi2_hat))
