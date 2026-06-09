
library(haven)
library(readr)
library(BART)
library(hal9001)
library(SuperLearner)
library(caret)
library(tidycensus)
library(here)
library(tidyverse)

theme_set(theme_bw(base_size = 16))

# import functions
source(here('code', 'utils.R'))

# all cases
# note that treatment initiation is always defined as ref_year + 1
# ignore cases with which_treat='raw'; prospectivity scores cannot
# be compared in absolute terms across shales
cases = expand.grid(ref_year = 2007:2008,
                    event_time = c(-4:-2, 0:4),
                    y = c('emp_tot', 'r_income_tot'),
                    which_treat = c('raw', 'within_shale'))


# Evaluate ----------------------------------------------------------------

# number of cross-fitting folds
cf_folds = 10

# nuisance function estimator
nuisance_estimator = 'sl'

# delta to examine
delta_vals = seq(from=-100, to=100, length.out=201)

# plot all cases
results.all = list()
for (case in 1:nrow(cases)) {

  # extract case values
  ref_year = cases$ref_year[case]
  event_time = cases$event_time[case]
  y_name = as.character(cases$y[case])
  which_treat = cases$which_treat[case]

  # read results
  estimates = readRDS(here('results', 'fracking',
                             paste0('treat', gsub('_shale', '', which_treat), '_',
                             nuisance_estimator, '_',
                             'cfK', cf_folds, '_',
                             gsub('_', '', y_name), '_', 
                             'refyear', ref_year,
                             '_lag', gsub('-', 'neg', event_time),
                             '.rds')))

  # format results
  results = data.frame(delta = delta_vals,
                       psihat = estimates$psihat,
                       ci_lower = estimates$psihat - qnorm(0.975)*sqrt(estimates$variance),
                       ci_upper = estimates$psihat + qnorm(0.975)*sqrt(estimates$variance),
                       ref_year = ref_year,
                       event_time = event_time,
                       y_name = y_name,
                       which_treat = which_treat)
  
  # store results in a list
  results.all[[case]] = results

}


# format combined results
results.all.df = bind_rows(results.all)
results.all.df = results.all.df %>%
  mutate(calendar_time = ref_year + 1 + event_time,
         exp_year = ref_year + 1,
         delta_char = ifelse(delta==-10, paste0('δ = ', -10), as.character(delta)),
         delta_char = factor(delta_char, levels = c(paste0('δ = ', -10), sort(unique(delta))[-1])),
         exp_year_char = ifelse(exp_year==min(exp_year), paste0('Cohort: ', min(exp_year)), as.character(exp_year)),
         exp_year_char = factor(exp_year_char, levels = c(paste0('Cohort: ', min(exp_year)), sort(unique(exp_year))[-1])))

# subset time period
results.all.df = results.all.df %>% filter(calendar_time <= 2012, event_time <= 3)

# plot function
myplot <- function(data, ylims=NULL) {
  
  gg = ggplot(data %>% filter(delta %in% c(-10, -5, 0, 5, 10)), 
              aes(x = calendar_time, y = psihat, ymin = ci_lower, ymax = ci_upper)) +
    geom_point() +
    geom_errorbar(width = 0.25) +
    geom_hline(yintercept=0, linetype='solid') +
    geom_vline(aes(xintercept=ref_year), linetype='dashed') +
    facet_grid(rows = vars(exp_year_char), cols = vars(delta_char)) +
    scale_x_continuous(breaks = seq(2005, 2011, by=3),
                       limits = c(2003.5, 2012.5)) +
    scale_y_continuous(limits = ylims) +
    labs(x = 'Time', y = 'Average stochastic dose effect among the treated') +
    theme(panel.grid.minor = element_blank())
  print(gg)
}


# plot for emp_tot, within shale treatment
results.tmp = results.all.df %>% 
  filter(y_name == 'emp_tot',
         which_treat == 'within_shale')

ylims = c(-0.07, 0.13)

myplot(results.tmp, ylims=ylims)
ggsave(here('results', 'fracking', 'main',
              paste0('treatwithin_',
              nuisance_estimator, '_',
              'cfK', cf_folds, '_',
              'emptot',
              '.png')),
       width=7, height=4.5, units='in', scale=1.4)



# plot for income, within shale treatment
results.tmp = results.all.df %>% 
  filter(y_name == 'r_income_tot',
         which_treat == 'within_shale')

myplot(results.tmp)
ggsave(here('results', 'fracking', 'main',
              paste0('treatwithin_',
              nuisance_estimator, '_',
              'cfK', cf_folds, '_',
              'rincometot',
              '.png')),
       width=7, height=4.5, units='in', scale=1.4)



# plot to highlight dose response
results.tmp = results.all.df %>% 
  filter(y_name == 'emp_tot',
         which_treat == 'within_shale',
         exp_year == 2008,
         event_time == 3)

ggplot(results.tmp %>% filter(abs(delta) <= 10),
       aes(x = delta, y = psihat, ymin = ci_lower, ymax = ci_upper)) +
  geom_line() +
  geom_ribbon(alpha = 0, color = 'black', linetype = 'dashed') +
  labs(x = 'δ', y = 'Average stochastic dose effect among the treated')
ggsave(here('results', 'fracking', 'main',
              paste0('treatwithin_',
              nuisance_estimator, '_',
              'cfK', cf_folds, '_',
              'emptot_',
              'g2008_',
              'e3',
              '.png')),
       width=7, height=4.5, units='in', scale=1.4)


# get exact values to put in manuscript
results.all.df %>% filter(delta==10, exp_year==2008, event_time==3, y_name=='emp_tot', which_treat=='within_shale')



# Gaussian concentration analysis -----------------------------------------

cases = cases %>% filter(which_treat == 'within_shale')

# number of cross-fitting folds
cf_folds = 10

# nuisance function estimator
nuisance_estimator = 'sl'

# concentration points to examine
delta_vals = seq(from=0.01, to=0.99, length.out=100)

# plot all cases
results.all = list()
for (case in 1:nrow(cases)) {
  
  # extract case values
  ref_year = cases$ref_year[case]
  event_time = cases$event_time[case]
  y_name = as.character(cases$y[case])
  which_treat = cases$which_treat[case]
  
  # read results
  estimates = readRDS(here('results', 'fracking',
                             paste0('treat', gsub('_shale', '', which_treat), '_',
                             nuisance_estimator, '_',
                             'cfK', cf_folds, '_',
                             gsub('_', '', y_name), '_', 
                             'refyear', ref_year,
                             '_lag', gsub('-', 'neg', event_time),
                             '_gauss.rds')))
  
  # format results
  results = data.frame(delta = delta_vals,
                       psihat = estimates$psihat,
                       ci_lower = estimates$psihat - qnorm(0.975)*sqrt(estimates$variance),
                       ci_upper = estimates$psihat + qnorm(0.975)*sqrt(estimates$variance),
                       ref_year = ref_year,
                       event_time = event_time,
                       y_name = y_name,
                       which_treat = which_treat)
  
  # store results in a list
  results.all[[case]] = results
  
}

# format combined results
results.all.df = bind_rows(results.all)
results.all.df = results.all.df %>%
  mutate(calendar_time = ref_year + 1 + event_time,
         exp_year = ref_year + 1,
         exp_year_char = ifelse(exp_year==min(exp_year), paste0('Cohort: ', min(exp_year)), as.character(exp_year)),
         exp_year_char = factor(exp_year_char, levels = c(paste0('Cohort: ', min(exp_year)), sort(unique(exp_year))[-1])))

# subset time period
results.all.df = results.all.df %>% filter(calendar_time <= 2012, event_time <= 3)

data.tmp = results.all.df %>%
  filter(y_name == 'emp_tot',
         which_treat == 'within_shale') %>%
  mutate(calendar_time = ref_year + 1 + event_time,
         exp_year = ref_year + 1)
ggplot(data.tmp %>% filter(delta > 0.5), 
       aes(x = delta, y = psihat, ymin = ci_lower, ymax = ci_upper)) +
  geom_line() +
  geom_ribbon(alpha = 0.1) +
  geom_hline(yintercept=0, linetype='solid') +
  facet_grid(rows = vars(exp_year), cols = vars(event_time)) +
  labs(x = "Concentration point d'", y = 'Average stochastic dose effect among the treated') +
  theme(panel.grid.minor = element_blank())


delta_vals_plot = delta_vals[c(50, 75, 100)]
data.tmp2 = data.tmp %>%
  filter(delta %in% delta_vals_plot) %>%
  mutate(delta_plot = paste0("d' = ", round(delta, 2)))

ggplot(data.tmp2, 
       aes(x = calendar_time, y = psihat, ymin = ci_lower, ymax = ci_upper)) +
  geom_point() +
  geom_errorbar(width = 0.25) +
  geom_hline(yintercept=0, linetype='solid') +
  geom_vline(aes(xintercept=ref_year), linetype='dashed') +
  facet_grid(rows = vars(exp_year_char), cols = vars(delta_plot)) +
  labs(x = "Time", y = 'Average stochastic dose effect among the treated') +
  theme(panel.grid.minor = element_blank())
ggsave(here('results', 'fracking', 'main', 
              paste0('treatwithin_',
              nuisance_estimator, '_',
              'cfK', cf_folds, '_',
              'emptot', 
              '_gauss.png')),
       width=7, height=4.5, units='in', scale=1.4)


ggplot(data.tmp %>% filter(delta > 0.5, exp_year==2008, event_time==3), 
       aes(x = delta, y = psihat, ymin = ci_lower, ymax = ci_upper)) +
  geom_line() +
  geom_ribbon(alpha = 0.1) +
  geom_hline(yintercept=0, linetype='dashed') +
  labs(x = "Concentration point d'", y = 'Average stochastic dose effect among the treated') +
  theme(panel.grid.minor = element_blank())
