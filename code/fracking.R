
library(haven)
library(readr)
library(BART)
library(hal9001)
library(SuperLearner)
library(caret)
library(tidycensus)
library(tigris)
library(viridis)
library(here)
library(tidyverse)

theme_set(theme_bw(base_size = 14))

# import functions
source(here('code', 'utils.R'))

# collect arguments passed in from job array
args = commandArgs(trailingOnly=TRUE)
case = as.numeric(args[1])

# all cases
# note that treatment initiation is always defined as ref_year + 1
# ignore cases with which_treat='raw'; prospectivity scores cannot
# be compared in absolute terms across shales
cases = expand.grid(ref_year = 2007:2008,
                    event_time = c(-4:-2, 0:4),
                    y = c('emp_tot', 'r_income_tot'),
                    which_treat = c('raw', 'within_shale'))

# extract case values
ref_year = cases$ref_year[case]
exp_year = ref_year + 1
event_time = cases$event_time[case]
y_name = as.character(cases$y[case])
which_treat = cases$which_treat[case]

set.seed(892385 + case)


# Data --------------------------------------------------------------------

# load cleaned data
data_analysis = readRDS(here('data', 'fracking', 'data_analysis.rds'))

# load demo data
demo = readRDS(here('data', 'fracking', 'demo.rds'))

# subset units to treated and not-yet-treated
# for negative event times (i.e., pre-tests), define treated and untreated groups at exp_year
data_analysis = data_analysis %>%
  filter(year %in% c(ref_year, exp_year+event_time), 
         (playyear1 == exp_year) | (playyear1 > max(exp_year, exp_year+event_time)), 
         !is.na(valScoreM1))


# more clean-up
data_analysis[, 'y'] = data_analysis[, y_name]  # rename outcome to generic "y"
data_analysis = data_analysis %>% filter(!is.na(y))  # subset to counties w/ non-missing outcomes
data_analysis$y = log(data_analysis$y)  # log scale

# compute relevant data
data_analysis = data_analysis %>%
  group_by(fips) %>%
  summarize(deltay = y[year==(exp_year+event_time)] - y[year==ref_year],
            treat = valScoreM1[year==ref_year],
            treat_shale = valScoreM1_2[year==ref_year],
            playyear1 = playyear1[year==ref_year],
            shale_play1 = shale_play1[year==ref_year]) %>%
  ungroup() %>%
  mutate(treat = ifelse(playyear1 > (ref_year+1), 0, treat),
         treat_shale = ifelse(playyear1 > (ref_year+1), 0, treat_shale))

# join demo data
# inner join bc left join results in one county with missing demo because it was created in 2001
data_analysis = inner_join(data_analysis,
                           demo %>% select(-name, -total_pop),
                           by = c('fips' = 'geoid'))

# format data
deltay = data_analysis$deltay
treat = data_analysis$treat
treat = treat / max(treat)
if (which_treat == 'within_shale') (treat = data_analysis$treat_shale)
x = model.matrix(~ pop_log + female + age + white + black + hispanic, data=data_analysis)
x = x[,-1]
x = scale(x)


# Analysis ----------------------------------------------------------------

# number of cross-fitting folds
cf_folds = 10

# nuisance function estimator
nuisance_estimator = 'sl'

# delta to examine
delta_vals = seq(from=-100, to=100, length.out=201)

# perform analysis
estimates = estimate_psi(deltay, treat, x,
                         delta_vals = delta_vals, K = cf_folds,
                         nuisance_estimator = nuisance_estimator)

# save results
saveRDS(estimates, here('results', 'fracking',
                        paste0('treat', gsub('_shale', '', which_treat), '_',
                               nuisance_estimator, '_',
                               'cfK', cf_folds, '_',
                               gsub('_', '', y_name), '_', 
                               'refyear', ref_year,
                               '_lag', gsub('-', 'neg', event_time),
                               '.rds')))

