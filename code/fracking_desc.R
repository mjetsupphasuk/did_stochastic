
setwd("/work/users/m/j/mjets/dissertation/paper2")
.libPaths("/nas/longleaf/home/mjets/RLibs")

library(haven)
library(readr)
library(BART)
library(hal9001)
library(SuperLearner)
library(caret)
library(tidycensus)
library(tigris)
library(viridis)
library(tidyverse)

theme_set(theme_bw(base_size = 14))

source('code/utils.R')


# Demographic data --------------------------------------------------------

# census api key
census_api_key("INSERT YOUR KEY HERE")

# extract data from ACS
demo = get_acs(geography = 'county',
               variables = c('total_pop' = 'B01003_001',
                             'female' = 'B01001_026',
                             'white' = 'B02001_002',
                             'black' = 'B02001_003',
                             'hispanic' = 'B03001_003',
                             'age' = 'B01002_001'),
               year = 2009,
               survey = 'acs5')

# format demo variable
colnames(demo) = tolower(colnames(demo))
demo = pivot_wider(demo %>% select(-moe), names_from = 'variable', values_from = 'estimate')

# create variables
demo = demo %>%
  mutate(pop_log = log(total_pop),
         female = female / total_pop,
         white = white / total_pop,
         black = black / total_pop,
         hispanic = hispanic / total_pop) %>%
  select(-hhincome_med)

# save
saveRDS(demo, 'data/fracking/demo.rds')


# Explore -----------------------------------------------------------------

# read in data
data = read_dta('data/fracking/county_clean_long.dta')
data_flat = read_dta('data/fracking/county_flat.dta')
hdpi = read_dta('data/fracking/hpdi_allwells_shale.dta')

# subset of data used in Callaway, Goodman-Bacon, Sant'Anna (CGS)
data_cgs = data %>% filter(year >= 2001, year <= 2014)

# valScoreM1: area-weighted max rystad score for the play with the earliest first frac date
# see 3C in 2_DataConstruction_ReadME and merge_shale_geo.do
# matches what was used in CGS
hist(data_cgs$valScoreM1)
range(data_cgs$valScoreM1[data_cgs$valScoreM1 > 0], na.rm=TRUE)
length(unique(data_cgs$valScoreM1))

hist(log(data_cgs$emp_tot))

# 38 never-treated
sum(data_flat$valScoreM1==0, na.rm=TRUE)
data_cgs %>%
  group_by(year) %>%
  summarize(sum(valScoreM1==0, na.rm=TRUE))

# not yet treated and treated
data_cgs %>%
  mutate(playyear1 = ifelse(valScoreM1==0, 2020, playyear1)) %>%
  filter(!is.na(valScoreM1)) %>%
  group_by(year) %>%
  summarize(not_yet = sum(playyear1 > year | valScoreM1==0, na.rm=TRUE),
            treated = sum(playyear1 == year & valScoreM1 > 0, na.rm=TRUE))

# same as above but also removing missing outcomes (specify)
sample_size = data_cgs %>%
  mutate(playyear1 = ifelse(valScoreM1==0, 2020, playyear1)) %>%
  filter(!is.na(valScoreM1), !is.na(emp_tot)) %>%
  group_by(year) %>%
  summarize(not_yet = sum(playyear1 > year | valScoreM1==0, na.rm=TRUE),
            treated = sum(playyear1 == year & valScoreM1 > 0, na.rm=TRUE))

# save sample size
write_csv(sample_size, 'results/fracking/sample_size.csv')


# shale plays
table(data_cgs$shale_play1[data_cgs$year==2008])

# verify that valScoreM1 is not time-varying in this dataset
data %>%
  select(fips, year, valScoreM1) %>%
  group_by(fips) %>%
  summarize(unique_treat = length(unique(valScoreM1))) %>%
  pull(unique_treat) %>%
  table(useNA='ifany')



# Process data ------------------------------------------------------------

# read in data
data = read_dta('data/fracking/county_clean_long.dta')

# divide prospectivity score by max within each shale for comparability
max_shale_prospect = data %>%
  filter(year==2007, !is.na(shale_play1)) %>%
  group_by(shale_play1) %>%
  summarize(max_shale_prospect = max(valScoreM1, na.rm=TRUE)) %>%
  ungroup()

data = data %>%
  left_join(max_shale_prospect, by = 'shale_play1') %>%
  mutate(valScoreM1_2 = valScoreM1 / max_shale_prospect)

# plot prospectivity score by shale play
# (recall that prosp scores non-time-varying in this dataset)
ggplot(data %>% filter(year==2007, !is.na(valScoreM1)), aes(valScoreM1_2)) +
  geom_histogram() +
  facet_grid(vars(shale_play1))

# if prospectivity score is 0, set playyear1 (exposure cohort) to large number
# i.e., consider never treated
data = data %>%
  mutate(playyear1 = ifelse(valScoreM1==0, 2020, playyear1))

# save data
saveRDS(data, 'data/fracking/data_analysis.rds')



# Maps of prospectivity scores --------------------------------------------


# plotting function
plot_map <- function(data) {
  gg = ggplot(data = data, aes(fill = treat_shale)) +
    geom_sf(color=gray(0.5), alpha=0.75) +
    coord_sf(xlim = c(-125, -68), ylim = c(26, 48)) +
    scale_fill_viridis('Prospectivity score  ') +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = 'bottom',
          panel.grid = element_blank(),
          panel.border = element_blank(),
          legend.key.width = unit(30, 'pt')) +
    labs(x='', y='')
  print(gg)
}


# load counties shapefile
counties_sf = counties(cb = TRUE, year = 2021)

# data to plot - any year
data.plot = data %>% select(fips, treat_shale=valScoreM1_2) %>% distinct()
data.plot = left_join(counties_sf, data.plot, by=c('GEOID'='fips'))

# plot and save
plot_map(data.plot)
ggsave('results/fracking/prospectivity_map_any.png')


# data to plot - 2008
data.plot = data %>% 
  filter(playyear1 >= 2008) %>% 
  mutate(treat_shale = ifelse(playyear1 > 2008, 0, valScoreM1_2)) %>%
  select(fips, treat_shale) %>% 
  distinct()
data.plot = left_join(counties_sf, data.plot, by=c('GEOID'='fips'))

# plot and save
plot_map(data.plot)
ggsave('results/fracking/prospectivity_map_2008.png')


# also plot histogram
ggplot(data.plot %>% filter(!is.na(treat_shale)), aes(x=treat_shale)) +
  geom_histogram(fill='grey', binwidth=0.1) +
  labs(x = 'Prospectivity score (rescaled)', y = 'Number of counties') +
  theme(panel.grid = element_blank())
ggsave('results/fracking/treatment_shale_hist.png', width=6, height=4.5, units='in')

# data to plot - 2008
data.plot = data %>% 
  filter(playyear1 >= 2009) %>% 
  mutate(treat_shale = ifelse(playyear1 > 2009, 0, valScoreM1_2)) %>%
  select(fips, treat_shale) %>% 
  distinct()
data.plot = left_join(counties_sf, data.plot, by=c('GEOID'='fips'))

# plot and save
plot_map(data.plot)
ggsave('results/fracking/prospectivity_map_2009.png')


# Compare pre-treatment characteristics -----------------------------------

# cases
ref_year = 2007
exp_year = ref_year + 1
event_time = 0
y_name = 'emp_tot'
which_treat = 'within_shale'

# subset units to treated and not-yet-treated
# for negative event times (i.e., pre-tests), define treated and untreated groups at exp_year
data_analysis = data %>%
  filter(year %in% c(ref_year, exp_year+event_time), 
         (playyear1 == exp_year) | (playyear1 > max(exp_year, exp_year+event_time)), 
         !is.na(valScoreM1))


# more clean-up
data_analysis[, 'y'] = data_analysis[, y_name]  # rename outcome to generic "y"
data_analysis = data_analysis %>% filter(!is.na(y))  # subset to counties w/ non-missing outcomes
data_analysis$y = log(data_analysis$y)  # log scale



datay = data_analysis %>% 
  filter(year==2008) %>%
  select(fips, y, valScoreM1_2, playyear1) %>%
  mutate(treat_shale = ifelse(playyear1 > (ref_year+1), 0, valScoreM1_2))

data.plot = left_join(counties_sf, datay, by=c('GEOID'='fips'))
# symmetric_limits = c(-max(abs(data.plot$deltay), na.rm=TRUE), max(abs(data.plot$deltay), na.rm=TRUE))
ggplot(data = data.plot, aes(fill = y)) +
  geom_sf(color=gray(0.5), alpha=0.75) +
  coord_sf(xlim = c(-125, -68), ylim = c(26, 48)) +
  scale_fill_viridis_c('Total employment (log)') +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'bottom',
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.key.width = unit(30, 'pt')) +
  labs(x='', y='')
ggsave('results/fracking/y_map_2008.png')

mean(datay$y)
sd(datay$y)


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

# plot map of delta Y
data.plot = left_join(counties_sf, data_analysis, by=c('GEOID'='fips'))
symmetric_limits = c(-max(abs(data.plot$deltay), na.rm=TRUE), max(abs(data.plot$deltay), na.rm=TRUE))
ggplot(data = data.plot, aes(fill = deltay)) +
  geom_sf(color=gray(0.5), alpha=0.75) +
  coord_sf(xlim = c(-125, -68), ylim = c(26, 48)) +
  scale_fill_distiller('Change in total employment', palette = 'RdBu', limits=symmetric_limits, direction = 1) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'bottom',
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.key.width = unit(30, 'pt')) +
  labs(x='', y='')
ggsave('results/fracking/deltay_map_2008.png')




# means
cov_means = data_analysis %>%
  mutate(treat_bin = treat > 0) %>%
  select(treat_bin, female:hhincome_med_log) %>%
  group_by(treat_bin) %>%
  summarize_all(mean)

# standard deviations
cov_sds = data_analysis %>%
  mutate(treat_bin = treat > 0) %>%
  select(treat_bin, female:hhincome_med_log) %>%
  group_by(treat_bin) %>%
  summarize_all(sd)

# t-tests
t.test(data_analysis$white[data_analysis$treat==0], data_analysis$white[data_analysis$treat > 0])
t.test(data_analysis$black[data_analysis$treat==0], data_analysis$black[data_analysis$treat > 0])
t.test(data_analysis$hispanic[data_analysis$treat==0], data_analysis$hispanic[data_analysis$treat > 0])
t.test(data_analysis$female[data_analysis$treat==0], data_analysis$female[data_analysis$treat > 0])
t.test(data_analysis$age[data_analysis$treat==0], data_analysis$age[data_analysis$treat > 0])
t.test(data_analysis$pop_log[data_analysis$treat==0], data_analysis$pop_log[data_analysis$treat > 0])
t.test(data_analysis$poverty[data_analysis$treat==0], data_analysis$poverty[data_analysis$treat > 0])
t.test(data_analysis$hhincome_med_log[data_analysis$treat==0], data_analysis$hhincome_med_log[data_analysis$treat > 0])

# save
write_csv(cov_means, paste0('results/fracking/cov_means_', exp_year, '.csv'))
write_csv(cov_sds, paste0('results/fracking/cov_sds_', exp_year, '.csv'))



