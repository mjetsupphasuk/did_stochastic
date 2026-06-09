
# load libraries
library(here)
library(tidyverse)

# set ggplot theme
theme_set(theme_bw(base_size = 18))

# import functions
source(here('code', 'utils.R'))

# data generating scenario
scenario = 2
source(here('code', 'sims_scenarios_mtp.R'))

# delta values want
delta_vals = seq(-100, 100, length.out = 201)

# cross-fitting
cf_folds = 5


# Disaggregated parameter results -----------------------------------------

num_sims = 1000
results_path = here('results', 'sims', 'interim')

# read in and aggregate results
results.all = list()
for (sim in 1:num_sims) {
  
  # read in results
  read_string = file.path(results_path, paste0('sim_res_mtp_sc', scenario,
                                               '_cfK', cf_folds,
                                               '_sim', sim, '.rds'))
  results.tmp = readRDS(read_string)

  # format and stash
  results.all[[sim]] = format_results_mtp(results.tmp, sim)
  
  # track progress
  if (sim %% 100 == 0) {
    print(sprintf('Completed sim %i at time %s', sim, Sys.time()))
  }
  
}
results.all = bind_rows(results.all)


# get truth
truth_file_path <- here("results", "sims", paste0("asdt_mtp_true_sc", scenario, ".rds"))
asdt_true = readRDS(truth_file_path)
results.all = results.all %>% left_join(asdt_true, by=c('g','t','e','delta'))

# average psi hat
results.all = results.all %>%
  group_by(nuisance, delta, g, t) %>%
  mutate(psihat_avg = mean(psihat)) %>%
  ungroup()

# nicer nuisance names and time indices
results.all = results.all %>%
  mutate(nuisance2 = case_when(nuisance=='bart' ~ 'BART',
                               nuisance=='glm' ~ 'GLMs',
                               nuisance=='oracle' ~ 'Oracle'),
         g_label = paste0("g=", g),
         e_label = paste0("e=", e))

# plot for each nuisance estimator
for (nuisance_estimator in unique(results.all$nuisance2)) {
  
  ggplot(results.all %>% filter(abs(delta) <= 10, nuisance2 == nuisance_estimator)) +
    geom_line(mapping = aes(x=delta, y=psihat, group=sim), color='grey', alpha=0.1) +
    geom_line(mapping = aes(x=delta, y=asdt), inherit.aes=FALSE, color='red') +
    geom_line(mapping = aes(x=delta, y=psihat_avg), inherit.aes=FALSE, color='blue', linetype='dashed') +
    facet_grid(rows = vars(e_label), cols = vars(g_label)) +
    labs(
      title = "Disaggregation Estimates",
      subtitle = paste("Scenario", scenario, "| Estimator:", nuisance_estimator),
      x = expression(delta), y = 'ASDT'
    ) +
    theme(panel.grid = element_blank())
  
  ggsave(here('results', 'sims', 
              paste0('result_plot_mtp_sc', scenario, '_cfK', cf_folds, '_', nuisance_estimator, '.png')),
         width=6, height=4.5, units='in')
}



# compute coverage
results.all <- results.all %>%
  mutate(cover = (asdt >= ci_lower) & (asdt <= ci_upper))

coverage_ged <- results.all %>%
  group_by(nuisance2, g_label, e_label, delta) %>%
  summarize(
    coverage = mean(cover, na.rm = TRUE),
    bias = mean(psihat - asdt, na.rm = TRUE),
    mse = mean((psihat - asdt)^2, na.rm = TRUE),
    .groups = "drop"
  )

# plot coverage
for (nuisance_estimator in unique(results.all$nuisance2)) {
  
  ggplot(coverage_ged %>% filter(abs(delta) <= 10, nuisance2 == nuisance_estimator), 
         aes(x=delta, y=coverage)) +
    geom_line(color="forestgreen") +
    geom_point(color="forestgreen", size=1) +
    geom_hline(yintercept=0.95, color='red', linetype='dashed') +
    facet_grid(rows = vars(e_label), cols = vars(g_label)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    labs(
      title = "Disaggregated Coverage",
      subtitle = paste("Scenario", scenario, "| Estimator:", nuisance_estimator),
      x = expression(delta), y = 'Coverage Rate'
    ) +
    theme(panel.grid = element_blank())
  
  ggsave(here('results', 'sims', 
              paste0('result_coverage_mtp_sc', scenario, '_cfK', cf_folds, '_', nuisance_estimator, '.png')),
         width=6, height=4.5, units='in')
}


# summary stats
sim_stats = results.all %>%
  group_by(nuisance, delta, g, e) %>%
  summarize(bias = mean(psihat - asdt),
            mse = mean((psihat - asdt)^2),
            emp_se = sd(psihat),
            avg_se = mean(sqrt(psihat_var)),
            coverage = mean(cover)) %>%
  mutate_all(round, 3)
sim_stats %>% filter(delta %in% c(-10,-5,0,5,10))



# Aggregated parameter results --------------------------------------------

num_sims = 1000
results_path = here('results', 'sims', 'interim')

# read in and aggregate results
agg_results.list <- list()
for (sim in 1:num_sims) {
  
  # read in results for this sim
  agg_file_name <- paste0(
    "sim_res_mtp-agg_sc", scenario,
    "_cfK", cf_folds,
    "_sim", sim,
    ".rds"
  )
  agg_read_string <- here(results_path, agg_file_name)
  agg_res = readRDS(agg_read_string)
  
  # format
  results.agg = lapply(1:length(agg_res), function(i) {
    df = agg_res[[i]]
    df$nuisance = names(agg_res)[i]
    df$sim = sim
    return(df)
  })
  agg_results.list[[sim]] = bind_rows(results.agg)
}
agg_results.all <- bind_rows(agg_results.list)


# get truth
agg_truth_file <- here("results", "sims", paste0("asdt_mtp-agg_true_sc", scenario, ".rds"))
agg.true = readRDS(agg_truth_file)
agg_results.all = agg_results.all %>% 
  left_join(agg.true, by=c("type", "target", "delta"))

# compute coverage
agg_results.all <- agg_results.all %>%
  mutate(
    ci_lower = estimate - qnorm(0.975) * se,
    ci_upper = estimate + qnorm(0.975) * se,
    cover = (true_estimate >= ci_lower) & (true_estimate <= ci_upper)
  )

# average estimates across sims
agg_results.all = agg_results.all %>%
  group_by(nuisance, delta, type, target) %>%
  mutate(est_mean = mean(estimate)) %>%
  ungroup()

# nicer nuisance names and type/target names
# also change "simple" to "overall"
# and "group" to "cohort"
agg_results.all = agg_results.all %>%
  mutate(nuisance2 = case_when(nuisance=='bart' ~ 'BART',
                               nuisance=='glm' ~ 'GLMs',
                               nuisance=='oracle' ~ 'Oracle'),
         type = ifelse(type=='simple', 'overall', type),
         type2 = gsub('_', ' ', type) %>% str_to_title(),
         type2 = ifelse(type2=='Group', 'Cohort', type2))


# cycle through and process each type of aggregation
agg_types <- unique(agg_results.all$type2)

for (atype in agg_types) {
  cat("\n -> Processing type:", toupper(atype), "\n")
  
  # subset to specified aggregation type
  df_type <- agg_results.all %>% filter(type2 == atype)
  if (atype != "Overall") {
    df_type <- df_type %>% filter(tolower(target) != "overall")
  }
  
  # Generate and Save Estimation Plot
  for (nuisance_estimator in unique(df_type$nuisance2)) {
    
    p_est <- ggplot(df_type %>% filter(abs(delta) <= 10, nuisance2 == nuisance_estimator)) +
      geom_line(aes(x = delta, y = estimate, group = sim), 
                color = "gray60", alpha = 0.15, linewidth = 0.3) +
      geom_line(aes(x = delta, y = true_estimate), 
                color = "red", linewidth = 1.2) +
      geom_line(aes(x = delta, y = est_mean), 
                color = "steelblue", linewidth = 1.2, linetype = "dashed") +
      facet_grid(. ~ target) +
      labs(
        title = paste("Aggregation:", atype),
        subtitle = paste("Scenario", scenario, "| Estimator:", nuisance_estimator),
        x = expression(delta), y = paste("ASDT", atype)
      ) +
      theme(strip.background = element_rect(fill = "gray90"))
    
    ggsave(here('results', 'sims', 
                paste0('result_plot_mtp-agg_', atype, '_sc', scenario, '_cfK', cf_folds, '_', nuisance_estimator, '.png')),
           plot = p_est,
           width=6, height=4.5, units='in')
  
  }
  
  
  # Generate and Save Coverage Plot
  
  # compute coverage
  df_summary <- df_type %>%
    group_by(target, delta, nuisance2) %>%
    summarise(
      coverage = mean(cover, na.rm = TRUE),
      .groups = "drop"
    )
  
  for (nuisance_estimator in unique(df_type$nuisance2)) {
    p_cov <- ggplot(data = df_summary %>% filter(abs(delta) <= 10, nuisance2 == nuisance_estimator), 
                    aes(x = delta, y = coverage)) +
      geom_line(color = "forestgreen", linewidth = 1) +
      geom_point(color = "forestgreen", size = 1.5) +
      geom_hline(yintercept = 0.95, color = "red", linetype = "dashed", linewidth = 0.8) +
      facet_grid(. ~ target) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
      labs(
        title = paste("Coverage:", atype),
        subtitle = paste("Scenario", scenario, "| Estimator:", nuisance_estimator),
        x = expression(delta), y = "Coverage Rate"
      ) +
      theme(strip.background = element_rect(fill = "gray90"))
    
    # print(p_cov)
    
    ggsave(here('results', 'sims',
                paste0('result_coverage_mtp-agg_', atype, '_sc', scenario, '_cfK', cf_folds, '_', nuisance_estimator, '.png')),
           plot = p_cov,
           width=6, height=4.5, units='in')
  }
}



# summary stats
sim_stats = agg_results.all %>%
  group_by(nuisance, delta, type, target) %>%
  summarize(bias = mean(estimate - true_estimate),
            mse = mean((estimate - true_estimate)^2),
            emp_se = sd(estimate),
            avg_se = mean(se),
            coverage = mean(cover)) %>%
  mutate_if(is.numeric, round, 3)
sim_stats %>% filter(delta %in% c(-10,-5,0,5,10))



