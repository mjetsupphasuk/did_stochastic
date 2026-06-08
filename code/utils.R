

# Simulations -------------------------------------------------------------


gen_data <- function(n, pi_bin, mu.0_f, shape1_alpha_f, shape2_beta_f, mu.d_f) {
  
  # covariates: 10 standard uniform
  X = matrix(runif(n*10, -1, 1), nrow=n, ncol=10)
  colnames(X) = paste0('X', 1:ncol(X))
  
  # treatment
  treat_bin = rbinom(n, size=1, prob=pi_bin(X))
  treat_d = rbeta(n, shape1 = shape1_alpha_f(X), shape2 = shape2_beta_f(X))
  treat = ifelse(treat_bin==1, treat_d, 0)
  
  # outcome
  epsilon = rnorm(n, 0, 1)
  deltaY = (1-treat_bin)*mu.0_f(X) + treat_bin*mu.d_f(X, treat) + epsilon
  
  # put in data frame
  data = data.frame(deltaY=deltaY, treat=treat, treat_bin=treat_bin, treat_d=treat_d, X)
  return(data)
  
}


# Simulate the panel data for MTP (Multiple Time Periods) with staggered adoption design
# and continuous treatment doses
gen_staggered_data <- function(n, shape1_alpha_f, shape2_beta_f,
                               T_periods, c_fun, theta_t,
                               cohort_predictor) {
  
  # covariates: 10 standard uniform variables ranging from -1 to 1
  X <- matrix(runif(n * 10, min = -1, max = 1), nrow = n, ncol = 10)
  colnames(X) <- paste0("X", 1:ncol(X))
  
  # X_index Multinomial Logit Assignment
  # Determines the latent score for cohort (adoption time) assignment
  if (!is.null(cohort_predictor)) {
    Xindex <- as.numeric(cohort_predictor(X)) 
  } else {
    # Default linear combination if no non-linear function is provided
    b <- seq(-0.1, 0.1, length.out = ncol(X))
    Xindex <- as.numeric(X %*% b)
  }
  
  # Probability of being in adoption time g (using a multinomial logistic link)
  # Divided by 2*T_periods to scale the exponent and satisfy positivity/overlap assumptions
  g_support <- 1:(T_periods + 1)
  gamG <- g_support / (2 * T_periods)
  prG <- exp(outer(Xindex, gamG)) / apply(exp(outer(Xindex, gamG)), 1, sum)
  
  # adoption time G: Sampled based on the calculated multinomial probabilities
  G <- apply(prG, 1, function(pvec) sample(g_support, size = 1, prob = pvec))
  
  # True probability of adopting treatment by time t
  PAgt0_true <- sapply(1:T_periods, function(t) {
    rowSums(prG[, 1:t, drop = FALSE])
  })
  
  # treatment dose  (D | X, G  ~ Beta(alpha(X,G), beta(X,G)))
  # Dose is drawn from a Beta distribution with covariate-dependent shape parameters
  treat_d <- sapply(1:n, function(i) {
    rbeta(1,
          shape1 = shape1_alpha_f(X[i, , drop = FALSE], G[i]),
          shape2 = shape2_beta_f(X[i, , drop = FALSE], G[i]))
  })
  treat_d <- as.numeric(treat_d)
  
  # Create t periods range from 0,...., T
  t_periods <- 0:T_periods
  
  # treatment path: when t = 0, receive no treatment. Otherwise follows d*1(t>=G)
  treat <- sapply(t_periods, function(t) {
    if (t == 0) return(rep(0, n))
    treat_d * (t >= G)
  })
  # Binary indicator of whether treatment has started
  treat_bin <- 1*(treat > 0)
  
  # Idiosyncratic error term drawn from standard normal distribution
  epsilon_it <- matrix(rnorm(n * length(t_periods), 0, 1),
                       nrow = n, ncol = length(t_periods))
  
  # Generate absolute potential outcomes for each time period
  Y_matrix <- sapply(seq_along(t_periods), function(j) {
    
    t <- t_periods[j]
    D_it <- treat[, j]
    Abin_it <- treat_bin[, j]
    base_trend <- theta_t[j] # Baseline time fixed effect
    
    # Non-parallel trend covariate effect: trend slope function * t
    if (!is.null(trend_slope_f)) {
      cov_effect <- as.numeric(trend_slope_f(X)) * t
    } else {
      cov_effect <- as.numeric(X %*% gamma_t[j, ])
    }
    
    # Polynomial Treatment Effect (Dynamic based on dose D and event time e = t - G)
    if (!is.null(trend_slope_f)) {
      coefs <- c_fun(t-G, G) 
      dose_effect <- D_it * coefs[,1] + (D_it^2) * coefs[,2] + (D_it^3) * coefs[,3]
      trt_effect <- Abin_it * as.numeric(dose_effect)
    } else {
      trt_effect <- Abin_it * as.numeric(0.5 * c_fun(t-G, G) * D_it)
    }
    
    # Final observed outcome combines base trend, covariate trend, treatment effect, and noise
    base_trend + cov_effect + trt_effect + epsilon_it[, j]
  })
  
  # put in wide data frame format initially
  data <- data.frame(id = 1:n, G = G, treat_d = treat_d, X)
  data[paste0("treat", t_periods)] <- treat
  data[paste0("treat_bin", t_periods)] <- treat_bin
  data[paste0("Y", t_periods)] <- Y_matrix
  
  # keep PAgt0_true only for t=1,...,T
  data[paste0("PAgt0_true",1:T_periods)] <- PAgt0_true
  
  # Transform to panel data (Long format for estimation)
  panel_data <- tidyr::pivot_longer(
    data,
    cols = c( dplyr::matches("^Y\\d+$"),
              dplyr::matches("^treat\\d+$"),
              dplyr::matches("^treat_bin\\d+$"),
              dplyr::matches("^PAgt0_true\\d+$")
    ),
    names_to = c(".value", "periods"),
    names_pattern = "^(.*?)(\\d+)$"
  )
  
  panel_data$periods <- as.numeric(panel_data$periods)
  # Calculate relative exposure time (event time) e = t - G
  panel_data$e <- panel_data$periods - panel_data$G
  panel_data <- panel_data[order(panel_data$id, panel_data$periods), ]
  
  return(panel_data)
}


# function to compute true ASDT
compute_true_asdt <- function(data, delta, shape1_alpha_f, shape2_beta_f) {
  
  # P(A > 0), needs to be estimated if logistic model used
  pAgt0_true = mean(data$treat_bin)
  
  # difference in mu_0(X) and mu_D(X)
  beta0 = alphad - alpha0
  
  # sequence of density values
  d_seq = seq(0.001, 0.999, length.out=500)
  
  # true prop scores
  X = data %>% select(starts_with('X'))
  prop_scores_mat = apply(X, 1, function(x) true_density_f(d_seq, t(x), shape1_alpha_f, shape2_beta_f)) %>% t()
  
  # computes asdt for a single delta
  compute_asdt_singledelta <- function(delta) {

    # compute tilts
    q_delta = exp_tilt(params=list('delta'=delta), pihats=prop_scores_mat, 
                       treat_points=d_seq, design_points=d_seq)
    
    # calculate dose function
    # hard-coded for now for a specific nonlinear dose effect (scenarios 1+2)
    d_coef = c_fun(0,1)
    dose_effect = d_seq*d_coef[1] + (d_seq^2)*d_coef[2] + (d_seq^3)*d_coef[3]
    
    # calculate \int f(d) q(d|X) d(d)
    int_d = rowMeans(q_delta %*% diag(dose_effect))
    
    # calculate ASDT
    asdt = beta0 + mean((data$treat_bin / pAgt0_true) * int_d)
    
    return(asdt)
  }
  
  # compute asdt for all delta
  asdt = sapply(delta, compute_asdt_singledelta)
  
  # return vector of asdt
  return(asdt)
  
}


# true ASDT Computation for a specific cohort g and time t
compute_true_stag_asdt <- function(panel_data, delta, g_spec, t_spec, 
                                   shape1_alpha_f, shape2_beta_f, c_fun, mu.0_f = NULL, mu.d_f = NULL,
                                   gamma_t = NULL, theta_t = NULL,
                                   T_periods = NULL) {
  
  d_seq <- seq(0.001, max(panel_data$treat_d, na.rm=TRUE), length.out = 500)
  e_spec = t_spec - g_spec
  
  unit <- panel_data %>% dplyr::distinct(id, .keep_all = TRUE) %>% arrange(id)
  X <- unit %>% select(starts_with('X')) %>% as.matrix()
  
  idx_g <- which(unit$G == g_spec)
  if (length(idx_g) == 0) return(rep(NA_real_, length(delta)))
  
  # Calculate true propensity scores matrix
  prop_scores_mat_g <- t(apply(X, 1, function(x) {
    true_density_f(d_seq, matrix(x, nrow=1), shape1_alpha_f, shape2_beta_f)
  }))
  
  # computes asdt for a single delta
  compute_stag_asdt_singledelta <- function(delta) {
    # Apply exponential tilt
    q_delta_g <- exp_tilt(params=list('delta'=delta), pihats=prop_scores_mat_g, 
                          treat_points=d_seq, design_points=d_seq)
    
    # Calculate true polynomial dose effect
    coefs <- c_fun(e_spec, g_spec)
    if (length(coefs) >= 3) {
      if(is.matrix(coefs)) coefs <- as.numeric(coefs[1,]) 
      dose_effect <- (t_spec >= g_spec) * (d_seq * coefs[1] + (d_seq^2) * coefs[2] + (d_seq^3) * coefs[3])
    } else {
      dose_effect <- (t_spec >= g_spec) * (0.5 * as.numeric(coefs) * d_seq)
    }
    
    # Integrate conditional mean over the tilted dose distribution
    int_mu1_g <- rowMeans(q_delta_g %*% diag(dose_effect))
    asdt = mean( ((unit$G==g_spec) / mean(unit$G==g_spec)) * int_mu1_g )
    
    return(asdt)
  }
  
  return(sapply(delta, compute_stag_asdt_singledelta))
}



# Estimators --------------------------------------------------------------


# `estimate_pi`, `aggregate_muhat` and `estimate_psi` adapted from
# Schindl et al. (2026), "Incremental effects for continuous exposures"
# https://arxiv.org/abs/2409.11967


# kernel transforming treatment density into mean estimation problem
# returns matrix with dim: (number of obs [rows]) x (number of design points [cols])
# each value in the matrix is equal to pi(a|x)
estimate_pi <- function(train_a, design_a, train_x, validate_x, nuisance_estimator) {
  
  #Initializing empty list
  pi_predictions <- list()
  
  #For each design point a_i we estimate E[1/hK((A-a_i)/h)|X]
  for (i in 1:length(design_a)) {
    
    # set bandwidth  
    h_opt = 0.05
    
    #Estimating E[1/hK((A-a_i)/h)|X] now that we have selected the optimal bandwidth.
    h_out = dnorm((train_a - design_a[i])/h_opt) / h_opt
    
    # just use glm for now
    pimod <- lm(h_out ~ ., data=data.frame(train_x))
    pi_predictions[[i]] = predict(pimod, newdata=data.frame(validate_x))
    
    # #Track progress
    # if (i %% 10 == 0) {
    #   print(sprintf('Completed kernel regression for design point # %i / %i; Time: %s',
    #                 i, length(design_a), Sys.time()))
    # }
    
  }
  
  # Creates a data frame of estimates of \hat{\pi} for each design point
  pihats = do.call(cbind, pi_predictions)
  
  # Enforce non-negativity for probabilities
  pihats[pihats < 0] = 0
  
  # Numerical integration to ensure density sums to 1
  pihats = pihats / rowMeans(pihats)
  
  return(pihats)
  
}


# takes in mu model and predicts for all doses (points in validate_a)
aggregate_muhat <- function(mumod, validate_a, validate_x, nuisance_estimator) {
  
  # For each design point d_1,...,d_m we obtain predicted values for mu_d(X)
  predicted_mu <- function(z) {
    
    if (nuisance_estimator == 'sl') {
      sl.pred = predict(mumod, newdata=cbind(validate_x, treat=z))
      return(sl.pred$pred[,1])
    }
    
    if (nuisance_estimator == 'hal') {
      return(predict(mumod, new_data=cbind(validate_x, treat=z), type='response'))
    }
    
    if (nuisance_estimator == 'bart') {
      return(predict(mumod, newdata = cbind(validate_x, treat=z)) %>% colMeans())
    }
    
    if (nuisance_estimator == 'glm') {
      return(predict(mumod, newdata=data.frame(validate_x, treat=z)))
    }
    
    if (nuisance_estimator == 'rf') {
      return(predict(mumod, data=data.frame(validate_x, treat=z))$predictions)
    }
  
  }
  
  #Creating dataframe of predicted values for each design point
  return(do.call(cbind, lapply(validate_a, predicted_mu)))
  
}


# estimator of ASDT
estimate_psi <- function(deltay, treat, x, stochastic_policy,
                         policy_params, nuisance_estimator, 
                         K=5, oracle_functions = NULL,
                         design_d = seq(0.001, 0.999, length.out = 500),
                         g=1, t=1, T_periods=1) {
  
  # sample size
  n = length(deltay)
  
  # create indicator for A=0 and A>0 groups
  treat_bin = as.integer(treat > 0)
  
  # split data into K folds
  folds = caret::createFolds(deltay, k=K, list=TRUE, returnTrain=FALSE)
  
  # store estimates for each k iteration
  psi_hat_k = list()
  psi_hat_var_k = list()
  psi1_hat_k = list()
  psi1_hat_var_k = list()
  psi2_hat_k = numeric(K)
  psi2_hat_var_k = numeric(K)
  covvals = list()
  psi_plugin_k = list()
  psi1_plugin_k = list()
  psi2_plugin_k = numeric()
  
  # initializing matrix to store influence functions
  # for variance estimation of aggregated parameters
  IF_matrix = matrix(0, nrow = n, ncol = length(delta_vals))
  colnames(IF_matrix) <- paste0("delta_", delta_vals)
  
  # analyze iteratively for each fold
  for (k in 1:K) {
    
    # nuisance (training) and estimation (validation) indices
    estimation_set = folds[[k]]
    nuisance_set = setdiff(1:n, estimation_set)
    if (K==1) (nuisance_set = estimation_set)
    
    # training sets are all folds except the validation set
    train_x = x[nuisance_set,]
    train_treat = treat[nuisance_set]
    train_treat_bin = treat_bin[nuisance_set]
    train_deltay = deltay[nuisance_set]
    
    # validation set
    validate_x = x[estimation_set,]
    validate_treat = treat[estimation_set]
    validate_treat_bin = treat_bin[estimation_set]
    validate_deltay = deltay[estimation_set]
    
    # Estimate Psi_2
    
    # nuisance function estimation
    if (nuisance_estimator == 'sl') {
      
      # outcome regression
      mu.model = SuperLearner(Y = train_deltay, X = cbind(train_x,treat=train_treat),
                              SL.library = c('SL.glm', 'SL.myBART', 'SL.myhal9001'))
      sl.mu.pred = predict(mu.model, newdata = cbind(validate_x,treat=0))
      mu_0_validate = sl.mu.pred$pred[,1]
      
      # binary treatment propensity score
      sl.prop = SuperLearner(Y = train_treat_bin, X = train_x,
                             family = 'binomial',
                             SL.library = c('SL.glm', 'SL.myBART', 'SL.myhal9001'))
      sl.prop.pred = predict(sl.prop, newdata=validate_x)
      pi_Agt0_validate = sl.prop.pred$pred[,1]
      
    }
    if (nuisance_estimator == 'hal') {
      
      # outcome regression
      mu.model = fit_hal(X=cbind(train_x,treat=train_treat), Y=train_deltay, max_degree=3, num_knots=c(25,10,5))
      mu_0_validate = predict(mu.model, new_data = cbind(validate_x,treat=0), type='response')
      
      # binary treatment propensity score
      hal.prop = fit_hal(X=train_x, Y=train_treat_bin, family = 'binomial', max_degree=3, num_knots=c(25,10,5))
      pi_Agt0_validate = predict(hal.prop, new_data=validate_x, type='response')
      
    }
    if (nuisance_estimator == 'bart') {
      
      # outcome regression
      mu.model = gbart(x.train=cbind(train_x,treat=train_treat), y.train=train_deltay,
                       ndpost=500, nskip=500, keepevery=1, printevery=1000)
      mu_0_validate = predict(mu.model, newdata = cbind(validate_x,treat=0)) %>% colMeans()
      
      # binary treatment propensity score
      bart.prop = gbart(x.train=train_x, y.train=train_treat_bin, type='pbart', 
                        ndpost=500, nskip=500, keepevery=1, printevery=1000)
      pi_Agt0_validate = predict(bart.prop, newdata=validate_x)$prob.test.mean
      
    }
    if (nuisance_estimator == 'glm') {
      
      # outcome regression
      mu.model = lm(deltay ~ ., data = data.frame(deltay=train_deltay, train_x, treat=train_treat))
      mu_0_validate = predict(mu.model, newdata = data.frame(validate_x,treat=0))
      
      # binary treatment propensity score
      pi_Agt0_train = glm(treat_bin ~ ., data = data.frame(treat_bin=train_treat_bin, train_x), family='binomial')
      pi_Agt0_validate = predict(pi_Agt0_train, newdata = data.frame(validate_x), type='response')
      
    }
    if (nuisance_estimator == 'oracle') {
      
      # outcome regression
      mu_0_validate = as.numeric(mu.0_f(as.matrix(validate_x), g=g, t=t))
      
      # binary treatment propensity score
      if (T_periods == 1) {
        pi_Agt0_validate = as.numeric(pi_bin(as.matrix(validate_x)))
      } else {
        Xindex <- as.numeric(cohort_predictor(validate_x))
        g_support <- 1:(T_periods + 1)
        gamG <- g_support / (2 * T_periods)
        PrG_val <- exp(outer(Xindex, gamG)) / apply(exp(outer(Xindex, gamG)), 1, sum)
        
        # compute Pr(G=g|X, I(G=g)+I(G>t)=1)
        num = PrG_val[, g]
        if (t < T_periods) {
          denom = PrG_val[, g] + rowSums(PrG_val[, (t + 1):(T_periods + 1), drop = FALSE])
        } else {
          denom = PrG_val[, g] + PrG_val[, (T_periods + 1)]
        }
        pi_Agt0_validate = num / denom
      }
      
    }
    
    # plug-in OR-based estimator
    pAgt0 = mean(validate_treat_bin)
    psi2_plugin = mean((validate_treat_bin/pAgt0) * mu_0_validate)
    psi2_plugin_k[k] = psi2_plugin
    
    # hajek weights
    psi2_weights = (1-validate_treat_bin) * pi_Agt0_validate / (1-pi_Agt0_validate)
    
    # one-step bias correction
    psi2_onestep = (psi2_weights / mean(psi2_weights)) * (validate_deltay - mu_0_validate)
    psi2_hat_k[k] = psi2_plugin + mean(psi2_onestep)
    
    # influence function
    phi2_i = psi2_onestep + (validate_treat_bin / pAgt0)*(mu_0_validate - psi2_hat_k[k])
    
    # variance estimate
    psi2_hat_var_k[k] = mean(phi2_i^2)
    
    
    # Estimate Psi_1
    
    # estimates of mu_d(X) = E[delta Y | X, A>0, D=d] across all design points
    if (nuisance_estimator != 'oracle') {
      muhats_validate = aggregate_muhat(mu.model, design_d, validate_x[validate_treat>0,], nuisance_estimator)
    }
    if (nuisance_estimator == 'oracle') {
      muhats_validate = sapply(design_d, function(d) mu.d_f(as.matrix(validate_x[validate_treat>0,]), d,
                                                            effect_coefs = c_fun(t-g, g),
                                                            g = g,
                                                            t = t))
    }
    
    # estimates of pi(d|X) for all design points
    # returns a matrix where ncol = n; nrow = number of design points
    if (nuisance_estimator != 'oracle') {
      pihats_validate = estimate_pi(train_treat[train_treat>0], design_d, train_x[train_treat>0,], validate_x[validate_treat>0,], nuisance_estimator)
    }
    if (nuisance_estimator == 'oracle') {
      pihats_validate = apply(validate_x[validate_treat>0,], 1, function(row) true_density_f(design_d, t(row), shape1_alpha_f, shape2_beta_f)) %>% t()
    }
    
    # init empty vectors for iterating over the first parameter given in policy_params
    # call it "delta_vals" here since main stochastic policy is exponential tilt
    delta_vals = policy_params[[1]]
    delta_psi1hat = numeric(length(delta_vals))
    delta_psi1var = numeric(length(delta_vals))
    delta_covvals = numeric(length(delta_vals))
    delta_plugin = numeric(length(delta_vals))
    
    for (j in 1:length(delta_vals)) {
      
      delta <- delta_vals[j]
      
      # stochastic policy
      if (stochastic_policy=='exp_tilt') {
        stoch_func = exp_tilt
        stoch_params = list('delta' = delta)
      }
      if (stochastic_policy=='gauss_conc') {
        stoch_func = gauss_conc
        stoch_params = list('dprime' = delta, 'l' = policy_params[['l']])
      }
      
      # estimate q(D|X) / pi(D|X)
      q_div_pi = stoch_func(params = stoch_params, 
                            pihats = pihats_validate, 
                            treat_points = validate_treat[validate_treat>0], 
                            design_points = design_d, 
                            divide_by_pi = TRUE)
      
      # estimate int_b mu(X) q(b|X) db
      q_hat = stoch_func(params = stoch_params, 
                         pihats = pihats_validate, 
                         treat_points = validate_treat[validate_treat>0], 
                         design_points = design_d, 
                         divide_by_pi = FALSE)
      
      muqhat_validate = rowMeans(muhats_validate * q_hat)
      
      # influence function values
      psi1_i = numeric(length(estimation_set))
      psi1_i[validate_treat > 0] <- q_div_pi * (validate_deltay[validate_treat>0] - muqhat_validate) + muqhat_validate
      psi1_i = psi1_i / pAgt0
      phi1_i = psi1_i - (validate_treat_bin/pAgt0)*mean(psi1_i)
      
      # estimators of Psi 1, variance, covariance, and plug-in estimator
      delta_psi1hat[j] = mean(psi1_i)
      delta_psi1var[j] = mean(phi1_i^2)
      delta_covvals[j] = mean(phi1_i * phi2_i)
      delta_plugin[j] = mean(muqhat_validate)
      
      # keep influence functions
      phi_i = phi1_i - phi2_i
      # IF_matrix[estimation_set, j] <- phi_i*(N_tot/n)
      IF_matrix[estimation_set, j] <- phi_i
      
    }
    
    # compiling parameter and variance estimates
    psi1_hat_k[[k]] <- delta_psi1hat
    psi1_hat_var_k[[k]] <- delta_psi1var
    covvals[[k]] = delta_covvals
    psi1_plugin_k[[k]] = delta_plugin
    
    # Final estimates
    psi_hat_k[[k]] = psi1_hat_k[[k]] - psi2_hat_k[[k]]
    psi_hat_var_k[[k]] = psi1_hat_var_k[[k]] + psi2_hat_var_k[[k]] - 2*covvals[[k]]
    psi_plugin_k[[k]] = psi1_plugin_k[[k]] - psi2_plugin_k[k]
    
    #Track progress
    print(sprintf('Completed fold %i out of %i; Time: %s', k, K, Sys.time()))
  }
  
  # psihat = as.numeric(t(do.call(rbind, psi_hat_k)) %*% (sapply(folds, length) / n))
  
  return(list('psihat' = colMeans(do.call(rbind, psi_hat_k)), 
              'variance' = colMeans(do.call(rbind, psi_hat_var_k)) / n,
              'psi1_hat' = colMeans(do.call(rbind, psi1_hat_k)),
              'psi2_hat' = mean(psi2_hat_k),
              'psi1_hat_var' = colMeans(do.call(rbind, psi1_hat_var_k)) / n,
              'psi2_hat_var' = mean(psi2_hat_var_k) / n,
              'psi_plugin' = colMeans(do.call(rbind, psi_plugin_k)),
              'psi1_plugin' = colMeans(do.call(rbind, psi1_plugin_k)),
              'psi2_plugin' = mean(psi2_plugin_k),
              'IF_matrix' = IF_matrix))
  
}

# estimator of ASDT for multiple time periods
estimate_psi_mtp <- function(Y_t, Y_g1, treat, G, x, stochastic_policy, 
                             policy_params, nuisance_estimator,
                             g, t, T_periods, K = 5,
                             design_d = seq(0.001, 0.999, length.out = 500)) {
  
  # only oracle/glm/bart supported for nuisance function estimators
  nuisance_estimator <- match.arg(nuisance_estimator, c("oracle", "glm", "bart"))
  
  # Calculate outcome difference (change from pre-treatment period g-1 to t)
  deltay = Y_t - Y_g1
  
  # Restrict sample to the treated cohort 'g' and available controls (not yet treated by time t)
  I_g = as.integer(G == g)
  I_ctrl = as.integer(G > t)
  keep = (I_g == 1) | (I_ctrl == 1)
  
  # entire sample size
  N_tot = length(G)
  
  # Store influence functions for variance estimation
  IF_matrix = matrix(0, nrow = N_tot, ncol = length(delta_vals))
  
  # Subset data to relevant observations
  deltay = deltay[keep]
  treat = treat[keep]
  x = x[keep,,drop = FALSE]
  G_keep = G[keep]
  treat_bin = I_g[keep]
  n = length(deltay)
  
  # estimate asdt (psi) using subsetted data
  estimates = estimate_psi(deltay = deltay,
                           treat = treat,
                           x = x,
                           stochastic_policy = stochastic_policy,
                           policy_params = policy_params,
                           nuisance_estimator = nuisance_estimator,
                           K = K,
                           design_d = design_d, 
                           g = g,
                           t = t,
                           T_periods = T_periods)
  
  # augment influence functions to account for P(G=g| I(G=g)+I(G>t)=1) 
  # being estimated above instead of P(G=g)
  estimates$IF_matrix = (N_tot / n)*estimates$IF_matrix
  
  # pad estimated IF_matrix with 0 for obs where G=g or G>t not satisfied
  IF_matrix[keep,] = estimates$IF_matrix
  
  # return list with main estimates and IF matrix as first layer in the list
  estimates$IF_matrix = NULL
  return(list('main' = estimates,
              'IF_matrix' = IF_matrix))
}


# Aggregation  and wif Functions 
# Referencing from Callaway's did reposistory
# https://github.com/bcallaway11/did/

# Helper function to compute the influence function of the weights
wif <- function(keepers, pg, weights.ind, G, group) {
  if1 <- sapply(keepers, function(k) {
    (weights.ind * 1 * (G == group[k]) - pg[k]) / sum(pg[keepers])
  })
  if2 <- base::rowSums(sapply(keepers, function(k) {
    weights.ind * 1 * (G == group[k]) - pg[k]
  })) %*% t(pg[keepers] / (sum(pg[keepers])^2))
  return(if1 - if2)
}

# Constructs the aggregated influence function to compute standard errors
get_agg_inf_func <- function(att, inffunc1, whichones, weights.agg, wif = NULL) {
  weights.agg <- as.matrix(weights.agg)
  thisinffunc <- inffunc1[, whichones, drop = FALSE] %*% weights.agg
  if (!is.null(wif)) {
    # Adjust for uncertainty in estimating the aggregation weights
    thisinffunc <- thisinffunc + wif %*% as.matrix(att[whichones])
  }
  return(thisinffunc)
}

# Core function to summarize specific cohort-time estimates into interpretable parameters
aggregate_continuous_did <- function(results_df, IF_list, G, T_periods, 
                                     type = "all", nuisance_estimator = 'oracle',
                                     min_e = -Inf, max_e = Inf, balance_e = 1) { 
  
  # subset results_df to specified nuisance estimator
  results_df = results_df %>% filter(nuisance == nuisance_estimator)
  
  # vector of time points
  tlist <- sort(unique(panel_data$periods))
  
  # total sample size before subsetting
  N_tot <- length(G)
  
  # weights (not used; set to 1)
  weights.ind <- rep(1, N_tot)
  
  # vector of cohorts
  originalglist <- 1:T_periods
  
  # Probability of being in cohort g
  pg_overall <- sapply(originalglist, function(g) mean(weights.ind * (G == g)))
  
  # stochastic intervention parameter
  delta_vals <- sort(unique(results_df$delta))
  
  # initialize result list for aggregated estimates
  res_list <- list()
  
  # Loop over each tilt parameter delta to aggregate estimates separately
  for (d_idx in 1:length(delta_vals)) {
    d <- delta_vals[d_idx]
    
    # Filter for post-treatment estimates corresponding to the current delta
    df_d <- results_df[results_df$delta == d & results_df$t >= results_df$g, ]
    if (nrow(df_d) == 0) next
    
    # extract estimates and indices
    att <- df_d$psihat
    group <- df_d$g
    t <- df_d$t
    e <- df_d$e
    
    # Extract influence functions from the list
    inffunc1 <- do.call(cbind, lapply(1:nrow(df_d), function(i) {
      IF_list[[paste(group[i], t[i], sep="_")]][[nuisance_estimator]][, d_idx]
    }))
    
    pg <- pg_overall[group]
    keepers <- which(group <= t & t <= (group + max_e))
    
    # Simple (Overall) Aggregation: Average effect across all post-treatment periods and cohorts
    if (type %in% c("simple", "all")) {
      simple.att <- sum(att[keepers] * pg[keepers]) / sum(pg[keepers])
      simple.wif <- wif(keepers, pg, weights.ind, G, group)
      simple.if <- get_agg_inf_func(att, inffunc1, keepers, pg[keepers] / sum(pg[keepers]), simple.wif)
      simple.se <- sqrt(mean(as.numeric(simple.if)^2) / N_tot)
      
      res_list[[length(res_list) + 1]] <- data.frame(
        type = "simple", target = "Overall", delta = d, estimate = simple.att, se = simple.se
      )
    }
    
    # Group-specific Aggregation: Average effect specific to each adoption cohort 'g'
    if (type %in% c("group", "all")) {
      glist_current <- sort(unique(group))
      selective.att.g <- sapply(glist_current, function(g_val) {
        whichg <- which((group == g_val) & (g_val <= t) & (t <= (g_val + max_e)))
        mean(att[whichg])
      })
      selective.inf.func.g <- do.call(cbind, lapply(glist_current, function(g_val) {
        whichg <- which((group == g_val) & (g_val <= t) & (t <= (g_val + max_e)))
        inf.func.g <- get_agg_inf_func(att, inffunc1, whichg, rep(1/length(whichg), length(whichg)), NULL)
        as.numeric(inf.func.g)
      }))
      pgg_current <- pg_overall[glist_current]
      selective.att <- sum(selective.att.g * pgg_current) / sum(pgg_current)
      selective.wif <- wif(1:length(glist_current), pgg_current, weights.ind, G, glist_current)
      selective.inf.func <- get_agg_inf_func(selective.att.g, selective.inf.func.g, 1:length(glist_current), pgg_current / sum(pgg_current), selective.wif)
      selective.se <- sqrt(mean(as.numeric(selective.inf.func)^2) / N_tot)
      
      res_list[[length(res_list) + 1]] <- data.frame(type = "group", target = "Overall", delta = d, estimate = selective.att, se = selective.se)
      for (i in 1:length(glist_current)) {
        se.g <- sqrt(mean(selective.inf.func.g[, i]^2) / N_tot)
        res_list[[length(res_list) + 1]] <- data.frame(type = "group", target = paste0("g=", glist_current[i]), delta = d, estimate = selective.att.g[i], se = se.g)
      }
    }
    
    # Dynamic Aggregation: Average effect by event time 'e' (Time since treatment)
    if (type %in% c("dynamic", "all")) {
      eseq <- sort(unique(e))
      eseq <- eseq[(eseq >= min_e) & (eseq <= max_e)]
      dynamic.att.e <- sapply(eseq, function(e_val) {
        whiche <- which(e == e_val)
        sum(att[whiche] * (pg[whiche] / sum(pg[whiche])))
      })
      dynamic.inf.func.e <- do.call(cbind, lapply(eseq, function(e_val) {
        whiche <- which(e == e_val)
        pge <- pg[whiche] / sum(pg[whiche])
        wif.e <- wif(whiche, pg, weights.ind, G, group)
        as.numeric(get_agg_inf_func(att, inffunc1, whiche, pge, wif.e))
      }))
      epos <- eseq >= 0
      if (any(epos)) {
        dynamic.att <- mean(dynamic.att.e[epos])
        dynamic.inf.func <- get_agg_inf_func(dynamic.att.e[epos], as.matrix(dynamic.inf.func.e[, epos, drop=FALSE]), 1:sum(epos), rep(1/sum(epos), sum(epos)), NULL)
        dynamic.se <- sqrt(mean(as.numeric(dynamic.inf.func)^2) / N_tot)
        res_list[[length(res_list) + 1]] <- data.frame(type = "dynamic", target = "Overall", delta = d, estimate = dynamic.att, se = dynamic.se)
      }
      for (i in 1:length(eseq)) {
        se.e <- sqrt(mean(dynamic.inf.func.e[, i]^2) / N_tot)
        res_list[[length(res_list) + 1]] <- data.frame(type = "dynamic", target = paste0("e=", eseq[i]), delta = d, estimate = dynamic.att.e[i], se = se.e)
      }
    }
    
    # Balanced Dynamic Aggregation: Dynamic aggregation restricted to cohorts observed up to balance_e periods
    if (type %in% c("dynamic_balanced", "all") && !is.null(balance_e)) {
      eseq_bal <- sort(unique(e))
      include_balanced <- (group + balance_e <= T_periods)
      eseq_bal <- eseq_bal[eseq_bal <= balance_e & eseq_bal >= min_e] 
      if (length(eseq_bal) > 0) {
        dynamic_bal.att.e <- sapply(eseq_bal, function(e_val) {
          whiche <- which((e == e_val) & include_balanced)
          if(length(whiche) == 0) return(NA)
          sum(att[whiche] * (pg[whiche] / sum(pg[whiche]))) 
        })
        dynamic_bal.inf.func.e <- do.call(cbind, lapply(eseq_bal, function(e_val) {
          whiche <- which((e == e_val) & include_balanced)
          if(length(whiche) == 0) return(rep(0, N_tot))
          pge <- pg[whiche] / sum(pg[whiche])
          wif.e <- wif(whiche, pg, weights.ind, G, group)
          as.numeric(get_agg_inf_func(att, inffunc1, whiche, pge, wif.e))
        }))
        epos_bal <- eseq_bal >= 0
        valid_idx <- !is.na(dynamic_bal.att.e)
        epos_valid <- epos_bal & valid_idx
        if (any(epos_valid)) {
          dynamic_bal.att <- mean(dynamic_bal.att.e[epos_valid])
          dynamic_bal.inf.func <- get_agg_inf_func(dynamic_bal.att.e[epos_valid], as.matrix(dynamic_bal.inf.func.e[, epos_valid, drop=FALSE]), 1:sum(epos_valid), rep(1/sum(epos_valid), sum(epos_valid)), NULL)
          dynamic_bal.se <- sqrt(mean(as.numeric(dynamic_bal.inf.func)^2) / N_tot)
          res_list[[length(res_list) + 1]] <- data.frame(type = "dynamic_balanced", target = "Overall", delta = d, estimate = dynamic_bal.att, se = dynamic_bal.se)
        }
        for (i in 1:length(eseq_bal)) {
          if (!valid_idx[i]) next
          se.e_bal <- sqrt(mean(dynamic_bal.inf.func.e[, i]^2) / N_tot)
          res_list[[length(res_list) + 1]] <- data.frame(type = "dynamic_balanced", target = paste0("e=", eseq_bal[i]), delta = d, estimate = dynamic_bal.att.e[i], se = se.e_bal)
        }
      }
    }
  } 
  
  # Bind all delta results and calculate standard 95% Wald-like Confidence Intervals
  final_res <- do.call(rbind, res_list)
  final_res$ci_lower <- final_res$estimate - qnorm(0.975) * final_res$se
  final_res$ci_upper <- final_res$estimate + qnorm(0.975) * final_res$se
  return(final_res)
}


# Stochastic policies -----------------------------------------------------

# Exponential tilting function to shift the dose distribution by increment parameter delta
# optionally divides by pi i.e., returns q(D|X) / pi(D|X)
exp_tilt <- function(params, pihats, treat_points, design_points, divide_by_pi=FALSE) {
  
  # =========================================
  # params: list of parameters with names:
  #   - delta: scalar, increment parameter
  # pihats: matrix of pi(d|X) where rows correspond to observations and cols correspond to design points
  #   - nrow = n of training set, ncol = number of design points
  # treat_points: vector of observed doses D
  # design_points: vector of design points for Monte Carlo integration
  # =========================================
  
  # unpack parameters
  delta = params[['delta']]
  
  # estimating denominator of q(D|X)
  # i.e. integral of pi(d|X) exp(delta d)
  # using Monte Carlo integration with uniform sampling of design points
  q_denom = rowMeans(pihats %*% diag(exp(delta * design_points)))
  
  # estimate numerator of q(D|X) or q(D|X) / pi(D|X)
  if (divide_by_pi) {
    q_numer = exp(delta * treat_points)
  } else {
    q_numer = pihats %*% diag(exp(delta * design_points))
  }
  
  return(q_numer / q_denom)
  
}


# Gaussian concentration to shift dose distribution around specified point d'
# optionally divides by pi i.e., returns q(D|X) / pi(D|X)
gauss_conc <- function(params, pihats, treat_points, design_points, divide_by_pi=FALSE) {
  
  # =========================================
  # params: list of parameters with names:
  #   - dprime: scalar between 0 and 1, concentration point
  #   - l: scalar length parameter, l > 0, 1/l^2 is concentration
  # pihats: matrix of pi(d|X) where rows correspond to observations and cols correspond to design points
  #   - nrow = n of training set, ncol = number of design points
  # treat_points: vector of observed doses D
  # design_points: vector of design points for Monte Carlo integration
  # =========================================
  
  # unpack parameters
  dprime = params[['dprime']]
  l = params[['l']]
  
  # estimating denominator of q(D|X)
  # i.e. integral of pi(d|X) exp(delta d)
  # using Monte Carlo integration with uniform sampling of design points
  q_denom = rowMeans(pihats %*% diag(exp(-(design_points - dprime)^2 / (2*l^2))))
  
  # estimate numerator of q(D|X) or q(D|X) / pi(D|X)
  if (divide_by_pi) {
    q_numer = exp(-(treat_points - dprime)^2 / (2*l^2))
  } else {
    q_numer = pihats %*% diag(exp(-(design_points - dprime)^2 / (2*l^2)))
  }
  
  return(q_numer / q_denom)
  
}


# Superlearner wrappers ---------------------------------------------------

SL.myhal9001 <- function(...) (SL.hal9001(..., max_degree=3, num_knots=c(25, 10, 5)))

SL.myBART <- function(Y, X, newX, family, obsWeights, id, ...) {
  
  # set type of bart depending on family
  if (family$family=='gaussian') {
    type='wbart'
    ndpost=500
    nskip=500
    keepevery=1
  }
  else if (family$family=='binomial') {
    type='pbart'
    ndpost=500
    nskip=500
    keepevery=10
  }
  
  # make sure X is a matrix
  if (!is.matrix(X)) (X = as.matrix(X))
  
  # model
  model = BART::gbart(x.train=X, y.train=Y, type=type, 
                      ndpost=ndpost, nskip=nskip, 
                      keepevery=keepevery,
                      ...)
  
  # predictions
  if (type=='wbart') {
    pred = colMeans(predict(model, newdata=newX))
  }
  else if (type=='pbart') {
    pred = predict(model, newdata=newX)$prob.test.mean
  }
  
  # return model and predictions
  fit = list(object = model)
  class(fit) = c('SL.myBART')
  out = list(pred = pred, fit = fit)
  return(out)
  
}
predict.SL.myBART <- function(object, newdata, ...) {
  if (!is.matrix(newdata)) (newdata = as.matrix(newdata))
  pred = predict(object = object$object, newdata = newdata)
  if (is.matrix(pred)) (pred = colMeans(pred))
  else (pred = pred$prob.test.mean)
  pred
}



# Misc --------------------------------------------------------------------

# takes multiple time periods dataset and specified
# cohort g and time t; outputs a wide dataset for
# with data relevant to times g and t
make_gt_cross_section <- function(panel_data, g, t) {
  
  # reference year for outcome
  base_t <- g - 1
  
  # subset data to reference time and time t
  df <- panel_data %>%
    filter(periods %in% c(t, base_t)) %>%
    select(id, G, periods, Y, treat, starts_with("X"))
  
  # convert to wide dataset
  wide <- df %>%
    pivot_wider(
      id_cols=c(id,G,starts_with("X")),
      names_from=periods,
      values_from=c(Y,treat),
      names_sep=""
    ) %>%
    arrange(id)
  
  y_t_col <- paste0("Y",t)
  y_g1_col <- paste0("Y",base_t)
  treat_col <- paste0("treat",t)
  
  list(
    Y_t  = wide[[y_t_col]],
    Y_g1 = wide[[y_g1_col]],
    treat= wide[[treat_col]],
    x    = wide %>% select(starts_with("X")) %>% as.matrix(),
    G    = wide$G,
    ids  = wide$id
  )
}


# used in sims_eval.R to format results of simulation.R
format_results <- function(results, sim_num) {
  
  # names of list tell what was used to estimate nuisance functions
  nuisance_names = names(results)
  
  # loop through list and format results
  formatted_results = list()
  for (i in 1:length(results)) {
    
    # extract estimates
    estimates = results[[i]]
    
    # format into a table
    results.df = data.frame(delta = delta_vals,
                            psihat = estimates[[1]],
                            psihat_var = estimates[[2]],
                            psi1_hat = estimates[[3]],
                            psi2_hat = estimates[[4]],
                            psi1_hat_var = estimates[[5]],
                            psi2_hat_var = estimates[[6]],
                            ci_lower = estimates[[1]] - qnorm(0.975)*sqrt(estimates[[2]]),
                            ci_upper = estimates[[1]] + qnorm(0.975)*sqrt(estimates[[2]]),
                            nuisance = nuisance_names[i],
                            sim = sim_num)
    
    # store
    formatted_results[[i]] = results.df
    
  }
  
  return(bind_rows(formatted_results))
  
}

# format results for multiple time periods
format_results_mtp <- function(results, sim_num) {
  
  # results is a list with names that specify the
  # cohort g and time t; each element of the list
  # is another list that is the same format as 
  # passed into format_results(...)
  
  results_list = list()
  
  for (gt in 1:length(results)) {
    
    # extract g and t from results
    gt_name = names(results)[gt] %>% str_split_1('_')
    g = as.integer(gt_name[1])
    t = as.integer(gt_name[2])
    
    # format results
    formatted.results = format_results(results[[gt]], sim_num)
    
    # add in g, t, e
    formatted.results = formatted.results %>%
      mutate(g = g, t = t, e = t-g)
    
    # store
    results_list[[gt]] = formatted.results
    
  }
  
  return(bind_rows(results_list))
  
}
