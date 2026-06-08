# defines parameters for different simulation scenarios
# multiple time periods simulations

# Scenario 1
if (scenario == 1) {
  
  # sample size
  n = 5000
  
  # number of time periods
  T_periods <- 3
  
  # covariate predictor for cohort assignment dgp
  cohort_predictor <- function(x) {
    if (is.null(dim(x))) x <- matrix(x, nrow = 1)
    3 * cos(x[,1]^2) - 0.5 * exp(x[,2]) - 1.5 * x[,3]^2 + 4 * (x[,4]<0) * abs(x[,5])
  }
  
  # parameters for conditional dose distribution
  shape1_alpha_f <- function(x, g=1) { 
    if (is.null(dim(x))) x <- matrix(x, nrow = 1)
    exp(as.numeric(x %*% seq(-2, 2, length.out=ncol(x)))/ncol(x)) 
  }
  
  shape2_beta_f  <- function(x, g=1) { 
    if (is.null(dim(x))) x <- matrix(x, nrow = 1)
    exp(as.numeric(x %*% seq(-2, 2, length.out=ncol(x)))/ncol(x)) 
  }
  
  # generalized dose propensity score pi(d|x)
  true_density_f <- function(d, x, shape1_alpha_f, shape2_beta_f) {
    dbeta(d, shape1 = shape1_alpha_f(x), shape2 = shape2_beta_f(x))
  }
  
  # time fixed effects (set to zero)
  make_theta_t <- function(T_periods, slope = 0) {
    slope * (0:T_periods)
  }
  theta_t <- make_theta_t(T_periods = T_periods, slope = 0)
  
  # function to generate coefficients for nonlinear dose effect
  c_fun <- function(e, g) { c(0.2 + 0.1*e,  0.5,  0.1) }
  
  # outcome regression for A=0
  mu.0_f <- function(x, g=1, t=1, theta_t=c(0,0,0,0)) {
    if (is.null(dim(x))) x <- matrix(x, nrow = 1)
    j_t <- t + 1; j_g1 <- (g - 1) + 1
    
    the_diff <- theta_t[j_t] - theta_t[j_g1]
    time_diff <- (t - g + 1)
    
    # covariate predictor for outcome regression dgp
    trend_slope_f <- function(x) {
      if (is.null(dim(x))) x <- matrix(x, nrow = 1)
      plogis(x[,1])^2 * x[,2]*x[,3] + x[,4] - 2*x[,4]^2 + 3*sin(x[,5]^2)
    }
    cov_diff <- time_diff * as.numeric(trend_slope_f(x))
    
    return(the_diff + cov_diff)
  }
  
  # outcome regression for A>0
  mu.d_f <- function(x, d, effect_coefs, g=1, t=1, theta_t=c(0,0,0,0)) {
    if (is.null(dim(x))) x <- matrix(x, nrow = 1)
    e <- t - g
    
    dose_effect <- d*effect_coefs[1] + (d^2)*effect_coefs[2] + (d^3)*effect_coefs[3]
    
    return(mu.0_f(x, g, t, theta_t) + (t >= g) * as.numeric(dose_effect))
  }
}

# Scenario 2
if (scenario == 2) {
  
  # sample size
  n = 5000
  
  # number of time periods
  T_periods <- 3
  
  # covariate predictor for cohort assignment dgp
  cohort_predictor <- function(x) {
    if (is.null(dim(x))) x <- matrix(x, nrow = 1)
    3 * cos(x[,1]^2) - 0.5 * exp(x[,2]) - 1.5 * x[,3]^2 + 4 * (x[,4]<0) * abs(x[,5])
  }
  
  # parameters for conditional dose distribution
  shape1_alpha_f <- function(x, g=1) { 
    if (is.null(dim(x))) x <- matrix(x, nrow = 1)
    (2) 
  }
  
  shape2_beta_f  <- function(x, g=1) { 
    if (is.null(dim(x))) x <- matrix(x, nrow = 1)
    (1 + 4*plogis(x %*% seq(0, 1, length.out=ncol(x)))) 
  }
  
  # generalized dose propensity score pi(d|x)
  true_density_f <- function(d, x, shape1_alpha_f, shape2_beta_f) {
    dbeta(d, shape1 = shape1_alpha_f(x), shape2 = shape2_beta_f(x))
  }
  
  # time fixed effects (set to zero)
  make_theta_t <- function(T_periods, slope = 0) {
    slope * (0:T_periods)
  }
  theta_t <- make_theta_t(T_periods = T_periods, slope = 0)
  
  # function to generate coefficients for nonlinear dose effect
  c_fun <- function(e, g) { c(0.2 + 0.1*e,  0.5,  0.1) }
  
  # outcome regression for A=0
  mu.0_f <- function(x, g=1, t=1, theta_t=c(0,0,0,0)) {
    if (is.null(dim(x))) x <- matrix(x, nrow = 1)
    j_t <- t + 1; j_g1 <- (g - 1) + 1
    
    the_diff <- theta_t[j_t] - theta_t[j_g1]
    time_diff <- (t - g + 1)
    
    # covariate predictor for outcome regression dgp
    trend_slope_f <- function(x) {
      if (is.null(dim(x))) x <- matrix(x, nrow = 1)
      plogis(x[,1])^2 * x[,2]*x[,3] + x[,4] - 2*x[,4]^2 + 3*sin(x[,5]^2)
    }
    cov_diff <- time_diff * as.numeric(trend_slope_f(x))
    
    return(the_diff + cov_diff)
  }
  
  # outcome regression for A>0
  mu.d_f <- function(x, d, effect_coefs, g=1, t=1, theta_t=c(0,0,0,0)) {
    if (is.null(dim(x))) x <- matrix(x, nrow = 1)
    e <- t - g
    
    dose_effect <- d*effect_coefs[1] + (d^2)*effect_coefs[2] + (d^3)*effect_coefs[3]
    
    return(mu.0_f(x, g, t, theta_t) + (t >= g) * as.numeric(dose_effect))
  }
}
