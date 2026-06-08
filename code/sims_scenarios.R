# defines parameters for different simulation scenarios

if (scenario == 1) {
  
  # sample size
  n = 5000
  
  # binary propensity score pi_{A>0}(X)
  pi_bin = function(x) {
    0.5*plogis(-0.5 + 2*cos(x[,1]^2) - 0.2*exp(x[,2]) - 0.5*x[,3]^2) + 0.3*(x[,4]<0)*abs(x[,5])
  }
  
  # generalized dose propensity score pi(d|x)
  shape1_alpha_f = function(x) exp(x %*% seq(-2, 2, length.out=ncol(x))/ncol(x))
  shape2_beta_f = function(x) exp(x %*% seq(-2, 2, length.out=ncol(x))/ncol(x))
  true_density_f <- function(d, x, shape1_alpha_f, shape2_beta_f) {
    dbeta(d, shape1 = shape1_alpha_f(x), shape2 = shape2_beta_f(x))
  }
  
  # outcome regression E[Delta Y | X, A]
  alpha0 = 1
  alphad = 1
  c_fun <- function(e, g) {c(2,-1, 0)}  # effect coefficients (structured this way to align w/ MTP)
  mu.0_f = function(x, g=1, t=1) {
    alpha0 + plogis(x[,1])^2 * x[,2]*x[,3] + x[,4] - 2*x[,4]^2 + 3*sin(x[,5]^2)
  }
  mu.d_f = function(x, d, effect_coefs, g=1, t=1) {
    alphad + d*effect_coefs[1] + (d^2)*effect_coefs[2] + (d^3)*effect_coefs[3] + 
      plogis(x[,1])^2 * x[,2]*x[,3] + x[,4] - 2*x[,4]^2 + 3*sin(x[,5]^2)
  }
  
}

if (scenario == 2) {
  
  # sample size
  n = 5000
  
  # binary propensity score pi_{A>0}(X)
  pi_bin = function(x) {
    0.5*plogis(-0.5 + 2*cos(x[,1]^2) - 0.2*exp(x[,2]) - 0.5*x[,3]^2) + 0.3*(x[,4]<0)*abs(x[,5])
  }
  
  # generalized dose propensity score pi(d|x)
  shape1_alpha_f = function(x) (2)
  shape2_beta_f = function(x) (1 + 4*plogis(x %*% seq(0, 1, length.out=ncol(x))))
  true_density_f <- function(d, x, shape1_alpha_f, shape2_beta_f) {
    dbeta(d, shape1 = shape1_alpha_f(x), shape2 = shape2_beta_f(x))
  }
  
  # outcome regression E[Delta Y | X, A]
  alpha0 = 1
  alphad = 1
  c_fun <- function(e, g) {c(2,-1, 0)}  # effect coefficients (structured this way to align w/ MTP)
  mu.0_f = function(x, g=1, t=1) {
    alpha0 + plogis(x[,1])^2 * x[,2]*x[,3] + x[,4] - 2*x[,4]^2 + 3*sin(x[,5]^2)
  }
  mu.d_f = function(x, d, effect_coefs, g=1, t=1) {
    alphad + d*effect_coefs[1] + (d^2)*effect_coefs[2] + (d^3)*effect_coefs[3] + 
      plogis(x[,1])^2 * x[,2]*x[,3] + x[,4] - 2*x[,4]^2 + 3*sin(x[,5]^2)
  }
  
}