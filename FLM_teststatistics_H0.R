# Load required libraries
library(ggplot2)

# Function to generate a single path of Brownian motion
generate_brownian_motion <- function(dt) {
  t <- seq(0, 1, by = dt)
  dW <- rnorm(length(t) - 1, mean = 0, sd = sqrt(dt))
  W <- c(0, cumsum(dW))
  return(W)
}

# Function to calculate the integral using the trapezoidal method 
trapezoidal_integral <- function(f, dt) {
  integral <- sum((f[-1] + f[-length(f)]) * dt / 2)
  return(integral)
}

# Function to compute the L^2 integral norm
l2_norm <- function(f, dt) {
  return(sqrt(trapezoidal_integral(f^2,dt)))
}

# Function to compute the empirical cross-covariance operator
empirical_cross_cov_operator <- function(X, Y, n, dt) {
  t <- seq(0,1, by = dt)
  
  delta_n <- numeric(length(t))
  for(i in 1:n) {
    delta_n <- delta_n + (1/n * (X[,i] * Y[i]))
  }
  return(delta_n)
}

# Function to compute the operator A_n 
operator_A_n <- function(V, x, p_n, lambda, dt) {
  t <- seq(0,1, by = dt)
  A_n <- numeric(length(t))
  
  for(i in 1:p_n) {
    V_i_x <- V[,i] * x
    scalar_product <- trapezoidal_integral(V_i_x, dt)
    A_n <- A_n + lambda[i]^{-1/2} * scalar_product * V[,i]
  }
  return(A_n)
}

# Function to generate a matrix of functions V_j
generate_V_functions <- function(p_n, dt) {
  t <- seq(0, 1, by = dt)
  V <- matrix(0, nrow = length(t), ncol = p_n)
  
  for (j in 1:p_n) {
    V[,j] <-  sqrt(2) * sin((j - 1/2) * pi * t)
  }
  
  return(V)
}

# Function to generate a vector lambda_inv of size n
generate_lambda <- function(p_n) {
  
  lambda <- numeric(p_n)
  
  for (j in 1:p_n) {
    lambda[j] <- 1/((j - 1/2) * pi)^2
  }
  
  return(lambda)
}


n <- 200        # Number of path         
dt <- 0.005       # Time step
t <- seq(0, 1, by = dt)

p_n <- 50             #p_n value 
phi <- numeric(length(t))
phi <- 0 * t             #regression function


run_experiment <- function(seed) {
  set.seed(seed)
  
  
  # Generate X data 
  X <- replicate(n, generate_brownian_motion(dt))
  
  sigma <- 1 
  # Transform X into output data Y
  Y <- numeric(n)
  for (i in 1:n) {
    epsilon <- rt(1, df = 3) * sqrt(1/3) * sigma
    Y[i] <- trapezoidal_integral(X[,i] * phi, dt) + epsilon
  }
  
  #Eigenvalues and Eigenfunctions of X
  lambda <- generate_lambda(p_n)
  V <- generate_V_functions(p_n, dt)
  
  
  # Compute empirical cross-covariance operator
  delta_n <- empirical_cross_cov_operator(X, Y, n,  dt)
  
  
  # Compute A_n(delta_n)
  A_n_delta_n <- operator_A_n(V, delta_n, p_n, lambda, dt)
  
  
  # Compute the norm of sqrt(n) * A_n(delta_n)
  l2_norm_A_n_delta_n <-  l2_norm(sqrt(n) * A_n_delta_n, dt)
  
  
  # Compute the test statistics  
  
  
  D_n <- (1 /sigma^2) * l2_norm_A_n_delta_n^2 
  T_n <- 1/(sqrt(2 * p_n)) * (D_n - p_n)
  
  
  return(list(value1 = D_n, value2 = T_n))
}

m <- 800   #number of experiment replications
seeds <- 1:m
level <- 0.95

D_n_list <- numeric(m) 
T_n_list <- numeric(m)

D_n_hyp <- numeric(m)
T_n_hyp <- numeric(m)
chi_squared_quantile <- qchisq(level, p_n)
gaussian_quantile <- qnorm(level)


for(j in 1:m) {
  results <- run_experiment(seeds[j])
  D_n_list[j] <- results$value1
  T_n_list[j] <- results$value2
  D_n_hyp[j] <- ifelse(D_n_list[j] < chi_squared_quantile, "accept", "reject")
  T_n_hyp[j] <- ifelse(T_n_list[j] < gaussian_quantile, "accept", "reject")
  
}

results_data <- data.frame( 
  D_n = D_n_list,
  Hypothesis = D_n_hyp,
  T_n = T_n_list,
  Hypothesis = T_n_hyp
)


print(results_data)
D_n_accept <- sum(D_n_hyp == "accept")
D_n_reject <- m - D_n_accept
T_n_accept <- sum(T_n_hyp == "accept")
T_n_reject <- m - T_n_accept
cat("Number of accept for D_n:", D_n_accept," Number of accept for T_n:", T_n_accept, "\n")
cat("Number of reject for D_n:", D_n_reject," Number of reject for T_n:", T_n_reject, "\n")
cat(level,"-quantile of chi-squared distribution with", p_n, "degrees of freedom:", chi_squared_quantile, "\n")
cat(level,"-quantile of standard gaussian distribution:", gaussian_quantile, "\n")

level_D_n <- mean(D_n_list > chi_squared_quantile)
level_T_n <- mean(T_n_list > gaussian_quantile)

cat("empirical level of D_n under", m, "Observations:", level_D_n, "\n")
cat("empirical level of T_n under", m, "Observations:", level_T_n, "\n")

# test1 <- runif(10, min = -1, max = 1)
# test2 <- rnorm(10, sd = 1)
# print(test1)
# print(mean(test1))
# print(test2)
# print(mean(test2))