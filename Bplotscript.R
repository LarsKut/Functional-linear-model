# Load required libraries
library(ggplot2)

# Function to generate a single path of Brownian motion
generate_brownian_motion <- function(T, dt) {
  t <- seq(0, T, by = dt)
  dW <- rnorm(length(t) - 1, mean = 0, sd = sqrt(dt))
  W <- c(0, cumsum(dW))
  return(W)
}

n <- 5
dt <- 0.005
T <- 1
t <- seq(0,T, by = dt)

paths <- replicate(n, generate_brownian_motion(T,dt))


#Prepare data for plotting
path_data <- data.frame(
    Time = rep(t, n),
   Value = as.vector(paths),
   Pfade = factor(rep(1:n, each = length(t)))
 )
 #
 # #Plot the generated Brownian motion paths
ggplot(path_data, aes(x = Time, y = Value, group = Pfade, color = Pfade)) +
geom_line(alpha = 0.8) +
theme_minimal() +
labs(x = "Zeit", y = "Werte")