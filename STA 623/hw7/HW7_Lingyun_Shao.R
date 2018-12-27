###################
# MCMC HW7        #
###################

# Model
# -----
#
# Modified : added c and d
# theta ~ beta(a,b)                                    : Level of chemical M
# X | theta ~ normal(logit(theta), sigma)              : Test X, measuring level of M
# pi | c, theta ~ Beta(1 + c*theta, 1 + c*(1 - theta)) : Probability that patient has sickness, high levels of M tend to give higher probability
# Sick | pi ~ Bernoulli(pi)                            : Patient having the sickness
# Y | d, Sick ~ Beta(1 + d*Sick, 1 + d(1 - Sick))      : Test Y, sick patients score tend to be much higher
# a, b ~ Gamma(2,1)                                    : Population parameters of levels of M
# c ~ Exp(1)
# d ~ Gamma(2,2)

library(rstan)


# Prior and known parameters
sigma <- 0.5
alpha <- 2
beta <- 1
lambda = 1 # Modified
alpha_d = 2 # Modified
beta_d = 2 # Modified

# Create data
set.seed(2112)
n <- 100

# Latent theta, probability of illness, and illness indicator.
th0 <- rbeta(n,2,2)
# Modified : The original code is commented, I'm not sure whether this code is wrong or not
# But judging from the formula in Lab 2 notes, I think the parameters should not be multiplied by 2
# And I also set 'c' here to be its expectation in a exponential distribution, that is 1.
# pi0 <- rbeta(n,2*(1 + 2*th0), 2*(1 + 2*(1-th0)))
pi0 = rbeta(n,(1 + 1*th0), (1 + 1*(1-th0)))
S0 <- rbinom(n,1,pi0)

# Test results
X <- rnorm(n, log(th0)-log(1-th0), sigma)
Y <- rbeta(n, 1 + S0, 1 + (1 - S0))

# A plot. (remember we do not know the health of the data, only X,Y) 
plot(X,Y,col=S0+1)

# Data format for Stan.
stan_data <- list(n = n, Y = Y, X = X, sigma = sigma, alpha = alpha, beta = beta,
                  lambda = lambda, alpha_d = alpha_d, beta_d = beta_d) 

# Execute this line and find the file (only once)
if(!exists('stan_file')){ stan_file <- file.choose() } # Modified : added quote

T <- 2000
B <- 500
# This is the Stan execution, may take a while.
fit <- stan(file = stan_file, data = stan_data, 
            iter = B+T, warmup = B, chains=1) # Modified : fixed warmup to B

# Extract MCMC chains, we will need pi and theta.
draws <- extract(fit, pars = c("a","b",'c', 'd', "Pi_s","theta_s"))  # Modified


# Some trace plots
plot(draws$a, type="l", ylab="a")
plot(draws$b, type="l", ylab="b")
plot(draws$Pi_s, type="l", ylab="pi")
plot(draws$theta_s, type="l", ylab="theta") # Modified : fixed label to theta
plot(draws$c, type="l", ylab="c") # Modified
plot(draws$d, type="l", ylab="d") # Modified

## EXERCISE 3: here

# Since we are using samples from stan, we don't need to sample in R


##

# Short names for mcmc chains. Erase these two lines if using EXERCISE 3 code.
theta <- draws$theta_s
pi <- draws$Pi_s
c <- draws$c # Modified
d <- draws$d # Modified

# Logit of theta's
logit_theta <- log(theta)-log(1-theta)

# Define the Monte Carlo approximation of gamma(x,y)
gamma <- function(x,y){
  
  numerator <- dbeta(y,1+d,1) * pi * dnorm(x, logit_theta, sigma) # Modified
  denominator <- (dbeta(y,1+d,1)* pi + dbeta(y,1,1+d) * (1-pi)) * dnorm(x, logit_theta, sigma) # Modified
  
  # Probability value
  mean(numerator)/mean(denominator)
  
}

# Now let's set a range of interest for both test values.
seq_x <- seq(-3,3,0.05)
seq_y <- seq(0.01,0.99,0.01)
gamma_matrix <- array(0,c(length(seq_x),length(seq_y)))

# Define a loss matrix
Loss <- matrix(c(0,1,4,0),2,2, byrow=TRUE, dimnames=list(c("Treat","!Treat"),c("Sick","!Sick")))

# Decision rule matrix
Delta <- array(0,c(length(seq_x),length(seq_y)))

# This loop evaluates the probability that the patient is sick,
# given the test values x, y. The idea is to integrate out both
# pi and theta.
for(s in 1:length(seq_x)){ for(t in 1:length(seq_y)){
  
  x <- seq_x[s]
  y <- seq_y[t]
  
  # Evaluate probability
  pr <- gamma(x,y)
  
  # Save for probability table
  gamma_matrix[s,t] <- pr
  
  # Bayesian expected losses
  B0 <- pr*Loss["!Treat","Sick"] + (1-pr)*Loss["!Treat","!Sick"]
  B1 <- pr*Loss["Treat","Sick"]  + (1-pr)*Loss["Treat","!Sick"]
  
  # Optimal decision for this particular (x,y) combination.
  Delta[s,t] <- 1*(B1 < B0)
  
}}

# Equivalent threshold for probabilities.
thr <- (Loss[2,2] - Loss[1,2])/( (Loss[2,2] - Loss[1,2]) + (Loss[1,1] - Loss[2,1]))

# Color plotting of probabilities. 
image(seq_x,seq_y,gamma_matrix, xlab= "Test X", ylab="Test Y")
contour(seq_x,seq_y,gamma_matrix, levels = thr,add=TRUE, lwd=2,col=4, labcex=1)

# Contour plot gives numerical values.
contour(seq_x,seq_y,gamma_matrix, xlab= "Test X", ylab="Test Y")

# Decision rule (only 1 or 0)
image(seq_x,seq_y,Delta, xlab= "Test X", ylab="Test Y")

