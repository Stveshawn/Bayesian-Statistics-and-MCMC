data{
  int<lower=1> n; 
  vector[n] Y;
  vector[n] X;
  real<lower=0> sigma;
  real<lower=0> alpha;
  real<lower=0> beta;
  real<lower=0> lambda; // Modified
  real<lower=0> alpha_d; // Modified
  real<lower=0> beta_d; // Modified
}

parameters{
  real<lower=0> a;
  real<lower=0> b;
  real<lower=0> c; // Modified
  real<lower=0> d; // Modified
  vector<lower=0, upper=1>[n] theta;
  vector<lower=0, upper=1>[n] Pi;  
  // Additional independent chain (pi,theta) 
  // to estimate integrals.
  
}

model{
  
  a ~ gamma(alpha,beta); // target += gamma_lpdf(a | alpha,beta);
  b ~ gamma(alpha,beta); // target += gamma_lpdf(b | alpha,beta);
  c ~ exponential(lambda); // Modified
  d ~ gamma(alpha_d,beta_d); // Modified
  
  for(i in 1:n){
    
    theta[i] ~ beta(a,b); // target += beta_lpdf(theta[i] | a,b);
    
    Pi[i] ~ beta(1 + c*theta[i], 1 + c*(1-theta[i])); // Modified
    
    X[i] ~ normal(logit(theta[i]), sigma);
    // target += normal_lpdf(logit(theta[i]), sigma);
    
    target += log_sum_exp(log(Pi[i]) + beta_lpdf(Y[i] | 1+d, 1), log(1-Pi[i]) + beta_lpdf(Y[i] | 1, 1+d)); // Modified
  }
  
}

generated quantities{
  real<lower=0, upper=1> Pi_s;
  real<lower=0, upper=1> theta_s;
  
  theta_s = beta_rng(a,b);
  Pi_s = beta_rng(1 + c*theta_s, 1 + c*(1-theta_s)); // Modified
}


