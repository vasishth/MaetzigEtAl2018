data {
  int<lower=1> N_obs;
  int<lower=1> N_subj;
//  int<lower=1> N_item;
  
  int<lower=0, upper=1> acc[N_obs];
  vector<lower=-1, upper=1>[N_obs] x;
}

parameters {
  real alpha;
  real beta;
//  real<lower=0> sigma;
//  real<lower=0> tau_u;
//  real<lower=0> tau_v;
// 
//  vector[N_subj] u;
//  vector[N_item] v;
}

model {
  // PRIORS HERE
//  sigma ~ normal(0, 1);
//  tau_u ~ normal(0, 1)
//  tau_v ~ normal(0, 1)

for (n in 1:N_obs) {
  acc[n] ~ bernoulli_logit(alpha + x[n] * beta);
}
}

generated quantities {
  real<lower=0> odds[N_;
  odds = exp(
}
  
//  for (i in 1:N_subj) {
//    u[i] ~ normal(0, tau_u)
//  }
//  
//  for (j in 1:N_item) {
//    v[j] ~ normal(0, tau_v)
//  }
  
//  acc ~ bernoulli_logit(alpha + beta * x);
