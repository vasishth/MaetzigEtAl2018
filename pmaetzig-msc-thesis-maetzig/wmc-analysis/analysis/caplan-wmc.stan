data {
  int<lower=1> n_obs;
  vector[n_obs] acc;
  vector<lower=-1, upper=1>[n_obs] condition;
}

parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}

model {
  alpha ~ uniform(0, 100);
  beta ~ uniform(0, 100);
  sigma ~ uniform(0, 100);
  
  for (i in 1:n_obs) {
    acc[i] ~ logit(alpha + beta * condition[i], sigma);
  }
}