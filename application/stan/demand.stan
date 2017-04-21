data{
  int<lower=1> N;
  vector[N] x1;
  vector[N] x2;
  vector[N] z;
}
parameters{
  real beta[4];
  real gamma[2];
  real<lower=0> delta[2];
  real<lower=0> sd_nu;
  real<lower=0> sd_eta[2];
}
transformed parameters{
  //real<lower=0> denom;
  real denom;
  denom = delta[1] * delta[2] - beta[4]^2;
}
model{
  //priors
  beta ~ normal(0, 1);
  gamma ~ normal(0, 1);
  delta ~ normal(0, 1);
  sd_nu ~ normal(0, 1);
  sd_eta ~ normal(0, 1);

  { vector[N] mu[2];
    //demand
    mu[1] = (delta[2] * beta[2] + beta[4] * beta[3])
              + (delta[2] * gamma[1] + beta[4] * gamma[2]) * z;
    mu[2] = (delta[1] * beta[3] + beta[4] * beta[2])
              + (delta[1] * gamma[2] + beta[4] * gamma[1]) * z;

    x1 ~ normal(mu[1]/denom, sd_eta[1]/denom);
    x2 ~ normal(mu[2]/denom, sd_eta[2]/denom);
  }
}
