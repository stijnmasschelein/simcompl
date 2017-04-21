data{
  int<lower=1> N;
  vector[N] y;
  vector[N] x1;
  vector[N] x2;
  vector[N] z;
}
transformed data{
  vector[N] x1_2;
  vector[N] x2_2;
  vector[N] x12;
  x1_2 = x1 .* x1;
  x2_2 = x2 .* x2;
  x12 = x1 .* x2;
}
parameters{
  real beta[4];
  real gamma[2];
  real<lower=0> delta[2];
  real<lower=0> sd_nu;
  real<lower=0> sd_eta[2];
}
transformed parameters{
  real<lower=0> denom;
  denom = delta[1] * delta[2] - beta[4]^2;
}
model{
  //priors
  beta ~ normal(0, 1);
  gamma ~ normal(0, 1);
  delta ~ normal(0, 1);
  sd_nu ~ normal(0, 1);
  sd_eta ~ normal(0, 1);

  { vector[N] mu[3];
    vector[N] eta[2];
    //demand
    mu[1] = (delta[2] * beta[2] + beta[4] * beta[3])
              + (delta[2] * gamma[1] + beta[4] * gamma[2]) * z;
    mu[2] = (delta[1] * beta[3] + beta[4] * beta[2])
              + (delta[1] * gamma[2] + beta[4] * gamma[1]) * z;

    x1 ~ normal(mu[1]/denom, sd_eta[1]/denom);
    eta[1] = x1 * denom - mu[1];
    x2 ~ normal(mu[2]/denom, sd_eta[2]/denom);
    eta[2] = x2 * denom - mu[2];
    //production
    mu[3] = beta[1] + (beta[2] + gamma[1] * z + eta[1]) .* x1
                    + (beta[3] + gamma[2] * z + eta[2]) .* x2
            + beta[4] * (x12)
            - .5 * (delta[1] * (x1_2) + delta[2] * (x2_2));
    y ~ normal(mu[3], sd_nu);
  }
}
