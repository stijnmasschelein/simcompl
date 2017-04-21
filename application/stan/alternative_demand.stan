data{
  int<lower=1> N;
  vector[N] y;
  vector[N] x1;
  vector[N] x2;
  vector[N] z;
}
transformed data{
  vector[2] zero2;
  vector[N] x1_2;
  vector[N] x2_2;
  vector[N] x12;
  zero2 = rep_vector(0, 2);
  x1_2 = x1 .* x1;
  x2_2 = x2 .* x2;
  x12 = x1 .* x2;
}
parameters{
  real beta[3];
  real gamma[2];
  real<lower=0> deltahat;
  real<lower=-1, upper=1> beta12hat;
  real<lower=-1, upper=1> rho;
  real<lower=0> sd_nu;
  vector<lower=0>[2] sd_eta;
}
transformed parameters{
  real<lower=0> denom;
  vector [2] sd_x;
  denom = 1 - beta12hat^2;
  // This is how the correlation is taken into account.
  sd_x[1] = sqrt((1/deltahat ^ 2 * sd_eta[1] ^ 2) +
                   (beta12hat ^ 2 * sd_eta[2] ^ 2)
                 + 2 * rho * deltahat * sd_eta[1] * beta12hat * sd_eta[2])
  / denom;
  sd_x[2] = sqrt((deltahat ^ 2 * sd_eta[2] ^ 2) +
                   (beta12hat ^ 2 * sd_eta[1] ^ 2)
                 + 2 * rho * deltahat * sd_eta[2] * beta12hat * sd_eta[1])
  / denom;
}
model{
  //priors
  beta ~ normal(0, 1);
  gamma ~ normal(0, 1);
  // rescaling
  beta12hat ~ normal(0, .5);
  deltahat ~ lognormal(0, 1);
  // rho and beta12hat with the same prior is equal probability to
  // both hypotheses?
  rho ~ normal(0, .5);
  sd_nu ~ lognormal(0, 1);
  sd_eta ~ lognormal(0, 1);
  {
    vector[N] mu[3];
    vector[N] eta[2];

    // demand
    mu[1] = (1/deltahat * (beta[2] + gamma[1] * z)
             + beta12hat * (beta[3] + gamma[2] * z)) / denom;
    mu[2] = (deltahat * (beta[3]  + gamma[2] * z)
             + beta12hat * (beta[2] + gamma[1] * z)) / denom;
    x1 ~ normal(mu[1], sd_x[1]);
    x2 ~ normal(mu[2], sd_x[2]);
 }
}
