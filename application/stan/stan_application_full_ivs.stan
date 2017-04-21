data{
  int<lower=1> N;
  vector[N] y;
  vector[N] x1;
  vector[N] x2;
  vector[N] z;
  vector[N] w;
  real prior;
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
  real<lower=0> scale;
  real<lower=-1, upper=1> rho;
  vector<lower=0>[2] sd_eta;
  real<lower=0> sd_nu;
}
transformed parameters{
  real<lower=0> denom;
  vector<lower=0>[2] sd_x;
  real<lower=0> delta[2];
  real beta12;

  delta[1] = deltahat * scale;
  delta[2] = 1/deltahat * scale;
  beta12 = scale * beta12hat;
  denom = delta[1] * delta[2] - beta12^2;
  // This is how the correlation is taken into account.
  sd_x[1] = sqrt((delta[2] ^ 2 * sd_eta[1] ^ 2)
                 + (beta12 ^ 2 * sd_eta[2] ^ 2)
                 + 2 * rho * delta[2] * sd_eta[1] * beta12 * sd_eta[2])
  / denom;
  sd_x[2] = sqrt(
    (delta[2] * sd_eta[1]) ^ 2
    + (beta12 / denom + 1 / delta[2]) ^ 2 * sd_eta[2] ^ 2
    + 2 * rho * (delta[2] / denom) * (beta12 / denom + 1 / delta[2])
      * sd_eta[1] * sd_eta[2]
    ); // conditional on environment

}
model{
  //priors
  beta ~ normal(0, prior);
  gamma ~ normal(0, prior);
  // rescaling
  beta12hat ~ normal(0, .5);
  deltahat ~ lognormal(0, 1);
  // rho and beta12hat with the same prior is equal probability to
  // both hypotheses?
  rho ~ normal(0, .5);
  scale ~ lognormal(0, 1);
  sd_eta ~ lognormal(0, 1);
  sd_nu ~ lognormal(0, 1);

  {
    // demand
    vector[N] mu[3];
    vector[N] envir[2];
    vector[N] E_out[2];
    vector[N] eta[2];

    envir[1] = gamma[1] * z;
    envir[2] = gamma[2] * w;
    mu[1] = (delta[2] * (beta[2] + envir[1])
             + beta12 *(beta[3] + envir[2])) / denom;
    mu[2] = (beta[3] + beta12 * mu[1] + envir[2])/delta[2]; //conditional on x_1
    x1 ~ normal(mu[1], sd_x[1]);
    x2 ~ normal(mu[2], sd_x[2]);

    // unobservables
    E_out[1] = x1 - mu[1];
    E_out[2] = x2 - mu[2];
    eta[1] = delta[1] * E_out[1] - beta12 * E_out[2];
    eta[2] = delta[2] * E_out[2] - beta12 * E_out[1];

    // production
    mu[3] = beta[1] + (beta[2] + envir[1] + eta[1]) .* x1
                     + (beta[3] + envir[2] + eta[2]) .* x2
             + beta12 * x12
             - .5 * (delta[1] * x1_2 + delta[2] * x2_2);
    y ~ normal(mu[3], sd_nu);
  }
}
generated quantities{
}
