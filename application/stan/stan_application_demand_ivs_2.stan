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
  vector[N] x1_2;
  vector[N] x2_2;
  vector[N] x12;
  vector[2] x[N];
  real<lower=0> scale;
  x1_2 = x1 .* x1;
  x2_2 = x2 .* x2;
  x12 = x1 .* x2;
  scale = 1;
  for (i in 1:N){
    x[i][1] = x1[i];
    x[i][2] = x2[i];
  }
}
parameters{
  real beta[2];
  real gamma[2];
  real<lower=0> deltahat;
  real<lower=-1, upper=1> beta12hat;
  real<lower=-1, upper=1> rho;
  vector<lower=0>[2] sd_eta;
}
transformed parameters{
  real<lower=0> denom;
  cov_matrix[2] var_x;
  real<lower=0> delta[2];
  real beta12;
  // vector<lower=0>[2] sd_eta;

  delta[1] = deltahat * scale;
  delta[2] = 1/deltahat * scale;
  // sd_eta[1] = deltahat * sd_eta_hat[1];
  // sd_eta[2] = 1/deltahat * sd_eta_hat[2];

  beta12 = scale * beta12hat;
  denom = delta[1] * delta[2] - beta12^2;
  // This is how the correlation is taken into account.
  var_x[1,1] = (delta[2]^2 * sd_eta[1]^2
                 + (beta12^2 * sd_eta[2]^2)
                 + 2 * rho * delta[2] * sd_eta[1] * beta12 * sd_eta[2])
  / denom^2;
  var_x[2,2] = (delta[1]^2 * sd_eta[2]^2
                  + (beta12^2 * sd_eta[1]^2)
                  + 2 * rho * delta[1] * beta12 * sd_eta[1] * sd_eta[2])
  / denom^2;
  var_x[1,2] = (delta[2] * beta12 * sd_eta[1]^2
            + delta[1] * beta12 * sd_eta[2]^2
            + rho * beta12 * (delta[1] + delta[2]) * sd_eta[1] * sd_eta[2])
  / denom^2;
  var_x[2,1] = var_x[1,2];
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
  //scale ~ lognormal(0, 1);
  sd_eta ~ lognormal(0, 1);

  // demand
  {
    vector[2] mu;
    for (i in 1:N){
      mu[1] = (delta[2] * (beta[1] + gamma[1] * z[i])
               + beta12 * (beta[2] + gamma[2] * w[i])) / denom;
      mu[2] = (delta[1] * (beta[2] + gamma[2] * w[i])
               + beta12 * (beta[1] + gamma[1] * z[i])) / denom;
      x[i] ~ multi_normal(mu, var_x);
    }
  }
}
generated quantities{
}
