data{
  int<lower=1> N;
  vector[N] x1;
  vector[N] x2;
  vector[N] z;
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
  vector[2] beta;
  vector[2] gamma;
  real<lower=0> deltahat;
  vector<lower=0>[2] sd_eta;
  vector[2] br;
  real<lower=0> lambda;
}
transformed parameters{
  real<lower=0> denom;
  real<lower=0> delta[2];
  real<lower=-1, upper=1> beta12hat;
  real<lower=-1, upper=1> rho;
  real beta12;

  beta12hat = 1 - 2 * inv_logit(br[1]);
  rho = 1 - 2 * inv_logit(br[2]);

  delta[1] = deltahat * scale;
  delta[2] = 1/deltahat * scale;

  beta12 = scale * beta12hat;
  denom = delta[1] * delta[2] - beta12^2;

}

model{
  matrix[2, 2] var_x;
  matrix[2, 2] L_x;
  // This is how the correlation is taken into account.
  var_x[1,1] = (delta[2]^2 * sd_eta[1]^2
                + beta12^2 * sd_eta[2]^2
                + 2 * rho * delta[2] * sd_eta[1] * beta12 * sd_eta[2])
  / denom^2;
  var_x[2,2] = (delta[1]^2 * sd_eta[2]^2
                + beta12^2 * sd_eta[1]^2
                + 2 * rho * delta[1] * beta12 * sd_eta[1] * sd_eta[2])
  / denom^2;
  var_x[1,2] = (beta12 * (delta[2] * sd_eta[1]^2 + delta[1] * sd_eta[2]^2)
                + rho * (beta12^2 + delta[1] * delta[2]) * sd_eta[1] * sd_eta[2])
  / denom^2;
  var_x[2,1] = var_x[1,2];
  L_x = cholesky_decompose(var_x);

  //priors
  beta ~ normal(0, prior);
  gamma ~ normal(0, prior);
  br ~ normal(0, 1.2);
  deltahat ~ normal(0, lambda);
  // scale ~ lognormal(0, 1);
  sd_eta[1] ~ normal(0, lambda);
  sd_eta[2] ~ normal(0, 1/lambda);
  lambda ~ lognormal(0, 1);

  // demand
  {
    vector[2] mu;
    for (i in 1:N){
      mu[1] = (delta[2] * (beta[1] + gamma[1] * z[i])
               + beta12 * (beta[2] + gamma[2] * z[i])) / denom;
      mu[2] = (delta[1] * (beta[2] + gamma[2] * z[i])
               + beta12 * (beta[1] + gamma[1] * z[i])) / denom;
      x[i] ~ multi_normal_cholesky(mu, L_x);
    }
  }
}
generated quantities{
}
