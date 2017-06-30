data{
  int<lower=1> N;
  vector[N] x1;
  vector[N] x2;
  vector[N] z;
  real prior;
  // vector<lower=0>[2] sd_br;
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
  real<lower=0> deltahat;
  vector<lower=0>[2] sd_eta;
  real<lower=0> lambda;
  real<lower=0, upper=1>pseudoR2;
  vector[2] gamma;
  vector[2] u_raw;
}
transformed parameters{
  real<lower=0> denom;
  real<lower=0> delta[2];
  real<lower=-1, upper=1> beta12hat;
  real<lower=-1, upper=1> rho;
  real beta12;
  vector<lower=-1, upper=1>[2] u;
  // real gammarho;

  u = u_raw / sqrt(dot_self(u_raw));
  beta12hat = (1 - sqrt(1 - (pseudoR2 * u[1] ^ 2)))
               / (sqrt(pseudoR2) * u[1]);
  rho = sqrt(pseudoR2) * u[2];

  // {
  //  gammarho = sqrt(pseudoR2) * u[3];
  //  gamma[1] = gamma1;
  //  gamma[2] = gammarho * gamma1 * gamma2;
  // }


  delta[1] = deltahat * scale;
  delta[2] = 1/deltahat * scale;
  // implies scale = sqrt(delta[1] * delta[2])

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
  deltahat ~ normal(0, lambda);
  // scale ~ lognormal(0, 1);
  sd_eta[1] ~ normal(0, lambda);
  sd_eta[2] ~ normal(0, 1/lambda);
  lambda ~ lognormal(0, 1);

  u_raw ~ normal(0, 1);
  pseudoR2 ~ beta(3/2.0, (1 - .25) / .25 * 3/2.0);
  gamma ~ normal(0, prior);

  // sd_br ~ normal(0, .6);

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

