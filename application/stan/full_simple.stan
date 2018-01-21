data{
  int<lower=1> N;
  vector[N] x1;
  vector[N] x2;
  vector[N] y;
  vector[N] z;
  real prior;
}
transformed data{
  vector[N] x1_2;
  vector[N] x2_2;
  vector[N] x12;
  vector[2] x[N];
  x1_2 = x1 .* x1;
  x2_2 = x2 .* x2;
  x12 = x1 .* x2;
  for (i in 1:N){
    x[i][1] = x1[i];
    x[i][2] = x2[i];
  }
}
parameters{
  vector[3] beta;
  real<lower=0> deltahat;
  vector<lower=0>[2] sd_eta;
  real<lower=0> sd_nu;
  vector<lower=0, upper=1>[2] corr; //change to fix rho
  vector[3] gamma;
  real<lower=0> scale;
  real<lower=0> scale_eqn;
}
transformed parameters{
  real<lower=0> denom;
  real<lower=0> delta[2];
  real beta12;
  real<lower=-1, upper=1> beta12hat;
  real<lower=-1, upper=1> rho;

  beta12hat = 2 * corr[1] - 1;
  rho = 2 * corr[2] - 1;
  delta[1] = deltahat * scale;
  delta[2] = 1/deltahat * scale;

  beta12 = scale * beta12hat;
  denom = scale * (1 - beta12hat ^ 2);
}

model{
  matrix[2, 2] var_x;
  matrix[2, 2] L_x;
  var_x[1,1] = (1/deltahat^2 * sd_eta[1]^2 + beta12hat^2 * sd_eta[2]^2
                + 2 * rho * 1/deltahat * sd_eta[1] * beta12hat * sd_eta[2]);
  var_x[2,2] = (deltahat^2 * sd_eta[2]^2 + beta12hat^2 * sd_eta[1]^2
                + 2 * rho * deltahat * beta12hat * sd_eta[1] * sd_eta[2]);
  var_x[1,2] = (beta12hat * (1/deltahat * sd_eta[1]^2
                             + deltahat * sd_eta[2]^2)
                + rho * (beta12hat^2 + 1) * sd_eta[1] * sd_eta[2]);
  var_x[2,1] = var_x[1,2];
  L_x = cholesky_decompose(var_x/denom^2);

  //priors
  // beta(4,4) 52% of the distribution in [-.25, .25]
  // beta(8,8) 99% of the distribution in [-.5, .5]
  corr[2] ~ beta(8, 8);
  beta ~ normal(0, prior);
  deltahat ~ lognormal(0, 1);
  sd_eta ~ normal(0, 1);

  scale ~ normal(0, 1);
  scale_eqn ~ normal(0, 1);
  sd_nu ~ normal(0, 1);

  gamma ~ normal(0, prior);

  {
    vector[2] mu;
    vector[2] eta;
    vector[N] muy;
    for (i in 1:N){
      mu[1] = (1/deltahat * (beta[1] + gamma[1] * z[i])
               + beta12hat * (beta[2] + gamma[2] * z[i]));
      mu[2] = (deltahat * (beta[2] + gamma[2] * z[i])
               + beta12hat * (beta[1] + gamma[1] * z[i]));
      x[i] ~ multi_normal_cholesky(mu/denom, L_x);

      eta[1] = delta[1] * x1[i] - beta12 * delta[2] * x2[i] - mu[1];
      eta[2] = delta[2] * x2[i] - beta12 * delta[2] * x1[i] - mu[2];
      muy[i] = beta[3] +
              scale_eqn * ((beta[1] + gamma[1] * z[i] + eta[1]) .* x1[i]
              + (beta[2] + gamma[2] * z[i] + eta[2]) .* x2[i]
              + beta12 * x12[i]
              - .5 * (delta[1] * x1_2[i] + delta[2] * x2_2[i]))
              + gamma[3] * z[i];
    }
    y ~ normal(muy, sd_nu);
  }
}
generated quantities{
}

