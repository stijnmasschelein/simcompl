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
  vector[3] xy[N];
  x1_2 = x1 .* x1;
  x2_2 = x2 .* x2;
  x12 = x1 .* x2;
  for (i in 1:N){
    xy[i][1] = x1[i];
    xy[i][2] = x2[i];
    xy[i][3] = y[i];
  }
}
parameters{
  real beta[3];
  real gamma[2];
  real<lower=0> deltahat;
  real<lower=-1, upper=1> beta12hat;
  vector<lower=0>[3] sds;
  cholesky_factor_corr[3] L;
  real<lower=0> scale;
}
transformed parameters{
  real<lower=0> denom;
  real<lower=0> delta[2];
  real beta12;

  delta[1] = deltahat * scale;
  delta[2] = 1/deltahat * scale;

  beta12 = scale * beta12hat;
  denom = delta[1] * delta[2] - beta12^2;

}
model{
  //priors
  beta ~ normal(0, prior);
  gamma ~ normal(0, prior);
  // rescaling
  beta12hat ~ normal(0, .5);
  deltahat ~ lognormal(0, 1);
  scale ~ lognormal(0, 1);
  sds ~ lognormal(0, 1);
  L ~ lkj_corr_cholesky(2.0);

  for (i in 1:N){
    vector[3] mu;
    vector[2] E_out;
    vector[2] eta;
    mu[1] = (delta[2] * (beta[2] + gamma[1] * z[i])
             + beta12 * (beta[2] + gamma[2] * w[i])) / denom;
    mu[2] = (delta[1] * (beta[3] + gamma[2] * w[i])
             + beta12 * (beta[1] + gamma[1] * z[i])) / denom;

    E_out[1] = x1[i] - mu[1];
    E_out[2] = x2[i] - mu[2];
    eta[1] = delta[1] * E_out[1] - beta12 * E_out[2];
    eta[2] = delta[2] * E_out[2] - beta12 * E_out[1];

    mu[3] = beta[1] + (beta[2] + gamma[1] * z[i] + eta[1]) .* x1[i]
                   + (beta[3] + gamma[2] * w[i] + eta[2]) .* x2[i]
           + beta12 * x12[i]
           - .5 * (delta[1] * x1_2[i] + delta[2] * x2_2[i]);
    xy[i] ~ multi_normal_cholesky(mu, diag_pre_multiply(sds, L));
  }
}
generated quantities{
  corr_matrix[3] corr_L;
  corr_L = multiply_lower_tri_self_transpose(L);
}
