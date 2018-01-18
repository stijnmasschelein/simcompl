data{
  int<lower=1> N;
  matrix[N, 2] x;
  vector[N] y;
  matrix[N, 2] z;
  real prior_gamma;
  real prior_lkj;
}
parameters{
  vector[2] beta0;
  corr_matrix[2] B;
  cholesky_factor_corr[2] L_rho;
  matrix[2, 2] Gamma;
  real<lower=0> delta;
  vector<lower=0>[2] sigma;
  real<lower=0> sigma_nu;
  vector<lower=0>[2] sigma_x;
  real<lower=0> lambda;
  matrix[N, 2] unit;
}
transformed parameters{
  row_vector[2] scaled_delta = [delta, 1/delta];
}
model{
  matrix[2, 2] C = quad_form_diag(B, scaled_delta);
  matrix[2, 2] inv_C = inverse_spd(C);
  matrix[N, 2] Mu;
  matrix[N, 2] Mu_temp;
  vector[N] muy;
  for(i in 1:N){
    // row_vector[2] eps;
    row_vector[2] zGamma = z[i] * Gamma;
    Mu[i] = zGamma * inv_C;
    // eps = (x[i] - Mu[i]) * C;
    // muy[i] = z[i] * beta0 + (zGamma + eps) * x[i]'
    //   - lambda * quad_form(C, x[i]') / 2;
  }

  Mu_temp = z * Gamma + diag_pre_multiply(sigma, L_rho) * unit;
  Mu = Mu_temp * inv_C;
  muy = z * beta0 + rows_dot_product(Mu_temp, x)
        - lambda/2 * rows_dot_product(x * C, x);

  // add an additional independent sigma. Similar to latent GP. This would
  // require a full rewrite of x, z into matrices
  to_vector(unit) ~ normal(0, 1);
  x[1:N, 1] ~ normal(Mu[1:N, 1], sigma_x[1]);
  x[1:N, 2] ~ normal(Mu[1:N, 2], sigma_x[2]);
  y ~ normal(muy, sigma_nu);

  B ~ lkj_corr(prior_lkj);
  beta0 ~ normal(0, 1);
  L_rho ~ lkj_corr_cholesky(prior_lkj);
  to_vector(Gamma) ~ normal(0, prior_gamma);
  delta ~ lognormal(0, 1);
  sigma ~ normal(0, 1);
  sigma_nu ~ normal(0, 1);
  sigma_x ~ normal(0, 1);
  lambda ~ normal(0, 1);
}
generated quantities{
  matrix[2, 2] Rho = L_rho * L_rho';
}
