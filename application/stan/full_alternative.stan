data{
  int<lower=1> N;
  matrix[N, 2] x;
  vector[N] y;
  matrix[N, 2] z;
  real prior_scale;
  real prior_gamma;
  real prior_lkj;
}
transformed data{
}
parameters{
  vector[2] beta0;
  corr_matrix[2] B;
  cholesky_factor_corr[2] L_rho;
  matrix[2, 2] Gamma;
  real<lower=0> delta;
  vector<lower=0>[2] sigma;
  real<lower=0> sigma_nu;
  real<lower=0> lambda;
  vector<lower=0>[2] sigma_x;
  real<lower=0> scale;
  matrix[N, 2] unit;
}
transformed parameters{
  row_vector[2] scaled_delta = [delta, 1/delta];
}
model{
  matrix[2, 2] DBD = quad_form_diag(B, sqrt(scaled_delta));
  matrix[2, 2] inv_DBD = inverse_spd(DBD);
  matrix[N, 2] Mu_temp = z * Gamma
                         + unit * diag_pre_multiply(sigma, L_rho');
  matrix[N, 2] Mu = Mu_temp * inv_DBD / lambda;
  vector[N] muy = z * beta0
                  + scale * rows_dot_product(Mu_temp - .5 * x * DBD, x);

  to_vector(unit) ~ normal(0, 1);
  x[1:N, 1] ~ normal(Mu[1:N, 1], sigma_x[1]);
  x[1:N, 2] ~ normal(Mu[1:N, 2], sigma_x[2]);
  y ~ normal(muy, sigma_nu);

  B ~ lkj_corr(prior_lkj);
  beta0 ~ normal(0, 1);
  L_rho ~ lkj_corr_cholesky(prior_lkj);
  to_vector(Gamma) ~ normal(0, prior_gamma);
  delta ~ normal(0, 1);
  sigma ~ normal(0, 1);
  sigma_nu ~ normal(0, 1);
  sigma_x ~ normal(0, prior_scale);
  lambda ~ normal(0, 1);
  scale ~ normal(0, 1);
}
generated quantities{
  matrix[2, 2] Rho = multiply_lower_tri_self_transpose(L_rho);
}
