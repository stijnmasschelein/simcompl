data{
  int<lower=1> N;
  matrix[N, 2] x;
  vector[N] y;
  matrix[N, 2] z;
  real prior_gamma;
  real prior_lkj;
}
transformed data{
  vector<lower=0>[2] sigma_x = [.001, .001]';
  row_vector[2] x_array[N];
  for (i in 1:N){
    x_array[i] = x[i, 1:2];
  }
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
  real<lower=0> scale;
  // matrix[N, 2] unit;
}
transformed parameters{
  row_vector[2] scaled_delta = [delta, 1/delta];
}
model{
  matrix[2, 2] DBD = quad_form_diag(B, sqrt(scaled_delta));
  matrix[2, 2] inv_DBD = inverse_spd(DBD);
  matrix[N, 2] Mu_temp = z * Gamma;
                         // + unit * diag_pre_multiply(sigma, L_rho');
                         // + unit * diag_matrix(sigma);
  matrix[N, 2] Mu = Mu_temp * inv_DBD / lambda;
  // simplification possible?
  matrix[N, 2] E = lambda * x * DBD - Mu_temp;
  vector[N] muy = z * beta0 + scale
                  * rows_dot_product(Mu_temp + E - .5 * lambda * x * DBD, x);
  matrix[2, 2] L_x = cholesky_decompose(
                      crossprod(
                        diag_pre_multiply(sigma, L_rho) * inv_DBD / lambda));
  row_vector[2] Mu_array[N];
  // to_vector(unit) ~ normal(0, 1);
  // x[1:N, 1] ~ normal(Mu[1:N, 1], sigma_x[1]);
  // x[1:N, 2] ~ normal(Mu[1:N, 2], sigma_x[2]);
  for (i in 1:N){
    Mu_array[i] = Mu[i, 1:2];
  }
  x_array ~ multi_normal_cholesky(Mu_array, L_x);
  y ~ normal(muy, sigma_nu);

  B ~ lkj_corr(prior_lkj);
  beta0 ~ normal(0, 1);
  L_rho ~ lkj_corr_cholesky(prior_lkj);
  to_vector(Gamma) ~ normal(0, prior_gamma);
  delta ~ normal(0, 1);
  sigma ~ normal(0, 1);
  sigma_nu ~ normal(0, 1);
  lambda ~ normal(0, 1);
  scale ~ normal(0, 1);
}
generated quantities{
  matrix[2, 2] Rho = multiply_lower_tri_self_transpose(L_rho);
}
