data{
  int<lower=1> N;
  row_vector[2] x[N];
  vector[N] y;
  row_vector[2] z[N];
  real prior_gamma;
  real prior_lkj;
}
parameters{
  corr_matrix[2] B;
  corr_matrix[2] Rho;
  matrix[2, 2] Gamma;
  real<lower=0> delta;
  vector<lower=0>[2] sigma;
}
transformed parameters{
  row_vector[2] scaled_delta;
  scaled_delta = [delta, 1/delta];
}
model{
  matrix[2, 2] inv_C = inverse_spd(quad_form_diag(B, scaled_delta));
  matrix[2, 2] Sigma = quad_form(quad_form_diag(Rho, sigma), inv_C);
  matrix[2, 2] L = cholesky_decompose(Sigma);
  row_vector[2] Mu[N];
  for(i in 1:N){
    Mu[i] = z[i] * Gamma * inv_C;
  }
  x ~ multi_normal_cholesky(Mu, L);

  B ~ lkj_corr(prior_lkj);
  Rho ~ lkj_corr(prior_lkj);
  to_vector(Gamma) ~ normal(0, prior_gamma);
  delta ~ lognormal(0, 1);
  sigma ~ normal(0, 1);
}
