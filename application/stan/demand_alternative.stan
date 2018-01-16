data{
  int<lower=1> N;
  matrix[N, 2] x;
  vector[N] y;
  matrix[N, 2] z;
  real prior_gamma;
  real prior_lkj;
}
transformed data{
  vector[2] zeros = rep_vector(0, 2);
}
parameters{
  corr_matrix[2] B;
  corr_matrix[2] Rho;
  matrix[2, 2] Gamma;
  real<lower=0> delta;
  vector<lower=0>[2] sigma;
}
transformed parameters{
  vector[2] scaled_delta;
  scaled_delta = [delta, 1/delta]';
}
model{
  // Steps to optimise
  // 1. Make Mu an array of vectors
  // 2. Adjust x and z
  // 3. Vectorise multi_normal
  // 4. Use cholesky version

  matrix[2, 2] inv_C = inverse_spd(quad_form_diag(B, scaled_delta));
  matrix[2, 2] Sigma = quad_form(quad_form_diag(Rho, sigma), inv_C);
  matrix[N, 2] Mu = z * Gamma * inv_C;
  for(i in 1:N){
    x[i, 1:2] ~ multi_normal(Mu[i,1:2], Sigma);
  }

  B ~ lkj_corr(prior_lkj);
  Rho ~ lkj_corr(prior_lkj);
  to_vector(Gamma) ~ normal(0, prior_gamma);
  delta ~ normal(0, 1);
  sigma ~ normal(0, 1);
}
