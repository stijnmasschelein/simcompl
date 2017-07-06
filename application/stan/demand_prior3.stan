functions{
  real beta_calc (real partial, real ratio){
    real result;
    real result1;
    real result2;
    real partial2;
    real ratio2;
    partial2 = partial ^ 2;
    ratio2 = ratio ^ 2;
    result1 = sqrt(ratio2 + partial2 * (4 - ratio2));
    result2 = ratio * sqrt(1 - partial2);
    result = (result1 - result2) / (2 * partial);
    if (is_nan(result)){
      print(partial, ratio)
    }
    return(result);
  }
  real ratio_calc (real d1, real d2, real sd1, real sd2){
    real result;
    real ratio;
    ratio = sqrt(d1 / d2) * (sd2/sd1);
    result = ratio + 1/ratio;
    if (is_nan(result)){
      print(d1, d2, sd1, sd2);
    }
    return(result);
  }
}
data{
  int<lower=1> N;
  vector[N] x1;
  vector[N] x2;
  vector[N] z;
  real prior;
  real<lower=0, upper=1> R2;
  // vector<lower=0>[2] sd_br;
}
transformed data{
  vector[N] x1_2;
  vector[N] x2_2;
  vector[N] x12;
  vector[2] x[N];

  real<lower=0> scale;
  int<lower=2> ncor;
  x1_2 = x1 .* x1;
  x2_2 = x2 .* x2;
  x12 = x1 .* x2;
  scale = 1;
  ncor = 2;
  for (i in 1:N){
    x[i][1] = x1[i];
    x[i][2] = x2[i];
  }
}
parameters{
  vector[2] beta;
  real<lower=0> deltahat;
  vector<lower=0>[2] sd_eta;
  real<lower=0, upper=1>pseudoR2;
  vector[ncor] u_raw;
  real<lower=0> lambda;
  // real<lower=-1, upper=1> gamma1;
  // real<lower=0> gamma2;
  vector[2] gamma;
}
transformed parameters{
  real<lower=0> denom;
  real<lower=0> delta[2];
  real<lower=-1, upper=1> beta12hat;
  real<lower=-1, upper=1> rho;
  real beta12;
  vector<lower=-1, upper=1>[ncor] u;
  // real<lower=-1, upper=1> grho;
  real avratio;

  delta[1] = deltahat * scale;
  delta[2] = 1/deltahat * scale;
  // implies scale = sqrt(delta[1] * delta[2])

  avratio = ratio_calc(delta[1], delta[2], sd_eta[1], sd_eta[2]);

  u = u_raw / sqrt(dot_self(u_raw));
  beta12hat = beta_calc(sqrt(pseudoR2) * u[1], avratio);
  rho = sqrt(pseudoR2) * u[2];

  // grho = sqrt(pseudoR2) * u[3];
  // gamma[1] = gamma1 * gamma2;
  // gamma[2] = grho / gamma1 * gamma2;

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
  lambda ~ normal(0, 1.5);
  deltahat ~ lognormal(0, lambda);
  // scale ~ lognormal(0, 1);
  sd_eta[1] ~ lognormal(0, lambda);
  sd_eta[2] ~ lognormal(0, 1/lambda);

  u_raw ~ normal(0, 1);
  pseudoR2 ~ beta(ncor/2.0, (1 - R2) / R2 * ncor/2.0);
  // gamma2 ~ normal(0, prior * sqrt(gamma1^2));
  gamma ~ normal(0, prior);

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

