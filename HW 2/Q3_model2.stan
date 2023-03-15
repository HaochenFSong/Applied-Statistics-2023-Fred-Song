// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector[N] dist;
  vector[N] arsenic;
  int<lower=0,upper=1> p[N];
}
parameters {
  real beta0;
  real beta1;
  real beta2;
  real beta3;
}
model {
  p ~ bernoulli_logit(beta0 + beta1 * (dist - mean(dist)) + beta2 * (log(arsenic) - mean(log(arsenic))) + beta3 * (dist - mean(dist)) .* (log(arsenic) - mean(log(arsenic))));
  
  // priors
  beta0 ~ normal(0,1);
  beta1 ~ normal(0,1);
  beta2 ~ normal(0,1);
  beta3 ~ normal(0,1);
}

generated quantities {
  vector[N] log_lik;    // pointwise log-likelihood for LOO
  vector[N] p_hat; // replications from posterior predictive dist

  for (n in 1:N) {
    real weight_n = beta0 + beta1 * (dist[n] - mean(dist)) + beta2 * (log(arsenic[n]) - mean(log(arsenic))) + beta3 * (dist[n] - mean(dist)) * (log(arsenic[n]) - mean(log(arsenic)));
    p_hat[n] = bernoulli_rng(inv_logit(weight_n));
    log_lik[n] = bernoulli_logit_lpmf( p[n] | weight_n);
  }
}

