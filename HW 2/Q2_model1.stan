
// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N; // number of age group 50 -100
  int y[N]; //number of deaths
  vector[N] x; //age group
  vector[N] log_p; //log of population
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real log_alpha;
  real beta;
}
transformed parameters{
  vector[N] log_lambda;
  
  log_lambda = log_alpha + beta * x;
}

model {
  y ~ poisson_log(log_lambda+log_p);
  log_alpha ~ normal(-12.28,3);
  beta ~ normal(0.1159, 0.05);
}

