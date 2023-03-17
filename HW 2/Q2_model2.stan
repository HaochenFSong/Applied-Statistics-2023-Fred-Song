
// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N; // number of age group 50 -100
  int<lower=0> Time; // number of years 1990 -2020
  int y[N,Time]; //number of deaths
  vector[N] x; //age group
  matrix[N, Time] log_p; //log of population
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[Time] log_alpha;
  vector[Time] beta;
}
transformed parameters{
  matrix[N,Time] log_lambda;
  for (t in 1:Time){
    log_lambda[,t] = log_alpha[t] + x * beta[t];
  }

}

model {
  for (n in 1:N){
    for (t in 1:Time){
      y[n,t] ~ poisson_log(log_lambda[n,t]+log_p[n,t]);
    }
  }
  
  log_alpha ~ normal(-12.28,10); //relax some standard deviation
  beta ~ normal(0.1159, 0.1);//relax some standard deviation
}

