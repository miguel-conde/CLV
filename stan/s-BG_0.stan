functions {
  // Retention Rate
  real rr(int i, real alpha, real beta) {
  
    return (beta + i - 1) / (alpha + beta + i - 1);
  }
  
  // log-Survivor function
  real log_survivor_fun(int t, real alpha, real beta) {
    return lbeta(alpha, beta + t) - lbeta(alpha, beta);
  }
  
  // Survivor function
  real survivor_fun(int t, real alpha, real beta) {
    return exp(log_survivor_fun(t, alpha, beta));
  }
  
  // Loglikelihood 
  real ll(int[] n, real alpha, real beta, int n_0) {
    int T = num_elements(n);
    vector[T] prob;
    vector[T] lprob;
    real loglikelihood;
    
    prob[1] = alpha / (alpha + beta);
    lprob[1] = n[1] * log(prob[1]);
    for(t in 2:T) {
      prob[t] = (beta + t - 2) / (alpha + beta + t - 1) * prob[t-1];
      lprob[t] = n[t] * log(prob[t]);
    }
    
    loglikelihood = sum(lprob) + (n_0 - sum(n)) * log_survivor_fun(T, alpha, beta);
    
    return loglikelihood;
  }
}

data {
  int<lower=0> N;    // Nº de muestras
  int<lower=0> n[N]; // "Muertes" en el periodo
  int<lower=0> n_0;  // población inicial
  
  real<lower=0> mean_alpha;
  real<lower=0> mean_beta;
}

parameters {
  real<lower=0> beta;
  real<lower=0> alpha;
}

transformed parameters {
}

model {
  
  alpha ~ exponential(mean_alpha);
  beta  ~ exponential(mean_beta);
  
  target += ll(n, alpha, beta, n_0);
}

generated quantities {
  real sim_r[N];
  real sim_S[N];
  
  for (t in 1:num_elements(n)) {
    sim_r[t] = rr(t, alpha, beta);
    sim_S[t] = survivor_fun(t, alpha, beta);
  }
}
