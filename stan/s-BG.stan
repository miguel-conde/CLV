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
  int<lower=0> N;      // Nº de muestras
  int<lower=0> G;      // Nº de grupos (cohortes)
  int<lower=0> N_g[G]; // Nº muestras por grupo
  int<lower=0> g[N];   // grupo
  int<lower=0> n[N];   // "Muertes" en el periodo
  int<lower=0> n_0[G]; // poblaciones iniciales
}

transformed data {
  int init_idx[G];
  int end_idx[G];
  
  // Índices de inicio y fin de cada grupo
  init_idx[1] = 1;
  end_idx[1]  = N_g[1];
  
  if (G > 1) {
    for (group in 2:G) {
      init_idx[group] = end_idx[group-1] + 1;
      end_idx[group]  = init_idx[group] + N_g[group] - 1;
    }
  }
}

parameters {
  // Hyper-parameters
  real<lower=0> avg_beta;
  real<lower=0> avg_alpha;
  
  // Parameters
  real<lower=0> beta[G];
  real<lower=0> alpha[G];
}

transformed parameters {
}

model {
  
  // Hyper-priors
  avg_alpha ~ exponential(1);
  avg_beta  ~  exponential(1);
  
  for (group in 1:G) {
    // Priors
    alpha[group] ~ exponential(avg_alpha);
    beta[group]  ~ exponential(avg_beta);
    
    // Likelihood
    target += ll(n[init_idx[group]:end_idx[group]], alpha[group], beta[group], n_0[group]);
  }
}

generated quantities {
  real sim_r[N]; // Retention rates
  real sim_S[N]; // Survival function
  
  for (i in 1:N) {
    sim_r[i] = rr(i - init_idx[g[i]] + 1, alpha[g[i]], beta[g[i]]);
    sim_S[i] = survivor_fun(i - init_idx[g[i]] + 1, alpha[g[i]], beta[g[i]]);
  }
}

