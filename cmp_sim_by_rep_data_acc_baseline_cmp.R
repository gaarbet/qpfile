rmse <- function(estimates, true_values) {
  sqrt(mean((estimates - true_values)^2))
}
bias <- function(estimated, true){
  mean(estimated - true)
}

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(cmdstanr)
library(loo)
library(COMPoissonReg)



baseline <- write_stan_file('
functions{

  real approximation(real log_lambda, real nu){
  
      real log_mu = log_lambda / nu;
      real nu_mu = nu * exp(log_mu);
      real nu2 = nu^2;
      // first 4 terms of the residual series
      real log_sum_resid = log1p(
        nu_mu^(-1) * (nu2 - 1) / 24 +
        nu_mu^(-2) * (nu2 - 1) / 1152 * (nu2 + 23) +
        nu_mu^(-3) * (nu2 - 1) / 414720 * (5 * nu2^2 - 298 * nu2 + 11237)+
        nu_mu^(-4) * (nu2 - 1) / 39813120 * (5 * nu2^3 - 1887*nu2^2 - 241041*nu^2 + 2482411)

      );
      return nu_mu + log_sum_resid  -
            ((log(2 * pi()) + log_mu) * (nu - 1) / 2 + log(nu) / 2);  
  
  }
  real summation(real log_lambda, real nu){
  
    real z = negative_infinity();
    real z_last = 0;
    real t = 0;
    
    for (j in 0:1000000) {
      t = t + 1;
      z_last = z;
      z = log_sum_exp(z, j * log_lambda  - nu * lgamma(j+1));

      if ((abs(z - z_last) < 1e-05)) {
        break;
  
      }
      if( t > 999998){
        reject("max terms hit, prob bad value");
      }

    }

    
    return z;
  }
  
  
  real partial_sum(array[] int slice_n, int start, int end,
                   array[] int y_slice,
                   array[] int jj_slice, array[] int ii_slice, vector theta,
                   vector beta, vector nu, vector alpha) {
    real partial_target = 0.0;
    
    for (n in start : end) {
      
      real log_lambda = nu[ii_slice[n]]*(alpha[ii_slice[n]]*theta[jj_slice[n]]+ beta[ii_slice[n]]);
       real log_prob = 0;
      
      
      if (log_lambda / nu[ii_slice[n]] > log(1.5) && log_lambda > log(1.5)) {
      
    //if((log_lambda > 8.25 && nu[ii_slice[n]] > 4 && log_lambda < 40) ||
    //    (log_lambda > 1 && log_lambda < 11 && nu[ii_slice[n]] > 1.5 && nu[ii_slice[n]] < 3.1))
   // {
      
      
        log_prob = y_slice[n] * log_lambda -
                   nu[ii_slice[n]] * lgamma(y_slice[n] + 1) -
                   approximation(log_lambda, nu[ii_slice[n]]);
      } else {
        log_prob = y_slice[n] * log_lambda -
                   nu[ii_slice[n]] * lgamma(y_slice[n] + 1) -
                   summation(log_lambda, nu[ii_slice[n]]);
      }
      
      partial_target += log_prob;
    }
    return partial_target;
  }
}

data {
  int<lower=0> I;
  int<lower=0> J;
  int<lower=1> N;
  int<lower=1> K;
  array[N] int<lower=1, upper=I> ii;
  array[N] int<lower=1, upper=J> jj;
  array[N] int<lower=1, upper=K> kk;
  array[I] int<lower=1,upper=K> item_type_for_beta;
  array[N] int<lower=0> y;
  array[N] int seq_N;
  int<lower=0> grainsize;
}
parameters {
  vector[J] theta;
  vector[I] beta;
  vector<lower=0.01>[I] nu;  //.2

  vector<lower=0>[I] alpha;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_nu;
  
  real mu_beta;
  real<lower=0> sigma_beta;
}
model {

  mu_beta ~ normal(0,5);
  sigma_beta ~ cauchy(0,5);
  beta ~ normal(mu_beta,sigma_beta);

  sigma_nu ~ cauchy(0,5);
  nu ~ lognormal(0, sigma_nu);
  
  theta ~ normal(0, .3); 
  
  alpha ~ lognormal(0,sigma_alpha);
  sigma_alpha ~ cauchy(0,5);

  target += reduce_sum(partial_sum, seq_N, grainsize, y, jj, ii, theta,
                       beta, nu, alpha);
}
generated quantities {
  array[N] real log_lik;

  for (n in 1:N) {
      real log_lambda = nu[ii[n]] * (alpha[ii[n]] * theta[jj[n]] + beta[ii[n]]);
      if (log_lambda / nu[ii[n]] > log(1.5) && log_lambda > log(1.5)) {
          log_lik[n] = y[n] * log_lambda - nu[ii[n]] * lgamma(y[n] + 1) - approximation(log_lambda, nu[ii[n]]);
      } else {
          log_lik[n] = y[n] * log_lambda - nu[ii[n]] * lgamma(y[n] + 1) - summation(log_lambda, nu[ii[n]]);
      }
  }
}


')


model_proposed <- write_stan_file('
functions{

  real approximation(real log_lambda, real nu){
  
      real log_mu = log_lambda / nu;
      real nu_mu = nu * exp(log_mu);
      real nu2 = nu^2;
      // first 4 terms of the residual series
      real log_sum_resid = log1p(
        nu_mu^(-1) * (nu2 - 1) / 24 +
        nu_mu^(-2) * (nu2 - 1) / 1152 * (nu2 + 23) +
        nu_mu^(-3) * (nu2 - 1) / 414720 * (5 * nu2^2 - 298 * nu2 + 11237)+
        nu_mu^(-4) * (nu2 - 1) / 39813120 * (5 * nu2^3 - 1887*nu2^2 - 241041*nu^2 + 2482411)

      );
      return nu_mu + log_sum_resid  -
            ((log(2 * pi()) + log_mu) * (nu - 1) / 2 + log(nu) / 2);  
  
  }
  real summation(real log_lambda, real nu){
  
    real z = negative_infinity();
    real z_last = 0;
    real t = 0;
    
    for (j in 0:1000000) {
      t = t + 1;
      z_last = z;
      z = log_sum_exp(z, j * log_lambda  - nu * lgamma(j+1));

      if ((abs(z - z_last) < 1e-05)) {
        break;
  
      }
      if( t > 999998){
        reject("max terms hit, prob bad value");
      }

    }

    
    return z;
  }
  
  
  real partial_sum(array[] int slice_n, int start, int end,
                   array[] int y_slice,
                   array[] int jj_slice, array[] int ii_slice, vector theta,
                   vector beta, vector nu, vector alpha, array[] int kk_slice, vector gamma) {
    real partial_target = 0.0;
    
    for (n in start : end) {
      
      real log_lambda = nu[ii_slice[n]]*(alpha[ii_slice[n]]*theta[jj_slice[n]]+ beta[ii_slice[n]] + gamma[kk_slice[n]]);
       real log_prob = 0;
      
      
      if (log_lambda / nu[ii_slice[n]] > log(1.5) && log_lambda > log(1.5)) {
      
      
        log_prob = y_slice[n] * log_lambda -
                   nu[ii_slice[n]] * lgamma(y_slice[n] + 1) -
                   approximation(log_lambda, nu[ii_slice[n]]);
      } else {
        log_prob = y_slice[n] * log_lambda -
                   nu[ii_slice[n]] * lgamma(y_slice[n] + 1) -
                   summation(log_lambda, nu[ii_slice[n]]);
      }
      
      partial_target += log_prob;
    }
    return partial_target;
  }
}

data {
  int<lower=0> I;
  int<lower=0> J;
  int<lower=1> N;
  int<lower=1> K;
  array[N] int<lower=1, upper=I> ii;
  array[N] int<lower=1, upper=J> jj;
  array[N] int<lower=1, upper=K> kk;
  array[I] int<lower=1,upper=K> item_type_for_beta;
  array[N] int<lower=0> y;
  array[N] int seq_N;
  int<lower=0> grainsize;
}
parameters {
  vector[J] theta;
  vector[I] beta;
  vector[K] gamma;
  vector<lower=0.01>[I] nu;  //.2
  vector<lower=0>[K] sigma_beta_k;
  
  vector<lower=0>[I] alpha;

  real<lower=0> sigma_alpha;
  real<lower=0> sigma_gamma;
  real mu_gamma;
  real<lower=0> sigma_nu;
}
model {

  sigma_nu ~ cauchy(0,5);
  nu ~ lognormal(0, sigma_nu);
  
  theta ~ normal(0, .3); 
  
  alpha ~ lognormal(0,sigma_alpha);
  sigma_alpha ~ cauchy(0,5);
  
  gamma ~ normal(mu_gamma,sigma_gamma);
  mu_gamma ~ normal(0,5);
  sigma_gamma ~ cauchy(0,5);

  sigma_beta_k ~ cauchy(0,5);
  
  for (i in 1:I) { 
    beta[i] ~ normal(gamma[item_type_for_beta[i]], sigma_beta_k[item_type_for_beta[i]]);
  }
  

  target += reduce_sum(partial_sum, seq_N, grainsize, y, jj, ii, theta,
                       beta, nu, alpha, kk, gamma);
}
generated quantities {
  array[N] real log_lik;

  for (n in 1:N) {
      real log_lambda = nu[ii[n]] * (alpha[ii[n]] * theta[jj[n]] + beta[ii[n]] + gamma[kk[n]]);
      if (log_lambda / nu[ii[n]] > log(1.5) && log_lambda > log(1.5)) {
          log_lik[n] = y[n] * log_lambda - nu[ii[n]] * lgamma(y[n] + 1) - approximation(log_lambda, nu[ii[n]]);
      } else {
          log_lik[n] = y[n] * log_lambda - nu[ii[n]] * lgamma(y[n] + 1) - summation(log_lambda, nu[ii[n]]);
      }
  }
}


')

#library(COMPoissonReg)
nrep <- 10
simulation_conditions <- list(c(100,20), c(100,40), c(100,60), c(200,20),c(200,40),
                              c(200,60), c(300,20), c(300,40), c(300,60))
nrep <-10
simulation_conditions <- list(c(200,60), c(300,20), c(300,40), c(300,60))
#simulation_conditions <- list(c(100,40))
mod2_base <- cmdstan_model(baseline, compile = T, cpp_options = list(stan_threads=T))
mod2_prop <- cmdstan_model(model_proposed, compile = T, cpp_options = list(stan_threads=T))

# Initialize a list to store results of all conditions
#all_results <- list() 

# Loop through each simulation condition
for(condition in simulation_conditions) {
  J <- condition[1]
  I <- condition[2]
  K <- 4
  
  condition_results <- list()
  
  # Loop through each replication
  for(x in 1:nrep) {
    results <- list()
    seed <- 34 * x
    set.seed(seed)
    
    print(paste("Working on CMP condition (J,I), data gen acc extended: (", J, ",", I, ") rep: ",
                x, " started at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    
    start_time <- Sys.time()
    
    N <- I * J
    ii <- rep(1:I, times = J)
    jj <- rep(1:J, each = I)
    kk <- rep(1:K, times = N / K)
    item_type_for_beta <- kk[1:I]
    beta <- numeric(length = I)
    
    seed <- 34*x
    theta <- rnorm(J, 0, .3)
    set.seed(34*x+1)
    beta <- rnorm(I, 3, .5)
    set.seed(34*x+2)
    alpha <- rlnorm(I, 0, 0.3)
    set.seed(34*x+3)
    nu <- rlnorm(I,0,.25)
    set.seed(Sys.time())
    y <- numeric(N)
  #nu[1] <- .3
    for(n in 1:N) {
      mu_n <- exp(nu[ii[n]] * (alpha[ii[n]] * theta[jj[n]] + beta[ii[n]]))
      y[n] <- rcmp(1, mu_n, nu[ii[n]])
      #print(tcmp(log(mu_n),nu[ii[n]]))
    }

    stan_data <- list(I = I,
                      J = J,
                      N = N,
                      K = K,
                      kk = kk,
                      item_type_for_beta = item_type_for_beta,
                      ii = ii,
                      jj = jj,
                      y = y,
                      grainsize = 1,
                      seq_N = 1:N)
    
    fit_baseline <- mod2_base$sample(data = stan_data, chains = 4, threads_per_chain = 13, iter_warmup = 2000, iter_sampling = 2000,
                                     refresh = 50, thin = 1, max_treedepth = 10, save_warmup = TRUE, adapt_delta = .7)
    
    
    log_lik_matrix <- fit_baseline$draws(variables = "log_lik", format = "matrix")
    loo_1 <- waic(log_lik_matrix)
    rm(log_lik_matrix)
    print(paste("Working on Proposed Model condition (J,I): (", J, ",", I, ") rep: ",
                x, " started at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    
    fit_proposed <- mod2_prop$sample(data = stan_data, chains =4, threads_per_chain = 13, iter_warmup = 2000, iter_sampling = 2000,
                                     refresh = 50, thin = 1, max_treedepth = 10, save_warmup = TRUE, adapt_delta = .7)
    
    log_lik_matrix <- fit_proposed$draws(variables = "log_lik", format = "matrix")
    loo_2 <- waic(log_lik_matrix)
    rm(log_lik_matrix)
    
    print(paste("Working on LOO condition (J,I): (", J, ",", I, ") rep: ",
                x, " started at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    
    #loo_1 <- fit_baseline$loo(cores = 8)
    #loo_2 <- fit_proposed$loo(cores = 8)
    
    # Assuming your fit object is named `fit_baseline`

    
    #loo_compare(loo_1,loo_2)

    beta_summary <- fit_baseline$summary("beta")
    theta_summary <- fit_baseline$summary("theta")
    alpha_summary <- fit_baseline$summary("alpha")
    nu_summary <- fit_baseline$summary("nu")
    
    beta_summary_p <- fit_proposed$summary("beta")
    theta_summary_p <- fit_proposed$summary("theta")
    alpha_summary_p <- fit_proposed$summary("alpha")
    nu_summary_p <- fit_proposed$summary("nu")
    
    replication_results <- list(
      seed = seed,
      I = I,
      J = J,
      K = K,
      alpha = alpha,
      theta = theta,
      beta = beta,
      gamma = gamma,
      nu = nu,
      comp <- loo_compare(loo_1, loo_2),
      elpd_diff = comp[, 1],
      se_diff = comp[, 2],
      rmse_beta_rpcm = rmse(beta_summary$mean, beta),
      rmse_beta_proposed = rmse(beta_summary_p$mean, beta),
      rmse_theta_rpcm = rmse(theta_summary$mean, theta),
      rmse_theta_proposed = rmse(theta_summary_p$mean, theta),
      rmse_nu_rpcm = rmse(nu_summary$mean, nu),
      rmse_nu_proposed = rmse(nu_summary_p$mean, nu),
      rmse_alpha_rpcm = rmse(alpha_summary$mean, alpha),
      rmse_alpha_proposed = rmse(alpha_summary_p$mean, alpha),
      bias_theta_rpcm = bias(theta_summary$mean, theta),
      bias_theta_proposed = bias(theta_summary_p$mean, theta),
      bias_nu_rpcm = bias(nu_summary$mean, nu),
      bias_nu_proposed = bias(nu_summary_p$mean, nu),
      bias_beta_rpcm = bias(beta_summary$mean, beta),
      bias_beta_proposed = bias(beta_summary_p$mean, beta),
      bias_alpha_rpcm = bias(alpha_summary$mean, alpha),
      bias_alpha_proposed = bias(alpha_summary_p$mean, alpha)
    )
    #Save the replication results
    # saveRDS(replication_results, file = paste0("C:/Users/halo2/Desktop/GRA Dropbox/Box Sync/Manuscript/Model Comparison Simulation/cmp/data_gen_according_to_baseline_CMP/",
    #                                            J, "_", I, "/results_", x, "_", J, "_", I, ".rds"))
    #saveRDS(replication_results, file = paste0(getwd(), "/data/results_", J, "_", I, "_", x, ".rds"))
    
    saveRDS(replication_results, file = paste0(getwd(), "/data/",
                                               J, "_", I, "/results_", x, "_", J, "_", I, ".rds"))
    
    #condition_results[[x]] <- replication_results
    rm(fit_baseline)
    rm(fit_proposed)
    rm(loo_1)
    rm(loo_2)
    gc()
    end_time <- Sys.time()
    elapsed_time <- end_time - start_time
    print(paste("Rep", x, "took in seconds: ", as.numeric(elapsed_time, units = "secs")))
  }
  #all_results[[length(all_results) + 1]] <- condition_results
}

gamma_draws <- fit_baseline$draws(variables = "nu", inc_warmup = TRUE)

# Extract the samples for gamma[2]
gamma_2_draws <- gamma_draws[, , "nu[3]"]
# Create the traceplot for gamma[2]
mcmc_trace(gamma_2_draws) +
  ggtitle("Traceplot of gamma[2]") +
  theme_minimal()
