library(tidyverse)
library(phylosamp)
library(doParallel)
source('utils_sens_spec.R')

##################################
## Epidemiological scenarios
##################################
## Load empirical distribution of the number of generations between infected
## individuals (from the phylosamp package) for R of 1.3
empirical_dist_g_gen_R_1.3 <- get_distribution_generations_func_R(1.3)$vec_gen_dist # Empirical distribution for R = 1.3
g_gen_max_R_1.3 <- get_distribution_generations_func_R(1.3)$n_gens # Maximum number of generations for R = 1.3

## Load dataframe with pathogens characteristics
df_pathogen_char <- read_csv('../input/characteristics_pathogens.csv') %>% 
  mutate(subs_rate_per_year = subs_rate_in_mut_per_site_per_year * genome_length,
         subs_rate_per_day = subs_rate_per_year / 365.25,
         mean_delay_between_mut = 1./subs_rate_per_day)

##################################
## Generation time characteristics used in the simulations
##################################
mean_GT <- 4.9 # Mean of the generation time distribution
sd_GT <- 4.8 # SD of the generation time distribution

## Corresponding parameters assuming a Gamma distribution
alpha_GT <- as.numeric(get_gamma_param_from_mean_sd(mean_gamma = mean_GT, sd_gamma = sd_GT)[1])
beta_GT <- as.numeric(get_gamma_param_from_mean_sd(mean_gamma = mean_GT, sd_gamma = sd_GT)[2])

##################################
## Define functions used to compute lambda (migration rate) from omega (assortativity parameter)
## and mu (mutation rate) from p (probability that transmission occurs before mutation)
##################################
get_mu_from_proba_trans_before_mut <- function(alpha_gen_time, beta_gen_time, proba_trans_before_mut){
  mu <- beta_gen_time * (exp(- log(proba_trans_before_mut) / alpha_gen_time) - 1.)
  return(as.numeric(mu))
}

get_lambda_from_omega <- function(alpha_gen_time, beta_gen_time, omega){
  lambda <- beta_gen_time * (exp(- log(omega) / alpha_gen_time) - 1.)
  return(as.numeric(lambda))
}

##################################
## Define the parameter space to be explored
##################################
df_params <- expand.grid(delta = 0:15, # Genetic distance threshold
                         omega = seq(0.01, 0.99, 0.01), # Assortativity parameter
                         proba_trans_before_mut = seq(0.01, 0.99, 0.01) # Probability that transmission occurs before mutation
                         ) %>% 
  group_by(omega, proba_trans_before_mut, delta) %>% 
  mutate(mu = get_mu_from_proba_trans_before_mut(alpha_GT, beta_GT, proba_trans_before_mut), # Corresponding mutation rate
         lambda = get_lambda_from_omega(alpha_GT, beta_GT, omega) # Corresponding migration rate
         )

##################################
## Compute PPV and F1 across parameter space
##################################
cl <- parallel::makeCluster(6)
registerDoParallel(cl)
t0 <- Sys.time()
## PPV 
df_params$ppv <- foreach(i_param = 1:nrow(df_params), .combine = 'c') %dopar% {
  curr_mu <- df_params$mu[i_param]
  curr_lambda <- df_params$lambda[i_param]
  curr_delta <- df_params$delta[i_param]
  
  empirical_dist_g_gen <- empirical_dist_g_gen_R_1.3
  g_gen_max <- g_gen_max_R_1.3
  
  get_ppv_delta(curr_delta, curr_mu, curr_lambda, alpha_GT, beta_GT, 
                empirical_dist_g_gen, g_gen_max)
}

## PPV corresponding to a threshold M = Delta instead of M <= Delta
df_params$ppv_exactly_delta <- foreach(i_param = 1:nrow(df_params), .combine = 'c') %dopar% {
  curr_mu <- df_params$mu[i_param]
  curr_lambda <- df_params$lambda[i_param]
  curr_delta <- df_params$delta[i_param]
  
  empirical_dist_g_gen <- empirical_dist_g_gen_R_1.3
  g_gen_max <- g_gen_max_R_1.3
  
  get_ppv_exactly_delta(curr_delta, curr_mu, curr_lambda, alpha_GT, beta_GT, 
                        empirical_dist_g_gen, g_gen_max)
}

## F1-score
df_params$f1_score <- foreach(i_param = 1:nrow(df_params), .combine = 'c') %dopar% {
  curr_mu <- df_params$mu[i_param]
  curr_lambda <- df_params$lambda[i_param]
  curr_delta <- df_params$delta[i_param]
  empirical_dist_g_gen <- empirical_dist_g_gen_R_1.3
  g_gen_max <- g_gen_max_R_1.3
  
  get_f1_score_delta(curr_delta, curr_mu, curr_lambda, alpha_GT, beta_GT, 
                     empirical_dist_g_gen, g_gen_max)
}

## Sensitivity
df_params$sensitivity <- foreach(i_param = 1:nrow(df_params), .combine = 'c') %dopar% {
  curr_mu <- df_params$mu[i_param]
  curr_lambda <- df_params$lambda[i_param]
  curr_delta <- df_params$delta[i_param]
  
  empirical_dist_g_gen <- empirical_dist_g_gen_R_1.3
  g_gen_max <- g_gen_max_R_1.3
  
  get_sensitivity_delta(curr_delta, curr_mu, curr_lambda, alpha_GT, beta_GT, 
                        empirical_dist_g_gen, g_gen_max)
}

## Specificity
df_params$specificity <- foreach(i_param = 1:nrow(df_params), .combine = 'c') %dopar% {
  curr_mu <- df_params$mu[i_param]
  curr_lambda <- df_params$lambda[i_param]
  curr_delta <- df_params$delta[i_param]
  
  empirical_dist_g_gen <- empirical_dist_g_gen_R_1.3
  g_gen_max <- g_gen_max_R_1.3
  
  get_specificity_delta(curr_delta, curr_mu, curr_lambda, alpha_GT, beta_GT, 
                        empirical_dist_g_gen, g_gen_max)
}

t1 <- Sys.time()
stopCluster(cl)

write_csv(df_params, '../results/df_sens_spec_PPV_F1_across_parameter_space_R_1.3.csv')