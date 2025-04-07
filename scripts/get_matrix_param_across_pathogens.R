library(tidyverse)
library(phylosamp)
source('utils_sens_spec.R')

##################################
## Evolutionary characteristics across pathogens
##################################

## Load dataframe with pathogens characteristics
df_pathogen_char <- read_csv('../input/characteristics_pathogens.csv') %>% 
  mutate(subs_rate_per_year = subs_rate_in_mut_per_site_per_year * genome_length,
         subs_rate_per_day = subs_rate_per_year / 365.25,
         mean_delay_between_mut = 1./subs_rate_per_day)

## Estimate parameters of the Gamma distribution for the generation time
## based on each pathogen's mean and SD of the generation time
df_pathogen_char <- bind_cols(df_pathogen_char, 
                              Reduce('bind_rows', lapply(1:nrow(df_pathogen_char), FUN = function(i_row){
                                vec_gamma_param <- get_gamma_param_from_mean_sd(df_pathogen_char$mean_GT[i_row], df_pathogen_char$sd_GT[i_row])
                                names(vec_gamma_param) <- c('alpha_GT', 'beta_GT')
                                vec_gamma_param
                              })))

## Compute the estimated probability that a transmission event occurs before a mutation one
df_pathogen_char$proba_trans_before_mut <- sapply(1:nrow(df_pathogen_char), FUN = function(i_path){
  get_proba_transm_before_mut(alpha_gen_time = df_pathogen_char$alpha_GT[i_path],
                              beta_gen_time = df_pathogen_char$beta_GT[i_path],
                              mu = df_pathogen_char$subs_rate_per_day[i_path])
})

##################################
## Epidemiological scenarios
##################################
## Load empirical distribution of the number of generations between infected
## individuals (from the phylosamp package) for R of 1.3, 1.5 and 1.7
### R = 1.3
empirical_dist_g_gen_R_1.3 <- get_distribution_generations_func_R(1.3)$vec_gen_dist # Empirical distribution for R = 1.3
g_gen_max_R_1.3 <- get_distribution_generations_func_R(1.3)$n_gens # Maximum number of generations for R = 1.3
empirical_dist_g_gen_R_1.5 <- get_distribution_generations_func_R(1.5)$vec_gen_dist # Empirical distribution for R = 1.5
g_gen_max_R_1.5 <- get_distribution_generations_func_R(1.5)$n_gens # Maximum number of generations for R = 1.5
empirical_dist_g_gen_R_1.7 <- get_distribution_generations_func_R(1.7)$vec_gen_dist # Empirical distribution for R = 1.7
g_gen_max_R_1.7 <- get_distribution_generations_func_R(1.7)$n_gens # Maximum number of generations for R = 1.7

##################################
## Parameter definition for the linkage criterion and the mixing process
##################################
## Vector of omega values (assortativity parameters)
vec_omega <- seq(0.5, 0.99, 0.01)
## Highest genetic distance threshold considered for the linkage criterion
n_mut_max <- 15
## Vector of delta values (genetic distance threshold used to define the linkage criterion)
vec_delta <- 0:n_mut_max

## Dataframe with the scenarios to be explored
df_params <- expand.grid(omega = vec_omega, delta = vec_delta) %>% 
  as_tibble()

##################################
## Compute coefficients associated with the confusion matrix across pathogens
##################################
t0 <- Sys.time()
df_sens_spec_FDR_across_pathogens_R_1.3 <- Reduce('bind_rows', lapply(1:nrow(df_pathogen_char), FUN = function(i_pathogen){
  curr_alpha_GT <- df_pathogen_char$alpha_GT[i_pathogen]
  curr_beta_GT <- df_pathogen_char$beta_GT[i_pathogen]
  curr_pathogen <- df_pathogen_char$pathogen[i_pathogen]
  curr_mu <- df_pathogen_char$subs_rate_per_day[i_pathogen]
  curr_p <- df_pathogen_char$proba_trans_before_mut[i_pathogen]
  
  Reduce('bind_rows', lapply(1:nrow(df_params), FUN = function(i_param){
    
    curr_omega <- df_params$omega[i_param]
    curr_delta <- df_params$delta[i_param]
    curr_lambda <- curr_beta_GT * (exp(- 1./curr_alpha_GT * log(curr_omega)) - 1)
    
    tibble(delta = curr_delta,
           omega = curr_omega,
           lambda = curr_lambda,
           pathogen = curr_pathogen,
           proba_trans_before_mut = curr_p,
           sensitivity = get_sensitivity_delta(curr_delta, curr_mu, curr_lambda, curr_alpha_GT, curr_beta_GT, 
                                               empirical_dist_g_gen_R_1.3, g_gen_max_R_1.3),
           specificity = get_specificity_delta(curr_delta, curr_mu, curr_lambda, curr_alpha_GT, curr_beta_GT, 
                                               empirical_dist_g_gen_R_1.3, g_gen_max_R_1.3),
           ppv = get_ppv_delta(curr_delta, curr_mu, curr_lambda, curr_alpha_GT, curr_beta_GT, 
                               empirical_dist_g_gen_R_1.3, g_gen_max_R_1.3),
           f1_score = get_f1_score_delta(curr_delta, curr_mu, curr_lambda, curr_alpha_GT, curr_beta_GT, 
                                         empirical_dist_g_gen_R_1.3, g_gen_max_R_1.3))
  }))
})) %>% 
  mutate(fdr = 1. - ppv, R = 1.3)
df_sens_spec_FDR_across_pathogens_R_1.5 <- Reduce('bind_rows', lapply(1:nrow(df_pathogen_char), FUN = function(i_pathogen){
  curr_alpha_GT <- df_pathogen_char$alpha_GT[i_pathogen]
  curr_beta_GT <- df_pathogen_char$beta_GT[i_pathogen]
  curr_pathogen <- df_pathogen_char$pathogen[i_pathogen]
  curr_mu <- df_pathogen_char$subs_rate_per_day[i_pathogen]
  curr_p <- df_pathogen_char$proba_trans_before_mut[i_pathogen]
  
  Reduce('bind_rows', lapply(1:nrow(df_params), FUN = function(i_param){
    
    curr_omega <- df_params$omega[i_param]
    curr_delta <- df_params$delta[i_param]
    curr_lambda <- curr_beta_GT * (exp(- 1./curr_alpha_GT * log(curr_omega)) - 1)
    
    tibble(delta = curr_delta,
           omega = curr_omega,
           lambda = curr_lambda,
           pathogen = curr_pathogen,
           proba_trans_before_mut = curr_p,
           sensitivity = get_sensitivity_delta(curr_delta, curr_mu, curr_lambda, curr_alpha_GT, curr_beta_GT, 
                                               empirical_dist_g_gen_R_1.5, g_gen_max_R_1.5),
           specificity = get_specificity_delta(curr_delta, curr_mu, curr_lambda, curr_alpha_GT, curr_beta_GT, 
                                               empirical_dist_g_gen_R_1.5, g_gen_max_R_1.5),
           ppv = get_ppv_delta(curr_delta, curr_mu, curr_lambda, curr_alpha_GT, curr_beta_GT, 
                               empirical_dist_g_gen_R_1.5, g_gen_max_R_1.5),
           f1_score = get_f1_score_delta(curr_delta, curr_mu, curr_lambda, curr_alpha_GT, curr_beta_GT, 
                                         empirical_dist_g_gen_R_1.5, g_gen_max_R_1.5))
  }))
})) %>% 
  mutate(fdr = 1. - ppv, R = 1.5)
df_sens_spec_FDR_across_pathogens_R_1.7 <- Reduce('bind_rows', lapply(1:nrow(df_pathogen_char), FUN = function(i_pathogen){
  curr_alpha_GT <- df_pathogen_char$alpha_GT[i_pathogen]
  curr_beta_GT <- df_pathogen_char$beta_GT[i_pathogen]
  curr_pathogen <- df_pathogen_char$pathogen[i_pathogen]
  curr_mu <- df_pathogen_char$subs_rate_per_day[i_pathogen]
  curr_p <- df_pathogen_char$proba_trans_before_mut[i_pathogen]
  
  Reduce('bind_rows', lapply(1:nrow(df_params), FUN = function(i_param){
    
    curr_omega <- df_params$omega[i_param]
    curr_delta <- df_params$delta[i_param]
    curr_lambda <- curr_beta_GT * (exp(- 1./curr_alpha_GT * log(curr_omega)) - 1)
    
    tibble(delta = curr_delta,
           omega = curr_omega,
           lambda = curr_lambda,
           pathogen = curr_pathogen,
           proba_trans_before_mut = curr_p,
           sensitivity = get_sensitivity_delta(curr_delta, curr_mu, curr_lambda, curr_alpha_GT, curr_beta_GT, 
                                               empirical_dist_g_gen_R_1.7, g_gen_max_R_1.7),
           specificity = get_specificity_delta(curr_delta, curr_mu, curr_lambda, curr_alpha_GT, curr_beta_GT, 
                                               empirical_dist_g_gen_R_1.7, g_gen_max_R_1.7),
           ppv = get_ppv_delta(curr_delta, curr_mu, curr_lambda, curr_alpha_GT, curr_beta_GT, 
                               empirical_dist_g_gen_R_1.7, g_gen_max_R_1.7),
           f1_score = get_f1_score_delta(curr_delta, curr_mu, curr_lambda, curr_alpha_GT, curr_beta_GT, 
                                         empirical_dist_g_gen_R_1.7, g_gen_max_R_1.7))
  }))
})) %>% 
  mutate(fdr = 1. - ppv, R = 1.7)
t1 <- Sys.time()
print(t1 - t0)

## Save results
write_csv(df_sens_spec_FDR_across_pathogens_R_1.3, '../results/df_sens_spec_FDR_F1_across_pathogens_R_1.3.csv')
write_csv(df_sens_spec_FDR_across_pathogens_R_1.5, '../results/df_sens_spec_FDR_F1_across_pathogens_R_1.5.csv')
write_csv(df_sens_spec_FDR_across_pathogens_R_1.7, '../results/df_sens_spec_FDR_F1_across_pathogens_R_1.7.csv')
