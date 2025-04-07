## Function to get Alpha and Beta parameters for a Gamma distribution
## from its mean and standard deviation
get_gamma_param_from_mean_sd <- function(mean_gamma, sd_gamma){
  alpha_gamma <- (mean_gamma^2) / (sd_gamma)^2
  beta_gamma <- mean_gamma / (sd_gamma)^2
  
  return(c(
    'alpha' = alpha_gamma, 
    'beta' = beta_gamma
  ))
}

## Function to compute the probability of the number of mutations conditional on the number of generations
get_proba_m_mutations_conditional_g_generations <- function(vec_nb_mutations_away_obs, 
                                                            mu, alpha_gen_time, beta_gen_time, 
                                                            g_gen_away){
  
  ## vec_nb_mutations_away_obs: Vector of number of mutations for which this probability is computed
  ## mu: pathogen mutation rate
  ## alpha_gen_time: shape of the generation time distribution (assumed Gamma distributed)
  ## beta_gen_time: scale of the generation time distribution (assumed Gamma distributed)
  ## g_gen_away: number of generations separating two individuals (on which we are conditioning)
  
  ## Ratio used in the parametrisation of the negative binomial distribution
  tmp_ratio <- beta_gen_time / (mu + beta_gen_time)
  
  ## Probability mass function of the number of mutation conditional on the number o generations
  proba_m_mut_cond_g_gen <- dnbinom(x = vec_nb_mutations_away_obs,  
                                    size = alpha_gen_time*g_gen_away, prob = tmp_ratio)
  
  return(proba_m_mut_cond_g_gen)
}
get_log_proba_m_mutations_conditional_g_generations <- function(vec_nb_mutations_away_obs, 
                                                                mu, alpha_gen_time, beta_gen_time, 
                                                                g_gen_away){
  
  ## vec_nb_mutations_away_obs: Vector of number of mutations for which this probability is computed
  ## mu: pathogen mutation rate
  ## alpha_gen_time: shape of the generation time distribution (assumed Gamma distributed)
  ## beta_gen_time: scale of the generation time distribution (assumed Gamma distributed)
  ## g_gen_away: number of generations separating two individuals (on which we are conditioning)
  
  ## Ratio used in the parametrisation of the negative binomial distribution
  tmp_ratio <- beta_gen_time / (mu + beta_gen_time)
  
  ## Probability mass function of the number of mutation conditional on the number o generations
  log_proba_m_mut_cond_g_gen <- dnbinom(x = vec_nb_mutations_away_obs,  
                                        size = alpha_gen_time*g_gen_away, 
                                        prob = tmp_ratio, 
                                        log = T)
  
  return(log_proba_m_mut_cond_g_gen)
}

## Return a matrix of the probability of the number of mutations conditional on the number of generations
get_mat_proba_m_mut_conditional_g_gen <- function(n_mut_max, g_gen_max, mu, alpha_gen_time, beta_gen_time){
  
  vec_nb_mutations_away_obs <- 0:n_mut_max
  vec_g_gen <- 1:g_gen_max
  
  ## Matrix with mutations in rows (0 to n_mut_max) and generations in columns (1 to g_gen_max)
  ## where coefficients indicate the probability that the number of mutations is equal to i the row id
  ## conditional on the number of generations being equal to j (column id)
  mat_proba_m_mut_conditional_g_gen <- sapply(vec_g_gen, FUN = function(curr_g_gen_away){
    get_proba_m_mutations_conditional_g_generations(vec_nb_mutations_away_obs = vec_nb_mutations_away_obs, 
                                                    mu = mu, 
                                                    alpha_gen_time = alpha_gen_time, 
                                                    beta_gen_time = beta_gen_time, 
                                                    g_gen_away = curr_g_gen_away)
    
    
  })
  
  
  return(mat_proba_m_mut_conditional_g_gen)
}

get_log_mat_proba_m_mut_conditional_g_gen <- function(n_mut_max, g_gen_max, mu, alpha_gen_time, beta_gen_time){
  
  vec_nb_mutations_away_obs <- 0:n_mut_max
  vec_g_gen <- 1:g_gen_max
  
  ## Matrix with mutations in rows (0 to n_mut_max) and generations in columns (1 to g_gen_max)
  ## where coefficients indicate the probability that the number of mutations is equal to i the row id
  ## conditional on the number of generations being equal to j (column id)
  mat_log_proba_m_mut_conditional_g_gen <- sapply(vec_g_gen, FUN = function(curr_g_gen_away){
    get_log_proba_m_mutations_conditional_g_generations(vec_nb_mutations_away_obs = vec_nb_mutations_away_obs, 
                                                        mu = mu, 
                                                        alpha_gen_time = alpha_gen_time, 
                                                        beta_gen_time = beta_gen_time, 
                                                        g_gen_away = curr_g_gen_away)
    
    
  })
  
  
  return(mat_log_proba_m_mut_conditional_g_gen)
}


## Function to compute the probability of the number of jumps conditional on the number of generations
get_proba_j_jumps_conditional_g_generations <- function(vec_nb_jumps_away_obs, 
                                                        lambda, alpha_gen_time, beta_gen_time, 
                                                        g_gen_away){
  
  ## vec_nb_jumps_away_obs: Vector of number of jumps for which this probability is computed
  ## lambda: pathogen jump rate (migration rate between population groups)
  ## alpha_gen_time: shape of the generation time distribution (assumed Gamma distributed)
  ## beta_gen_time: scale of the generation time distribution (assumed Gamma distributed)
  ## g_gen_away: number of generations separating two individuals (on which we are conditioning)
  
  ## Ratio used in the parametrisation of the negative binomial distribution
  tmp_ratio <- beta_gen_time / (lambda + beta_gen_time)
  
  ## Probability mass function of the number of mutation conditional on the number o generations
  proba_j_jumps_cond_g_gen <- dnbinom(x = vec_nb_jumps_away_obs,  
                                      size = alpha_gen_time*g_gen_away, prob = tmp_ratio)
  
  return(proba_j_jumps_cond_g_gen)
}

get_log_proba_j_jumps_conditional_g_generations <- function(vec_nb_jumps_away_obs, 
                                                            lambda, alpha_gen_time, beta_gen_time, 
                                                            g_gen_away){
  
  ## vec_nb_jumps_away_obs: Vector of number of jumps for which this probability is computed
  ## lambda: pathogen jump rate (migration rate between population groups)
  ## alpha_gen_time: shape of the generation time distribution (assumed Gamma distributed)
  ## beta_gen_time: scale of the generation time distribution (assumed Gamma distributed)
  ## g_gen_away: number of generations separating two individuals (on which we are conditioning)
  
  ## Ratio used in the parametrisation of the negative binomial distribution
  tmp_ratio <- beta_gen_time / (lambda + beta_gen_time)
  
  ## Probability mass function of the number of mutation conditional on the number o generations
  log_proba_j_jumps_cond_g_gen <- dnbinom(x = vec_nb_jumps_away_obs,  
                                          size = alpha_gen_time*g_gen_away, 
                                          prob = tmp_ratio,
                                          log = T)
  
  return(log_proba_j_jumps_cond_g_gen)
}
## Return a matrix
get_mat_proba_j_jumps_conditional_g_gen <- function(n_jump_max, g_gen_max, lambda, alpha_gen_time, beta_gen_time){
  
  vec_nb_jumps_away_obs <- 0:n_jump_max
  vec_g_gen <- 1:g_gen_max
  
  ## Matrix with mutations in rows (0 to n_jump_max) and generations in columns (1 to g_gen_max)
  ## where coefficients indicate the probability that the number of jumps is equal to i the row id
  ## conditional on the number of generations being equal to j (column id)
  mat_proba_j_jumps_conditional_g_gen <- sapply(vec_g_gen, FUN = function(curr_g_gen_away){
    get_proba_j_jumps_conditional_g_generations(vec_nb_jumps_away_obs = vec_nb_jumps_away_obs, 
                                                lambda = lambda, 
                                                alpha_gen_time = alpha_gen_time, 
                                                beta_gen_time = beta_gen_time, 
                                                g_gen_away = curr_g_gen_away)
  })
  
  return(mat_proba_j_jumps_conditional_g_gen)
}

get_mat_log_proba_j_jumps_conditional_g_gen <- function(n_jump_max, g_gen_max, lambda, alpha_gen_time, beta_gen_time){
  
  vec_nb_jumps_away_obs <- 0:n_jump_max
  vec_g_gen <- 1:g_gen_max
  
  ## Matrix with mutations in rows (0 to n_jump_max) and generations in columns (1 to g_gen_max)
  ## where coefficients indicate the probability that the number of jumps is equal to i the row id
  ## conditional on the number of generations being equal to j (column id)
  mat_log_proba_j_jumps_conditional_g_gen <- sapply(vec_g_gen, FUN = function(curr_g_gen_away){
    get_log_proba_j_jumps_conditional_g_generations(vec_nb_jumps_away_obs = vec_nb_jumps_away_obs, 
                                                    lambda = lambda, 
                                                    alpha_gen_time = alpha_gen_time, 
                                                    beta_gen_time = beta_gen_time, 
                                                    g_gen_away = curr_g_gen_away)
  })
  
  return(mat_log_proba_j_jumps_conditional_g_gen)
}


## Function to compute the probability of the number of jumps conditional on the number of mutations
get_proba_j_jumps_conditional_m_mutations <- function(nb_mutations_away_obs, nb_jumps_away, 
                                                      mu, lambda, alpha_gen_time, beta_gen_time, 
                                                      empirical_dist_g_gen, g_gen_max){
  
  vec_nb_mutations_away_obs <- nb_mutations_away_obs:nb_mutations_away_obs
  vec_nb_jumps_away_obs <- nb_jumps_away:nb_jumps_away
  vec_g_gen <- 1:g_gen_max
  
  ## Vector where coefficients indicate the probability that the number of mutations is equal to nb_mutations_away_obs
  ## conditional on the number of generation being equal to j (vector id)
  vec_log_proba_given_mut_conditional_g_gen <- sapply(vec_g_gen, FUN = function(curr_g_gen_away){
    get_log_proba_m_mutations_conditional_g_generations(vec_nb_mutations_away_obs = vec_nb_mutations_away_obs, 
                                                        mu = mu, 
                                                        alpha_gen_time = alpha_gen_time, 
                                                        beta_gen_time = beta_gen_time, 
                                                        g_gen_away = curr_g_gen_away)  
  })
  
  ## Vector where coefficients indicate the probability that the number of jumps is equal to nb_jumps_away
  ## conditional on the number of generation being equal to j (vector id)
  vec_log_proba_given_jump_conditional_g_gen <- sapply(vec_g_gen, FUN = function(curr_g_gen_away){
    get_log_proba_j_jumps_conditional_g_generations(vec_nb_jumps_away_obs = vec_nb_jumps_away_obs, 
                                                    lambda = lambda, 
                                                    alpha_gen_time = alpha_gen_time, 
                                                    beta_gen_time = beta_gen_time, 
                                                    g_gen_away = curr_g_gen_away)
  })
  
  vec_log_numerator <- sapply(vec_g_gen, FUN = function(curr_g_gen_away){
    vec_log_proba_given_jump_conditional_g_gen[curr_g_gen_away] +
      vec_log_proba_given_mut_conditional_g_gen[curr_g_gen_away] + 
      log(empirical_dist_g_gen[curr_g_gen_away])
  })
  
  vec_log_denominator <- sapply(vec_g_gen, FUN = function(curr_g_gen_away){
    vec_log_proba_given_mut_conditional_g_gen[curr_g_gen_away] + 
      log(empirical_dist_g_gen[curr_g_gen_away])
  })
  
  return(sum(exp(vec_log_numerator)) / sum(exp(vec_log_denominator)))
}

## Function to get the unconditional probability that M = m
get_proba_m_mutations <- function(nb_mutations_away_obs, mu,
                                  alpha_gen_time, beta_gen_time,
                                  empirical_dist_g_gen, g_gen_max){
  
  vec_nb_mutations_away_obs <- nb_mutations_away_obs:nb_mutations_away_obs
  vec_g_gen <- 1:g_gen_max
  
  
  ## Vector where coefficients indicate the probability that the number of mutations is equal to nb_mutations_away_obs
  ## conditional on the number of generation being equal to j (vector id)
  vec_log_proba_given_mut_conditional_g_gen <- sapply(vec_g_gen, FUN = function(curr_g_gen_away){
    get_log_proba_m_mutations_conditional_g_generations(vec_nb_mutations_away_obs = vec_nb_mutations_away_obs, 
                                                        mu = mu, 
                                                        alpha_gen_time = alpha_gen_time, 
                                                        beta_gen_time = beta_gen_time, 
                                                        g_gen_away = curr_g_gen_away)  
  })
  
  vec_log_for_sum <- sapply(vec_g_gen, FUN = function(curr_g_gen_away){
    vec_log_proba_given_mut_conditional_g_gen[curr_g_gen_away] + 
      log(empirical_dist_g_gen[curr_g_gen_away])
  })
  return(sum(exp(vec_log_for_sum)))
}

## Function to get the unconditional probability that J = j
get_proba_j_jumps <- function(nb_jumps_away, lambda,
                              alpha_gen_time, beta_gen_time,
                              empirical_dist_g_gen, g_gen_max){
  
  vec_nb_jumps_away_obs <- nb_jumps_away:nb_jumps_away
  vec_g_gen <- 1:g_gen_max
  
  ## Vector where coefficients indicate the probability that the number of jumps is equal to nb_jumps_away
  ## conditional on the number of generation being equal to j (vector id)
  vec_log_proba_given_jump_conditional_g_gen <- sapply(vec_g_gen, FUN = function(curr_g_gen_away){
    get_log_proba_j_jumps_conditional_g_generations(vec_nb_jumps_away_obs = vec_nb_jumps_away_obs, 
                                                    lambda = lambda, 
                                                    alpha_gen_time = alpha_gen_time, 
                                                    beta_gen_time = beta_gen_time, 
                                                    g_gen_away = curr_g_gen_away)
  })
  
  vec_log_for_sum <- sapply(vec_g_gen, FUN = function(curr_g_gen_away){
    vec_log_proba_given_jump_conditional_g_gen[curr_g_gen_away] + 
      log(empirical_dist_g_gen[curr_g_gen_away])
  })
  return(sum(exp(vec_log_for_sum)))
}

#### Functions to compute sensitivity
get_eta_prime_delta <- function(delta, mu, lambda, alpha_gen_time, beta_gen_time, 
                                empirical_dist_g_gen, g_gen_max){

  proba_0_jumps_conditional_delta_mut <- get_proba_j_jumps_conditional_m_mutations(nb_mutations_away_obs = delta,
                                                                                   nb_jumps_away = 0, 
                                                                                   mu, lambda, alpha_gen_time, beta_gen_time, 
                                                                                   empirical_dist_g_gen, g_gen_max)
  
  proba_1_jump_conditional_delta_mut <- get_proba_j_jumps_conditional_m_mutations(nb_mutations_away_obs = delta,
                                                                                  nb_jumps_away = 1, 
                                                                                  mu, lambda, alpha_gen_time, beta_gen_time, 
                                                                                  empirical_dist_g_gen, g_gen_max)
  
  proba_delta_mutations <- get_proba_m_mutations(nb_mutations_away_obs = delta,
                                                 mu,
                                                 alpha_gen_time, beta_gen_time,
                                                 empirical_dist_g_gen, g_gen_max)
  
  proba_0_jump <- get_proba_j_jumps(nb_jumps_away = 0,
                                    lambda,
                                    alpha_gen_time, beta_gen_time,
                                    empirical_dist_g_gen, g_gen_max)
  
  proba_1_jump <- get_proba_j_jumps(nb_jumps_away = 1,
                                    lambda,
                                    alpha_gen_time, beta_gen_time,
                                    empirical_dist_g_gen, g_gen_max)
  
  log_eta_prime <- 
    log(proba_0_jumps_conditional_delta_mut + proba_1_jump_conditional_delta_mut) +
    log(proba_delta_mutations) -
    log(proba_0_jump + proba_1_jump)
  
  return(exp(log_eta_prime))
}

## Sensitivity corresponding to a threshold M <= Delta
get_sensitivity_delta <- function(delta, mu, lambda, alpha_gen_time, beta_gen_time, 
                                  empirical_dist_g_gen, g_gen_max){
  vec_delta_to_sum <- 0:delta
  vec_eta_primes <- sapply(vec_delta_to_sum, FUN = function(curr_delta){
    get_eta_prime_delta(curr_delta, mu, lambda, alpha_gen_time, beta_gen_time, 
                        empirical_dist_g_gen, g_gen_max)
  })
  return(sum(vec_eta_primes))
}

## Sensitivity corresponding to a threshold M = Delta
get_sensitivity_exactly_delta <- function(delta, mu, lambda, alpha_gen_time, beta_gen_time, 
                                  empirical_dist_g_gen, g_gen_max){
  vec_delta_to_sum <- delta:delta
  vec_eta_primes <- sapply(vec_delta_to_sum, FUN = function(curr_delta){
    get_eta_prime_delta(curr_delta, mu, lambda, alpha_gen_time, beta_gen_time, 
                        empirical_dist_g_gen, g_gen_max)
  })
  return(sum(vec_eta_primes))
}
  
#### Functions to compute specificity
get_chi_prime_delta <- function(delta, mu, lambda, alpha_gen_time, beta_gen_time, 
                                empirical_dist_g_gen, g_gen_max){
  proba_0_jumps_conditional_delta_mut <- get_proba_j_jumps_conditional_m_mutations(nb_mutations_away_obs = delta,
                                                                                   nb_jumps_away = 0, 
                                                                                   mu, lambda, alpha_gen_time, beta_gen_time, 
                                                                                   empirical_dist_g_gen, g_gen_max)
  
  proba_1_jump_conditional_delta_mut <- get_proba_j_jumps_conditional_m_mutations(nb_mutations_away_obs = delta,
                                                                                  nb_jumps_away = 1, 
                                                                                  mu, lambda, alpha_gen_time, beta_gen_time, 
                                                                                  empirical_dist_g_gen, g_gen_max)
  
  proba_delta_mutations <- get_proba_m_mutations(nb_mutations_away_obs = delta,
                                                 mu,
                                                 alpha_gen_time, beta_gen_time,
                                                 empirical_dist_g_gen, g_gen_max)
  
  proba_0_jump <- get_proba_j_jumps(nb_jumps_away = 0,
                                    lambda,
                                    alpha_gen_time, beta_gen_time,
                                    empirical_dist_g_gen, g_gen_max)
  
  proba_1_jump <- get_proba_j_jumps(nb_jumps_away = 1,
                                    lambda,
                                    alpha_gen_time, beta_gen_time,
                                    empirical_dist_g_gen, g_gen_max)
  
  
  log_chi_prime <- 
    log(1 - proba_0_jumps_conditional_delta_mut - proba_1_jump_conditional_delta_mut) +
    log(proba_delta_mutations) -
    log(1 - proba_0_jump - proba_1_jump)
  
  return(exp(log_chi_prime))
}

## Specificity corresponding to a threshold M <= Delta
get_specificity_delta <- function(delta, mu, lambda, alpha_gen_time, beta_gen_time, 
                                  empirical_dist_g_gen, g_gen_max){
  vec_delta_to_sum <- 0:delta
  vec_chi_primes <- sapply(vec_delta_to_sum, FUN = function(curr_delta){
    get_chi_prime_delta(curr_delta, mu, lambda, alpha_gen_time, beta_gen_time, 
                        empirical_dist_g_gen, g_gen_max)
  })
  return(1 - sum(vec_chi_primes))
}

#### Function to compute PPV
## PPV corresponding to a threshold M <= Delta
get_ppv_delta <- function(delta, mu, lambda, alpha_gen_time, beta_gen_time, 
                          empirical_dist_g_gen, g_gen_max){
  log_sensitivity_delta <- log(
    get_sensitivity_delta(delta, mu, lambda, alpha_gen_time, beta_gen_time, 
                          empirical_dist_g_gen, g_gen_max)
    )
  
  proba_0_jump <- get_proba_j_jumps(nb_jumps_away = 0,
                                     lambda,
                                     alpha_gen_time, beta_gen_time,
                                     empirical_dist_g_gen, g_gen_max)
  
  proba_1_jump <- get_proba_j_jumps(nb_jumps_away = 1,
                                    lambda,
                                    alpha_gen_time, beta_gen_time,
                                    empirical_dist_g_gen, g_gen_max)
  
  
  proba_below_delta_mutations <- sum(sapply(0:delta, FUN = function(curr_mut){
    get_proba_m_mutations(nb_mutations_away_obs = curr_mut,
                          mu,
                          alpha_gen_time, beta_gen_time,
                          empirical_dist_g_gen, g_gen_max)
  }))
  
  vec_log_ppv <- log_sensitivity_delta + log(proba_0_jump + proba_1_jump) - log(proba_below_delta_mutations)
  return(exp(vec_log_ppv))
}

## PPV corresponding to a threshold M <= Delta
get_ppv_exactly_delta <- function(delta, mu, lambda, alpha_gen_time, beta_gen_time, 
                          empirical_dist_g_gen, g_gen_max){
  log_sensitivity_delta <- log(
    get_sensitivity_exactly_delta(delta, mu, lambda, alpha_gen_time, beta_gen_time, 
                                  empirical_dist_g_gen, g_gen_max)
  )
  
  proba_0_jump <- get_proba_j_jumps(nb_jumps_away = 0,
                                    lambda,
                                    alpha_gen_time, beta_gen_time,
                                    empirical_dist_g_gen, g_gen_max)
  
  proba_1_jump <- get_proba_j_jumps(nb_jumps_away = 1,
                                    lambda,
                                    alpha_gen_time, beta_gen_time,
                                    empirical_dist_g_gen, g_gen_max)
  
  
  proba_below_delta_mutations <- sum(sapply(delta:delta, FUN = function(curr_mut){
    get_proba_m_mutations(nb_mutations_away_obs = curr_mut,
                          mu,
                          alpha_gen_time, beta_gen_time,
                          empirical_dist_g_gen, g_gen_max)
  }))
  
  vec_log_ppv <- log_sensitivity_delta + log(proba_0_jump + proba_1_jump) - log(proba_below_delta_mutations)
  return(exp(vec_log_ppv))
}

## PPV that would be associated with two random sequences from the dataset
get_ppv_random <- function(mu, lambda, alpha_gen_time, beta_gen_time, 
                           empirical_dist_g_gen, g_gen_max){
  
  proba_0_jump <- get_proba_j_jumps(nb_jumps_away = 0,
                                    lambda,
                                    alpha_gen_time, beta_gen_time,
                                    empirical_dist_g_gen, g_gen_max)
  
  proba_1_jump <- get_proba_j_jumps(nb_jumps_away = 1,
                                    lambda,
                                    alpha_gen_time, beta_gen_time,
                                    empirical_dist_g_gen, g_gen_max)
  return(proba_0_jump + proba_1_jump)
}

#### Function to compute F1 score corresponding to a threshold M <= Delta
get_f1_score_delta <- function(delta, mu, lambda, alpha_gen_time, beta_gen_time, 
                               empirical_dist_g_gen, g_gen_max){
  
  ppv <- get_ppv_delta(delta, mu, lambda, alpha_gen_time, beta_gen_time, 
                       empirical_dist_g_gen, g_gen_max)
  
  sensitivity <- get_sensitivity_delta(delta, mu, lambda, alpha_gen_time, beta_gen_time, 
                                       empirical_dist_g_gen, g_gen_max)
  
  log_f1_score <- log(2) + log(sensitivity) + log(ppv) - log(sensitivity + ppv)
  return(exp(log_f1_score))
}

### Function to get the probability that transmission occurs before mutation
get_proba_transm_before_mut <- function(alpha_gen_time, beta_gen_time, mu){
  return(get_proba_m_mutations_conditional_g_generations(0, mu, alpha_gen_time, beta_gen_time, g_gen_away = 1))
}
#### Function to get the probability that transmission occurs before migration
get_proba_stays_within_group <- function(alpha_gen_time, beta_gen_time, lambda){
  return(get_proba_j_jumps_conditional_g_generations(0, lambda, alpha_gen_time, beta_gen_time, g_gen_away = 1))
}

## Function to get the distribution of the number of generation separating two individuals
## from the phylosamp package (as a function of the reproduction number R)
get_distribution_generations_func_R <- function(curr_R){
  
  gen_dist <- genDistSim %>% filter(R == curr_R)
  n_gens <- gen_dist[2] %>% as.numeric()
  vec_gen_dist <- gen_dist[3:(n_gens*2 + 2)] %>% as.numeric()
  
  return(
    list(n_gens = 2 * n_gens, vec_gen_dist = vec_gen_dist)
  )
}

