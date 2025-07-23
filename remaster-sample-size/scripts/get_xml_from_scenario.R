source('scripts/utils_remaster.R')

args <- commandArgs(trailingOnly = T)
proba_trans_before_mut <- as.numeric(args[1])
assortativity <- as.numeric(args[2])
p_seq <- as.numeric(args[3])
xml_file_path <- as.character(args[4])
dir_results <- as.character(args[5])
output_migration_rate_mat <- as.character(args[6])

## Define generation time parameters
rate_out_of_E <- 0.33 # Rate out which individuals exit the E compartment
rate_out_of_I <- rate_out_of_E # Rate out which individuals exit the I compartment

## Corresponding parameters for a gamma distribution
alpha_GT <- 2
beta_GT <- 1. / rate_out_of_E

get_mu_from_proba_trans_before_mut <- function(alpha_GT, beta_GT, proba_trans_before_mut){
  return(
    beta_GT * (exp(- 1/ alpha_GT * log(proba_trans_before_mut)) - 1)
  )
}

## Compute mutation rate corresponding to probability that transmission occurs before mutation
## This is a per genome mutation rate
mu_per_genome <- get_mu_from_proba_trans_before_mut(alpha_GT, beta_GT, proba_trans_before_mut) 
sequence_length <- 3000 # Length of simulated genome
mu_per_site <- mu_per_genome / sequence_length

## Demographic parameters
n_demes <- 4
pop_per_deme <- 50000
vec_S_init_per_deme <- rep(pop_per_deme, n_demes)

## Transmission parameters
R0 <- 2.0 # Basic reproduction number
beta_rate_SEIR <- R0 * rate_out_of_I
beta_rate = R0 * rate_out_of_I / (pop_per_deme) # Remaster rate for the force of infection (needs to be scaled by population size)

## Define mixing matrix
curr_mat_proba_migration <- (1. - assortativity) / 12 * 
  matrix(c(0., 1, 4, 7, 1, 0, 7, 4, 4, 7, 0, 1, 7, 4, 1, 0), nrow = 4, ncol = 4, byrow = T)
diag(curr_mat_proba_migration) <- rep(assortativity, n_demes)

## Compute corresponding migration rate matrix
get_migration_rate <- function(alpha_GT, beta_GT, assortativity, curr_mat_proba_migration){
  
  # Get matrix with proportion of infectors going in each group (conditional on transmission occurring outside group)
  mat_rates <- curr_mat_proba_migration
  diag(mat_rates) <- rep(0, nrow(mat_rates))
  mat_rates <- mat_rates/ (1. - assortativity)
  
  # Compute corresponding total migration rate from assortativity
  total_mig_rate <- beta_GT * (exp(- 1/ alpha_GT * log(assortativity)) - 1)
  
  mat_rates <- mat_rates * total_mig_rate
  return(mat_rates)
}

migration_rate_mat <- get_migration_rate(alpha_GT, beta_GT, assortativity, curr_mat_proba_migration)
write.csv(migration_rate_mat, output_migration_rate_mat)

## Write XML
write_xml(output_file = xml_file_path, dir_results = dir_results,
          mat_proba_migration = curr_mat_proba_migration,
          sequence_length = sequence_length, mutation_rate = mu_per_site, 
          beta_rate = beta_rate, 
          rate_out_of_E = rate_out_of_E, 
          rate_out_of_I = rate_out_of_I, 
          p_sample = p_seq,
          n_demes = n_demes, min_sample_per_deme = -1 , 
          vec_S_init_per_deme = vec_S_init_per_deme)