source('scripts/utils_remaster.R')

args <- commandArgs(trailingOnly = T)
mutation_scenario <- as.character(args[1])
migration_scenario <- as.character(args[2])
xml_file_path <- as.character(args[3])
dir_results <- as.character(args[4])

# Define parameters to run scenarios
## Model parameters
rate_out_of_E <- 0.33 # Rate out which individuals exit the E compartment
rate_out_of_I <- 0.33 # Rate out which individuals exit the I compartment

n_demes <- 2

sequencing_proba <- 0.1 # Sequencing probability
pop_per_deme <- 1000 # Population size per deme
seed <- 1996 # Seed for the random number generator
sequence_length <- 3000 # Length of simulated genome
vec_S_init_per_deme <- rep(pop_per_deme, n_demes)

R0 <- 1.5 # Basic reproduction number
beta_rate_SEIR <- R0 * rate_out_of_I
beta_rate = R0 * rate_out_of_I / (pop_per_deme) # Remaster rate for the force of infection (needs to be scaled by population size)

##

mutation_rate_low <- 0.00002 # Pathogen substitution rate (low scenario)
mutation_rate_high <- mutation_rate_low * 4 # Pathogen substitution rate (high scenario)


if(mutation_scenario == "Low"){
  mutation_rate <- mutation_rate_low
} else{
  mutation_rate <- mutation_rate_high
}

assortativity_high <- 0.98
assortativity_low <- 0.5

if(migration_scenario == "Low"){
  assortativity <- assortativity_high
} else{
  assortativity <- assortativity_low
}

## Define migration matrix
curr_mat_proba_migration <- matrix(rep((1. - assortativity)/(n_demes - 1), n_demes * n_demes),
                                   nrow = n_demes, ncol = n_demes)
diag(curr_mat_proba_migration) <- rep(assortativity, n_demes)

## Write the XML
write_xml(output_file = xml_file_path, 
          dir_results = dir_results,
          mat_proba_migration = curr_mat_proba_migration,
          sequence_length = sequence_length, mutation_rate = mutation_rate, 
          beta_rate = beta_rate, 
          rate_out_of_E = rate_out_of_E, 
          rate_out_of_I = rate_out_of_I, 
          p_sample = sequencing_proba,
          n_demes = n_demes, min_sample_per_deme = -1 , 
          vec_S_init_per_deme = vec_S_init_per_deme)
