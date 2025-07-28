args <- commandArgs(trailingOnly = T)
assortativity <- as.numeric(args[1])
output_migration_proba <- as.character(args[2])

mat_proba <- function(assortativity){
  n_demes <- 4
  curr_mat_proba_migration <- (1. - assortativity) / 12 * 
    matrix(c(0., 1, 4, 7, 1, 0, 7, 4, 4, 7, 0, 1, 7, 4, 1, 0), nrow = 4, ncol = 4, byrow = T)
  diag(curr_mat_proba_migration) <- rep(assortativity, n_demes)
  return(curr_mat_proba_migration)
}
mat_migration_rate <- function(assortativity){
  mat_proba_migration <- mat_proba(assortativity)
  mat_migration_rate <- - log(1. - mat_proba_migration)
  diag(mat_migration_rate) <- rep(0., nrow(mat_proba_migration))
  return(mat_migration_rate)
}

mat_migration_proba <- mat_proba(assortativity)

write.csv(mat_migration_proba, output_migration_proba)
