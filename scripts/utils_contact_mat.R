# Function to compute the total number of contacts occurring between
# two age groups of the population (across all individuals)
get_df_contacts_from_mat <- function(contact_mat){
  # contact_mat denotes the contact matrix
  df_contact <- contact_mat %>% 
    mutate(age_i = 0:(n() - 1)) %>% 
    pivot_longer(cols = -'age_i', names_prefix = 'V', names_to = 'age_j', values_to = 'contact_rate') %>% 
    mutate(age_i = as.character(age_i), age_j = as.character(as.numeric(age_j) - 1))
  
  return(df_contact)
}

get_contact_mat_decade_from_df_contacts <- function(df_contacts, age_dist, name_age_group = 'age_decade'){
  
  tmp_age_dist <- age_dist %>% 
    rename(age_aggregated = name_age_group) %>% 
    select(age, n_indiv, age_aggregated)
  
  tmp_age_dist_aggregated <- tmp_age_dist %>% 
    group_by(age_aggregated) %>% summarise(n_indiv = sum(n_indiv))
  
  contact_mat_aggregated <- df_contacts %>% 
    left_join(tmp_age_dist, by = c('age_i' = 'age')) %>% 
    rename(age_aggregated_i = age_aggregated, n_indiv_i = n_indiv) %>% 
    left_join(tmp_age_dist, by = c('age_j' = 'age')) %>% 
    rename(age_aggregated_j = age_aggregated, n_indiv_j = n_indiv) %>% 
    mutate(n_tot_contacts_ij = contact_rate * n_indiv_i) %>% 
    group_by(age_aggregated_i, age_aggregated_j) %>% 
    summarise(n_tot_contacts_ij = sum(n_tot_contacts_ij)) %>% 
    left_join(tmp_age_dist_aggregated, by = c('age_aggregated_i' = 'age_aggregated')) %>% 
    rename(n_indiv_i = n_indiv) %>% 
    left_join(tmp_age_dist_aggregated, by = c('age_aggregated_j' = 'age_aggregated')) %>% 
    rename(n_indiv_j = n_indiv) %>% 
    mutate(n_contact_per_day_ij = n_tot_contacts_ij/n_indiv_i) %>% 
    select(age_decade_i, age_decade_j, n_contact_per_day_ij) %>% 
    ungroup()
  
  return(contact_mat_aggregated)
}

# Function to compute the contact matrix in decades from the
# number of contacts occurring between singe-year age groups
get_prop_contacts_within_group <- function(df_contacts, age_dist, name_age_group = 'age_decade'){

  tmp_age_dist <- age_dist %>% 
    rename(age_aggregated = name_age_group) %>% 
    select(age, n_indiv, age_aggregated)
  
  tmp_age_dist_aggregated <- tmp_age_dist %>% 
    group_by(age_aggregated) %>% summarise(n_indiv = sum(n_indiv))
  
  prop_contacts_within_group <- df_contacts %>% 
    left_join(tmp_age_dist, by = c('age_i' = 'age')) %>% 
    rename(age_aggregated_i = age_aggregated, n_indiv_i = n_indiv) %>% 
    left_join(tmp_age_dist, by = c('age_j' = 'age')) %>% 
    rename(age_aggregated_j = age_aggregated, n_indiv_j = n_indiv) %>% 
    mutate(n_tot_contacts_ij = contact_rate * n_indiv_i) %>% 
    group_by(age_aggregated_i, age_aggregated_j) %>% 
    summarise(n_tot_contacts_ij = sum(n_tot_contacts_ij)) %>% 
    left_join(tmp_age_dist_aggregated, by = c('age_aggregated_i' = 'age_aggregated')) %>% 
    rename(n_indiv_i = n_indiv) %>% 
    left_join(tmp_age_dist_aggregated, by = c('age_aggregated_j' = 'age_aggregated')) %>% 
    rename(n_indiv_j = n_indiv) %>% 
    filter(age_aggregated_i >= age_aggregated_j) %>% ungroup() %>% 
    summarise(n_tot_contacts = sum(n_tot_contacts_ij), n_tot_contacts_within_group = sum(n_tot_contacts_ij[age_aggregated_i == age_aggregated_j])) %>% 
    mutate(prop_within_group = n_tot_contacts_within_group/n_tot_contacts)
  
  
  return(prop_contacts_within_group$prop_within_group)
}

# Function to compute the probability for a contact taken from a random
# individual in the population of being within the same age group 
get_prop_contacts_within_group_across_individuals <- function(df_contacts, age_dist, name_age_group = 'age_decade'){
  
  tmp_age_dist <- age_dist %>% 
    rename(age_aggregated = name_age_group) %>% 
    select(age, n_indiv, age_aggregated)
  
  tmp_age_dist_aggregated <- tmp_age_dist %>% 
    group_by(age_aggregated) %>% summarise(n_indiv = sum(n_indiv))
  
  prop_contacts_within_group <- df_contacts %>% 
    left_join(tmp_age_dist, by = c('age_i' = 'age')) %>% 
    rename(age_aggregated_i = age_aggregated, n_indiv_i = n_indiv) %>% 
    left_join(tmp_age_dist, by = c('age_j' = 'age')) %>% 
    rename(age_aggregated_j = age_aggregated, n_indiv_j = n_indiv) %>% 
    mutate(n_tot_contacts_ij = contact_rate * n_indiv_i) %>% 
    group_by(age_aggregated_i, age_aggregated_j) %>% 
    summarise(n_tot_contacts_ij = sum(n_tot_contacts_ij)) %>% 
    left_join(tmp_age_dist_aggregated, by = c('age_aggregated_i' = 'age_aggregated')) %>% 
    rename(n_indiv_i = n_indiv) %>% 
    left_join(tmp_age_dist_aggregated, by = c('age_aggregated_j' = 'age_aggregated')) %>% 
    rename(n_indiv_j = n_indiv) %>% 
    mutate(n_daily_contacts_ij = n_tot_contacts_ij/n_indiv_i) %>% 
    group_by(age_aggregated_i, n_indiv_i) %>% 
    summarise(prop_contacts_within_i = n_daily_contacts_ij[age_aggregated_i == age_aggregated_j]/sum(n_daily_contacts_ij)) %>% 
    ungroup() %>% 
    summarise(prop_contact_within_group_across_indiv = sum(prop_contacts_within_i*n_indiv_i)/sum(n_indiv_i))
  
  return(prop_contacts_within_group$prop_contact_within_group_across_indiv)
}
get_prop_contacts_within_group_across_individuals_by_age_group <- function(df_contacts, age_dist, name_age_group = 'age_decade'){
  
  tmp_age_dist <- age_dist %>% 
    rename(age_aggregated = name_age_group) %>% 
    select(age, n_indiv, age_aggregated)
  
  tmp_age_dist_aggregated <- tmp_age_dist %>% 
    group_by(age_aggregated) %>% summarise(n_indiv = sum(n_indiv))
  
  prop_contacts_within_group <- df_contacts %>% 
    left_join(tmp_age_dist, by = c('age_i' = 'age')) %>% 
    rename(age_aggregated_i = age_aggregated, n_indiv_i = n_indiv) %>% 
    left_join(tmp_age_dist, by = c('age_j' = 'age')) %>% 
    rename(age_aggregated_j = age_aggregated, n_indiv_j = n_indiv) %>% 
    mutate(n_tot_contacts_ij = contact_rate * n_indiv_i) %>% 
    group_by(age_aggregated_i, age_aggregated_j) %>% 
    summarise(n_tot_contacts_ij = sum(n_tot_contacts_ij)) %>% 
    left_join(tmp_age_dist_aggregated, by = c('age_aggregated_i' = 'age_aggregated')) %>% 
    rename(n_indiv_i = n_indiv) %>% 
    left_join(tmp_age_dist_aggregated, by = c('age_aggregated_j' = 'age_aggregated')) %>% 
    rename(n_indiv_j = n_indiv) %>% 
    mutate(n_daily_contacts_ij = n_tot_contacts_ij/n_indiv_i) %>% 
    group_by(age_aggregated_i, n_indiv_i) %>% 
    summarise(prop_contacts_within_i = n_daily_contacts_ij[age_aggregated_i == age_aggregated_j]/sum(n_daily_contacts_ij)) 
  
  return(prop_contacts_within_group)
}


## Function to allocate age groups based on an age group binwidth
allocate_age_group <- function(age_group_width = 10){
  breaks_age_group <- c(seq(0., 80, age_group_width), 130)
  vec_age_group_to_allocate <- 0:85
  vec_upper_bound_age_group <- sapply(vec_age_group_to_allocate, FUN = function(curr_age){
    breaks_age_group[(which(breaks_age_group > curr_age)[1])] - 1
  })
  vec_lower_bound_age_group <- sapply(vec_age_group_to_allocate, FUN = function(curr_age){
    breaks_age_group[(which(breaks_age_group > curr_age)[1]) - 1]
  })
  vec_name_age_groups <- sapply(1:length(vec_age_group_to_allocate), FUN = function(i_age){
    paste0(vec_lower_bound_age_group[i_age], '-', vec_upper_bound_age_group[i_age], 'y')
  })
  return(tibble(
    age = as.character(vec_age_group_to_allocate),
    age_aggregated = vec_name_age_groups
  ))
}