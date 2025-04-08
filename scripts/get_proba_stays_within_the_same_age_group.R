library(tidyverse)
source('utils_contact_mat.R')

##  Load WA age distribution from Mistry et al. 10.1038/s41467-020-20544-y (1-year age bin)
age_dist_WA <- read.csv('../input/contact_data/United_States_subnational_Washington_age_distribution_85.csv', header = F) %>% 
  as_tibble() %>% rename(age = V1, n_indiv = V2) %>% 
  mutate(age = as.character(age))

## Get contact matrices by setting for WA state for our age groups from Mistry et al. 10.1038/s41467-020-20544-y
# The contact rate represents the per capita probability of contact for an individual of age i (row) with individuals of age j (column)
contact_mat_overall_WA <- read.csv('../input/contact_data/United_States_subnational_Washington_M_overall_contact_matrix_85.csv', header = F)

## Define the binning windows used to defined age groups
vec_age_group_width <- c(1, 2, 5, 10, 15, 20)

## Function to compute the probability for contacts of occurring within the same age group
## for a given age group width
get_prop_contacts_given_age_group_width <- function(age_dist_WA, contact_mat_overall_WA, age_group_width = 10){
  
  ## Allocate age group with width of interest
  tmp_age_dist <- age_dist_WA %>% 
    left_join(allocate_age_group(age_group_width))
  
  ## Dataframe with number of contacts occurring between two age groups (1-year age bin)
  df_contact_overall <- get_df_contacts_from_mat(contact_mat_overall_WA)
  
  return(get_prop_contacts_within_group(df_contact_overall, tmp_age_dist, name_age_group = 'age_aggregated'))
}

## Function to compute the probability for contacts of occurring within the same age group
## for a given age group width by age group
get_prop_contacts_by_age_group_given_age_group_width <- function(age_dist_WA, contact_mat_overall_WA, age_group_width = 10){
  
  ## Allocate age group with width of interest
  tmp_age_dist <- age_dist_WA %>% 
    left_join(allocate_age_group(age_group_width))
  
  ## Dataframe with number of contacts occurring between two age groups (1-year age bin)
  df_contact_overall <- get_df_contacts_from_mat(contact_mat_overall_WA)
  
  return(get_prop_contacts_within_group_across_individuals_by_age_group(df_contact_overall, tmp_age_dist, name_age_group = 'age_aggregated'))
}

## Compute probability that a randomly chosen contact is within the same group
vec_proba_stays_within_group <- sapply(vec_age_group_width, FUN = function(curr_width){
  get_prop_contacts_given_age_group_width(age_dist_WA, contact_mat_overall_WA, curr_width)
})

## Save results
df_proba_stays_within_age_group <- 
  tibble(name_group = paste0('Age groups (', vec_age_group_width, '-year width)'),
         binning_width = vec_age_group_width,
         proba_stays_within_group = vec_proba_stays_within_group)

write_csv(df_proba_stays_within_age_group, '../results/df_proba_stays_within_age_group_WA.csv')

## Compute probability that a contact occurs within the same age group by age group
df_proba_stays_within_age_group_by_age_group <- get_prop_contacts_by_age_group_given_age_group_width(age_dist_WA, contact_mat_overall_WA, 10) %>% 
  mutate(age_aggregated_i = ifelse(age_aggregated_i == '80-129y', '80y+', age_aggregated_i)) %>% 
  ungroup()

write_csv(df_proba_stays_within_age_group_by_age_group, '../results/df_proba_stays_within_age_group_by_age_decade_WA.csv')
