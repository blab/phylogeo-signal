library(tidyverse)

args <- commandArgs(trailingOnly = T)
file_path_df_distance <- as.character(args[1])
file_path_metadata <- as.character(args[2])
threshold <- as.numeric(args[3])
output_path_df_RR <- as.character(args[4])

# Load distance matrix
df_distance <- read_tsv(file_path_df_distance) %>% filter(n_mutations <= threshold)

# Load metadata
metadata <- read_csv(file_path_metadata) %>% mutate(subgroup = as.factor(subgroup))

## Compute dataframe with number of pairs between groups
compute_n_pairs <- function(df_distance, metadata){
  
  df_distance_sym <- df_distance %>% rename(strain_for_1 = strain_2) %>% 
    rename(strain_2 = strain_1) %>% rename(strain_1 = strain_for_1)
  
  tmp_df <- df_distance %>% bind_rows(df_distance_sym) %>% 
    left_join(metadata %>% select(- time), by = c('strain_1' = 'strain')) %>% 
    rename(subgroup_1 = subgroup) %>% 
    left_join(metadata %>% select(- time), by = c('strain_2' = 'strain')) %>% 
    rename(subgroup_2 = subgroup)
  
  ## Compute number of pairs
  df_n_pairs <- tmp_df %>% 
    group_by(subgroup_1, subgroup_2) %>% 
    summarise(n_pairs = n())
  
  ## Add 0 when necessary
  group_id <- unique(metadata$subgroup) %>% as.character() %>% sort()
  
  df_n_pairs_with_0 <- expand.grid(subgroup_1 = group_id, subgroup_2 = group_id) %>% 
    as_tibble() %>% 
    left_join(df_n_pairs) %>% 
    mutate(n_pairs = replace_na(n_pairs, 0))
  
  return(df_n_pairs_with_0)
}

df_n_pairs <- compute_n_pairs(df_distance, metadata) %>% mutate(n_pairs = as.numeric(n_pairs))

## Compute RR
df_RR <- df_n_pairs %>% 
  group_by(subgroup_1) %>% mutate(n_pairs_1_x = sum(n_pairs)) %>% 
  group_by(subgroup_2) %>% mutate(n_pairs_x_2 = sum(n_pairs)) %>% 
  ungroup() %>% 
  mutate(n_pairs_x_x = sum(n_pairs), 
         RR = n_pairs / n_pairs_1_x / n_pairs_x_2 * n_pairs_x_x)
  
## Save results
write_tsv(df_RR, output_path_df_RR)