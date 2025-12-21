library(tidyverse)
library(tidycensus)
library(covidcast)

## Function to load and process connectivity data from Pullano et al.
load_US_connectivity_by_month <- function(month_id){
  file_path <- paste0('../input/mobility_data/pullano_connectivity/monthly_matrix_sd_county_month', month_id, '..txt')
  connectivity_data <- read_table(file_path, col_names = c('origin', 'destination', 'month', 'proba')) %>% 
    mutate(origin = as.character(origin), destination = as.character(destination))
  return(connectivity_data)
}

## Function to load population estimate and match them to the GEOID used in Pullano et al.,
load_US_county_population_estimates_from_pullano <- function(){
  index_counties <- read_csv('../input/mobility_data/pullano_connectivity/counties_fips_index.csv', col_types = 'cc')
  index_counties %>% 
    left_join(county_census %>% mutate(FIPS = as.character(as.numeric(FIPS))), by = c('GEO_ID' = 'FIPS')) %>% 
    return()
}

## Get US county population size corresponding to the indices used in Pullano et al.
pop_size_US_counties <- load_US_county_population_estimates_from_pullano()

df_adjusted_movement <- load_US_connectivity_by_month(1) %>% 
  left_join(pop_size_US_counties %>% select(index, STNAME, POPESTIMATE2019) %>% rename(state = STNAME),
            by = c('origin' = 'index')) %>% 
  rename(pop_origin = POPESTIMATE2019, state_origin = state) %>% 
  mutate(n_movement_origin_destination = proba * pop_origin) %>% 
  left_join(pop_size_US_counties %>% select(index, STNAME) %>% rename(state = STNAME),
            by = c('destination' = 'index')) %>% 
  rename(state_destination = state) %>% 
  group_by(state_origin, state_destination) %>% 
  summarise(n_movement_origin_destination = sum(n_movement_origin_destination)) 

write_tsv(df_adjusted_movement, '~/Documents/projects/ncov-usa-mig/safegraph_between_states_adjusted_pop_size_pullano.tsv')
