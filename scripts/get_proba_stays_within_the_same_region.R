library(tidyverse)
library(tidycensus)
library(covidcast)

source('utils_connectivity.R')

#####################
## Assortativity estimates from US county mobility data
## made publicly available by Pullano et al.
#####################
## Get US county population size corresponding to the indices used in Pullano et al.
pop_size_US_counties <- load_US_county_population_estimates_from_pullano()

## Estimate assortativity parameter at the county level in the US
assortativity_county_US <- load_US_connectivity_by_month(1) %>% 
  left_join(pop_size_US_counties %>% select(index, POPESTIMATE2019),
            by = c('origin' = 'index')) %>% 
  rename(pop_origin = POPESTIMATE2019) %>% 
  filter(origin == destination) %>% arrange(as.numeric(origin)) %>% 
  mutate(n_pop_stay_within = proba * pop_origin) %>% 
  summarise(prop_stay_within_county = sum(n_pop_stay_within)/sum(pop_origin)) %>% 
  unlist() %>% as.numeric()

assortativity_state_US <- load_US_connectivity_by_month(1) %>% 
  left_join(pop_size_US_counties %>% select(index, STNAME, POPESTIMATE2019) %>% rename(state = STNAME),
            by = c('origin' = 'index')) %>% 
  rename(pop_origin = POPESTIMATE2019, state_origin = state) %>% 
  mutate(n_movement_origin_destination = proba * pop_origin) %>% 
  left_join(pop_size_US_counties %>% select(index, STNAME) %>% rename(state = STNAME),
            by = c('destination' = 'index')) %>% 
  rename(state_destination = state) %>% 
  group_by(state_origin, state_destination) %>% 
  summarise(n_movement_origin_destination = sum(n_movement_origin_destination)) %>% 
  ungroup() %>% 
  summarise(prop_within_state = sum(n_movement_origin_destination[state_origin == state_destination])/sum(n_movement_origin_destination)) %>% 
  unlist() %>% as.numeric()

#####################
## Assortativity estimates from WA mobile phone
## visits between census block groups
#####################
df_assortativity_WA <- read_csv('../input/mobility_data/safegraph_washington/df_proba_stay_within_WA_geographies.csv')

## Aggregate estimates
df_all_estimates <- df_assortativity_WA %>% 
  rename(name_group = geography,
         proba_stays_within_group = proportion_within_unit) %>% 
  mutate(name_group = case_when(name_group == 'cbg' ~ 'Census block group',
                                name_group == 'puma' ~ 'PUMA',
                                name_group == 'tract' ~ 'Census tract',
                                name_group == 'county' ~ 'County (WA)')) %>% 
  mutate(name_for_plot = case_when(name_group == 'Census block group' ~ 'Census block group (~1,400 inhab.)',
                                   name_group == 'PUMA' ~ 'PUMA (~125,000 inhab.)',
                                   name_group == 'Census tract' ~ 'Census tract (~4,000 inhab.)',
                                   name_group == 'County (WA)' ~ 'County')) %>% 
  bind_rows(
    tibble(
      name_group = c('County (US)', 'State'),
      proba_stays_within_group = c(assortativity_county_US, assortativity_state_US),
      name_for_plot = c('County', 'State')
    )
  )

write_csv(df_all_estimates, '../results/df_proba_stay_within_geography.csv')
