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