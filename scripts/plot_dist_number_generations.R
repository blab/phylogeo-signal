library(tidyverse)
library(phylosamp)
source('utils_sens_spec.R')

empirical_dist_g_gen_R_1.3 <- get_distribution_generations_func_R(1.3)$vec_gen_dist
g_gen_max_R_1.3 <- get_distribution_generations_func_R(1.3)$n_gens
empirical_dist_g_gen_R_1.5 <- get_distribution_generations_func_R(1.5)$vec_gen_dist
g_gen_max_R_1.5 <- get_distribution_generations_func_R(1.5)$n_gens
empirical_dist_g_gen_R_1.8 <- get_distribution_generations_func_R(1.7)$vec_gen_dist
g_gen_max_R_1.8 <- get_distribution_generations_func_R(1.7)$n_gens

tibble(id_gen = 1:g_gen_max_R_1.3, proba = empirical_dist_g_gen_R_1.3) %>% 
  mutate(R = 1.3) %>% 
  bind_rows(tibble(id_gen = 1:g_gen_max_R_1.5, proba = empirical_dist_g_gen_R_1.5) %>% 
              mutate(R = 1.5)) %>% 
  bind_rows(tibble(id_gen = 1:g_gen_max_R_1.8, proba = empirical_dist_g_gen_R_1.8) %>% 
              mutate(R = 1.8)) %>% 
  mutate(R = paste0('R = ', R)) %>% 
  ggplot(aes(x = id_gen, y = proba)) +
  geom_bar(stat = 'identity') +
  facet_wrap(. ~ R)
