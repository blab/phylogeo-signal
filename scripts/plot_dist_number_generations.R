library(tidyverse)
library(phylosamp)
source('utils_sens_spec.R')

empirical_dist_g_gen_R_1.3 <- get_distribution_generations_func_R(1.3)$vec_gen_dist
g_gen_max_R_1.3 <- get_distribution_generations_func_R(1.3)$n_gens
empirical_dist_g_gen_R_1.5 <- get_distribution_generations_func_R(1.5)$vec_gen_dist
g_gen_max_R_1.5 <- get_distribution_generations_func_R(1.5)$n_gens
empirical_dist_g_gen_R_1.7 <- get_distribution_generations_func_R(1.7)$vec_gen_dist
g_gen_max_R_1.7 <- get_distribution_generations_func_R(1.7)$n_gens

plt_dist_gen <- tibble(id_gen = 1:g_gen_max_R_1.3, proba = empirical_dist_g_gen_R_1.3) %>% 
  mutate(R = 1.3) %>% 
  bind_rows(tibble(id_gen = 1:g_gen_max_R_1.5, proba = empirical_dist_g_gen_R_1.5) %>% 
              mutate(R = 1.5)) %>% 
  bind_rows(tibble(id_gen = 1:g_gen_max_R_1.7, proba = empirical_dist_g_gen_R_1.7) %>% 
              mutate(R = 1.7)) %>% 
  mutate(R = paste0('R = ', R)) %>% 
  ggplot(aes(x = id_gen, y = proba)) +
  geom_bar(stat = 'identity', fill = 'firebrick') +
  facet_wrap(. ~ R) +
  theme_bw() +
  scale_x_continuous(name = 'Number of generations') +
  scale_y_continuous(name = 'Probability') +
  theme(strip.background = element_rect(fill = 'gray22'),
        strip.text = element_text(colour = 'white'))

png('../figures/supplementary-figures/distribution_number_generations.png',
    height = 2.5, width = 6, res = 350, units = 'in')
plot(plt_dist_gen)
dev.off()

