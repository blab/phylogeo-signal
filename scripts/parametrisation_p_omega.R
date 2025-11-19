library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(ggrepel)
source('utils_sens_spec.R')

## Load dataframe with pathogens characteristics
df_pathogen_char <- read_csv('../input/characteristics_pathogens.csv') %>% 
  mutate(subs_rate_per_year = subs_rate_in_mut_per_site_per_year * genome_length,
         subs_rate_per_day = subs_rate_per_year / 365.25,
         mean_delay_between_mut = 1./subs_rate_per_day)

## Estimate parameters of the Gamma distribution for the generation time
## based on each pathogen's mean and SD of the generation time
df_pathogen_char <- bind_cols(df_pathogen_char, 
                              Reduce('bind_rows', lapply(1:nrow(df_pathogen_char), FUN = function(i_row){
                                vec_gamma_param <- get_gamma_param_from_mean_sd(df_pathogen_char$mean_GT[i_row], df_pathogen_char$sd_GT[i_row])
                                names(vec_gamma_param) <- c('alpha_GT', 'beta_GT')
                                vec_gamma_param
                              })))

## Compute the estimated probability that a transmission event occurs before a mutation one
df_pathogen_char$proba_trans_before_mut <- sapply(1:nrow(df_pathogen_char), FUN = function(i_path){
  get_proba_transm_before_mut(alpha_gen_time = df_pathogen_char$alpha_GT[i_path],
                              beta_gen_time = df_pathogen_char$beta_GT[i_path],
                              mu = df_pathogen_char$subs_rate_per_day[i_path])
})

# Plot substitution rate as a function of the expected number of substitution during a transm gen
plt_sub_rate <- df_pathogen_char %>% 
  ggplot(aes(y = subs_rate_per_day, 
             x = subs_rate_per_day * mean_GT)) +
  geom_point(aes(fill = mean_GT),
             shape = 21, color = 'black', size = 3) +
  geom_text_repel(aes(label = pathogen), 
                  size = 3, nudge_y = -0.002) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0., 0.05)),
                     name = expression(atop('Mutation rate '*mu, '(in events per genome per day)'))) +
  scale_x_continuous(expand = expansion(mult = c(0., 0.05)),
                     limits = c(0., NA),
                     breaks = seq(0., 2., 0.2),
                     name = 'Expected number of substitutions\nduring a transmission generation') +
  scale_fill_stepsn(
    colours = colorRampPalette(colors = brewer.pal(9, 'BuPu'))(11), 
    breaks = seq(0., 20, 2),
    limits = c(0., 20),
    name = 'Mean generation time (in days)') +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.key.width = unit(1.4, "cm"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

# 
plt_p <- df_pathogen_char %>% 
  ggplot(aes(x = subs_rate_per_day * mean_GT, y = proba_trans_before_mut, fill = mean_GT)) +
  geom_point(shape = 21, color = 'black', size = 3) +
  scale_y_continuous(name = expression(atop('Probability that transmission', 'occurs before mutation '*p)),
                     limits = c(0., 1.), expand = expansion(mult = c(0., 0.05)),
                     breaks = seq(0., 1., 0.2)) +
  scale_x_continuous(expand = expansion(mult = c(0., 0.05)),
                     limits = c(0., NA),
                     breaks = seq(0., 2., 0.2),
                     name = 'Expected number of substitutions\nduring a transmission generation')  +
  scale_fill_stepsn(
    colours = colorRampPalette(colors = brewer.pal(9, 'BuPu'))(11), 
    breaks = seq(0., 20, 2),
    limits = c(0., 20),
    name = 'Mean generation time (in days)') +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.key.width = unit(1.4, "cm"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))


## Get the mean GT and sd GT of Omicron
## This will be used as a baseline parametrization of the GT for plot C and D
mean_GT <- 4.9 # Mean of the generation time distribution
sd_GT <- 4.8 # SD of the generation time distribution

## Corresponding parameters assuming a Gamma distribution
alpha_GT <- as.numeric(get_gamma_param_from_mean_sd(mean_gamma = mean_GT, sd_gamma = sd_GT)[1])
beta_GT <- as.numeric(get_gamma_param_from_mean_sd(mean_gamma = mean_GT, sd_gamma = sd_GT)[2])


label_GT_plot <- round(alpha_GT / beta_GT * c(1, 2, 4), 2)
vec_mean_GT_plot <- alpha_GT / beta_GT * c(1, 2, 4)
palette_GT <- colorRampPalette(colors = brewer.pal(9, 'BuPu'))(11)[sapply(vec_mean_GT_plot, FUN = function(curr_GT){
  min(which(seq(0, 20, 2) > curr_GT))
})]


plt_p_mu <- expand.grid(mu = seq(1e-4, 0.2, 1e-4),
            mean_GT = vec_mean_GT_plot) %>% 
  as_tibble() %>% 
  mutate(p = get_proba_transm_before_mut(mean_GT * beta_GT, beta_GT, mu)) %>% 
  ggplot(aes(x = mu)) +
  geom_line(aes(y = p, colour = as.factor(mean_GT), group = as.factor(mean_GT)), size = 1) +
  scale_x_continuous(name = expression(paste('Mutation rate '*mu*' (in events per genome per day)')),
                     expand = expansion(mult = c(0., 0.)),
                     limits = c(0., NA)) +
  scale_y_continuous(name = expression(atop('Probability that transmission', 'occurs before mutation '*p)),
                     limits = c(0., 1.),
                     breaks = seq(0., 1., 0.2),
                     expand = expansion(mult = c(0., 0.))) +
  scale_colour_manual(values = palette_GT, breaks = vec_mean_GT_plot, labels = label_GT_plot,
                      name = 'Mean generation time (in days)') +
  theme_classic() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))


plt_omega_lambda <- expand.grid(lambda = seq(1e-3, 1., 1e-3),
                                mean_GT = vec_mean_GT_plot) %>% 
  as_tibble() %>% 
  mutate(omega = get_proba_stays_within_group(mean_GT * beta_GT, beta_GT, lambda)) %>% 
  ggplot(aes(x = lambda)) +
  geom_line(aes(y = omega, colour = as.factor(mean_GT), group = as.factor(mean_GT)), size = 1) +
  scale_x_continuous(name = expression(paste('Migration rate '*lambda*' (in events per day)')),
                     expand = expansion(mult = c(0., 0.)),
                     limits = c(0., NA),
                     breaks = seq(0., 1., 0.2)) +
  scale_y_continuous(name = expression(atop('Within-group transmission', 'probability '*omega)),
                     limits = c(0., 1.),
                     breaks = seq(0., 1., 0.2),
                     expand = expansion(mult = c(0., 0.))) +
  scale_colour_manual(values = palette_GT, 
                      breaks = vec_mean_GT_plot, labels = label_GT_plot,
                      name = 'Mean generation time (in days)') +
  theme_classic() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

panel_parametrisation_p_omega <- 
  ggarrange(ggarrange(plt_sub_rate, plt_p, ncol = 2, 
                      common.legend = T, legend = 'bottom', 
                      labels = c('A', 'B')),
            ggarrange(plt_p_mu, plt_omega_lambda, 
                      ncol = 2, common.legend = T, 
                      legend = 'bottom', labels = c('C', 'D')),
            nrow = 2)

pdf('../figures/figure_parametrisation.pdf', height = 8., width = 9.)
plot(panel_parametrisation_p_omega)
dev.off()
