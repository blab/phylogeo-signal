library(tidyverse)
library(RColorBrewer)
library(ggpubr)
source('utils_sens_spec.R')

## Load confusion matrix parameters across parameter space
df_ppv_conditional_1_jump <- read_csv('../results/df_sens_spec_PPV_F1_across_parameter_space_R_1.3_conditional_1_jump.csv')
df_ppv <- read_csv('../results/df_sens_spec_PPV_F1_across_parameter_space_R_1.3.csv')

## Load pathogen evolutionary characteristics
df_pathogen_char <- read_csv('../input/characteristics_pathogens.csv') %>% 
  mutate(subs_rate_per_year = subs_rate_in_mut_per_site_per_year * genome_length,
         subs_rate_per_day = subs_rate_per_year / 365.25,
         mean_delay_between_mut = 1./subs_rate_per_day)

## Estimate parameters of the Gamma distribution for the generation time
df_pathogen_char <- bind_cols(df_pathogen_char, 
                              Reduce('bind_rows', lapply(1:nrow(df_pathogen_char), FUN = function(i_row){
                                vec_gamma_param <- get_gamma_param_from_mean_sd(df_pathogen_char$mean_GT[i_row], df_pathogen_char$sd_GT[i_row])
                                names(vec_gamma_param) <- c('alpha_GT', 'beta_GT')
                                vec_gamma_param
                              })))

## Compute probability that transmission occurs before mutation
df_pathogen_char$proba_trans_before_mut <- sapply(1:nrow(df_pathogen_char), FUN = function(i_path){
  get_proba_transm_before_mut(alpha_gen_time = df_pathogen_char$alpha_GT[i_path],
                              beta_gen_time = df_pathogen_char$beta_GT[i_path],
                              mu = df_pathogen_char$subs_rate_per_day[i_path])
})

## Only select a subset of pathogens to be plotted
df_pathogen_char_to_plot <- df_pathogen_char[c(1, 2, 3, 5, 6, 7, 9, 10), ] %>% 
  mutate(pathogen = ifelse(pathogen == 'Influenza A (H3N2) - HA only', 'Influenza A (H3N2) - HA', pathogen))

## Display optimal Delta value maximising F1 score across the parameter space

### Plot characteristics
crop_delta_axis <- 15 # Maximim delta value used for plotting
my_pal <- colorRampPalette(colors = brewer.pal(11, 'PiYG'))(crop_delta_axis + 1) # Color palette for heatmap

## Define theme for figures
my_theme_classic <- function(){
  (theme_classic() +
     theme(axis.text = element_text(size = 12), 
           axis.title = element_text(size = 12),
           legend.text = element_text(size = 12), 
           legend.title = element_text(size = 12),
           strip.background = element_rect(fill = 'gray22'), 
           strip.text = element_text(size = 12, colour = 'white'))) %>% 
    return()
}
my_theme_bw <- function(){
  (theme_bw() + 
     theme(panel.grid = element_blank(),
           axis.text = element_text(size = 12), 
           axis.title = element_text(size = 12),
           legend.text = element_text(size = 12), 
           legend.title = element_text(size = 12),
           strip.background = element_rect(fill = 'gray22'), 
           strip.text = element_text(size = 12, colour = 'white'))) %>% 
    return()
}

panel_updated_ppv <- df_ppv %>% filter(delta == 0) %>% 
  left_join(df_ppv_conditional_1_jump) %>% 
  ggplot(aes(x = proba_trans_before_mut, y = omega)) +
  geom_tile(aes(fill = ppv_conditional_1_jump)) +
  scale_x_continuous(name = expression(paste('Probability that transmission occurs before mutation '*p)),
                     limits = c(0., 1.),
                     breaks = seq(0., 1., 0.1),
                     expand = expansion(mult = c(0., 0.0))) +
  scale_y_continuous(name = expression(paste('Within-group transmission probability ', omega)),
                     limits = c(0., 1.),
                     breaks = seq(0., 1., 0.1),
                     expand = expansion(mult = c(0., 0.0))) +
  scale_fill_stepsn(
    colours = colorRampPalette(colors = brewer.pal(9, 'Reds'))(11),
    breaks = seq(0., 1., 0.1),
    limits = c(0., 1.),
    name = 'Modified PPV') +
  coord_fixed() +
  my_theme_classic() +
  theme(legend.key.width = unit(1.8, 'cm'),
        legend.position = 'bottom',
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12)
  )

panel_decrease_ppv <- df_ppv %>% filter(delta == 0) %>% 
  left_join(df_ppv_conditional_1_jump) %>% 
  ggplot(aes(x = proba_trans_before_mut, y = omega)) +
  geom_tile(aes(fill = ppv - ppv_conditional_1_jump)) +
  scale_x_continuous(name = expression(paste('Probability that transmission occurs before mutation '*p)),
                     limits = c(0., 1.),
                     breaks = seq(0., 1., 0.1),
                     expand = expansion(mult = c(0., 0.0))) +
  scale_y_continuous(name = expression(paste('Within-group transmission probability ', omega)),
                     limits = c(0., 1.),
                     breaks = seq(0., 1., 0.1),
                     expand = expansion(mult = c(0., 0.0))) +
  scale_fill_stepsn(
    colours = colorRampPalette(colors = brewer.pal(9, 'Oranges'))(11),
    breaks = seq(0., 0.27, 0.03),
    limits = c(0., NA),
    name = 'Decrease in PPV') +
  coord_fixed() +
  my_theme_classic() +
  theme(legend.key.width = unit(1.8, 'cm'),
        legend.position = 'bottom',
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12)
  )


panel_conditional_1_jump <- ggarrange(panel_updated_ppv, panel_decrease_ppv, ncol = 2, labels = 'AUTO')

png('../figures/supplementary-figures/ppv_1_jump.png', height = 5.5, width = 11, res = 350, units = 'in')
plot(panel_conditional_1_jump)
dev.off()
