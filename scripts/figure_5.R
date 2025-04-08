library(tidyverse)
library(RColorBrewer)
source('utils_sens_spec.R')

## Load F1 score across parameter space
df_f1_score_R_1.3 <- read_csv('../results/df_sens_spec_PPV_F1_across_parameter_space_R_1.3.csv')

## Load age assortativity estimates
df_assortativity_age <- read_csv('../results/df_proba_stays_within_age_group_WA.csv') %>% 
  mutate(name_for_plot = paste0(binning_width, 'y bin'))

## Load spatial assortativity estimates
df_assortativity_space <- read_csv('../results/df_proba_stay_within_geography.csv')

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

col_age_groups <- 'dodgerblue3'
light_col_age_groups <- colorspace::lighten(col_age_groups, 0.8)
col_regions <- 'darkorange3'
light_col_regions <- colorspace::lighten(col_regions, 0.8)

col_pathogens <- 'gray10'
light_col_pathogens <- colorspace::lighten(col_pathogens, 0.8)

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

### Make heatmap
plt_with_annotations <- df_f1_score_R_1.3 %>% 
  group_by(omega, proba_trans_before_mut) %>% 
  filter(f1_score == max(f1_score)) %>% 
  mutate(delta_crop = ifelse(delta >= crop_delta_axis, crop_delta_axis, delta)) %>% 
  ggplot() +
  geom_tile(aes(x = proba_trans_before_mut, y = omega, fill = as.factor(delta_crop))) +
  annotate(geom = 'text', x = 0.5, y = 1.5, label = 'PATHOGENS', color = col_pathogens, 
           fontface = 'bold') +
  annotate(geom = 'rect', xmin = 0., xmax = 1.0, ymin = 1.03, ymax = 1.48, 
           fill = light_col_pathogens) +
  geom_point(data = df_pathogen_char_to_plot, 
             aes(x = proba_trans_before_mut, y = 1.01), col = col_pathogens) +
  geom_text(data = df_pathogen_char_to_plot,
            aes(x = proba_trans_before_mut, y = 1.04, label = pathogen),
            angle = 90, hjust = 0., col = col_pathogens) +
  geom_segment(data = df_pathogen_char_to_plot, 
               aes(x = proba_trans_before_mut, xend = proba_trans_before_mut, 
                   y = 1.01, yend = 1.03), col = col_pathogens) +
  
  annotate(geom = 'text', x = 1.17, y = 1.02, label = 'AGES', color = col_age_groups, 
           fontface = 'bold') +
  annotate(geom = 'rect', xmin = 1.1, xmax = 1.25, ymin = 0, ymax = 1., fill = light_col_age_groups) +
  geom_point(data = df_assortativity_age, aes(x = 1.01, y = proba_stays_within_group),
             color = col_age_groups) +
  geom_text(data = df_assortativity_age,
            aes(x = 1.17, y = proba_stays_within_group, label = name_for_plot), 
            col = col_age_groups) +
  geom_segment(data = df_assortativity_age, 
               aes(x = 1.01, xend = 1.1, y = proba_stays_within_group, yend = proba_stays_within_group),
               color = col_age_groups, linetype = 'dashed') +
  
  annotate(geom = 'text', x = 1.505, y = 1.02, label = 'REGIONS', color = col_regions, 
           fontface = 'bold') +
  annotate(geom = 'rect', xmin = 1.26, xmax = 1.75, ymin = 0, ymax = 1., fill = light_col_regions) +
  geom_point(data = df_assortativity_space, aes(x = 1.01, y = proba_stays_within_group),
             color = col_regions) +
  geom_text(data = df_assortativity_space,
            aes(x = 1.505, y = proba_stays_within_group, label = name_for_plot), 
            col = col_regions) +
  geom_segment(data = df_assortativity_space, 
               aes(x = 1.01, xend = 1.28, y = proba_stays_within_group, yend = proba_stays_within_group),
               color = col_regions, linetype = 'dashed') +
  scale_x_continuous(name = expression(atop('Probability that transmission', 'occurs before mutation '*p)),
                     limits = c(0., NA),
                     breaks = seq(0., 1., 0.1),
                     expand = expansion(mult = c(0., 0.05))) +
  scale_y_continuous(name = expression(atop('Probability that transmission', 'occurs before migration '*omega)),
                     limits = c(0., NA),
                     breaks = seq(0., 1., 0.1),
                     expand = expansion(mult = c(0., 0.05))) +
  scale_fill_manual(breaks = seq(0, crop_delta_axis, 1), 
                    values = my_pal,
                    labels = c(0:(crop_delta_axis -1), paste0(crop_delta_axis, '+')),
                    name = expression(paste('Optimal ', Delta, ' threshold'))) +
  coord_fixed() +
  my_theme_classic() +
  theme(legend.key.height = unit(0.8, 'cm')) 

#pdf('../figures/figure_5.pdf', height = 9, width = 11)
plot(plt_with_annotations)
#dev.off()
