library(tidyverse)
library(ggpubr)

## Load confusion matrix parameters across pathogens
df_sens_spec_FDR_across_pathogens_R_1.3 <- read_csv('../results/df_sens_spec_FDR_F1_across_pathogens_R_1.3.csv')

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

## Define colors for pathogens
name_pathogens <- unique(df_sens_spec_FDR_across_pathogens_R_1.3$pathogen)
vec_pathogens_to_keep <- c('Influenza A (H3N2) - HA only', 'SARS-CoV-2 (Omicron)', 'SARS-CoV')
vec_names_pathogens_to_keep <- c('Influenza A (H3N2) - HA only\nHigh p', 
                                 'SARS-CoV-2 (Omicron)\nMedium p',
                                 'SARS-CoV\nLow p')
vec_colors_pathogens_to_keep <- c('firebrick3', 'orange2', 'forestgreen')
vec_colors_all_pathogens <- sapply(name_pathogens, FUN = function(pathogen){
  ifelse(pathogen %in% vec_pathogens_to_keep, 
         vec_colors_pathogens_to_keep[which(vec_pathogens_to_keep == pathogen)], 
         'gray10')
})


## Subset the results table to only retain 3 focal pathogens
df_for_f1_score <- df_sens_spec_FDR_across_pathogens_R_1.3 %>% 
  filter(pathogen %in% vec_pathogens_to_keep)
## For each pathogen and mixing process (characterized by omega), 
## only retain the Delta with the highest F1 score
df_arrow_f1 <- df_for_f1_score %>% group_by(pathogen, omega) %>% 
  filter(f1_score == max(f1_score))

## Figure 4A - F1 score as a function of Delta
### Define characteristics of arrows that we plot
length_arrow_f1 <- 0.05
y_position_arrow_f1 <- 0.5

plt_f1_score_func_delta <- df_for_f1_score %>% filter(omega == 0.7) %>% 
  ggplot(aes(x = as.factor(delta), colour = pathogen, group = pathogen)) +
  geom_point(aes(y = f1_score)) +
  geom_line(aes(y = f1_score)) +
  geom_segment(data = df_arrow_f1 %>% filter(omega == 0.7), 
               aes(x = as.factor(delta), xend = as.factor(delta), 
                   y = y_position_arrow_f1 + length_arrow_f1, yend = y_position_arrow_f1),
               arrow = arrow(length = unit(0.03, "npc"))) +
  scale_x_discrete(breaks = 0:15, 
                   name = expression(atop('Genetic distance threshold', 'for linkage '*Delta))) +
  scale_y_continuous(limits = c(0., NA), breaks = seq(0., 1., 0.1),
                     name = expression(paste(F[1], ' score')),
                     expand = expansion(mult = c(0., 0.01))) +
  scale_colour_manual(name = '', 
                      breaks = vec_pathogens_to_keep, 
                      labels = vec_names_pathogens_to_keep, 
                      values = vec_colors_pathogens_to_keep) +
  my_theme_classic() +
  theme(legend.position = 'bottom')

## Figure 4B - Delta value maximizing F1 score across omega values
plt_optimal_threshold_func_omega <- df_sens_spec_FDR_across_pathogens_R_1.3 %>% 
  group_by(omega, pathogen) %>% 
  filter(f1_score == max(f1_score)) %>% 
  filter(pathogen %in% vec_pathogens_to_keep) %>% 
  ggplot(aes(x = omega, y = delta, colour = pathogen)) +
  scale_x_continuous(name = expression(atop('Probability that transmission', 'occurs within the same group '*omega))) +
  scale_y_continuous(breaks = 0:15, 
                     name = expression(atop('Genetic distance threshold '*Delta, ' maximizing '*F[1]*' score')),
                     limits = c(0, 14)) +
  scale_colour_manual(name = '', 
                      breaks = vec_pathogens_to_keep, 
                      labels = vec_names_pathogens_to_keep, 
                      values = vec_colors_pathogens_to_keep) +
  geom_line() +
  my_theme_classic() +
  theme(legend.position = 'bottom')


panel_f1_score <- 
  ggarrange(plt_f1_score_func_delta, plt_optimal_threshold_func_omega, 
          nrow = 1, ncol = 2, widths = c(1.1, 1.0),
          common.legend = T, legend = 'bottom', labels = 'AUTO')

#pdf('../figures/figure_4_f1_score.pdf', height = 4., width = 8)
plot(panel_f1_score)
#dev.off()

# png('../figures/figure_4_f1_score.png', height = 4., width = 8,
#     res = 350, units = 'in')
# plot(panel_f1_score)
# dev.off()

