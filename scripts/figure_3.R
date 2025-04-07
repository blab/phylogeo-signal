library(tidyverse)
library(ggpubr)
library(ggrepel)

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


## Figure 3A (Sensitivity across pathogens)
plt_sensitivity_with_colours <- df_sens_spec_FDR_across_pathogens_R_1.3 %>% 
  filter(delta == 0, omega == 0.7) %>% 
  ggplot(aes(x = proba_trans_before_mut, y = sensitivity, colour = pathogen)) +
  geom_point() +
  geom_text_repel(data = df_sens_spec_FDR_across_pathogens_R_1.3 %>% 
                    filter(delta == 0, omega == 0.7, pathogen %in% vec_pathogens_to_keep), aes(label = pathogen),
                  box.padding = 0.2, max.overlaps = Inf) +
  scale_x_continuous(name = expression(atop('Probability that transmission', 'occurs before mutation '*p)),
                     breaks = seq(0., 1., 0.2),
                     limits = c(0., 1.)) +
  scale_y_continuous(name = 'Sensitivity', limits = c(0., NA),
                     expand = expansion(mult = c(0.,.05))) +
  scale_colour_manual(values = vec_colors_all_pathogens, breaks = name_pathogens) +
  my_theme_classic() +
  theme(legend.position = 'none')

## Figure 3B (Specificity across pathogens)
plt_specificity_with_colours <- df_sens_spec_FDR_across_pathogens_R_1.3 %>% 
  filter(delta == 0, omega == 0.7) %>% 
  ggplot(aes(x = proba_trans_before_mut, y = specificity, colour = pathogen)) +
  geom_point() +
  scale_x_continuous(name = expression(atop('Probability that transmission', 'occurs before mutation '*p)),
                     breaks = seq(0., 1., 0.2),
                     limits = c(0., 1.)) +
  scale_y_continuous(name = 'Specificity', limits = c(0., 1.),
                     expand = expansion(mult = c(0.0,.05))) +
  scale_colour_manual(values = vec_colors_all_pathogens, breaks = name_pathogens) +
  my_theme_classic() +
  theme(legend.position = 'none')

## Figure 3C (PPV across pathogens)
plt_ppv_with_colours <- df_sens_spec_FDR_across_pathogens_R_1.3 %>% 
  filter(delta == 0, omega == 0.7) %>% 
  ggplot(aes(x = proba_trans_before_mut, y = ppv, colour = pathogen)) +
  geom_point() +
  scale_x_continuous(name = expression(atop('Probability that transmission', 'occurs before mutation '*p)),
                     breaks = seq(0., 1., 0.2),
                     limits = c(0., 1.)) +
  scale_y_continuous(name = 'PPV', limits = c(0., 1.),
                     expand = expansion(mult = c(0.,.05))) +
  scale_colour_manual(values = vec_colors_all_pathogens, breaks = name_pathogens) +
  my_theme_classic() +
  theme(legend.position = 'none')

## Figure 3D (Sensitivity as a function of Delta)
plt_sensitivity_func_delta <- df_sens_spec_FDR_across_pathogens_R_1.3 %>% 
  filter(omega == 0.7, delta <= 10,
         pathogen %in% vec_pathogens_to_keep) %>% 
  ggplot(aes(x = as.factor(delta), colour = pathogen, group = pathogen)) +
  geom_point(aes(y = sensitivity)) +
  geom_line(aes(y = sensitivity)) +
  scale_x_discrete(breaks = 0:10, 
                   name = expression(atop('Genetic distance threshold', 'for linkage '*Delta))) +
  scale_y_continuous(limits = c(0., NA), breaks = seq(0., 1., 0.2),
                     name = 'Sensitivity',
                     expand = expansion(mult = c(0., 0.05))) +
  scale_colour_manual(name = '', 
                      breaks = vec_pathogens_to_keep, 
                      labels = vec_names_pathogens_to_keep, 
                      values = vec_colors_pathogens_to_keep) +
  my_theme_classic() +
  theme(legend.position = 'bottom')

## Figure 3E (Specificity as a function of Delta)
plt_specificity_func_delta <- df_sens_spec_FDR_across_pathogens_R_1.3 %>% 
  filter(omega == 0.7, delta <= 10,
         pathogen %in% vec_pathogens_to_keep) %>% 
  ggplot(aes(x = as.factor(delta), colour = pathogen, group = pathogen)) +
  geom_point(aes(y = specificity)) +
  geom_line(aes(y = specificity)) +
  scale_x_discrete(breaks = 0:10, 
                   name = expression(atop('Genetic distance threshold', 'for linkage '*Delta))) +
  scale_y_continuous(limits = c(0., 1.), 
    name = 'Specificity',
    expand = expansion(mult = c(0., 0.05))) +
  scale_colour_manual(name = '', 
                      breaks = vec_pathogens_to_keep, 
                      labels = vec_names_pathogens_to_keep, 
                      values = vec_colors_pathogens_to_keep) +
  my_theme_classic() +
  theme(legend.position = 'bottom')

## Figure 3F (PPV as a function of Delta)
plt_ppv_func_delta <- df_sens_spec_FDR_across_pathogens_R_1.3 %>% 
  filter(omega == 0.7, delta <= 10,
         pathogen %in% vec_pathogens_to_keep) %>% 
  ggplot(aes(x = as.factor(delta), colour = pathogen, group = pathogen)) +
  geom_point(aes(y = ppv)) +
  geom_line(aes(y = ppv)) +
  scale_x_discrete(breaks = 0:10, 
                   name = expression(atop('Genetic distance threshold', 'for linkage '*Delta))) +
  scale_y_continuous(limits = c(0., NA), breaks = seq(0., 1., 0.2),
                     name = 'PPV',
                     expand = expansion(mult = c(0., 0.05))) +
  scale_colour_manual(name = '', 
                      breaks = vec_pathogens_to_keep, 
                      labels = vec_names_pathogens_to_keep, 
                      values = vec_colors_pathogens_to_keep) +
  my_theme_classic() +
  theme(legend.position = 'bottom')

## Figure 3G (Sensitivity as a function of omega)
plt_sensitivity_func_omega <- df_sens_spec_FDR_across_pathogens_R_1.3 %>% 
  filter(delta == 0,
         pathogen %in% vec_pathogens_to_keep) %>% 
  ggplot(aes(x = omega, colour = pathogen, 
             group = interaction(pathogen, delta))) +
  geom_line(aes(y = sensitivity)) +
  scale_x_continuous(name = expression(atop('Probability that transmission', 'occurs within the same group '*omega))) +
  scale_y_continuous(limits = c(0., 1.), name = 'Sensitivity',
                     expand = expansion(mult = c(0., 0.01))) +
  scale_colour_manual(name = '', 
                      breaks = vec_pathogens_to_keep, 
                      labels = vec_names_pathogens_to_keep, 
                      values = vec_colors_pathogens_to_keep) +
  my_theme_classic() +
  theme(legend.position = 'bottom')

## Figure 3H (Specificity as a function of omega)
plt_specificity_func_omega <- df_sens_spec_FDR_across_pathogens_R_1.3 %>% 
  filter(delta == 0,
         pathogen %in% vec_pathogens_to_keep) %>% 
  ggplot(aes(x = omega, colour = pathogen, 
             group = interaction(pathogen, delta))) +
  geom_line(aes(y = specificity)) +
  scale_x_continuous(name = expression(atop('Probability that transmission', 'occurs within the same group '*omega))) +
  scale_y_continuous(limits = c(0., 1.), name = 'Specificity',
                     expand = expansion(mult = c(0.0, 0.05))) +
  scale_colour_manual(name = '', 
                      breaks = vec_pathogens_to_keep, 
                      labels = vec_names_pathogens_to_keep,
                      values = vec_colors_pathogens_to_keep) +
  my_theme_classic() +
  theme(legend.position = 'bottom')

## Figure 3I (PPV as a function of omega)
plt_ppv_func_omega <- df_sens_spec_FDR_across_pathogens_R_1.3 %>% 
  filter(delta == 0,
         pathogen %in% vec_pathogens_to_keep) %>% 
  ggplot(aes(x = omega, colour = pathogen, 
             group = interaction(pathogen, delta))) +
  geom_line(aes(y = ppv)) +
  scale_x_continuous(name = expression(atop('Probability that transmission', 'occurs within the same group '*omega))) +
  scale_y_continuous(limits = c(0., 1.), name = 'PPV',
                     expand = expansion(mult = c(0., 0.01))) +
  scale_colour_manual(name = '', 
                      breaks = vec_pathogens_to_keep, 
                      labels = vec_names_pathogens_to_keep, 
                      values = vec_colors_pathogens_to_keep) +
  my_theme_classic() +
  theme(legend.position = 'bottom')


## Panel for Figure 3
panel_final <- 
  ggarrange(ggarrange(plt_sensitivity_with_colours, plt_specificity_with_colours, plt_ppv_with_colours,
                      labels = c('A', 'B', 'C'),
                      ncol = 3, nrow = 1),
            ggarrange(plt_sensitivity_func_delta, plt_specificity_func_delta, plt_ppv_func_delta,
                      plt_sensitivity_func_omega, plt_specificity_func_omega, plt_ppv_func_omega,
                      common.legend = T, legend = 'bottom', 
                      labels = c('D', 'E', 'F', 'G', 'H', 'I'),
                      ncol = 3, nrow = 2),
            nrow = 2, ncol = 1,
            heights = c(1., 2.2))


#pdf('../figures/figure_3.pdf', height = 9, width = 9)
plot(panel_final)
#dev.off()

######### 
## F1-score
df_for_f1_score <- df_sens_spec_FDR_across_pathogens_R_1.3 %>% 
  filter(pathogen %in% vec_pathogens_to_keep)
df_arrow_f1 <- df_for_f1_score %>% group_by(pathogen, omega) %>% 
  filter(f1_score == max(f1_score))

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

pdf('../figures/figure_4_f1_score.pdf', 
    height = 4., width = 8)
plot(panel_f1_score)
dev.off()

png('../figures/figure_4_f1_score.png', height = 4., width = 8,
    res = 350, units = 'in')
plot(panel_f1_score)
dev.off()

