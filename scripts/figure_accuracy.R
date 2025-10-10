library(tidyverse)
library(ggpubr)

## Load confusion matrix parameters across pathogens
df_sens_spec_FDR_across_pathogens_R_1.3 <- read_csv('../results/df_sens_spec_FDR_F1_across_pathogens_R_1.3.csv')

## Add accuracy calculations
df_sens_spec_FDR_across_pathogens_R_1.3 <- df_sens_spec_FDR_across_pathogens_R_1.3 %>% 
  mutate(accuracy = specificity + proba_0_1_jumps * (sensitivity - specificity),
         contribution_TP_accuracy = specificity * proba_0_1_jumps,
         contribution_TN_accuracy = sensitivity * (1. - proba_0_1_jumps))

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

## Accuracy as a function of p
plt_accuracy_func_p <- df_sens_spec_FDR_across_pathogens_R_1.3 %>% 
  filter(delta == 0, omega == 0.7) %>% 
  ggplot(aes(x = proba_trans_before_mut, y = accuracy, colour = pathogen)) +
  geom_point() +
  scale_x_continuous(name = expression(atop('Probability that transmission', 'occurs before mutation '*p)),
                     breaks = seq(0., 1., 0.2),
                     limits = c(0., 1.)) +
  scale_y_continuous(name = 'Accuracy', limits = c(0., NA),
                     expand = expansion(mult = c(0.,.05))) +
  scale_colour_manual(values = vec_colors_all_pathogens, breaks = name_pathogens,
                      guide = 'none') +
  my_theme_classic() +
  theme(legend.position = 'none')


plt_accuracy_func_p_omega_0.99 <- df_sens_spec_FDR_across_pathogens_R_1.3 %>% 
  filter(delta == 0, omega == 0.99) %>% 
  ggplot(aes(x = proba_trans_before_mut, y = accuracy, colour = pathogen)) +
  geom_point() +
  scale_x_continuous(name = expression(atop('Probability that transmission', 'occurs before mutation '*p)),
                     breaks = seq(0., 1., 0.2),
                     limits = c(0., 1.)) +
  scale_y_continuous(name = 'Accuracy', limits = c(0., NA),
                     expand = expansion(mult = c(0.,.05))) +
  scale_colour_manual(values = vec_colors_all_pathogens, breaks = name_pathogens,
                      guide = 'none') +
  my_theme_classic() +
  theme(legend.position = 'none')

## Accuracy as a function of Delta
plt_accuracy_func_delta <- df_sens_spec_FDR_across_pathogens_R_1.3 %>% 
  filter(omega == 0.7, delta <= 10,
         pathogen %in% vec_pathogens_to_keep) %>% 
  ggplot(aes(x = as.factor(delta), colour = pathogen, group = pathogen)) +
  geom_point(aes(y = accuracy)) +
  geom_line(aes(y = accuracy)) +
  scale_x_discrete(breaks = 0:10, 
                   name = expression(atop('Genetic distance threshold', 'for linkage '*Delta))) +
  scale_y_continuous(limits = c(0., NA), breaks = seq(0., 1., 0.2),
                     name = 'Accuracy',
                     expand = expansion(mult = c(0., 0.05))) +
  scale_colour_manual(name = '', 
                      breaks = vec_pathogens_to_keep, 
                      labels = vec_names_pathogens_to_keep, 
                      values = vec_colors_pathogens_to_keep) +
  my_theme_classic() +
  theme(legend.position = 'bottom')

plt_accuracy_func_delta_omega_0.99 <- df_sens_spec_FDR_across_pathogens_R_1.3 %>% 
  filter(omega == 0.99, delta <= 10,
         pathogen %in% vec_pathogens_to_keep) %>% 
  ggplot(aes(x = as.factor(delta), colour = pathogen, group = pathogen)) +
  geom_point(aes(y = accuracy)) +
  geom_line(aes(y = accuracy)) +
  scale_x_discrete(breaks = 0:10, 
                   name = expression(atop('Genetic distance threshold', 'for linkage '*Delta))) +
  scale_y_continuous(limits = c(0., NA), breaks = seq(0., 1., 0.2),
                     name = 'Accuracy',
                     expand = expansion(mult = c(0., 0.05))) +
  scale_colour_manual(name = '', 
                      breaks = vec_pathogens_to_keep, 
                      labels = vec_names_pathogens_to_keep, 
                      values = vec_colors_pathogens_to_keep) +
  my_theme_classic() +
  theme(legend.position = 'bottom')

## Accuracy as a function of omega
plt_accuracy_func_omega <- df_sens_spec_FDR_across_pathogens_R_1.3 %>% 
  filter(delta == 0,
         pathogen %in% vec_pathogens_to_keep) %>% 
  ggplot(aes(x = omega, colour = pathogen, 
             group = interaction(pathogen, delta))) +
  geom_line(aes(y = accuracy)) +
  scale_x_continuous(name = expression(atop('Within-group transmission', 'probability '*omega))) +
  scale_y_continuous(limits = c(0., 1.), name = 'Accuracy',
                     expand = expansion(mult = c(0., 0.01))) +
  scale_colour_manual(name = '', 
                      breaks = vec_pathogens_to_keep, 
                      labels = vec_names_pathogens_to_keep, 
                      values = vec_colors_pathogens_to_keep) +
  my_theme_classic() +
  theme(legend.position = 'bottom')



accuracy_panel <- ggarrange(plt_accuracy_func_p + 
                              theme(legend.position = 'none') +
                              ggtitle(label = expression(paste(omega, ' = 0.7'))),
                            plt_accuracy_func_p_omega_0.99 + theme(legend.position = 'none') +
                              ggtitle(label = expression(paste(omega, ' = 0.99'))),
                            plt_accuracy_func_delta +
                              ggtitle(label = expression(paste(omega, ' = 0.7'))),
                            plt_accuracy_func_delta_omega_0.99 +
                              ggtitle(label = expression(paste(omega, ' = 0.99'))),
                            plt_accuracy_func_omega,
                            nrow = 3, ncol = 2, 
                            common.legend = T, legend = 'bottom',
                            labels = 'AUTO')

plot(accuracy_panel)
png('../figures/supplementary-figures/accuracy_panel.png', height = 9, width = 6,
    res = 350, units = 'in')
plot(accuracy_panel)
dev.off()
pdf('../figures/supplementary-figures/accuracy_panel.pdf', height = 9, width = 6)
plot(accuracy_panel)
dev.off()

