library(tidyverse)
library(ggpubr)
library(ggrepel)

## Load confusion matrix parameters across pathogens
df_sens_spec_FDR_across_pathogens_R_1.3 <- read_csv('../results/df_sens_spec_FDR_F1_across_pathogens_R_1.3.csv') %>% mutate(R = 1.3)
df_sens_spec_FDR_across_pathogens_R_1.5 <- read_csv('../results/df_sens_spec_FDR_F1_across_pathogens_R_1.5.csv') %>% mutate(R = 1.5)
df_sens_spec_FDR_across_pathogens_R_1.8 <- read_csv('../results/df_sens_spec_FDR_F1_across_pathogens_R_1.8.csv') %>% mutate(R = 1.8)

df_sens_spec_FDR_across_pathogens <- df_sens_spec_FDR_across_pathogens_R_1.3 %>% 
  bind_rows(df_sens_spec_FDR_across_pathogens_R_1.5) %>% 
  bind_rows(df_sens_spec_FDR_across_pathogens_R_1.8)

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
name_pathogens <- unique(df_sens_spec_FDR_across_pathogens$pathogen)
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


plt_sensitivity_func_delta <- df_sens_spec_FDR_across_pathogens %>% 
  filter(omega == 0.7, delta <= 10, pathogen %in% vec_pathogens_to_keep) %>% 
  mutate(pathogen = factor(pathogen, levels = vec_pathogens_to_keep)) %>% 
  ggplot(aes(x = as.factor(delta), colour = as.factor(pathogen))) +
  geom_point(aes(y = sensitivity)) +
  geom_line(aes(y = sensitivity, linetype = as.factor(R), group = interaction(as.factor(R), pathogen))) +
  scale_x_discrete(breaks = 0:10, 
                   name = expression(atop('Genetic distance threshold', 'for linkage '*Delta))) +
  scale_y_continuous(limits = c(0., NA), breaks = seq(0., 1., 0.2),
                     name = 'Sensitivity',
                     expand = expansion(mult = c(0., 0.05))) +
  scale_colour_manual(name = '', 
                      breaks = vec_pathogens_to_keep, 
                      labels = vec_names_pathogens_to_keep, 
                      values = vec_colors_pathogens_to_keep) +
  scale_linetype_manual(name = 'R', values = 1:3) +
  my_theme_classic() +
  theme(legend.position = 'bottom')


plt_specificity_func_delta <- df_sens_spec_FDR_across_pathogens %>% 
  filter(omega == 0.7, delta <= 10, pathogen %in% vec_pathogens_to_keep) %>% 
  mutate(pathogen = factor(pathogen, levels = vec_pathogens_to_keep)) %>% 
  ggplot(aes(x = as.factor(delta), colour = as.factor(pathogen))) +
  geom_point(aes(y = specificity)) +
  geom_line(aes(y = specificity, linetype = as.factor(R), group = interaction(as.factor(R), pathogen))) +
  scale_x_discrete(breaks = 0:10, 
                   name = expression(atop('Genetic distance threshold', 'for linkage '*Delta))) +
  scale_y_continuous(limits = c(0., NA), breaks = seq(0., 1., 0.2),
                     name = 'Specificity',
                     expand = expansion(mult = c(0.03, 0.05))) +
  scale_colour_manual(name = '', 
                      breaks = vec_pathogens_to_keep, 
                      labels = vec_names_pathogens_to_keep, 
                      values = vec_colors_pathogens_to_keep) +
  scale_linetype_manual(name = 'R', values = 1:3) +
  my_theme_classic() +
  theme(legend.position = 'bottom')

plt_ppv_func_delta <- df_sens_spec_FDR_across_pathogens %>% 
  filter(omega == 0.7, delta <= 10, pathogen %in% vec_pathogens_to_keep) %>% 
  mutate(pathogen = factor(pathogen, levels = vec_pathogens_to_keep)) %>% 
  ggplot(aes(x = as.factor(delta), colour = as.factor(pathogen))) +
  geom_point(aes(y = ppv)) +
  geom_line(aes(y = ppv, linetype = as.factor(R), group = interaction(as.factor(R), pathogen))) +
  scale_x_discrete(breaks = 0:10, 
                   name = expression(atop('Genetic distance threshold', 'for linkage '*Delta))) +
  scale_y_continuous(limits = c(0., NA), breaks = seq(0., 1., 0.2),
                     name = 'PPV',
                     expand = expansion(mult = c(0.03, 0.05))) +
  scale_colour_manual(name = '', 
                      breaks = vec_pathogens_to_keep, 
                      labels = vec_names_pathogens_to_keep, 
                      values = vec_colors_pathogens_to_keep) +
  scale_linetype_manual(name = 'R', values = 1:3) +
  my_theme_classic() +
  theme(legend.position = 'bottom')

plt_sensitivity_func_omega <- df_sens_spec_FDR_across_pathogens %>% 
  filter(delta == 0, pathogen %in% vec_pathogens_to_keep) %>% 
  mutate(pathogen = factor(pathogen, levels = vec_pathogens_to_keep)) %>% 
  ggplot(aes(x = omega, colour = pathogen, group = interaction(pathogen, R))) +
  geom_line(aes(y = sensitivity, linetype = as.factor(R))) +
  scale_x_continuous(name = expression(atop('Probability that transmission', 'occurs within the same group '*omega))) +
  scale_y_continuous(limits = c(0., 1.), name = 'Sensitivity',
                     expand = expansion(mult = c(0., 0.01))) +
  scale_colour_manual(name = '', 
                      breaks = vec_pathogens_to_keep, 
                      labels = vec_names_pathogens_to_keep, 
                      values = vec_colors_pathogens_to_keep) +
  scale_linetype_manual(name = 'R', values = 1:3) +
  my_theme_classic() +
  theme(legend.position = 'bottom')

plt_specificity_func_omega <- df_sens_spec_FDR_across_pathogens %>% 
  filter(delta == 0,
         pathogen %in% vec_pathogens_to_keep) %>% 
  ggplot(aes(x = omega, colour = pathogen, 
             group = interaction(pathogen, R))) +
  geom_line(aes(y = specificity, linetype = as.factor(R))) +
  scale_x_continuous(name = expression(atop('Probability that transmission', 'occurs within the same group '*omega))) +
  scale_y_continuous(limits = c(0., 1.), name = 'Specificity',
                     expand = expansion(mult = c(0.0, 0.05))) +
  scale_colour_manual(name = '', 
                      breaks = vec_pathogens_to_keep, 
                      labels = vec_names_pathogens_to_keep,
                      values = vec_colors_pathogens_to_keep) +
  scale_linetype_manual(name = 'R', values = 1:3) +
  my_theme_classic() +
  theme(legend.position = 'bottom')

plt_ppv_func_omega <- df_sens_spec_FDR_across_pathogens %>% 
  filter(delta == 0, pathogen %in% vec_pathogens_to_keep) %>% 
  ggplot(aes(x = omega, colour = pathogen, 
             group = interaction(pathogen, as.factor(R)))) +
  geom_line(aes(y = ppv, linetype = as.factor(R))) +
  scale_x_continuous(name = expression(atop('Probability that transmission', 'occurs within the same group '*omega))) +
  scale_y_continuous(limits = c(0., 1.), name = 'PPV',
                     expand = expansion(mult = c(0., 0.01))) +
  scale_colour_manual(name = '', 
                      breaks = vec_pathogens_to_keep, 
                      labels = vec_names_pathogens_to_keep, 
                      values = vec_colors_pathogens_to_keep) +
  scale_linetype_manual(name = 'R', values = 1:3) +
  my_theme_classic() +
  theme(legend.position = 'bottom')


panel_final <- 
  ggarrange(plt_sensitivity_func_delta, plt_specificity_func_delta, plt_ppv_func_delta,
            plt_sensitivity_func_omega, plt_specificity_func_omega, plt_ppv_func_omega, 
            labels = 'AUTO', common.legend = T, legend = 'bottom', ncol = 3, nrow = 2)

plot(panel_final)


pdf('../figures/figure_sensitivity_R.pdf', height = 7, width = 9)
plot(panel_final)
dev.off()
png('../figures/figure_sensitivity_R.png', height = 7, width = 9, res = 350, units = 'in')
plot(panel_final)
dev.off()
