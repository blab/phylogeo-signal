library(tidyverse)
library(RColorBrewer)
source('utils_sens_spec.R')

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

## Figure depicting 
list_plots <- lapply(1:length(vec_names_pathogens_to_keep), FUN = function(i_pathogen){
  tmp_df_pathogen <- df_sens_spec_FDR_across_pathogens_R_1.3 %>% 
    filter(pathogen == vec_pathogens_to_keep[i_pathogen], omega == 0.7, delta <= 10) 
  pathogen_factor <- max(tmp_df_pathogen$proba_less_delta_mutations)/max(tmp_df_pathogen$ppv)
  plt <- tmp_df_pathogen %>% 
    ggplot(aes(x = delta)) + 
    scale_x_continuous(breaks = 0:10, 
                       name = expression(atop('Genetic distance threshold', 'for linkage '*Delta))) +
    scale_y_continuous(name = 'PPV', limits = c(0., NA), 
                       expand = expansion(mult = c(0., 0.05)),
                       sec.axis = sec_axis(~ . * pathogen_factor,
                                           name = expression("Linkage probability "~ P*"["*M <= Delta*"]"))) +
    geom_line(aes(y = ppv, group = pathogen)) +
    geom_line(aes(y = proba_less_delta_mutations/pathogen_factor, group = pathogen), linetype = 'dashed') +
    my_theme_classic() +
    theme(strip.background = element_rect(fill = vec_colors_pathogens_to_keep[i_pathogen], colour = vec_colors_pathogens_to_keep[i_pathogen])) +
    facet_wrap(. ~ pathogen, scales = 'free')
  plt
})

panel_trade_off <- ggarrange(plotlist = list_plots, ncol = 3)

pdf('../figures/trade_off_ppv_linkage_proba.pdf', height = 3, width = 11.5)
plot(panel_trade_off)
dev.off()


png('../figures/trade_off_ppv_linkage_proba.png', height = 3, width = 11.5, res = 350, units = 'in')
plot(panel_trade_off)
dev.off()



