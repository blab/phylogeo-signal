library(tidyverse)
library(ggtree)
library(ggpubr)

plot_tree <- function(mut_scenario, mig_scenario){
  scenario_name <- paste0('mut_', mut_scenario, '_mig_', mig_scenario)
  title_plot <- paste0('Migration ', mig_scenario, ' / Mutation ', mut_scenario)
  
  # Load tree and metadata
  divergence_tree <- read.tree(paste0('../remaster/results/augur/tree_', scenario_name, '.nwk'))
  metadata <- read_csv(paste0('../remaster/results/metadata/', scenario_name, '.csv'))
  
  metadata_for_tree <- as_tibble(divergence_tree) %>% left_join(metadata, by = c('label' = 'strain')) %>% 
    rename(tip.label = label)
  
  # Plot tree
  ggtree(divergence_tree) %<+%metadata_for_tree  +
    geom_tree(color = 'gray22', size = 0.05) +
    theme_tree() +
    geom_tippoint(aes(colour = as.factor(subgroup))) +
    scale_colour_manual(breaks = 0:1, values = c('darkgoldenrod1', 'deepskyblue2'), name = '') +
    theme(legend.position = 'none') + 
    scale_x_continuous(limits = c(0, 3e-2)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)))  %>% 
    return()
}

tree_mut_low_mig_low <- plot_tree('Low', 'Low')
tree_mut_low_mig_high <- plot_tree('Low', 'High')
tree_mut_high_mig_low <- plot_tree('High', 'Low')
tree_mut_high_mig_high <- plot_tree('High', 'High')

panel <- 
  ggarrange(ggarrange(tree_mut_low_mig_low %>% annotate_figure(top = 'Low', left = 'Low'),
                      tree_mut_low_mig_high %>% annotate_figure(top = 'High', left = '')),
            ggarrange(tree_mut_high_mig_low %>% annotate_figure(top = '', left = 'High'), 
                      tree_mut_high_mig_high %>% annotate_figure(top = '', left = '')),
            nrow = 2) %>% 
  annotate_figure(top = 'Migration rate', left = 'Mutation rate')

pdf('../figures/figure_1.pdf', height = 6, width = 10)
plot(panel)
dev.off()
