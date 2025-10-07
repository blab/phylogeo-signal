library(tidyverse)
library(RColorBrewer)
library(viridis)

## Load aggregated RR
df_for_comp <- Reduce('bind_rows', lapply(1:50, FUN = function(curr_seed){
  read_csv(paste0('../remaster-sample-size/results/aggregated_RR/df_comp_RR_seed_', curr_seed, '.csv'))
}))

## Define labels used for facet plots
delta_labels <- paste0('Delta == ', 0:15)
names(delta_labels) <- as.character(0:15)
p_seq_labels <- paste0('p[seq] == ', c(0.001, 0.005, 0.01, 0.05))
names(p_seq_labels) <- as.character(c(0.001, 0.005, 0.01, 0.05))

## Compute median correlation across replicates for each scenario
df_median_cor <- df_for_comp %>% 
  mutate(adj_RR = (n_pairs + 1) / (n_pairs_1_x + 1) / (n_pairs_x_2 + 1) * (n_pairs_x_x + 1)) %>% 
  filter(subgroup_1 >= subgroup_2) %>% 
  group_by(p_seq, assortativity, proba_trans_before_mut, seed, threshold) %>% 
  summarise(cor = cor(adj_RR, true_migration_rate, method = 'spearman')) %>%
  mutate(cor = replace_na(cor, 0)) %>% 
  group_by(p_seq, assortativity, proba_trans_before_mut, threshold) %>% 
  summarise(median_cor = median(cor)) %>% 
  mutate(median_cor_crop = ifelse(median_cor < 0., 0., median_cor))

## Function to get heatmaps with the median correlation
get_heatmap_median_cor <- function(df_median_cor){
  heatmap_median_cor <- df_median_cor %>% 
    filter(threshold %in% c(0, 5, 10, 15)) %>% 
    ggplot(aes(x = proba_trans_before_mut, y = assortativity, fill = median_cor_crop)) +
    geom_tile() +
    scale_x_continuous(breaks = seq(0.1, 0.9, 0.1), 
                       labels = c(0.1, '', 0.3, '', 0.5, '', 0.7, '', 0.9),
                       expand = expansion(mult = c(0., 0.)),
                       expression('Probability that transmission occurs before mutation '*p)) +
    scale_y_continuous(breaks = seq(0.1, 0.9, 0.1), 
                       labels = c(0.1, '', 0.3, '', 0.5, '', 0.7, '', 0.9),
                       expand = expansion(mult = c(0., 0.)),
                       expression('Within-group transmission probability '*omega)) +
    viridis::scale_fill_viridis(option = 'magma', limits = c(0., 1.),
                                name = 'Median correlation   ',
                                breaks = seq(0., 1., 0.2)) +
    facet_grid(threshold ~ p_seq, 
               labeller = as_labeller(c(delta_labels, p_seq_labels),
                                      default = label_parsed)) +
    coord_fixed() +
    theme(strip.background = element_rect(fill = 'gray82'),
          strip.text = element_text(colour = 'gray22', size = 14),
          axis.title = element_text(size = 14),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          axis.text = element_text(size = 14))
  return(heatmap_median_cor)
}

heatmap_median_cor <- get_heatmap_median_cor(df_median_cor) +
  theme(legend.key.width = unit(0.5, "in"),
        legend.position = 'top')

# pdf('~/Documents/GitHub/phylogeo-signal/figures/figure_6.pdf', height = 10.5, width = 10)
# plot(heatmap_median_cor)
# dev.off()
# png('~/Documents/GitHub/phylogeo-signal/figures/figure_6.png', height = 10.5, width = 10,
#     res = 350, units = 'in')
# plot(heatmap_median_cor)
# dev.off()
