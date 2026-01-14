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

## Heatmap depicting the maximum correlation obtained for a given
## p_seq value (across Delta thresholds)
heatmap_max_corr_across_delta_by_pseq <- df_median_cor %>% 
  group_by(proba_trans_before_mut, assortativity, p_seq) %>% 
  summarise(max_median_cor = max(median_cor)) %>% 
  ggplot(aes(x = proba_trans_before_mut, y = assortativity, fill = max_median_cor)) +
  geom_tile() +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.1), 
                     labels = c(0.1, '', 0.3, '', 0.5, '', 0.7, '', 0.9), 
                     expand = expansion(mult = c(0., 0.)),
                     expression('Probability that transmission occurs before mutation '*p)) +
  scale_y_continuous(breaks = seq(0.1, 0.9, 0.1), 
                     labels = c(0.1, '', 0.3, '', 0.5, '', 0.7, '', 0.9),
                     expand = expansion(mult = c(0., 0.)),
                     name = expression(atop("Within-group transmission", 
                                            paste("probability ", omega)))) +
  viridis::scale_fill_viridis(option = 'magma', limits = c(0., 1.),
                              name = 'Maximum\nmedian\ncorrelation\nacross\nthresholds\n',
                              breaks = seq(0., 1., 0.2)) +
  facet_grid(. ~ p_seq, 
             labeller = as_labeller(c(delta_labels, p_seq_labels),
                                    default = label_parsed)) +
  coord_fixed() +
  theme(strip.background = element_rect(fill = 'gray82'),
        strip.text = element_text(colour = 'gray22', size = 14),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14))


## Heatmap depicting the minimum sequencing effort required 
## to reach correlation of  50% and 90%
heatmap_min_seq_rate_correlation_threshold <- df_median_cor %>% 
  group_by(proba_trans_before_mut, assortativity, p_seq) %>% 
  summarise(max_median_cor = max(median_cor)) %>% 
  mutate(is_above_cor_50 = max_median_cor > 0.5,
         is_above_cor_90 = max_median_cor > 0.9) %>% 
  group_by(proba_trans_before_mut, assortativity) %>% 
  summarise(min_p_seq_90 = min(p_seq[is_above_cor_90]),
            min_p_seq_50 = min(p_seq[is_above_cor_50]),
  ) %>% 
  group_by(assortativity, proba_trans_before_mut) %>% 
  pivot_longer(cols = 3:4, values_to = 'min_p_seq', names_to = 'cor_threshold', names_prefix = 'min_p_seq_') %>% 
  mutate(cor_threshold = paste0('Correlation threshold: ', cor_threshold, '%')) %>% 
  ggplot(aes(x = proba_trans_before_mut, y = assortativity, fill = as.factor(min_p_seq))) +
  geom_tile() +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.1), 
                     labels = c(0.1, '', 0.3, '', 0.5, '', 0.7, '', 0.9),
                     expand = expansion(mult = c(0., 0.)),
                     expression('Probability that transmission occurs before mutation '*p)) +
  scale_y_continuous(breaks = seq(0.1, 0.9, 0.1), 
                     labels = c(0.1, '', 0.3, '', 0.5, '', 0.7, '', 0.9),
                     expand = expansion(mult = c(0., 0.)),
                     name = expression(atop("Within-group transmission", 
                                            paste("probability ", omega)))) +
  scale_fill_manual(name = 'Minimal sequencing\nrate to reach\ncorrelation\nthreshold\n',
                    values = brewer.pal(5, 'Blues'),
                    breaks = c(0.001, 0.005, 0.01, 0.05, Inf),
                    labels = c(0.001, 0.005, 0.01, 0.05, '0.05 is insufficient')) +
  coord_fixed() +
  facet_wrap(. ~ cor_threshold) +
  theme(strip.background = element_rect(fill = 'gray82'),
        strip.text = element_text(colour = 'gray22', size = 14),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14))

panel_figure_7 <- ggpubr::ggarrange(heatmap_max_corr_across_delta_by_pseq,
                                    heatmap_min_seq_rate_correlation_threshold, nrow = 2,
                                    labels = 'AUTO')

pdf('../figures/figure_7.pdf',
    height = 7, width = 11)
plot(panel_figure_7)
dev.off()
png('../figures/figure_7.png',
    height = 7, width = 11, res = 350, units = 'in')
plot(panel_figure_7)
dev.off()
