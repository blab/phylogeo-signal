library(tidyverse)

## Load dataframe with probability of a contact of occurring within the same age group
df_proba_stays_within_age_group <- read_csv('../results/df_proba_stays_within_age_group_WA.csv')

## Load dataframe with probability of a contact of occurring within the same age group 
## (by age group in decade)
df_proba_stays_within_age_group_by_age <- read_csv('../results/df_proba_stays_within_age_group_by_age_decade_WA.csv')

## Plot assortativity estimate at the population level as a function of binning width
plt_assortativity_estimate_age_population_level <- df_proba_stays_within_age_group %>% 
  ggplot(aes(x = as.factor(binning_width), y = proba_stays_within_group)) +
  geom_bar(stat = 'identity', fill = 'darkslateblue') +
  geom_text(aes(label = round(proba_stays_within_group, 2)), vjust = -0.2) +
  scale_x_discrete(name = 'Binning width') +
  scale_y_continuous(name = 'Population assortativity estimate', 
                     expand = expansion(mult = c(0., 0.08))) +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))

#png('../figures/supplementary-figures/assortativity_estimates_age.png', 
#    height = 3, width = 4, res = 350, units = 'in')
plot(plt_assortativity_estimate_age_population_level)
#dev.off()
#pdf('../figures/supplementary-figures/assortativity_estimates_age.pdf', 
#    height = 3, width = 4)
#plot(plt_assortativity_estimate_age_population_level)
#dev.off()

## Plot assortativity estimate at the population level as the age group in decade
plt_prop_contacts_within_age_group_across_age_group_with_text <- df_proba_stays_within_age_group_by_age %>% 
  ggplot(aes(x = age_aggregated_i, y = prop_contacts_within_i)) +
  geom_bar(stat = 'identity', fill = 'darkslateblue') +
  geom_text(aes(label = round(prop_contacts_within_i, 2)), vjust = -0.2) +
  scale_x_discrete(name = 'Age group') +
  scale_y_continuous(name = 'Proportion', expand = expansion(mult = c(0., 0.05))) +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1., hjust = 1.))

# png('../figures/supplementary-figures/prop_contacts_within_age_group_across_age_group.png', 
#     height = 3, width = 4, res = 350, units = 'in')
plot(plt_prop_contacts_within_age_group_across_age_group_with_text)
# dev.off()
# pdf('../figures/supplementary-figures/prop_contacts_within_age_group_across_age_group.pdf', 
#     height = 3, width = 4)
# plot(plt_prop_contacts_within_age_group_across_age_group_with_text)
# dev.off()



plot(plt_prop_contacts_within_age_group_across_age_group_with_text)
