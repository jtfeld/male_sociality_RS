library(dplyr)
library(lmerTest)
library(ggplot2)

yearly_data = readRDS(file = "./data/yearly_strong_ties_coal_betw_min_8_coal_dyads.rds")

# various ways to inspect the relationship between 
# yearly measures of strong ties and coalitionary betweenness:

summary(lm(coal_betw ~ n_high_dai_arr, data = yearly_data))

summary(lmer(coal_betw ~ n_high_dai_arr + (1|id), data = yearly_data))

summary(lmer(coal_betw ~ 
               n_high_dai_arr + 
               (1|id), 
             data = yearly_data[!yearly_data$alpha, ]))

# plotting the relationship:

alphastat = as_labeller(c(`FALSE` = "Subordinate males", `TRUE` = "Alpha males"))

yearly_data %>%
  ggplot(aes(x = n_high_dai_arr, y = coal_betw)) + geom_point() + 
  geom_smooth(method = "lm") +
  # or you could just use a different smooth, e.g. replace above with:
  # geom_smooth(se = F) +
  facet_wrap(~ alpha, labeller = alphastat) + 
  xlab("Count strong assn. ties") + 
  ylab("Coalitionary betweenness")
