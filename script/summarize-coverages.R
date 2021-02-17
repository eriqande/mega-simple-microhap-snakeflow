

library(tidyverse)


cov <- read_tsv("coverage/all_coverages.tsv") %>%
  select(
    genome_condition,
    sample,
    target,
    num_reads
    ) %>%
  pivot_wider(
    names_from = genome_condition,
    values_from = num_reads
  ) %>%
  mutate(
    fract_retained = full / thin
  ) %>%
  group_by(sample) %>%
  mutate(tot_mapped_indiv = sum(full)) %>%
  ungroup() %>%
  arrange(
    tot_mapped_indiv,
    sample
  ) %>%
  mutate(indiv_depth_order = as.integer(factor(sample, levels = unique(sample))))

ggplot(cov, aes(x = indiv_depth_order, y = full)) +
  geom_col(colour = "blue", fill = "blue") +
  facet_wrap(~ target)
