# redirect output and messages/errors to the log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)

input <- snakemake@input[[1]]
output <- snakemake@output[[1]]


cov <- read_tsv(input) %>%
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

g <- ggplot(cov, aes(x = indiv_depth_order, y = full)) +
  geom_col(colour = "blue", fill = "blue") +
  facet_wrap(~ target) +
  xlab("Indivs ordered by total read depth")


# get number of targets
NT <- n_distinct(cov$target)

# set the width and height according to NT
nt_width <- 5 + 2 * sqrt(NT)
nt_height <- 5 + 2 * sqrt(NT)

ggsave(
  g,
  filename = output,
  width = nt_width,
  height = nt_height
)
