# redirect output and messages/errors to the log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)

idxfile <- snakemake@input[[1]]

output <- snakemake@output



#### First arrange the data and spit out a table ###

idx <- read_tsv(
  idxfile,
  col_names = c("sample", "chrom", "length", "num_mapped", "num_unmapped")
)

# if snakemake@wildcards$map_type is "fullg" or "fullg-extracted",
# then we must toss any of the chrom's that have zero reads across
# all individuals.
if(snakemake@wildcards$map_type %in% c("fullg", "fullg-extracted")) {
  idx <- idx %>%
    group_by(chrom) %>%
    filter(sum(num_mapped) > 0) %>%
    ungroup()
}

# get samples and markers in order of most to least reads
ord_samp <- idx %>%
  group_by(sample) %>%
  summarise(tot_reads = sum(num_mapped)) %>%
  arrange(desc(tot_reads)) %>%
  pull(sample)
ord_markers <- idx %>%
  group_by(chrom) %>%
  summarise(tot_reads = sum(num_mapped)) %>%
  arrange(desc(tot_reads)) %>%
  pull(chrom)

idx2 <- idx %>%
  mutate(
    sample_f = factor(sample, levels = ord_samp),
    chrom_f = factor(chrom, levels = ord_markers)
  ) %>%
  arrange(sample_f, chrom_f)

# now, make a table with samples in rows and markers in columns
idx_table <- idx2 %>%
  select(sample_f, chrom_f, num_mapped) %>%
  pivot_wider(id_cols = sample_f, names_from = chrom_f, values_from = num_mapped) %>%
  rename(sample = sample_f) %>%
  mutate(sample = as.character(sample)) # so we can join more data back to it

# add some sample meta data onto that
samples_file <- file.path(
  snakemake@wildcards$run_dir,
  "samples.csv"
)
samples_tib <- read_csv(samples_file) %>%
  mutate(
    Marker_Set = snakemake@wildcards$marker_set,
    Map_Type = snakemake@wildcards$map_type,
    Genome_or_Target_Fasta = snakemake@wildcards$tf_or_gen
  ) %>%
  select(Marker_Set, Map_Type, Genome_or_Target_Fasta, NMFS_DNA_ID, Sample_ID, Sample_Project, sample)

idx_table2 <- left_join(idx_table, samples_tib, by = "sample") %>%
  select(Marker_Set, Map_Type, Genome_or_Target_Fasta, NMFS_DNA_ID, Sample_ID, Sample_Project, sample, everything())

# write that out
out_table <- output$csv
write_csv(idx_table2, file = out_table)


#### Now, use idx2 to make a heatmap figure ####

# to make the figure we want to factorize sample in reverse.
# So that the ones with the most reads are at the top.
idx3 <- idx2 %>%
  mutate(sample_f = factor(sample, levels = rev(ord_samp)))

# get a title string for it:
title <- paste(
  "Species: ",
  snakemake@wildcards$species_dir,
  "\nRun Dir.: ",
  snakemake@wildcards$run_dir,
  "\nMarker Set: ",
  snakemake@wildcards$marker_set,
  "\nMap Type: ",
  snakemake@wildcards$map_type,
  "\nGenome/Target_Fasta: ",
  snakemake@wildcards$tf_or_gen,
  "\nNum Markers: ",
  n_distinct(idx3$chrom),
  "\nNum Samples: ",
  n_distinct(idx3$sample),
  sep = "",
  collapse = ""
)
# now, ggplot that dude
g <- ggplot(idx3, aes(x = chrom_f, y = sample_f, fill = num_mapped)) +
  geom_tile() +
  scale_fill_viridis_c(trans = 'log', option = "plasma") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle(title) +
  theme(title = element_text(size = 35)) +
  xlab("Marker Name, Amplicon Name, or Chromosome") +
  ylab("Sample")

ggsave(g, filename = snakemake@output$heatmap, width = 24, height = 24)


#### Make a samples bar plot ####

g_samp <- idx3 %>%
  group_by(sample_f) %>%
  summarise(tot_reads = sum(num_mapped)) %>%
  ggplot(., aes(x=sample_f, y=tot_reads)) +
  geom_col(fill = "lightblue", colour = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle(title) +
  theme(title = element_text(size = 20)) +
  xlab("Sample") +
  ylab("Total reads mapped")

ggsave(g_samp, filename = snakemake@output$samp_bars, width = 24, height = 10)


### Make a markers reads-mapped plot


g_mark <- idx2 %>%
  group_by(chrom) %>%
  summarise(tot_reads = sum(num_mapped)) %>%
  arrange(tot_reads) %>%
  mutate(chrom_f = factor(chrom, levels = chrom)) %>%
  ggplot(., aes(x=chrom_f, y=tot_reads)) +
  geom_col(fill = "goldenrod", colour = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle(title) +
  theme(title = element_text(size = 20)) +
  xlab("Marker Name, Amplicon Name, or Chromosome") +
  ylab("Total reads mapped")


ggsave(g_mark, filename = snakemake@output$marker_bars, width = 24, height = 10)
