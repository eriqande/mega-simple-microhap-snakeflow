library(tidyverse)

#' Turn Sample_Sheet.csv into samples.csv and units.csv
#' in the same directory
#' @param path  The path to the Sample_Sheet.csv file.  This is a file
#' that comes out of the MiSeq.
#' @param NMFS_DNA_ID_from_Sample_ID  If this is FALSE (the default)
#' then the NMFS_DNA_ID is found as the part of the Sample_Plate column
#' before the first underscore.  (This is how our group tends to
#' record things).  For other uses, the Sample_Plate might be the same
#' for all individuals, and their unique identifiers will be stored
#' in the Sample_ID column.  Setting this option to TRUE will use the
#' Sample_ID as the NMFS_DNA_ID.
create_samples_and_units <- function(
  path,
  NMFS_DNA_ID_from_Sample_ID = FALSE
) {

  # first, find where to start reading the file
  lines <- read_lines(path, n_max = 100)
  skip <- min(which(str_detect(lines, "^\\[Data\\]")))

  # first make and retain only the columns we will need going forward
  T1 <- read_csv(path, skip = skip) %>%
    filter(!is.na(Sample_ID)) # if there were empty lines included in the file, toss them this way

  if(NMFS_DNA_ID_from_Sample_ID == TRUE) {
    T2 <- T1 %>%
      mutate(NMFS_DNA_ID = Sample_ID)
  } else {
    T2 <- T1 %>%
      mutate(NMFS_DNA_ID = str_split_fixed(Sample_Plate, "_", n = 3)[,1])
  }

  S <- T2 %>%
    mutate(
      sample = str_c("s", sprintf("%04d", 1:n()))
    ) %>%
    rename(
      Marker_Sets = Description
    ) %>%
    select(sample, Marker_Sets, NMFS_DNA_ID, Sample_ID, Sample_Name, Sample_Project)


  # now, look in the raw directory to get the fastq paths
  rawdir <- file.path(dirname(path), "raw")
  rawfiles <- dir(rawdir)
  r1 <- rawfiles[str_detect(rawfiles, "_L00[0-9]_R1_00[0-9]\\.fastq.gz")]
  r2 <- rawfiles[str_detect(rawfiles, "_L00[0-9]_R2_00[0-9]\\.fastq.gz")]
  r1tib <- tibble(file = r1) %>%
    mutate(
      Sample_Name = str_split_fixed(file, "_", n = 2)[,1],
      fq1 = file
    ) %>%
    select(-file)
  r2tib <- tibble(file = r2) %>%
    mutate(
      Sample_Name = str_split_fixed(file, "_", n = 2)[,1],
      fq2 = file
    ) %>%
    select(-file)

  samples <- S %>%
    left_join(r1tib, by = "Sample_Name") %>%
    left_join(r2tib, by = "Sample_Name")

  #### some light error checking ####
  missers <- samples %>%
    filter(is.na(fq1) | is.na(fq2))

  if(nrow(missers) > 0) {
    miss_str <- paste(missers$Sample_ID, collapse = ", ")
    stop("FASTQs not found for: ", miss_str)
  }

  both <- c(samples$fq1, samples$fq2)
  dupies <- both[duplicated(both)]
  if(length(dupies) < 0) {
    dup_str <- paste(dupies, collapse = ", ")
    stop("Bad News! FASTQ paths appearing more than once: ", dup_str)
  }


  #### Now, make units  ####
  units <- samples %>%
    mutate(Markers = str_split(Marker_Sets, "\\s*,\\s*")) %>%
    unnest(cols = Markers) %>%
    select(-Marker_Sets) %>%
    select(sample, Markers)

  #### Finally, write out samples and units ####
  write_csv(samples, file = file.path(dirname(path), "samples.csv"))
  write_csv(units, file = file.path(dirname(path), "units.csv"))
}
