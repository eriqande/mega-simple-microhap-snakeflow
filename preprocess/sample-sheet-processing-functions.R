library(tidyverse)

#' Turn Sample_Sheet.csv into samples.csv and units.csv
#' in the same directory
#' @param path  The path to the Sample_Sheet.csv file.  This is a file
#' that comes out of the MiSeq.
create_samples_and_units <- function(path) {

  # first, find where to start reading the file
  lines <- read_lines(path, n_max = 100)
  skip <- min(which(str_detect(lines, "^\\[Data\\]")))

  S <- read_csv(path, skip = skip)

}
