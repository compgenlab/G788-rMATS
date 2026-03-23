#!/usr/bin/env Rscript
#
# 05-find-interesting.R
# Filter rMATS output for high-confidence, biologically interesting splicing events.
#
# Usage: Rscript 05-find-interesting.R
#

library(dplyr)

# --- Configuration -----------------------------------------------------------

rmats_dir   <- "rmats_out/"
output_dir  <- "rmats_filtered/"

min_delta_psi  <- 0.2    # minimum |ΔPSI|
max_fdr        <- 0.05   # FDR cutoff
min_reads      <- 20     # minimum IJC + SJC per group
max_sd         <- 0.15   # maximum within-group PSI standard deviation

event_types <- c("SE", "A5SS", "A3SS", "MXE", "RI")

# -----------------------------------------------------------------------------

dir.create(output_dir, showWarnings = FALSE)

for (event in event_types) {
  infile <- file.path(rmats_dir, paste0(event, ".MATS.JCEC.txt"))
  if (!file.exists(infile)) {
    message("Skipping ", event, " — file not found: ", infile)
    next
  }

  df <- read.table(infile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  # Within-group PSI standard deviation
  df$sd1 <- sapply(strsplit(df$IncLevel1, ","), function(x) sd(as.numeric(x), na.rm = TRUE))
  df$sd2 <- sapply(strsplit(df$IncLevel2, ","), function(x) sd(as.numeric(x), na.rm = TRUE))

  interesting <- df %>%
    filter(
      FDR < max_fdr,
      abs(IncLevelDifference) > min_delta_psi,
      IJC_SAMPLE_1 + SJC_SAMPLE_1 >= min_reads,
      IJC_SAMPLE_2 + SJC_SAMPLE_2 >= min_reads,
      sd1 < max_sd,
      sd2 < max_sd
    ) %>%
    arrange(FDR)

  outfile <- file.path(output_dir, paste0(event, "_interesting.txt"))
  write.table(interesting, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
  message(event, ": ", nrow(interesting), " interesting events -> ", outfile)
}
