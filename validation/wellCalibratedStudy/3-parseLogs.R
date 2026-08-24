### INSTALLATION TraceR
library("devtools")
remove.packages("TraceR")
devtools::install_github("walterxie/TraceR")
library("TraceR")

library(tools)
library(TraceR)
library(ape)
require(phytools)
require(tidyverse)

mcmc_path <- "./xmls"
output_path <- "./xmls_parsed_results"

if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)  # Creates the directory along with any necessary parent directories
}

burnin <- 0.1
files <- dir(path=mcmc_path, pattern=".log")

for (f in files) {
  mcmc_log <- readMCMCLog(file.path(mcmc_path, f))
  traces <- getTraces(mcmc_log, burn.in = burnin)
  stats <- analyseTraces(traces)
  out <- paste0(file_path_sans_ext(basename(f)), "_stats.log")
  write_tsv(stats, file.path(output_path, out))
}

