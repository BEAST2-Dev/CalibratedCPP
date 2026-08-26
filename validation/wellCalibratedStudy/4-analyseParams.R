#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

# -----------------------------
# Paths
# -----------------------------
trueLogFolder <- "./lphyScripts"
resultFolder  <- "./xmls_parsed_results"
output_dir    <- "./xmls_summary"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# -----------------------------
# Extraction functions
# -----------------------------
extract_true_values <- function(filepath) {
  lines <- readLines(filepath, warn = FALSE)
  lines <- lines[grepl("\\S", lines) & !grepl("^#", lines)]

  if (length(lines) < 2) {
    warning(paste("File format error:", filepath))
    return(list(diversification = NA_real_, turnover = NA_real_))
  }

  delimiter <- if (grepl("\t", lines[1])) "\t" else "[[:space:]]+"
  params <- unlist(strsplit(trimws(lines[1]), delimiter))
  values <- as.numeric(unlist(strsplit(trimws(lines[2]), delimiter)))

  if (length(params) != length(values)) {
    warning(paste("Mismatch in", filepath))
    return(list(diversification = NA_real_, turnover = NA_real_))
  }

  get_param_value <- function(name) {
    idx <- which(params == name)
    if (length(idx) > 0) values[idx[1]] else NA_real_
  }

  list(diversification = get_param_value("diversification"),
       turnover     = get_param_value("turnover"))
}

extract_result_values <- function(filepath) {
  result_data <- read.table(filepath, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

  if (!("trace" %in% colnames(result_data))) {
    warning(paste("Missing 'trace' column in:", filepath))
    return(list(
      diversification = c(median=NA, HPD95.lower=NA, HPD95.upper=NA, stdev=NA),
      turnover     = c(median=NA, HPD95.lower=NA, HPD95.upper=NA, stdev=NA)
    ))
  }

  rownames(result_data) <- result_data$trace

  stats <- c("median", "HPD95.lower", "HPD95.upper", "stdev")
  parameters <- c("diversification", "turnover")

  get_stat <- function(param, stat) {
    if (param %in% colnames(result_data) && stat %in% rownames(result_data)) {
      as.numeric(result_data[stat, param])
    } else NA_real_
  }

  res <- lapply(parameters, function(p) {
    out <- sapply(stats, function(s) get_stat(p, s))
    names(out) <- stats
    out
  })
  names(res) <- parameters
  res
}

# -----------------------------
# Run one prefix (fixRootInference / fixStemInference)
# -----------------------------
run_prefix <- function(prefix) {
  # true logs: prefix_true-27.log
  true_logs <- list.files(trueLogFolder,
                          pattern = paste0("^", prefix, "-[0-9]+\\_true.log$"),
                          full.names = TRUE)
  if (length(true_logs) == 0) {
    warning(paste("No true logs found for", prefix, "in", trueLogFolder))
    return(NULL)
  }

  rs <- sub(paste0("^.*", prefix, "-([0-9]+)\\_true.log$"), "\\1", true_logs)
  rs <- as.character(sort(as.integer(rs)))

  process_r <- function(r_num) {
    true_log_file <- file.path(trueLogFolder, paste0(prefix, "-", r_num, "_true.log"))
    result_file   <- file.path(resultFolder,  paste0(prefix, "-", r_num, "_stats.log"))

    if (file.exists(true_log_file) && file.exists(result_file)) {
      true_values   <- extract_true_values(true_log_file)
      result_values <- extract_result_values(result_file)

      data.frame(
        r = as.integer(r_num),
        true_diversification = true_values$diversification,
        true_turnover = true_values$turnover,
        diversification_median = result_values$diversification["median"],
        diversification_HPD95_lower = result_values$diversification["HPD95.lower"],
        diversification_HPD95_upper = result_values$diversification["HPD95.upper"],
        diversification_stdev = result_values$diversification["stdev"],
        turnover_median = result_values$turnover["median"],
        turnover_HPD95_lower = result_values$turnover["HPD95.lower"],
        turnover_HPD95_upper = result_values$turnover["HPD95.upper"],
        turnover_stdev = result_values$turnover["stdev"]
      )
    } else {
      warning(paste("Missing files for", prefix, "r", r_num))
      NULL
    }
  }

  n_cores <- max(1, detectCores(logical = FALSE))
  all_results <- mclapply(rs, process_r, mc.cores = n_cores)
  all_results <- Filter(Negate(is.null), all_results)

  if (length(all_results) == 0) return(NULL)
  rbindlist(all_results, fill = TRUE)
}

# -----------------------------
# Do both prefixes, save separately
# -----------------------------
prefixes <- c("fixRootInference", "fixStemInference")

for (p in prefixes) {
  res <- run_prefix(p)
  if (!is.null(res) && nrow(res) > 0) {
    out_file <- file.path(output_dir, paste0(p, "_diversification_turnover_results.csv"))
    fwrite(res, out_file)
    cat("Saved:", out_file, "\n", sep = "")
  } else {
    cat("No results for:", p, "\n", sep = "")
  }
}