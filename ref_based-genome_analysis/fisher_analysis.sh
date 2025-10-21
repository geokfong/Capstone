suppressMessages(library(data.table))
library(parallel)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript run_fisher_analysis.R <group_name> <BT_summed.txt> <UT_summed.txt>")
}

group_name <- args[1]
bt_file <- args[2]
ut_file <- args[3]

bt_data <- fread(bt_file)
ut_data <- fread(ut_file)

setkey(bt_data, chr, pos)
setkey(ut_data, chr, pos)

combined <- merge(bt_data, ut_data, by = c("chr", "pos"), suffixes = c("_BT", "_UT"))

# Precompute coverage values (edit rates)
combined[, `:=`(
  edit_rate_BT = edit_count_BT / coverage_BT,
  edit_rate_UT = edit_count_UT / coverage_UT
)]

# Fisher's exact test (parallel)
fisher_p_value <- function(ref_BT, edit_BT, ref_UT, edit_UT) {
  mat <- matrix(
    c(ref_BT, edit_BT,
      ref_UT, edit_UT),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("BT", "UT"), c("Ref", "Edit"))
  )
  tryCatch(fisher.test(mat, alternative = "smaller")$p.value, error = function(e) 1)
}

num_cores <- detectCores() - 1
fisher_results <- mclapply(1:nrow(combined), function(i) {
  fisher_p_value(
    combined$ref_count_BT[i],
    combined$edit_count_BT[i],
    combined$ref_count_UT[i],
    combined$edit_count_UT[i]
  )
}, mc.cores = num_cores)

combined[, fisher_p := unlist(fisher_results)]
combined[, fdr := p.adjust(fisher_p, method = "fdr")]

output_dir <- "output"  # Make sure this directory exists or create it in your shell
sig_file <- file.path(output_dir, paste0(group_name, "_significant_sites_v3.txt"))
all_file <- file.path(output_dir, paste0(group_name, "_fisher_all_results_v3.txt"))

fwrite(combined[fisher_p < 0.05], file = sig_file, sep = "\t", quote = FALSE)
fwrite(combined, file = all_file, sep = "\t", quote = FALSE)
cat("Fisher's test complete.\n")
cat("Significant results saved to:", sig_file, "\n")
