############################################################
# Script: compile_source_data.R
# Purpose: Compile every CSV under source_data/{main_figs,supp_figs}/
#          into one Excel workbook (source_data/SourceData.xlsx),
#          one tab per CSV/panel, main figures first then
#          supplementary, in manuscript figure-number order.
#
# Run this AFTER all figure scripts (and CIcontact_sourceData.R)
# have populated source_data/.
############################################################

packages <- c("writexl")
installed <- rownames(installed.packages())
for (p in packages) {
  if (!(p %in% installed)) install.packages(p, dependencies = TRUE)
}
lapply(packages, library, character.only = TRUE)

sanitize_sheet_name <- function(x) {
  x <- gsub("[\\[\\]:*?/\\\\]", "_", x)
  substr(x, 1, 31)
}

# Explicit figure-number ordering; anything not matched falls at the
# end, alphabetically. Update these prefixes if figure numbers change.
main_order <- c("fig2", "secondOrderIxn_fig3", "AMMI_fig4", "CIcontact_fig5")
supp_order <- c("buildCorr_figS1", "phenoCorr_figS2", "phenoPCA_figS3", "set9Corr_figS4",
                "phyloSignal_Sfig2", "Sfig1", "mothersCurse_mch", "ixn_subsample",
                "weightDensity", "survPlot")

order_files <- function(files, order_prefixes) {
  key <- sapply(files, function(f) {
    idx <- which(sapply(order_prefixes, function(p) startsWith(basename(f), p)))
    if (length(idx) == 0) length(order_prefixes) + 1 else idx[1]
  })
  files[order(key, basename(files))]
}

main_files <- order_files(
  list.files("source_data/main_figs", full.names = TRUE, pattern = "\\.csv$"),
  main_order
)
supp_files <- order_files(
  list.files("source_data/supp_figs", full.names = TRUE, pattern = "\\.csv$"),
  supp_order
)

all_files <- c(main_files, supp_files)

if (length(all_files) == 0) {
  stop("No source_data CSVs found. Run the figure scripts first to populate source_data/.")
}

sheets <- lapply(all_files, read.csv)
names(sheets) <- sanitize_sheet_name(tools::file_path_sans_ext(basename(all_files)))
names(sheets) <- make.unique(names(sheets))  # guard against truncation collisions

write_xlsx(sheets, "source_data/SourceData.xlsx")
cat("Wrote source_data/SourceData.xlsx with", length(sheets), "sheets:\n")
cat(paste(" -", names(sheets)), sep = "\n")
