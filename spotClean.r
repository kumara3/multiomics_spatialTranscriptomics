#!/usr/bin/env Rscript

# Libraries
library(SpotClean)
library(Seurat)
library(ggplot2)
library(tibble)
library(readr)
library(purrr)
library(cowplot)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript process_visium_sample.R <sample_id> <matrix_dir> <tissue_csv> <image_path> <scale_factor_json>")
}

sample_id <- args[1]
matrix_dir <- args[2]
tissue_csv <- args[3]
image_path <- args[4]
scale_factor_json <- args[5]

dir.create("results", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)

log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", sample_id, msg))
  flush.console()
}

tryCatch({
  log_msg("Starting")
  raw_mat <- read10xRaw(matrix_dir)
  slide_info <- read10xSlide(tissue_csv, image_path, scale_factor_json)
  slide_obj <- createSlide(count_mat = raw_mat, slide_info = slide_info)

  gene_to_plot <- rownames(raw_mat)[1]
  raw_plot <- visualizeHeatmap(slide_obj, gene_to_plot) + ggtitle(paste("Raw:", gene_to_plot))
  ggsave(sprintf("plots/%s_raw_heatmap.png", sample_id), raw_plot, width = 6, height = 5)

  decont_obj <- spotclean(slide_obj)

  contam_plot <- visualizeHeatmap(decont_obj, metadata(decont_obj)$contamination_rate,
                                  logged = FALSE,
                                  legend_title = "Contamination Rate",
                                  legend_range = c(0, 1)) +
    ggtitle("Estimated contamination")
  ggsave(sprintf("plots/%s_contam_rate.png", sample_id), contam_plot, width = 6, height = 5)

  clean_plot <- visualizeHeatmap(decont_obj, gene_to_plot) + ggtitle(paste("Clean:", gene_to_plot))
  ggsave(sprintf("plots/%s_clean_heatmap.png", sample_id), clean_plot, width = 6, height = 5)

  # You can toggle between saving as Seurat or SpotClean
  saveRDS(decont_obj, file = file.path("results", paste0(sample_id, "_spotclean.rds")))
  log_msg("Saved SpotClean object")

}, error = function(e) {
  log_msg(paste("ERROR:", e$message))
})
