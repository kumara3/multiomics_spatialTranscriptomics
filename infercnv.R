#!/usr/bin/env Rscript

################################################################################
#' InferCNV Analysis for Spatial Transcriptomics
#'
#' @description
#' This script performs copy number variation (CNV) inference on spatial
#' transcriptomics data using inferCNV. It identifies reference (normal) cells
#' based on marker gene expression and compares tumor cells against them.
#'
#' @param sample_name Name of the sample to analyze
#' @param output_dir Full path to output directory for results
#'
#' @details
#' The pipeline:
#'   1. Loads integrated Seurat object with all samples
#'   2. Identifies reference cells based on neuronal/oligodendrocyte markers
#'   3. Extracts expression matrix for the target sample
#'   4. Runs inferCNV to detect copy number variations
#'   5. Optionally adds CNV results back to Seurat object
#'
#' @references
#'   - inferCNV: https://github.com/broadinstitute/infercnv
#'   - Based on glioma spatial transcriptomics analysis
#'
#' @note
#'   Reference cells are identified as late-stage neurons/oligodendrocytes
#'   expressing ZIC1, GABRA1, and MBP (GCL_N population)
#'
#' @examples
#'   Rscript infercnv.R DMG1 /path/to/output
#'
#' @author Your Name
#' @date 2024
################################################################################

################################################################################
# Parse Command Line Arguments
################################################################################

args <- commandArgs(trailingOnly = TRUE)

usage <- "
Usage:
  Rscript infercnv.R <sample_name> <output_dir>

Arguments:
  sample_name   Name of the sample to analyze (must match 'sample' metadata)
  output_dir    Full path to output directory

Example:
  Rscript infercnv.R DMG1 /scratch/infercnv_results/DMG1
"

if (length(args) != 2) {
  cat(usage)
  stop("\nError: Incorrect number of arguments.\n", call. = FALSE)
}

################################################################################
# Configuration Parameters
################################################################################

#' @section Input/Output Parameters:
SAMPLE_NAME <- args[1]        #' Sample to analyze
OUTPUT_DIR <- args[2]          #' Output directory

#' @section File Paths:
FILE_PATHS <- list(
  integrated_rds = "/scratch/Gliome_spatial/github_scripts/multiomics_spatialTranscriptomics/banksy_harmony_integration/horizontal_integration_allSamples_banksy_harmony.rds",
  gene_order = "/scratch/Gliome_spatial/github_scripts/multiomics_spatialTranscriptomics/inferCNV/hg38_gencode_v27.txt"
)

#' @section Reference Cell Identification Parameters:
#' Reference cells are non-tumor cells (GCL_N) expressing neuronal and
#' oligodendrocyte markers. These represent normal tissue for comparison.
REFERENCE_PARAMS <- list(
  markers = c("ZIC1", "GABRA1", "MBP"),  #' Marker genes for reference cells
  expression_threshold = 0.5,             #' Minimum expression level
  clusters = c(1, 3, 18),                 #' Seurat clusters containing GCL_N
  annotation_label = "GCL_N"              #' Label for reference cells
)

#' @section InferCNV Parameters:
INFERCNV_PARAMS <- list(
  cutoff = 0.1,                           #' Expression cutoff threshold
  cluster_by_groups = TRUE,               #' Cluster cells by annotation groups
  analysis_mode = "subclusters",          #' Analysis mode (subclusters, samples, cells)
  denoise = TRUE,                         #' Apply denoising
  hmm = TRUE,                             #' Use Hidden Markov Model
  num_threads = 8,                        #' Number of parallel threads
  smooth_method = "runmeans",             #' Smoothing method
  output_format = "pdf"                   #' Output format (pdf, png)
)

#' @section CNV Detection Thresholds:
CNV_THRESHOLDS <- list(
  scale_min = 0.7,                        #' Lower bound for CNV scaling
  scale_max = 1.3,                        #' Upper bound for CNV scaling
  outlier_lower = 0.7,                    #' Lower outlier threshold (deletions)
  outlier_upper = 1.3                     #' Upper outlier threshold (amplifications)
)

#' @section Tumor Subclustering:
SUBCLUSTER_PARAMS <- list(
  method = "leiden",                      #' Subclustering method (leiden, louvain, random_trees)
  top_n = 10                              #' Top N genes per CNV region to add to Seurat
)

#' @section Cell Type Annotations:
#' Based on published glioma spatial transcriptomics signatures:
#' - GCL_N (late peak): Neurons expressing ZIC1, GABRA1, MBP
#' - GCL_TI (one peak): Tumor-invasive cells expressing TNC, HOPX, PTPRZ1,
#'   VIM, CLU, LGALS1, CD44, SPARC, GAP43
CELL_TYPE_INFO <- list(
  reference = "GCL_N: Late-stage neurons/oligodendrocytes",
  tumor = "Tumor cells including GCL_TI invasive population"
)

################################################################################
# Setup Environment
################################################################################

#' @section System Configuration:
system("ulimit -n 10000")
options(echo = TRUE)
options(future.globals.maxSize = 4 * 1024^3)  # 4GB limit

#' @section Required Packages:
required_packages <- c(
  "Seurat",           # Single-cell analysis
  "SeuratData",       # Data utilities
  "dplyr",            # Data manipulation
  "SeuratWrappers",   # Additional wrappers
  "ggplot2",          # Plotting
  "gridExtra",        # Multi-panel plots
  "infercnv"          # CNV inference
)

invisible(lapply(required_packages, library, character.only = TRUE))

# Create output directory
dir.create(OUTPUT_DIR, showWarnings = TRUE, recursive = TRUE)

################################################################################
# Load Integrated Data
################################################################################

cat("\n", rep("=", 80), "\n", sep = "")
cat("InferCNV Analysis Pipeline\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("Step 1/6: Loading integrated Seurat object...\n")

if (!file.exists(FILE_PATHS$integrated_rds)) {
  stop(sprintf("Error: Integrated RDS file not found at:\n  %s\n", 
               FILE_PATHS$integrated_rds))
}

all_samples_obj <- readRDS(FILE_PATHS$integrated_rds)
DefaultAssay(all_samples_obj) <- "SCT"

cat(sprintf("  Loaded object with %d cells across %d samples\n",
            ncol(all_samples_obj),
            length(unique(all_samples_obj$sample))))

# Verify sample exists
if (!SAMPLE_NAME %in% unique(all_samples_obj$sample)) {
  stop(sprintf("\nError: Sample '%s' not found in integrated object.\n", 
               SAMPLE_NAME))
}

cat(sprintf("  Target sample: %s\n\n", SAMPLE_NAME))

################################################################################
# Identify Reference and Tumor Cells
################################################################################

cat("Step 2/6: Identifying reference (non-tumor) cells...\n")

# Set identities to clusters
Idents(all_samples_obj) <- "seurat_clusters"

# Get all cells from target sample
sample_cells <- WhichCells(
  all_samples_obj,
  expression = sample == SAMPLE_NAME
)

cat(sprintf("  Total cells in %s: %d\n", SAMPLE_NAME, length(sample_cells)))

# Identify reference cells (GCL_N population)
# These are mature neurons/oligodendrocytes expressing ZIC1, GABRA1, and MBP
ref_expression <- sprintf(
  "%s > %s & %s > %s & %s > %s & sample == '%s'",
  REFERENCE_PARAMS$markers[1], REFERENCE_PARAMS$expression_threshold,
  REFERENCE_PARAMS$markers[2], REFERENCE_PARAMS$expression_threshold,
  REFERENCE_PARAMS$markers[3], REFERENCE_PARAMS$expression_threshold,
  SAMPLE_NAME
)

ref_cells <- WhichCells(
  all_samples_obj,
  expression = ref_expression,
  idents = REFERENCE_PARAMS$clusters,
  slot = "data"
)

cat(sprintf("  Reference cells (GCL_N): %d (%.1f%%)\n",
            length(ref_cells),
            100 * length(ref_cells) / length(sample_cells)))

if (length(ref_cells) == 0) {
  stop("\nError: No reference cells identified. Check marker expression and clusters.\n")
}

# Annotate cells
all_samples_obj$annotation <- "tumor"
all_samples_obj$annotation[ref_cells] <- REFERENCE_PARAMS$annotation_label

cat("  Cell annotations:\n")
cat(sprintf("    - Tumor cells: %d\n", 
            sum(all_samples_obj$annotation[sample_cells] == "tumor")))
cat(sprintf("    - Reference cells (%s): %d\n\n", 
            REFERENCE_PARAMS$annotation_label,
            sum(all_samples_obj$annotation[sample_cells] == 
                REFERENCE_PARAMS$annotation_label)))

################################################################################
# Prepare InferCNV Inputs
################################################################################

cat("Step 3/6: Preparing inferCNV inputs...\n")

# Extract expression matrix for target sample
expr_matrix <- all_samples_obj@assays$SCT@counts[, sample_cells] %>%
  as.matrix() %>%
  as.data.frame()

cat(sprintf("  Expression matrix: %d genes x %d cells\n",
            nrow(expr_matrix), ncol(expr_matrix)))

# Extract cell annotations
cell_annotations <- all_samples_obj@meta.data[sample_cells, ] %>%
  dplyr::select(annotation)

cat(sprintf("  Cell annotations prepared\n"))

# Verify gene order file exists
if (!file.exists(FILE_PATHS$gene_order)) {
  stop(sprintf("Error: Gene order file not found at:\n  %s\n", 
               FILE_PATHS$gene_order))
}

cat(sprintf("  Gene order file: %s\n\n", basename(FILE_PATHS$gene_order)))

################################################################################
# Create InferCNV Object
################################################################################

cat("Step 4/6: Creating inferCNV object...\n")

infercnv_obj <- infercnv::CreateInfercnvObject(
  raw_counts_matrix = expr_matrix,
  gene_order_file = FILE_PATHS$gene_order,
  annotations_file = cell_annotations,
  ref_group_names = c(REFERENCE_PARAMS$annotation_label)
)

cat("  InferCNV object created successfully\n\n")

################################################################################
# Run InferCNV Analysis
################################################################################

cat("Step 5/6: Running inferCNV analysis...\n")
cat("  This may take some time depending on dataset size...\n\n")

infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = INFERCNV_PARAMS$cutoff,
  out_dir = OUTPUT_DIR,
  cluster_by_groups = INFERCNV_PARAMS$cluster_by_groups,
  final_scale_limits = c(CNV_THRESHOLDS$scale_min, CNV_THRESHOLDS$scale_max),
  outlier_lower_bound = CNV_THRESHOLDS$outlier_lower,
  outlier_upper_bound = CNV_THRESHOLDS$outlier_upper,
  analysis_mode = INFERCNV_PARAMS$analysis_mode,
  tumor_subcluster_partition_method = SUBCLUSTER_PARAMS$method,
  output_format = INFERCNV_PARAMS$output_format,
  denoise = INFERCNV_PARAMS$denoise,
  HMM = INFERCNV_PARAMS$hmm,
  num_threads = INFERCNV_PARAMS$num_threads,
  smooth_method = INFERCNV_PARAMS$smooth_method
)

cat("  InferCNV analysis complete\n\n")

################################################################################
# Save Results
################################################################################

cat("Step 6/6: Saving results...\n")

# Save inferCNV object
infercnv_rds <- file.path(OUTPUT_DIR, sprintf("%s_infercnv_obj.rds", SAMPLE_NAME))
saveRDS(infercnv_obj, infercnv_rds)
cat(sprintf("  Saved inferCNV object: %s\n", basename(infercnv_rds)))

# Optional: Add CNV results back to Seurat object
# Uncomment and modify the path below if you want to integrate CNV data
# cat("\nAdding CNV results to Seurat object...\n")
# seurat_obj_updated <- infercnv::add_to_seurat(
#   infercnv_output_path = OUTPUT_DIR,
#   seurat_obj = all_samples_obj,
#   top_n = SUBCLUSTER_PARAMS$top_n
# )
# 
# seurat_rds <- file.path(OUTPUT_DIR, sprintf("%s_seurat_with_cnv.rds", SAMPLE_NAME))
# saveRDS(seurat_obj_updated, seurat_rds)
# cat(sprintf("  Saved updated Seurat object: %s\n", basename(seurat_rds)))

################################################################################
# Analysis Summary
################################################################################

cat("\n", rep("=", 80), "\n", sep = "")
cat("Analysis Complete\n")
cat(rep("=", 80), "\n", sep = "")

cat(sprintf("Sample: %s\n", SAMPLE_NAME))
cat(sprintf("Output directory: %s\n\n", OUTPUT_DIR))

cat("Cell populations:\n")
cat(sprintf("  - Reference (%s): %d cells\n",
            REFERENCE_PARAMS$annotation_label,
            sum(all_samples_obj$annotation[sample_cells] == 
                REFERENCE_PARAMS$annotation_label)))
cat(sprintf("  - Tumor: %d cells\n",
            sum(all_samples_obj$annotation[sample_cells] == "tumor")))

cat("\nReference markers used:\n")
for (marker in REFERENCE_PARAMS$markers) {
  cat(sprintf("  - %s (threshold: %.1f)\n", 
              marker, REFERENCE_PARAMS$expression_threshold))
}

cat("\nOutput files:\n")
cat("  - infercnv.pdf (main CNV heatmap)\n")
cat("  - infercnv.observations.txt (CNV residuals for tumor cells)\n")
cat("  - infercnv.references.txt (CNV residuals for reference cells)\n")
cat(sprintf("  - %s (saved inferCNV object)\n", basename(infercnv_rds)))

cat("\nNote: CNV estimations are based on residual expression found in:\n")
cat("  - infercnv_obj@expr.data slot\n")
cat("  - infercnv.observations.txt and infercnv.references.txt files\n")

cat("\nFor more information:\n")
cat("  - https://github.com/broadinstitute/infercnv\n")
cat("  - GitHub issues: #609, #208, #338, #434\n")

cat(rep("=", 80), "\n", sep = "")

################################################################################
# Additional Notes
################################################################################

#' @section Biological Context:
#' 
#' Cell Type Signatures (from literature):
#' 
#' GCL_N (Late Peak) - Reference Population:
#'   - Neuronal markers: ZIC1, GABRA1
#'   - Mature oligodendrocyte marker: MBP
#'   - Represents normal, non-malignant cells
#' 
#' GCL_TI (One Peak) - Tumor-Invasive Population:
#'   - Radial glia markers: TNC, HOPX, PTPRZ1, VIM
#'   - Astrocyte markers: CLU, LGALS1
#'   - Migration/network formation: CD44, SPARC, GAP43
#'   - Represents invasive glioma cells
#' 
#' CNV Interpretation:
#'   - Values > 1.3: Amplifications (gains)
#'   - Values < 0.7: Deletions (losses)
#'   - Values ~1.0: Normal copy number
#'   - Tumor cells show CNV patterns distinct from reference cells