#!/usr/bin/env Rscript

################################################################################
#' Spatial Transcriptomics Analysis with BANKSY
#'
#' @description
#' This script performs spatial transcriptomics analysis on multiple samples
#' using BANKSY for spatial feature extraction. Each sample is processed
#' independently with clustering, visualization, and marker identification.
#'
#' @param matrix_list_file Full path to text file containing sample information
#'   Format: sample_name,path_to_spatial_data (one per line)
#' @param output_dir Full path to output directory for results
#'
#' @details
#' The script performs the following steps for each sample:
#'   1. Load 10X Visium spatial data
#'   2. Quality control filtering
#'   3. SCTransform normalization
#'   4. BANKSY spatial analysis
#'   5. PCA, clustering, and UMAP
#'   6. Visualization and marker gene identification
#'
#' @references
#'   BANKSY: https://prabhakarlab.github.io/Banksy/
#'
#' @examples
#'   Rscript seuratSpat.R samples.txt /path/to/output
#'
#' @author ashwani kumar
#' @date 2025-07-01
################################################################################

################################################################################
# Parse Command Line Arguments
################################################################################

args <- commandArgs(trailingOnly = TRUE)

usage <- "
Usage:
  Rscript seuratSpat.R <matrix_list_file> <output_dir>

Arguments:
  matrix_list_file  Text file with sample info (format: sample_name,data_path)
  output_dir        Output directory for results

Example:
  Rscript seuratSpat.R samples.txt /path/to/output
"

if (length(args) != 2) {
  cat(usage)
  stop("\nError: Incorrect number of arguments.\n", call. = FALSE)
}

################################################################################
# Configuration Parameters
################################################################################

#' @section Input/Output Parameters:
MATRIX_LIST_FILE <- args[1]  #' Path to sample list file
OUTPUT_DIR <- args[2]         #' Output directory

#' @section Quality Control Parameters:
QC_PARAMS <- list(
  min_features = 200,         #' Minimum features per spot
  max_mt_percent = 25         #' Maximum mitochondrial percentage
)

#' @section BANKSY Parameters:
BANKSY_PARAMS <- list(
  lambda = 0.2,               #' Spatial weight parameter (0=nonspatial, 1=spatial)
  k_geom = 15,                #' Number of spatial neighbors
  assay = "SCT",              #' Assay to use
  slot = "data",              #' Data slot to use
  features = "variable"       #' Features to use
)

#' @section Dimensionality Reduction Parameters:
REDUCTION_PARAMS <- list(
  n_pcs = 30,                 #' Number of principal components
  dims = 1:30                 #' Dimensions to use for clustering/UMAP
)

#' @section Clustering Parameters:
CLUSTER_PARAMS <- list(
  resolution = 0.5,           #' Clustering resolution
  cluster_name = "banksy_cluster"  #' Name for cluster metadata
)

#' @section Marker Gene Parameters:
MARKER_PARAMS <- list(
  only_pos = TRUE,            #' Only positive markers
  assay = "SCT",              #' Assay for marker identification
  min_pct = 0.25,             #' Minimum percentage of cells
  logfc_threshold = 0.25      #' Log fold-change threshold
)

#' @section Visualization Parameters:
PLOT_PARAMS <- list(
  height = 780,               #' Plot height in pixels
  width = 980,                #' Plot width in pixels
  pt_size = 0.25,             #' Point size in UMAP
  label_size = 3,             #' Label size
  alpha = 0.5,                #' Transparency for spatial plots
  pt_size_factor = 2          #' Point size factor for spatial plots
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
  "Banksy",           # Spatial analysis
  "SeuratData",       # Data loading
  "dplyr",            # Data manipulation
  "pryr",             # Memory utilities
  "SeuratWrappers",   # Additional Seurat functions
  "ggplot2",          # Plotting
  "gridExtra"         # Multi-panel plots
)

invisible(lapply(required_packages, library, character.only = TRUE))

# Create main output directory
dir.create(OUTPUT_DIR, showWarnings = TRUE, recursive = TRUE)

################################################################################
#' Load Sample Information
#'
#' @description Reads the sample list file and extracts paths
#' @return Character vector of sample information lines
################################################################################

cat("\n=== Loading sample information ===\n")
matrix_lines <- readLines(MATRIX_LIST_FILE)
cat(sprintf("Found %d samples to process\n\n", length(matrix_lines)))

################################################################################
#' Process Single Sample
#'
#' @description Main processing function for individual spatial samples
#' @param sample_line String containing "sample_name,data_path"
#' @return Processed Seurat object
################################################################################

process_sample <- function(sample_line) {
  
  # Parse sample information
  split_info <- strsplit(sample_line, ",")[[1]]
  sample_name <- split_info[1]
  data_path <- split_info[2]
  
  cat(sprintf("=== Processing Sample: %s ===\n", sample_name))
  cat(sprintf("Data path: %s\n\n", data_path))
  
  # ============================================================================
  # Step 1: Load Spatial Data
  # ============================================================================
  cat("Step 1/8: Loading 10X Visium data...\n")
  seobj <- Load10X_Spatial(data.dir = data_path)
  cat(sprintf("  Loaded %d spots\n\n", ncol(seobj)))
  
  # ============================================================================
  # Step 2: Quality Control
  # ============================================================================
  cat("Step 2/8: Quality control filtering...\n")
  
  # Calculate mitochondrial percentage
  seobj[["percent.mt"]] <- PercentageFeatureSet(seobj, pattern = "^MT-")
  
  # Filter based on QC parameters
  seobj <- subset(
    seobj,
    subset = nFeature_Spatial > QC_PARAMS$min_features & 
             percent.mt < QC_PARAMS$max_mt_percent
  )
  cat(sprintf("  Spots after QC: %d\n\n", ncol(seobj)))
  
  # ============================================================================
  # Step 3: Normalization
  # ============================================================================
  cat("Step 3/8: SCTransform normalization...\n")
  seobj <- SCTransform(seobj, assay = "Spatial", verbose = FALSE)
  cat("  Normalization complete\n\n")
  
  # ============================================================================
  # Step 4: BANKSY Spatial Analysis
  # ============================================================================
  cat("Step 4/8: Running BANKSY spatial analysis...\n")
  seobj <- RunBanksy(
    seobj,
    lambda = BANKSY_PARAMS$lambda,
    verbose = TRUE,
    assay = BANKSY_PARAMS$assay,
    slot = BANKSY_PARAMS$slot,
    features = BANKSY_PARAMS$features,
    k_geom = BANKSY_PARAMS$k_geom
  )
  cat("  BANKSY complete\n\n")
  
  # ============================================================================
  # Step 5: PCA
  # ============================================================================
  cat("Step 5/8: Running PCA...\n")
  seobj <- RunPCA(
    seobj,
    assay = "BANKSY",
    reduction.name = "pca.banksy",
    features = rownames(seobj),
    npcs = REDUCTION_PARAMS$n_pcs,
    verbose = FALSE
  )
  cat(sprintf("  Computed %d PCs\n\n", REDUCTION_PARAMS$n_pcs))
  
  # ============================================================================
  # Step 6: Clustering
  # ============================================================================
  cat("Step 6/8: Clustering...\n")
  
  # Find neighbors
  seobj <- FindNeighbors(
    seobj,
    dims = REDUCTION_PARAMS$dims,
    reduction = "pca.banksy",
    verbose = FALSE
  )
  
  # Find clusters
  seobj <- FindClusters(
    seobj,
    cluster.name = CLUSTER_PARAMS$cluster_name,
    resolution = CLUSTER_PARAMS$resolution,
    verbose = FALSE
  )
  
  # Set active identity
  Idents(seobj) <- CLUSTER_PARAMS$cluster_name
  
  n_clusters <- length(unique(Idents(seobj)))
  cat(sprintf("  Identified %d clusters\n\n", n_clusters))
  
  # ============================================================================
  # Step 7: UMAP
  # ============================================================================
  cat("Step 7/8: Running UMAP...\n")
  seobj <- RunUMAP(
    seobj,
    dims = REDUCTION_PARAMS$dims,
    reduction = "pca.banksy",
    verbose = FALSE
  )
  cat("  UMAP complete\n\n")
  
  # ============================================================================
  # Step 8: Visualization
  # ============================================================================
  cat("Step 8/8: Generating visualizations...\n")
  
  # Create sample output directory
  sample_dir <- file.path(OUTPUT_DIR, sample_name)
  dir.create(sample_dir, showWarnings = TRUE, recursive = TRUE)
  
  # Generate plots
  p1 <- DimPlot(
    seobj,
    pt.size = PLOT_PARAMS$pt_size,
    label = TRUE,
    label.size = PLOT_PARAMS$label_size,
    repel = TRUE
  ) + ggtitle("UMAP Clusters")
  
  p2 <- SpatialDimPlot(
    seobj,
    group.by = CLUSTER_PARAMS$cluster_name,
    stroke = NA,
    label = TRUE,
    label.size = PLOT_PARAMS$label_size,
    repel = TRUE,
    alpha = PLOT_PARAMS$alpha,
    pt.size.factor = PLOT_PARAMS$pt_size_factor
  ) + ggtitle("Spatial Clusters")
  
  # Save plot
  png_file <- file.path(sample_dir, 
                        sprintf("%s.Spatial_dimPlot.png", sample_name))
  png(png_file, 
      height = PLOT_PARAMS$height, 
      width = PLOT_PARAMS$width)
  grid.arrange(p1, p2, ncol = 2)
  dev.off()
  
  cat(sprintf("  Saved plot: %s\n", basename(png_file)))
  
  # ============================================================================
  # Marker Gene Identification
  # ============================================================================
  cat("\nIdentifying marker genes...\n")
  
  # Find markers for all clusters
  markers <- FindAllMarkers(
    seobj,
    only.pos = MARKER_PARAMS$only_pos,
    assay = MARKER_PARAMS$assay,
    min.pct = MARKER_PARAMS$min_pct,
    logfc.threshold = MARKER_PARAMS$logfc_threshold,
    verbose = FALSE
  )
  
  cat(sprintf("  Found %d marker genes\n", nrow(markers)))
  
  # Save results
  rds_file <- file.path(sample_dir, sprintf("%s.rds", sample_name))
  saveRDS(seobj, file = rds_file)
  cat(sprintf("  Saved RDS: %s\n", basename(rds_file)))
  
  marker_file <- file.path(sample_dir, sprintf("%s.diff.txt", sample_name))
  write.table(
    data.frame(GENE_NAME = rownames(markers), markers),
    file = marker_file,
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
  )
  cat(sprintf("  Saved markers: %s\n\n", basename(marker_file)))
  
  cat(sprintf("=== Completed: %s ===\n\n", sample_name))
  
  return(seobj)
}

################################################################################
# Process All Samples
################################################################################

cat("\n" , rep("=", 80), "\n", sep = "")
cat("Starting analysis pipeline\n")
cat(rep("=", 80), "\n\n", sep = "")

# Process each sample
object_list <- lapply(matrix_lines, process_sample)

################################################################################
# Summary
################################################################################

cat("\n" , rep("=", 80), "\n", sep = "")
cat("Analysis Complete\n")
cat(rep("=", 80), "\n", sep = "")
cat(sprintf("Processed %d samples\n", length(object_list)))
cat(sprintf("Output directory: %s\n", OUTPUT_DIR))
cat("\nGenerated files per sample:\n")
cat("  - {sample_name}.rds (Seurat object)\n")
cat("  - {sample_name}.Spatial_dimPlot.png (Visualization)\n")
cat("  - {sample_name}.diff.txt (Marker genes)\n")
cat(rep("=", 80), "\n", sep = "")