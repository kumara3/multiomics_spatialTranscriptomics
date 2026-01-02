#!/usr/bin/env Rscript

################################################################################
# Spatial Transcriptomics Integration with BANKSY and Harmony
################################################################################
# Description: This script integrates multiple spatial transcriptomics samples
# using BANKSY for spatial feature extraction and Harmony for batch correction.
# It performs clustering, visualization, and marker gene identification.
#
# References:
# - BANKSY: https://prabhakarlab.github.io/Banksy/articles/batch-correction.html
# - GitHub: https://github.com/prabhakarlab/Banksy/issues/27
################################################################################

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Define usage message
usage <- "
Usage:
  Rscript seuratSpat.R <matrix_list_file> <output_dir>

Arguments:
  matrix_list_file  Full path to text file containing sample information
                    Format: sample_name,path_to_gene_matrix (one per line)
  output_dir        Full path to output directory

Example:
  Rscript seuratSpat.R samples.txt /path/to/output
"

# Validate arguments
if (length(args) != 2) {
  cat(usage)
  stop("\nError: Incorrect number of arguments provided.\n", call. = FALSE)
}

################################################################################
# Configuration Parameters
################################################################################

# Input/Output paths
MATRIX_LIST_FILE <- args[1]
OUTPUT_DIR <- args[2]

# Quality control thresholds
MIN_FEATURES <- 200           # Minimum number of features per spot
MAX_MT_PERCENT <- 25          # Maximum mitochondrial percentage

# Integration parameters
N_INTEGRATION_FEATURES <- 2000  # Number of features for integration
BANKSY_LAMBDA <- 0.2           # BANKSY spatial weight parameter
BANKSY_K_GEOM <- 15            # Number of spatial neighbors
N_PCS <- 30                    # Number of principal components

# Clustering parameters
CLUSTER_RESOLUTION <- 0.8      # Resolution for FindClusters

# Visualization parameters
PLOT_HEIGHT <- 980
PLOT_WIDTH <- 1080
POINT_SIZE <- 0.25
LABEL_SIZE <- 3

# Marker gene identification parameters
MIN_PCT <- 0.25               # Minimum percentage of cells expressing marker
LOGFC_THRESHOLD <- 0.25       # Log fold-change threshold for markers

# Color palettes
COLOR_CLUSTER <- c("#00ADFA", "#005e61", "#BF40BF", "#655b6d", "#FF0000", 
                   "#E88526", "#D39200", "#B79F00", "#93AA00", "#5EB300", 
                   "#6c8f9d", "#00BA38", "#762a83", "#4292c6", "#807dba", 
                   "#08519c", "#c6dbef", "#deebf7", "#f8766d", "#7CAE00", 
                   "#00BFC4", "#C77CFF")

COLOR_SAMPLES <- c("#00ADFA", "#005e61", "#BF40BF", "#655b6d", "#FF0000", 
                   "#E88526", "#D39200", "#B79F00", "#93AA00", "#5EB300", 
                   "#6c8f9d")

################################################################################
# Setup Environment
################################################################################

# Set system limits
system("ulimit -n 10000")

# Set options
options(echo = TRUE)
options(future.globals.maxSize = 4 * 1024^3)  # 4GB max object size

# Load required packages
required_packages <- c("Seurat", "Banksy", "SeuratData", "dplyr", "pryr",
                       "SeuratWrappers", "ggplot2", "gridExtra", "harmony", 
                       "cowplot")

invisible(lapply(required_packages, library, character.only = TRUE))

# Create output directory
dir.create(OUTPUT_DIR, showWarnings = TRUE, recursive = TRUE)

################################################################################
# Load Sample Information
################################################################################

cat("\n=== Loading sample information ===\n")

# Read sample matrix file paths
matrix_lines <- readLines(MATRIX_LIST_FILE)
cat(sprintf("Found %d samples to process\n", length(matrix_lines)))

################################################################################
# Process Individual Samples
################################################################################

cat("\n=== Processing individual samples ===\n")

object_list <- lapply(matrix_lines, function(line) {
  # Parse sample information
  split_info <- strsplit(line, ",")[[1]]
  sample_name <- split_info[1]
  matrix_path <- split_info[2]
  
  cat(sprintf("\nProcessing sample: %s\n", sample_name))
  cat(sprintf("  Matrix path: %s\n", matrix_path))
  
  # Load 10X Visium data
  seobj <- Load10X_Spatial(data.dir = matrix_path)
  
  # Calculate mitochondrial percentage
  seobj[["percent.mt"]] <- PercentageFeatureSet(seobj, pattern = "^MT-")
  
  # Add sample metadata
  seobj[["sample"]] <- sample_name
  
  # Quality control filtering
  seobj <- subset(seobj, 
                  subset = nFeature_Spatial > MIN_FEATURES & 
                           percent.mt < MAX_MT_PERCENT)
  
  cat(sprintf("  Spots after QC: %d\n", ncol(seobj)))
  
  # Normalize and scale using SCTransform
  seobj <- SCTransform(seobj, assay = "Spatial", verbose = FALSE)
  
  return(seobj)
})

################################################################################
# Merge and Integrate Samples
################################################################################

cat("\n=== Merging samples ===\n")

# Merge all Seurat objects
seob <- merge(x = object_list[[1]], 
              y = object_list[2:length(object_list)])

cat(sprintf("Total spots after merging: %d\n", ncol(seob)))

# Select integration features
cat("\n=== Selecting integration features ===\n")
ranked_features <- SelectIntegrationFeatures(object_list, 
                                             nfeatures = N_INTEGRATION_FEATURES)
VariableFeatures(seob) <- ranked_features

################################################################################
# BANKSY Spatial Analysis
################################################################################

cat("\n=== Running BANKSY spatial analysis ===\n")

seob <- RunBanksy(seob, 
                  lambda = BANKSY_LAMBDA,
                  assay = "SCT",
                  slot = "data",
                  k_geom = BANKSY_K_GEOM,
                  split.scale = FALSE)

################################################################################
# Dimensionality Reduction and Batch Correction
################################################################################

cat("\n=== Running PCA ===\n")
seob <- RunPCA(seob, 
               assay = "BANKSY",
               features = rownames(seob),
               npcs = N_PCS,
               verbose = FALSE)

cat("\n=== Running Harmony batch correction ===\n")
seob <- RunHarmony(seob, 
                   group.by.vars = "sample",
                   verbose = FALSE)

cat("\n=== Running UMAP ===\n")
seob <- RunUMAP(seob, 
                reduction = "harmony",
                reduction.name = "umap_harmony",
                dims = 1:N_PCS,
                verbose = FALSE)

################################################################################
# Clustering
################################################################################

cat("\n=== Clustering ===\n")
seob <- FindNeighbors(seob, reduction = "harmony", verbose = FALSE)
seob <- FindClusters(seob, resolution = CLUSTER_RESOLUTION, verbose = FALSE)

cat(sprintf("Identified %d clusters\n", 
            length(unique(seob@meta.data[[paste0("BANKSY_snn_res.", 
                                                  CLUSTER_RESOLUTION)]]))))

################################################################################
# Save Integrated Object
################################################################################

cat("\n=== Saving results ===\n")

rds_file <- file.path(OUTPUT_DIR, "horizontal_integration_allSamples_banksy_harmony.rds")
saveRDS(seob, rds_file)
cat(sprintf("Saved RDS file: %s\n", rds_file))

################################################################################
# Visualization
################################################################################

cat("\n=== Generating visualizations ===\n")

# Main integration plot
plot_file <- file.path(OUTPUT_DIR, "horizontal_integration_allSamples.png")
png(plot_file, height = PLOT_HEIGHT, width = PLOT_WIDTH)

grid.arrange(
  DimPlot(seob, 
          pt.size = POINT_SIZE,
          label = TRUE,
          label.size = LABEL_SIZE,
          repel = TRUE,
          group.by = "sample",
          cols = COLOR_SAMPLES) + 
    ggtitle("UMAP by Sample"),
  
  DimPlot(seob,
          pt.size = POINT_SIZE,
          label = TRUE,
          label.size = LABEL_SIZE,
          repel = TRUE,
          group.by = paste0("BANKSY_snn_res.", CLUSTER_RESOLUTION),
          cols = COLOR_CLUSTER) +
    ggtitle("UMAP by Cluster"),
  
  SpatialDimPlot(seob,
                 stroke = NA,
                 label = TRUE,
                 label.size = LABEL_SIZE,
                 repel = TRUE,
                 alpha = 0.5,
                 pt.size.factor = 2) +
    ggtitle("Spatial Clusters"),
  
  nrow = 3
)

dev.off()
cat(sprintf("Saved plot: %s\n", plot_file))

# Split by sample plot
split_plot_file <- file.path(OUTPUT_DIR, "horizontal_integration_splitBySample.png")
png(split_plot_file, height = PLOT_HEIGHT * 2, width = PLOT_WIDTH)

grid.arrange(
  DimPlot(seob,
          pt.size = POINT_SIZE,
          label = TRUE,
          label.size = LABEL_SIZE,
          repel = TRUE,
          split.by = "sample",
          cols = COLOR_CLUSTER),
  nrow = 4
)

dev.off()
cat(sprintf("Saved split plot: %s\n", split_plot_file))

################################################################################
# Identify Marker Genes
################################################################################

cat("\n=== Identifying marker genes ===\n")

# Prepare SCT assay for differential expression
seob <- PrepSCTFindMarkers(seob, assay = "SCT", verbose = FALSE)

# Find markers for all clusters
markers <- FindAllMarkers(seob,
                         only.pos = TRUE,
                         assay = "SCT",
                         min.pct = MIN_PCT,
                         logfc.threshold = LOGFC_THRESHOLD,
                         verbose = FALSE)

# Save marker genes
marker_file <- file.path(OUTPUT_DIR, "Integrated.FindMarkers.diff.txt")
write.table(data.frame(GENE_NAME = rownames(markers), markers),
            file = marker_file,
            row.names = FALSE,
            sep = "\t",
            quote = FALSE)

cat(sprintf("Saved marker genes: %s\n", marker_file))
cat(sprintf("Total markers identified: %d\n", nrow(markers)))

################################################################################
# Completion
################################################################################

cat("\n=== Analysis complete ===\n")
cat(sprintf("Output directory: %s\n", OUTPUT_DIR))
cat("\nGenerated files:\n")
cat(sprintf("  - %s\n", basename(rds_file)))
cat(sprintf("  - %s\n", basename(plot_file)))
cat(sprintf("  - %s\n", basename(split_plot_file)))
cat(sprintf("  - %s\n", basename(marker_file)))