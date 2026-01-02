# usage:
usage1 = "

Usage:
Rscript RCTD.R (followed by options below)
    1. <FULL PATH OF REFERNC GENEMATRIX>
    2. <FULL PATH OF METADATA FILE WITH CELL TYPE ANNOTATION>
    3. <FULL PATH OF OUTDIR>
    4. <PREFIX FOR OUTPUT FILES>
    5. <COLUMN NAME IN METADATA FILE WHICH HAS CELL TYPE ANNOTATION>
    6. <TRUE/FALSE - WHETHER TO REQUIRE INTEGER COUNTS IN REFERENCE DATA>
   "
args = commandArgs(trailingOnly = TRUE)
if(length(args)!=8) {
  cat(usage1)
  stop("\n Wrong parameters. See usage above.\n")
}

options(echo=TRUE)
system("ulimit -n 10000")
packages <- c("Seurat","spacexr")
invisible(lapply(packages, library, character.only=TRUE))
options(future.globals.maxSize = 4 * 1024^3)

matrix_ref = args[1]
metadata_ref = args[2]
sample_count = args[3]
sample_image = args[4]
output_dir = args[5]
prefix = args[6]
celltype_colname = args[7]  # column name in metadata file which has cell type annotation
require_int_val = args[8]  #(TRUE/FALSE)
dir.create(file.path(output_dir), showWarnings = TRUE)


counts <- read.table(matrix_ref, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
metadata <- read.table(metadata_ref, sep = "\t", header = TRUE, row.names = 1)

# Ensure rownames of metadata match column names of counts
all(colnames(counts) %in% rownames(metadata))  # should be TRUE
metadata <- metadata[colnames(counts), , drop = FALSE]  # re-order metadata to match columns

sanitize_labels <- function(x) {
  x <- as.character(x)
  x <- gsub("/", "-", x)                 # replace slashes
  x <- gsub("[[:space:]]+", "_", x)      # spaces -> underscore
  x <- gsub("[^A-Za-z0-9._-]", "_", x)   # keep only safe chars
  trimws(x)
}

# Create Seurat object
ref <- CreateSeuratObject(
  counts = counts,
  assay = "RNA",
  project = prefix,
  meta.data = metadata
)

Idents(ref) <- celltype_colname

# Remove Type with low cell counts min cell count > 25
min_cells <- 25                        # choose your cutoff
tab <- table(Idents(ref))
keep_ids <- names(tab)[tab >= min_cells]

# keep only sufficiently large identities
ref_filt <- subset(ref, idents = keep_ids)
Idents(ref_filt) <- droplevels(Idents(ref_filt))
# downsample each cell type to max 500 cells
ref_500 <- subset(ref_filt, downsample = 500)

# extract information to pass to the RCTD Reference function
counts <- ref_500[["RNA"]]$counts
cluster <- as.factor(ref_500@meta.data[[celltype_colname]])
cluster <- factor(sanitize_labels(cluster))
stopifnot(!any(grepl("/", levels(cluster))))

names(cluster) <- colnames(ref_500)
nUMI_count <- ref_500$nCount_RNA
names(nUMI_count) <- colnames(ref_500)
reference <- Reference(counts, cluster, min_UMI = 10,nUMI_count,require_int=require_int_val)


# set up query with the RCTD function SpatialRNA
# loop in all the samples and create one single spatialRNA object

sp_obj = readRDS("/scratch/Gliome_spatial/github_scripts/multiomics_spatialTranscriptomics/banksy_harmony_integration/horizontal_integration_allSamples_banksy_harmony.rds")
#sample_name = unique(sp_obj@meta.data$sample)
#names(sp_obj@images) <- sample_name
#image_slot_name <- seq_along(sample_name)
#names(sp_obj@images) - gives the image argument for tissue coordinate
#Layers(sp_obj[["Spatial"]]) - gives the count slot names

counts <- sp_obj[["Spatial"]][sample_count]
coords <- GetTissueCoordinates(sp_obj, image=sample_image)
coords <- coords[,c(1,2)]
colnames(coords) <- c("x", "y")
coords[is.na(colnames(coords))] <- NULL
query <- SpatialRNA(coords, counts, colSums(counts))

# run RCTD
RCTD <- create.RCTD(query, reference, max_cores = 32)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
tmp = paste(output_dir,"/","RCTD_",prefix,".rds",sep="")
saveRDS(RCTD, file = tmp)


# sp_obj <- AddMetaData(sp_obj, metadata = RCTD@results$results_df)
# sp_obj$first_type[is.na(sp_obj$first_type)] <- "Unknown"
# Idents(sp_obj) <- "first_type"
# cells <- CellsByIdentities(sp_obj)
##attributes(spatial_RNA_list[["DMG1"]])
#[1] "coords" "counts" "nUMI"
#spatial_RNA_list[["DMG1"]]@coords


#plot cells:
# plot_cell_types <- function(data, label,n) {
#   p <- ggplot(data, aes(x = get(label), y = n, fill = first_type)) +
#     geom_bar(stat = "identity", position = "stack") +
#     geom_text(aes(label = ifelse(n >= min_count_to_show_label, first_type, "")), position = position_stack(vjust = 0.5), size = 2) +
#     xlab(label) +
#     ylab("# of Spots") +
#     ggtitle(paste0("Distribution of Cell Types across ", label)) +
#     theme_minimal()
# }

## for specific cell types
#cell_type_banksy_counts <- query_obj[[]] %>%
#  dplyr::filter(first_type %in% excitatory_names) %>%
#  dplyr::count(full_first_type, banksy_cluster)

# min_count_to_show_label <- 20
# n = length(Cells(sp_obj))

# p <- plot_cell_types(query_obj@meta.data, "banksy_cluster",n)
# pdf("RCTD_Bhaduri_ref_500_celltype_banksy.pdf", height=10, width=10)
# print(p)
# dev.off()

# References for scRNA-seq data:
# Bhaduri : https://www.ebi.ac.uk/ena/browser/view/PRJNA579593; https://cells.ucsc.edu/?ds=gbm ; counts data available; cell ranger noramalization
# Nowakowsi : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102130; https://cells.ucsc.edu/?ds=cortex-dev ; No raw counts data available; https://cells.ucsc.edu/?ds=cortex-dev cpm counts
# Filbin : https://singlecell.broadinstitute.org/single_cell/study/SCP147/single-cell-analysis-in-pediatric-midline-gliomas-with-histone-h3k27m-mutation#/ ; No raw counts data available TPM counts
# Aldinger : https://cells.ucsc.edu/?ds=cbl-dev ; raw counts data available; cell ranger normalization
# https://github.com/dmcable/spacexr/issues/49
#https://github.com/dmcable/spacexr/issues/80

# matrix_ref = "/scratch//Gliome_spatial/github_scripts/multiomics_spatialTranscriptomics/RCTD/Bhaduri/exprMatrix.tsv"
# metadata_ref = "/scratch//Gliome_spatial/github_scripts/multiomics_spatialTranscriptomics/RCTD/Bhaduri/meta.tsv"
# sample_count = "counts.1"
# sample_image = "slice1"
# output_dir = "Bhaduri_RCTD"
# prefix = "DMG1"
# celltype_colname = "Cell.Type.Assignment" # column name in metadata file which has cell type annotation
# require_int_val = FALSE

# Create a new column in metadata to combine sample and cluster info
sp_obj$module <- paste(sp_obj$sample,sp_obj$BANKSY_snn_res.0.8,sep="_")

library(purrr)
library(dplyr)

rds_dir <- "/scratch/Gliome_spatial/github_scripts/multiomics_spatialTranscriptomics/Bhaduri_RCTD"
files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)

stopifnot(length(files) == 11)  # optional safety
df_list <- lapply(files, FUN=function(x){
  obj <- readRDS(x)
  df <- obj@results$results_df
  })

concatenated_dfs <- bind_rows(df_list, .id = "source")
sp_obj_concatenated <- AddMetaData(sp_obj, metadata = concatenated_dfs)

filename = "/scratch/Gliome_spatial/github_scripts/multiomics_spatialTranscriptomics/AddModuleScore_TranscriptionalModule/with_mapped_group.tsv"
df <- read.table('filename', header = TRUE,  sep = '\t')
df_distinct <- df %>% distinct(across(c(1,3,4)))
sp_meta_joined_df = sp_obj_concatenated@meta.data %>% rownames_to_column("cell") %>% left_join(df_distinct,by='module') %>% column_to_rownames("cell")
sp_obj_concatenated <- AddMetaData(sp_obj_concatenated, metadata = sp_meta_joined_df)
#write.table(all_df, file="all.concatenated.df.txt",sep="\t", col.names=NA)


cells_keep <- rownames(sp_obj_concatenated@meta.data)[!is.na(sp_obj_concatenated@meta.data$Cluster) & sp_obj_concatenated@meta.data$Cluster != "NA"]
sp_obj_concatenated_subset <- subset(sp_obj_concatenated, cells=cells_keep)

library(dplyr)
library(ggplot2)
library(scales)

# color map from earlier (edit hex codes as needed)
cols <- c(
  "Radial_Glia" = "#E86F2D",
  "Protoplasmic_Astrocyte" = "#66C2A4",
  "OPC" = "#E34A33",
  "Oligodendrocyte" = "#386CB0",
  "Neuron" = "#A6D854",
  "Immature_Astrocyte" = "#7FBF7B",
  "Dividing_OPC" = "#1B9E77",
  "Microglia" = "#B2B2B2",
  "Pericyte" = "#8DA0CB",
  "B_Cells" = "#000000"
)
cols_extra <- c(
  "Mature_IPC-Newborn_Neuron"   = "#66A61E",  # neuronal progenitor (darker green than Neuron)
  "Mixed_Progenitor-Neuron"     = "#FFD92F",  # progenitor → warm yellow
  "Glycolytic_Progenitor"       = "#FEC44F",  # progenitor (metabolic) → amber
  "Dividing_Progenitor"         = "#F0E442",  # dividing progenitor → bright yellow-green
  "Dividing_Neuron"             = "#4DAF4A",  # neuron but dividing → darker green
  "Dividing_B_Cells"            = "#4D4D4D",  # near B_Cells (black) but distinguishable
  "Unknown"                     = "#DDDDDD",  # neutral gray
  "CGE_iN"                      = "#7570B3",  # interneuron-like → purple
  "Endothelial"                 = "#80B1D3",  # close to Pericyte but distinct blue
  "Tumor_Associated_Macrophage" = "#8C510A",  # myeloid but distinct (brown)
  "Red_blood_cells"             = "#B2182B"   # RBC → deep red
)

cols <- c(cols, cols_extra)
meta <- sp_obj_concatenated_subset@meta.data %>%
  as.data.frame() %>%
  select(Cluster, first_type) %>%
  filter(!is.na(Cluster)) %>%                    # keep valid clusters
  mutate(first_type = trimws(first_type)) %>%    # tidy labels
  filter(!is.na(first_type), first_type != "Unknown") %>%  # <- exclude "Unknown"
  mutate(
    Cluster   = factor(Cluster),
    first_type = factor(first_type, levels = names(cols))
  )
pdf("Bhaduri_RCTD_celltype_by_clusterv2.pdf", height=6, width=8)
ggplot(meta, aes(x = Cluster, fill = first_type)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = cols, breaks = names(cols), name = "first_type") +
  labs(x = "Cluster", y = "% of cells", title = "Cell type composition by cluster") +
  theme_minimal()
dev.off()



df_Filbin <- read.table("/scratch/Gliome_spatial/github_scripts/multiomics_spatialTranscriptomics/RCTD/Filbin/PortalK27M_Metadata.vh20180223.txt", header = TRUE, sep = "\t",
                 quote = "", comment.char = "", check.names = FALSE,
                 stringsAsFactors = FALSE)

# Drop the metadata row that starts with "TYPE"
df_Filbin <- df_Filbin[df_Filbin$NAME != "TYPE", , drop = FALSE]

# Identify the last 5 columns by position
last5_idx   <- (ncol(df_Filbin) - 4):ncol(df_Filbin)
last5_names <- names(df_Filbin)[last5_idx]

# Coerce the last 5 columns to numeric (handles "NA" strings cleanly)
for (j in last5_idx) df_Filbin[[j]] <- suppressWarnings(as.numeric(df_Filbin[[j]]))

# For each row, grab the column name of the maximum among the last 5
df_Filbin$MaxCol <- apply(df_Filbin[, last5_idx, drop = FALSE], 1, function(v) {
  if (all(is.na(v))) NA_character_ else last5_names[which.max(v)]
})

df_Filbin$celltype <- ifelse(is.na(df_Filbin$MaxCol), df_Filbin$Type, paste(df_Filbin$Type, df_Filbin$MaxCol, sep = "_"))

tmp = "/scratch/Gliome_spatial/github_scripts/multiomics_spatialTranscriptomics/RCTD/Filbin/PortalK27M_Metadata.vh20180223v2.txt"
write.table(df_Filbin, file=tmp, sep = "\t", row.names = FALSE, quote = FALSE)
