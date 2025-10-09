# function 1: find the top 50 markers from each cluster filtered by log2FC > 0.25 and Padj < 0.05. Keep the markers shared between cluster
# output : getsets to be used for co-expression module in seurat
library(dplyr)
<<<<<<< HEAD
#library(tidyverse)
=======
>>>>>>> commit and sync
library(tibble)
library(gplots)
library(tools)
library(Seurat)
<<<<<<< HEAD
output_dir <- "top50markers"
=======
library(corrplot) 

output_dir <- "AddModuleScore_TranscriptionalModule_09082025"
>>>>>>> commit and sync
dir.create(file.path(output_dir), showWarnings = TRUE)
filter_markers <- function(A){
    base_filename <- tools::file_path_sans_ext(basename(A))
    marker_df <- read.table(A, header=TRUE, sep="\t")
    markers <- marker_df %>% group_by(cluster) %>% filter(avg_log2FC > 0.25 & p_val_adj < 0.05 ) %>% top_n(n=100, wt=avg_log2FC)

    # remove clustesr with fewer than 50 genes
    cluster_counts <- markers %>% count(cluster)
    valid_clusters <- cluster_counts %>% filter(n >= 50) %>% pull(cluster)
    markers <- markers %>% filter(cluster %in% valid_clusters)

     # Create list of gene vectors per cluster
    marker_list <- markers %>%
        group_by(cluster) %>%
        summarise(genes = list(unique(GENE_NAME))) %>%
        deframe()

    common_genes <- markers %>% count(GENE_NAME) %>% filter(n > 1) %>% pull(GENE_NAME)
    gene_cluster <- markers %>% filter(GENE_NAME %in% common_genes) %>% select(GENE_NAME, cluster) %>% arrange(GENE_NAME, cluster)
    tmp = paste(output_dir,"/",base_filename,".Top50.markers.txt", sep="")
    tmp_common = paste(output_dir,"/",base_filename,".common.shared.between.clusters.markers.txt", sep="")
    write.table(markers, file=tmp, sep="\t", col.name=NA)
    write.table(common_genes, file=tmp_common, sep="\t", col.name=NA)
    return (marker_list)
    
}

genematrixList="/scratch/Gliome_spatial/github_scripts/multiomics_spatialTranscriptomics/input_topmarkers_diff_gene.txt"
diffgene_matrix_vector = character(0)
gene_matrix <- file(genematrixList,'r')
linn <- readLines(gene_matrix)

## gene matrix to character vector ####
for(i in 1:length(linn))
	{
        diffgene_matrix_vector <- c(diffgene_matrix_vector,linn[i])
	}
class(diffgene_matrix_vector)
gene_set_list = lapply(diffgene_matrix_vector, filter_markers)
filenames <- sub("\\.diff\\.txt$", "", basename(diffgene_matrix_vector))
names(gene_set_list) <- filenames

feature <- unlist(
  lapply(names(gene_set_list), function(set_name) {
    lapply(names(gene_set_list[[set_name]]), function(cluster_id) {
      # Extract the gene vector
      genes <- gene_set_list[[set_name]][[cluster_id]]
<<<<<<< HEAD
      # Create a named list entry like "GBM5_1_8"
=======
>>>>>>> commit and sync
      setNames(list(genes), paste0(set_name, "_", cluster_id))
    })
  }),
  recursive = FALSE
)
feature_name = as.character(lapply(feature, function(x) names(x)))
unflattened_list <- unlist(feature, recursive = FALSE, use.names = TRUE)
names(unflattened_list)

<<<<<<< HEAD
# Total gene set
#sum(sapply(gene_set_list, function(x) length(x))) # 91 


# function 2
=======
obj = readRDS("/scratch/Gliome_spatial/github_scripts/multiomics_spatialTranscriptomics/banksy_harmony_integration/horizontal_integration_allSamples_banksy_harmony.rds")
# modify the gene list to cater to add module. which assay to take either 'SCT' or 'BANKSY'

DefaultAssay(obj) <- 'SCT'
# downsamples the seurat object to include subset of cells
Idents(obj) <- "sample"
#obj <- subset(obj, downsample = 1000)
available_genes <- rownames(obj[['SCT']])

# Filter feature lists
features_filtered <- lapply(unflattened_list, function(gene_set) {
  intersect(gene_set, available_genes)
})
features_filtered <- features_filtered[sapply(features_filtered, length) > 29]


rds_file = paste(output_dir,"/","addmodulescore_seuratobj.rds", sep="")
obj = AddModuleScore(obj, features = features_filtered, name=paste0(names(features_filtered),"_score"), assay = 'SCT')
saveRDS(obj, file=rds_file) # save R object

# function to clean the column names added by addmodule score seurat
clean_module_names <- function(names_vector) {
  # Remove _scoreXX where XX is one or more digits
  cleaned <- gsub("_score\\d+$", "", names_vector)
  return(cleaned)
}

# Get all metadata column names that match the score pattern
score_cols <- names(obj@meta.data)[grepl("_score\\d+$", names(obj@meta.data))]
cleaned_names <- clean_module_names(score_cols)

# Rename columns
names(obj@meta.data)[names(obj@meta.data) %in% score_cols] <- cleaned_names

#convert feature list to data frame
features_df <- do.call(rbind, lapply(names(features_filtered), function(set_name) {
  data.frame(
    GeneSet = set_name,
    Gene = features_filtered[[set_name]],
    stringsAsFactors = FALSE
  )
}))

write.table(features_df, file = paste(output_dir,"/","module_features.txt", sep = ""), row.names = FALSE, sep="\t")

mod_name = names(obj@meta.data)[10:66]
mod_name_df  = obj@meta.data[,mod_name]
#obj[["mod_score"]] <- CreateAssayObject(data = t(x = FetchData(object = obj, vars = mod_name)))
#mod_name = gsub("_", "-", names(obj@meta.data)[10:62])

# complex heatmap
MyMatrix_gm = as.matrix(t(mod_name_df)) ## transpose the matrix cells are columns and modules are genes
row_scale <- t(scale(t(MyMatrix_gm))) # row gene modules

# include the hierarchical clustering of modules on correlation of modules
hr <- hclust(as.dist(1-cor(t(MyMatrix_gm), method = "pearson")), method="complete") 

# cuttree
mycl <- cutree(hr, k=4)
clusterCols <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]


# gene module to gene and cluster mapping:
gene_cluster_mapping <- features_df %>%
  left_join(
    data.frame(
      GeneSet = names(mycl),
      Cluster = mycl,
      stringsAsFactors = FALSE
    ),
    by = "GeneSet"
  )
head(gene_cluster_mapping)
tmp <- paste(output_dir,"/","gene_cluster_mapping.txt", sep = "")
write.table(gene_cluster_mapping, file = tmp, row.names = FALSE, sep="\t", quote=FALSE)


# or simply add the cluster ID to your data
MyMatrixC <- cbind(MyMatrix_gm, clusterID=mycl,colr_name=myClusterSideBar)
MyMatrixC <- apply(MyMatrixC[hr$order,], 2,rev)
write.table(MyMatrixC, file = paste(output_dir,"/","transcriptional_module",".hc_clust.txt", sep = ""), sep = "\t",col.names = NA)

# plot heatmap using ComplexHeatmap
library(circlize)
library(ComplexHeatmap)

row_anno <- rowAnnotation(
  Cluster = as.factor(mycl),
  col = list(Cluster = setNames(clusterCols, unique(mycl))),
  width = unit(0.5, "cm")
)
#plot heatmap using ComplexHeatmap

heatmap_save <- paste(output_dir,"/","transcriptional_module",".heatmap.pdf", sep = "")
heatmap_save <- paste("transcriptional_moduleT",".heatmap.pdf", sep = "")
pdf(heatmap_save, height = 8, width = 8)
ht <- Heatmap(row_scale, col=colorRamp2(c(-2,-1,0,1,2), c("blue", "lightblue3", "floralwhite", "gold", "gold4")),
        show_column_names = FALSE, 
        width = unit(8, "cm"),
        show_row_dend = FALSE,
        heatmap_legend_param = list(at=c(-2,-1,0,1,2),color_bar="continuous",title = "Scaled Expression"),
        cluster_columns = FALSE,
        cluster_rows = hr,
        left_annotation = row_anno,
        row_names_gp = gpar(fontsize = 3), 
        show_column_dend = FALSE)
draw(ht)
dev.off()

## Plot correlation matrix of module scores
# calculate module score for correlation matrix
mod_score <- obj@meta.data[,grep("GBM|DMG",colnames(obj@meta.data))]
colnames(mod_score) <- names(features_filtered)
cor_matrix <- cor(mod_score, method = "pearson")
#cor_matrix[is.na(cor_matrix)] <- 0

corr_file = paste(output_dir,"/","addmodulescore_corrPlot_color.png", sep="")
png(corr_file, height=980, width=980)
corrplot(cor_matrix, method = "color", order = "hclust", hclust.method = "complete",addrect = 4, tl.cex = 0.6)
dev.off()

corr_file = paste(output_dir,"/","addmodulescore_corrPlot_square.png", sep="")
png(corr_file, height=980, width=980)
corrplot(cor_matrix, method = "square", order = "hclust", hclust.method = "complete",ddCoef.col = 'black', tl.pos = 'd', cl.pos = 'n', col = COL2('BrBG'))
dev.off()


#dmg_gene_lists <- lapply(dmg_vars, get)
#unflattened_list[dmg1_genes]
# Total gene set
#sum(sapply(features_filtered, length))
#sum(sapply(gene_set_list, function(x) length(x))) # 91 

# function 2: Transcription module
>>>>>>> commit and sync
# Input : get the unique list of genes expressed across samples and clusters from above function. This will be the geneset
# output: module score
#https://github.com/satijalab/seurat/issues/5549

<<<<<<< HEAD
# Transcription module
dir.create(file.path("AddModuleScore"), showWarnings = TRUE)
obj = readRDS("/scratch/s163196/Gliome_spatial/github_scripts/multiomics_spatialTranscriptomics/banksy_harmony_integration_excluding_GBM5_2/horizontal_integration_allSamples_banksy_harmony.rds")
DefaultAssay(obj) <- 'SCT'
obj = AddModuleScore(obj, features = unflattened_list, name=names(unflattened_list))
#rds_file = paste("AddModuleScore","/","addmodulescore_seuratobj.rds", sep="")
#saveRDS(obj, file=rds_file)
#gene list input to add module score

#unflattened_list$DMG1_5
=======


# Find the percentage of samples in each module
# install.packages(c("readr","dplyr"))
setwd("/scratch/Gliome_spatial/github_scripts/multiomics_spatialTranscriptomics/AddModuleScore_TranscriptionalModule")
library(readr)
library(dplyr)
>>>>>>> commit and sync

# Read your table (tab- or space-delimited with header)
# Replace "input.tsv" with your path
df <- read_delim("/scratch/Gliome_spatial/github_scripts/multiomics_spatialTranscriptomics/AddModuleScore_TranscriptionalModule/gene_cluster_mapping.txt", delim = "\t", col_types = cols())  # use delim=" " if space-separated

# Helper tokens:
first_one <- sub("_.*", "", df$GeneSet)
first_two <- sub("^([^_]+_[^_]+).*", "\\1", df$GeneSet)

# Map to new group per your rules
df_mapped <- df %>%
  mutate(
    MappedGroup = case_when(
      first_one %in% c("GBM1", "GBM4")                 ~ "GBM_IDH_MUT",
      first_one %in% c("GBM2", "GBM3")                 ~ "GBM_IDH_WT",
      first_two == "GBM5_1"                            ~ "GBM_IDH_WT",
      first_two == "GBM5_2"                            ~ "GBM5_2",
      grepl("^DMG[1-5]$", first_one)                   ~ "DMG",
    )
  )
write.table(df_mapped, file = "/scratch/Gliome_spatial/github_scripts/multiomics_spatialTranscriptomics/AddModuleScore_TranscriptionalModule/with_mapped_group.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Write out with the new column

library(readr)
library(dplyr)
library(ggplot2)

# Read the input file (tab-delimited, no header assumed)
df <- read.table("/scratch/Gliome_spatial/github_scripts/multiomics_spatialTranscriptomics/AddModuleScore_TranscriptionalModule/with_mapped_group.tsv", header=TRUE, sep="\t")

# Assign column names
colnames(df) <- c("Sample", "Gene", "Module", "Group")

# Summarize counts by Value (3rd col) and Group (4th col)
plot_data <- df %>%
  count(Module, Group) %>%
  group_by(Module) %>%
  mutate(Percent = 100 * n / sum(n))

# Plot stacked bar chart
pdf("percentage_of_samples_per_module.pdf", height=6, width=8)
ggplot(plot_data, aes(x = factor(Module), y = Percent, fill = Group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = c(
      "DMG" = "red",
      "GBM_IDH_MUT" = "lightgreen",
      "GBM_IDH_WT"  = "darkgreen",
      "GBM5_2" = "blue"
    )
  )+
  labs(x = "Module",
       y = "Percentage of GBM samples",
       fill = "Group") +
  theme_minimal()
dev.off()
