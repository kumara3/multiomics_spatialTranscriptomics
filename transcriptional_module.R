# function 1: find the top 50 markers from each cluster filtered by log2FC > 0.25 and Padj < 0.05. Keep the markers shared between cluster
# output : getsets to be used for co-expression module in seurat
library(dplyr)
#library(tidyverse)
library(gplots)
library(tools)
library(seurat)
output_dir <- "top50markers"
dir.create(file.path(output_dir), showWarnings = TRUE)
filter_markers <- function(A){
    base_filename <- tools::file_path_sans_ext(basename(A))
    marker_df <- read.table(A, header=TRUE, sep="\t")
    markers <- marker_df %>% group_by(cluster) %>% filter(avg_log2FC > 0.25 & p_val_adj < 0.05 ) %>% top_n(n=50, wt=avg_log2FC)

    # remove clustesr with fewer than 50 genes
    cluster_counts <- markers %>% count(cluster)
    valid_clusters <- cluster_counts %>% filter(n >= 50) %>% pull(cluster)
    markers <- markers %>% filter(cluster %in% valid_clusters)

    common_genes <- markers %>% count(GENE_NAME) %>% filter(n > 1) %>% pull(GENE_NAME)
    gene_cluster <- markers %>% filter(GENE_NAME %in% common_genes) %>% select(GENE_NAME, cluster) %>% arrange(GENE_NAME, cluster)
    tmp = paste(output_dir,"/",base_filename,".Top50.markers.txt", sep="")
    tmp_common = paste(output_dir,"/",base_filename,".common.shared.between.clusters.markers.txt", sep="")
    write.table(markers, file=tmp, sep="\t", col.name=NA)
    write.table(common_genes, file=tmp_common, sep="\t", col.name=NA)
    return (unique(markers$GENE_NAME))
    
}

genematrixList="input_topmarkers_diff_gene.txt"
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
# 
# function 2
# Input : get the unique list of genes expressed across samples and clusters from above function. This will be the geneset
# output: module score
#https://github.com/satijalab/seurat/issues/5549

# Transcription module
dir.create(file.path("AddModuleScore"), showWarnings = TRUE)
obj = readRDS("/scratch/s163196/Gliome_spatial/github_scripts/multiomics_spatialTranscriptomics/banksy_harmony_integration_excluding_GBM5_1_2_samples/horizontal_integration_allSamples_banksy_harmony.rds")
obj = AddModuleScore(obj, features = gene_set_list, name=filenames)
rds_file = paste("AddModuleScore","/","addmodulescore_seuratobj.rds", sep="")
saveRDS(obj, file=rds_file)



    







    



}