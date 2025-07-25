# https://prabhakarlab.github.io/Banksy/articles/batch-correction.html
# Spatial data integration with Harmony
# https://github.com/prabhakarlab/Banksy/issues/27

# usage:
usage1 = "

Usage:
Rscript seuratSpat.R (followed by options below)
    1. <FULL PATH OF A TEXT FILE WHICH HAS GENEMATRIX FOR ALL THE SAMPLES WITH PREFIX AS SAMPLE NAME>
    2. <FULL PATH OF OUTDIR>
   "
args = commandArgs(trailingOnly = TRUE)
if(length(args)!=2) {
  cat(usage1)
  stop("\n Wrong parameters. See usage above.\n")
}

options(echo=TRUE)
system("ulimit -n 10000")
packages <- c("Seurat","Banksy","SeuratData", "dplyr","pryr","SeuratWrappers","ggplot2","gridExtra","harmony","cowplot")
invisible(lapply(packages, library, character.only=TRUE))
options(future.globals.maxSize = 4 * 1024^3)

matrixList = args[1]
output_dir = args[2]
dir.create(file.path(output_dir), showWarnings = TRUE)

matrix_vector = character(0)
read_input_gene_matrix <- file(matrixList,'r')
linn <- readLines(read_input_gene_matrix)

## gene matrix to character vector ####
for(i in 1:length(linn))
	{
        matrix_vector <- c(matrix_vector,linn[i])
	}

object_list <- lapply(matrix_vector, FUN=function(i) {
    split_list <- strsplit(i,",")
    sample_name <- split_list[[1]][1]
    filtered_matrix_infile <- split_list[[1]][2]
    print("Processing:")
    print(sample_name)
    print(filtered_matrix_infile)
    seobj <- Load10X_Spatial(data.dir=filtered_matrix_infile)
    seobj[['percent.mt']] <- PercentageFeatureSet(seobj, pattern = '^MT-')
    seobj[['sample']] <- sample_name
    seobj <- subset(seobj, subset = nFeature_Spatial > 200 & percent.mt < 25)
    seobj <- SCTransform(seobj, assay = "Spatial")
    
    #return(list(sample_name = sample_name, seobj = seobj))
})

# Create a named list using the sample names
# object_list <- setNames(
#     lapply(object_list, function(x) x$seobj),
#     sapply(object_list, function(x) x$sample_name)
# )

seob = merge(x=object_list[[1]], y=object_list[2:length(object_list)]) # merge objects
ranked_features = SelectIntegrationFeatures(object_list, nfeatures = 2000)
VariableFeatures(seob) = ranked_features
seob = RunBanksy(seob, lambda=0.2,assay='SCT',slot='data',k_geom = 15, split.scale = FALSE)
seob = RunPCA(seob, assay = 'BANKSY', features = rownames(seob), npcs=30)
seob = RunHarmony(seob, group.by.vars='sample')
seob = RunUMAP(seob, reduction = 'harmony', reduction.name = 'umap_harmony', dims=1:30)
seob = FindNeighbors(seob, reduction = 'harmony')
seob = FindClusters(seob)

color_cluster=c("#00ADFA","#005e61","#BF40BF","#655b6d","#FF0000","#E88526","#D39200","#B79F00","#93AA00","#5EB300","6c8f9d","#00BA38","#762a83","#4292c6","#807dba","#08519c","#c6dbef","#deebf7","#f8766d","#7CAE00","#00BFC4","#C77CFF")
color_samples=c("#00ADFA","#005e61","#BF40BF","#655b6d","#FF0000","#E88526","#D39200","#B79F00","#93AA00","#5EB300")
#"#C77CFF","#FF69B4"
file = paste(output_dir,"/","horizontal_integration_allSamples.png", sep="")
rds_file = paste(output_dir,"/","horizontal_integration_allSamples_banksy_harmony.rds", sep="")
saveRDS(seob, rds_file)
png(file, height= 980, width=1080)
grid.arrange(
    DimPlot(seob, pt.size = 0.25, label = TRUE, label.size = 3, repel = TRUE,group.by = c('sample'), cols=color_samples),
    DimPlot(seob, pt.size = 0.25, label = TRUE, label.size = 3, repel = TRUE,group.by = c('BANKSY_snn_res.0.8'), cols=color_cluster),
    SpatialDimPlot(seob, stroke = NA, label = TRUE, label.size = 3,repel = TRUE, alpha = 0.5, pt.size.factor = 2),
    nrow = 3
)
dev.off()
