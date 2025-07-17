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
packages <- c("Seurat","Banksy","SeuratData", "dplyr","pryr","SeuratWrappers","ggplot2","gridExtra")
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
#spacranger.sample.list <- Load10X_Spatial("")
object_list <- lapply(matrix_vector, FUN=function(i) {
    split_list <- strsplit(i,",")
    sample_name <- split_list[[1]][1]
    filtered_matrix_infile <- split_list[[1]][2]
    print("Processing:")
    print(sample_name)
    print(filtered_matrix_infile)
    seobj <- Load10X_Spatial(data.dir=filtered_matrix_infile)
    seobj[['percent.mt']] <- PercentageFeatureSet(seobj, pattern = '^MT-')
    seobj <- subset(seobj, subset = nFeature_Spatial > 200 & percent.mt < 25)
    seobj <- SCTransform(seobj, assay = "Spatial")

    seobj <- RunBanksy(seobj, lambda = 0.2, verbose=TRUE,assay = 'SCT', slot = 'data', features = 'variable',k_geom = 15)

    seobj <- RunPCA(seobj, assay = 'BANKSY', reduction.name = "pca.banksy", features = rownames(seobj), npcs = 30)
    seobj <- FindNeighbors(seobj, dims = 1:30,reduction = "pca.banksy")
    seobj <- FindClusters(seobj, cluster.name = "banksy_cluster",resolution = 0.5)
    seobj <- RunUMAP(seobj, dims = 1:30, reduction = "pca.banksy")
    Idents(seobj) <- "banksy_cluster"
    p1 <- DimPlot(seobj,pt.size = 0.25, label = TRUE, label.size = 3, repel = TRUE)
    p2 <- SpatialDimPlot(seobj, group.by = "banksy_cluster", stroke = NA, label = TRUE, label.size = 3,repel = TRUE, alpha = 0.5, pt.size.factor = 2)
    dir.create(file.path(output_dir,sample_name, fsep = "/"), showWarnings = TRUE)
    png_file = paste(output_dir,"/",sample_name,"/",sample_name,".Spatial_dimPlot",".png", sep="")
    png(png_file, height=780, width=980)
    grid.arrange(p1,p2,ncol=2)
    dev.off()

    # Find markers
    DefaultAssay(seobj) <- 'Spatial'
    markers <- FindAllMarkers(seobj, only.pos = TRUE,assay = 'SCT', min.pct = 0.25, logfc.threshold = 0.25)
    #markers <- markers[markers$p_val_adj < 0.05]
    rds_file = paste(output_dir,"/",sample_name,"/",sample_name,".rds", sep="")
    saveRDS(seobj, file=rds_file)
    write.table(markers, file = paste(output_dir,"/",sample_name,"/",sample_name,".diff.txt", sep=""), col.names = NA)

})

