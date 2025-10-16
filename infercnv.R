usage1 = "

Usage:
Rscript infercnv.R (followed by options below)
    1. <SAMPLE NAME>
    2. <FULL PATH OF OUTDIR>
   "
args = commandArgs(trailingOnly = TRUE)
if(length(args)!=2) {
  cat(usage1)
  stop("\n Wrong parameters. See usage above.\n")
}

options(echo=TRUE)
system("ulimit -n 10000")
packages <- c("Seurat","SeuratData", "dplyr","SeuratWrappers","ggplot2","gridExtra","infercnv")
invisible(lapply(packages, library, character.only=TRUE))
options(future.globals.maxSize = 4 * 1024^3)

sample_name = args[1]
output_dir = args[2]
dir.create(file.path(output_dir), showWarnings = TRUE)

all_samples_obj <- readRDS("/scratch/Gliome_spatial/github_scripts/multiomics_spatialTranscriptomics/banksy_harmony_integration/horizontal_integration_allSamples_banksy_harmony.rds")
DefaultAssay(all_samples_obj) <- "SCT"


#"MBP4"     5.83438892198847e-17   1.10429019890168   0.641  0.56   1.43321763868647e-12   "11"  "MBP"
# "CLU3"     4.82055413261313e-85   0.502312550649358  0.998  0.985  1.18416912267642e-80   "14"  "CLU"
# "GABRA1"   0                      2.77922642194792   0.605  0.096  0                      "14"  "GABRA1"
# "GAP431"   4.66015585380818e-44   0.42419622163564   0.901  0.679  1.14476728548798e-39   "14"  "GAP43"
# "HOPX4"    1.50024135668177e-51   0.344367646926433  0.774  0.483  3.68534289268876e-47   "14"  "HOPX"
# "MBP5"     3.42621832255961e-147  0.258870280857774  0.933  0.55   8.41650530936769e-143  "14"  "MBP"

# "GABRA11"  0                      3.38621234451194   0.786  0.104  0                      "18"  "GABRA1"
# "HOPX5"    5.94989571435464e-24   0.64673162129909   0.742  0.49   1.46159188223122e-19   "18"  "HOPX"
# "MBP6"     3.33101097247902e-92   1.03311962050385   0.954  0.558  8.1826284538947e-88    "18"  "MBP"
# "ZIC13"    5.26781856048144e-20   1.66969139751524   0.458  0.3    1.29403962938227e-15   "18"  "ZIC1"


# "CLU"      5.86546325799879e-116  0.339404178144588  1      0.983  1.4408510493274e-111   "1"   "CLU"
# "MBP"      4.32258756738069e-80   1.00999772961915   0.693  0.548  1.06184363592707e-75   "1"   "MBP"
# "ZIC1"     9.18494640495332e-65   0.972652882455789  0.418  0.289  2.25628208437678e-60   "1"   "ZIC1"


# "MBP1"     2.18448453200359e-47   0.407425305121649  0.688  0.552  5.36618625286682e-43   "3"   "MBP"
# "ZIC11"    5.23285299280329e-60   1.23638327394478   0.424  0.292  1.28545033768213e-55   "3"   "ZIC1"


Idents(all_samples_obj) <- "seurat_clusters"
sample_cells <- WhichCells(all_samples_obj, expression = sample == sample_name)
ref_cells <- WhichCells(all_samples_obj, expression = "ZIC1" > 0.5 & "GABRA1" > 0.5 & "MBP" > 0.5 & sample %in% c(sample_name), idents = c(1,3,18), slot = "data")

#annotate cells as tumor and non-tumor
all_samples_obj$annotation <- "tumor"
all_samples_obj$annotation[ref_cells] <- "GCL_N"    # non tumor cells


#prepare for infercnv
sample_cells <- WhichCells(all_samples_obj, expression = sample == sample_name)
exprMatrix <- all_samples_obj@assays$SCT@counts[, sample_cells] %>% as.matrix() %>% as.data.frame()
cellAnnota <- all_samples_obj@meta.data[sample_cells,] %>% dplyr::select(annotation)
gene_order_file <- "/scratch/Gliome_spatial/github_scripts/multiomics_spatialTranscriptomics/inferCNV/hg38_gencode_v27.txt"

# create infer CNV object
infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix=exprMatrix,
                                               gene_order_file=gene_order_file,
                                               annotations_file=cellAnnota,
                                               ref_group_names=c("GCL_N"))




infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, 
                             out_dir=output_dir, 
                             cluster_by_groups=TRUE, 
                             final_scale_limits=c(0.7,1.3),
                             outlier_lower_bound=0.7,
                             outlier_upper_bound=1.3,
                             analysis_mode='subclusters',
                             tumor_subcluster_partition_method = c("leiden"),
                             output_format = "pdf",
                             denoise=TRUE,
                             HMM=T, 
                             num_threads=8,
                             smooth_method = "runmeans")

#while neuronal markers ZIC1, GABRA1, and mature oligodendrocyte marker MBP are highest in the GCL_N (late peak). 
#GCL_TI specifically upregulates RG markers TNC, HOPX, PTPRZ1, and VIM, astrocyte markers CLU, LGALS1, as well as genes involved in migration and glioma network formation such as CD44, SPARC, and GAP43 (One peak) (Fig. 5d and Supplementary Data 9). 
#These genes are collectively termed GCL_TI signature.





# Fetch all the cells corresponding to expression of MBP (cluster 12,9), ZIC1(cluster 1,4), and GABRA1(cluster 21) from sample DMG1
# Assign these cells as reference
# Assign all other cells as tumor
#reference <- FetchData(all_samples_obj, vars = c("sample","ZIC1","GABRA1","MBP","seurat_clusters"))
#ref_tmp <- subset(x = all_samples_obj, subset = ZIC1 > 0.5 & GABRA1 > 0.5 & MBP > 0.5, idents = c(1,4,9,12,21))





#The CNV estimations mentioned in the paragraphs you sent are based on the residual expression. This you can find those either directly in the infercnv object in the infercnv_obj@expr.data slot, or can read them from the infercnv.observations.txt and infercnv.references.txt tables
#https://github.com/broadinstitute/infercnv/issues/609
#https://github.com/broadinstitute/infercnv/issues/208
#https://github.com/broadinstitute/infercnv/issues/338
#https://github.com/broadinstitute/infercnv/issues/434

seurat_obj = infercnv::add_to_seurat(infercnv_output_path=path_to_your_infercnv_run_folder,
                                     seurat_obj=your_seurat_obj, # optional
                                     top_n=10
                                     )