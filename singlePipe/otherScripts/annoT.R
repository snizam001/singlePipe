#!/usr/bin/env Rscript
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
if(!require(SummarizedExperiment)){
       BiocManager::install("SummarizedExperiment")}
library(SummarizedExperiment)

if(!require(optparse)){
        install.packages("optparse")}
library(optparse)

if(!require(Matrix)){
        install.packages("Matrix")}
library(Matrix)

if(!require(Seurat)){
        install.packages("Seurat")}
library(Seurat)

if(!require(data.table)){
       install.packages("data.table")}
library(data.table)

if(!require(tidyverse)){
       install.packages("tidyverse")}
library(tidyverse)

if(!require(HGNChelper)){
       install.packages("HGNChelper")}
library(HGNChelper)

#--- 
option_list = list(

make_option(c("-i", "--inputFile"), 
    type="character", 
    default=NA,
    help="Input file", 
    metavar="character"),

make_option(c("-o", "--outputDir"), 
    type="character", 
    default=NA,
    help="output directory",
    metavar="character"
           ),
           
make_option(c("-g", "--gene_sets_prepare"), 
    type="character", 
    default=NA,
    help="R script provided by scType pakcage for preparation of gene set",
    metavar="character"),
        
make_option(c("-s", "--sctype_score"), 
    type="character", 
    default=NA,
    help="R script provided by scType pakcage for scoring cells",
    metavar="character"
           ),

make_option(c("-d", "--ScTypeDB_full"), 
    type="character", 
    default=NA,
    help="data set of scType pakcage for scoring cells",
    metavar="character"
           ),

make_option(c("-c", "--cellAnnotateTissue"), 
    type="character", 
    default=NA,
    help="cells of this tissue should be annotated",
    metavar="character"
           ),
           
make_option(c("-n", "--clusterPlotN1"), 
    type="integer", 
    default=1,
    help="Dimension number of UMAP/tsne/PCA to generate graph. Default is 1"),

make_option(c("-m", "--clusterPlotN2"), 
    type="integer", 
    default=2,
    help="Dimension number of UMAP/tsne/PCA to generate graph. Default is 2"),
    
make_option(c("-t", "--clustMethod"), 
    type="character", 
    default="Seurat",
    help="Should this pipeline use default method (Seurat) or method of Granja et al (2019). Default is Seurat.")                             
)

parser<- OptionParser(option_list=option_list)
opt <- parse_args(parser)
inputFile=opt$inputFile
outputDir=opt$outputDir
gene_sets_prepare=opt$gene_sets_prepare
sctype_score=opt$sctype_score
ScTypeDB_full=opt$ScTypeDB_full
cellAnnotateTissue=opt$cellAnnotateTissue
clusterPlotN1=opt$clusterPlotN1
clusterPlotN2=opt$clusterPlotN2
clustMethod=opt$clustMethod

#set.seed(1)
sobj = readRDS(inputFile)


#--- Find cell types using sc-Type package: https://github.com/IanevskiAleksandr/sc-type 

source(gene_sets_prepare)
source(sctype_score)

gs_list = gene_sets_prepare(ScTypeDB_full, cellAnnotateTissue) ############# download it into data folder of package and also set if anyone want to provide his own list of markers

# get cell-type by cell matrix

es.max = sctype_score(scRNAseqData = sobj[["RNA"]]@scale.data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# merge by cluster

cL_resutls = do.call("rbind", lapply(unique(sobj@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(sobj@meta.data[sobj@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(sobj@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"

sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

# output the scores and number of cells

fwrite(sctype_scores,
	paste(outputDir, '/', 'sctype_scores','.txt', sep = ""),
	sep = '\t',
	quote=F
	)

# generate figure/s
 
sobj@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  sobj@meta.data$customclassif[sobj@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

if (clustMethod == "Seurat") {
	jpeg(paste(outputDir, '/', clustMethod, 'DimPlotOverlayingAnnotation-UMAP.jpeg', sep = ""),
		unit = "in",
		res = 300,
		height = 6,
		width =  8)
		
	DimPlot(sobj, reduction = "umap", dims = c(clusterPlotN1,clusterPlotN2), label = TRUE, repel = TRUE, group.by = 'customclassif') 

	dev.off()

	jpeg(paste(outputDir, '/', clustMethod, 'DimPlotOverlayingAnnotation-tSNE.jpeg', sep = ""),
		unit = "in",
		res = 300,
		height = 6,
		width =  8)
		
	DimPlot(sobj, reduction = "tsne", dims = c(clusterPlotN1,clusterPlotN2), label = TRUE, repel = TRUE, group.by = 'customclassif') 

	dev.off()

	jpeg(paste(outputDir, '/', clustMethod, 'DimPlotOverlayingAnnotation-PCA.jpeg', sep = ""),
		unit = "in",
		res = 300,
		height = 6,
		width =  8)
		
	DimPlot(sobj, reduction = "pca", dims = c(clusterPlotN1,clusterPlotN2), label = TRUE, repel = TRUE, group.by = 'customclassif')  
	dev.off()
} else if (clustMethod == "GranjaEtal")  {

	print("Dimplot for the Granja et al methods after annotation is ongoing. Please see table MetadataInformation.txt")

}


# scType vs. Clusters output in table 

fwrite(table(sobj@meta.data$seurat_clusters,sobj@meta.data$customclassif),
	paste(outputDir, '/', 'sctypeVsCluster','.txt', sep = ""),
	sep = '\t',
	quote=F
	)

fwrite(data.frame(sobj@meta.data),
	paste(outputDir, '/', 'MetadataInformation','.txt', sep = ""),
	sep = '\t',
	quote=F
	)
	
saveRDS(sobj, file = paste(outputDir, '/', 'sobj.rds', sep = ""))

