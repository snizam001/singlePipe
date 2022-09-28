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
           
make_option(c("-l", "--findMarkersLogFC"), 
    type="double", 
    default=0.25,
    help="Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is 0.25"),
        
make_option(c("-p", "--findMarkersMinPCT"), 
    type="double", 
    default=0.25,
    help="Only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. Default is 0.25"
           )  ,

make_option(c("-c", "--clustMethod"), 
    type="character", 
    default="Seurat",
    help="Should this pipeline use default method (Seurat) or method of Granja et al (2019). Default is Seurat.")           
)

parser<- OptionParser(option_list=option_list)
opt <- parse_args(parser)
inputFile=opt$inputFile
outputDir=opt$outputDir
findMarkersLogFC=opt$findMarkersLogFC
findMarkersMinPCT=opt$findMarkersMinPCT
clustMethod=opt$clustMethod

set.seed(1)
sobj = readRDS(inputFile)


if (clustMethod=="GranjaEtal") {
pident=as.factor(sobj@meta.data$seurat_clusters)
names(pident)=rownames(sobj@meta.data)    
Idents(sobj)=pident
}


sobj.markers <- FindAllMarkers(sobj, only.pos = T, min.pct = findMarkersMinPCT, logfc.threshold = findMarkersLogFC)

Alltable = sobj.markers %>% group_by(cluster)

fwrite(Alltable,
	paste(outputDir, '/', clustMethod, 'ClusterAssociatedGenesAll','.txt', sep = ""),
	sep = '\t',
	quote=F
	)

top10 = sobj.markers %>% group_by(cluster) %>% slice_max(n = 10, order_by = avg_log2FC)

fwrite(top10,
	paste(outputDir, '/', clustMethod, 'ClusterAssociatedGenestop10','.txt', sep = ""),
	sep = '\t',
	quote=F
	)
	
top1 = sobj.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC) 
jpeg(paste(outputDir, '/', clustMethod, 'Top1GeneAssociatedCluster.jpeg', sep = ""),
	unit = "in",
	res = 300,
	height = 0.5 * length(levels(as.factor(sobj@meta.data$seurat_clusters))),
	width =  8+4+4)
VlnPlot(sobj, features = top1$gene )
dev.off()

jpeg(paste(outputDir, '/', clustMethod, 'Top10GeneAssociatedHeatmap.jpeg', sep = ""),
	unit = "in",
	res = 300,
	height = 6,
	width =  6)	
DoHeatmap(sobj, features = top10$gene) + NoLegend()
dev.off()



