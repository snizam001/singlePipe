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

#--- 
option_list = list(

make_option(c("-i", "--inputFile"), 
    type="character", 
    default=NA,
    help="Input file", 
    metavar="character"),

make_option(c("-t", "--inType"), 
    type="character", 
    default=NA,
    help="Input file (csv,tsv,rds,....)", 
    metavar="character"),
        
make_option(c("-o", "--outputDir"), 
    type="character", 
    default=NA,
    help="output directory",
    metavar="character"
           ),

make_option(c("-n", "--nFeatureMin"), 
    type="character", 
    default="200",
    help="Cells less than this feature will be filtered. Default value is 200",
    metavar="character"
           ),

make_option(c("-x", "--nFeatureMax"), 
    type="character", 
    default="2500",
    help="Cells mote than this feature will be filtered. Default value is 2500",
    metavar="character"
           ),

make_option(c("-p", "--percMT"), 
    type="character", 
    default="5",
    help="Cells more than this percentage of mitochondria will be filtered. Default value is 5 percent",
    metavar="character"
           )
);

parser<- OptionParser(option_list=option_list)
opt <- parse_args(parser)
inType=opt$inType
inputFile=opt$inputFile
outputDir=opt$outputDir
nFeatureMin=opt$nFeatureMin
nFeatureMax=opt$nFeatureMax
percMT=opt$percMT

set.seed(1)
#--- reading and creating seurat object

if(inType=='rds') {
	print('File type is RDS.')
	data = readRDS(inputFile)
	if ( class(data) == "RangedSummarizedExperiment"){
		counts=assays(data)$counts
		CreateSeuratObject(
		counts = counts, 
		project = "my"
			) -> sobj
	}  else if (class(data) == "Seurat") {
		sobj = data
	}
	
} else if (inType=='csv') {
	print('File type is CSV')
	data = fread(inputFile, header  = T , sep = ",")
	counts=data
	CreateSeuratObject(
		counts = counts, 
		project = "my"
	) -> sobj

} else if (inType=='tsv') {
	print('File type is TSV')
	data = fread(inputFile, header  = T , sep = "\t")
	counts=data
	CreateSeuratObject(
		counts = counts, 
		project = "my"
	) -> sobj
} 

print ('QC: read file')
#--- percentage of mitochondrial genes

sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")

#---
jpeg(paste(outputDir, '/', 'QCMatrix-beforeFiltering.jpeg', sep = ""),
	unit = "in",
	res = 300,
	height = 3.5,
	width = 6 )
 
VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

dev.off()

#--- filtering of cells 

sobj <- subset(sobj, subset = nFeature_RNA > nFeatureMin & nFeature_RNA < nFeatureMax & percent.mt < percMT)
#---
jpeg(paste(outputDir, '/', 'QCMatrix-afterFiltering.jpeg', sep = ""),
	unit = "in",
	res = 300,
	height = 3.5,
	width = 6 )
 
VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

dev.off()

saveRDS(sobj, file = paste(outputDir, '/', 'sobj.rds', sep = ""))
print ('QC: Finished')

