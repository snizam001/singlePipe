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
       BiocManager::install("data.table")}
library(data.table)

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

make_option(c("-s", "--selMethod"), 
    type="character", 
    default='vst',
    help="selection method for variables. Options are vst, dispersion and mean.var.plot. For detail see: https://satijalab.org/seurat/reference/findvariablefeatures. Default is vst.",
    metavar="character"
           ),
                      
make_option(c("-v", "--VariableNumberFeatures"), 
    type="integer", 
    default=2000,
    help="Number of the features which are variable and used for dimension reduction. Default is 2000."
           )
)

parser<- OptionParser(option_list=option_list)
opt <- parse_args(parser)
inputFile=opt$inputFile
outputDir=opt$outputDir
VariableNumberFeatures=opt$VariableNumberFeatures
selMethod=opt$selMethod

set.seed(1)
sobj = readRDS(inputFile)
sobj <- FindVariableFeatures(sobj, selection.method = selMethod, nfeatures = VariableNumberFeatures)
top10 <- head(VariableFeatures(sobj), 10)
fwrite(data.frame(gene = VariableFeatures(sobj)),
	paste(outputDir, '/', 'VariableFeatures-',VariableNumberFeatures,'.txt', sep = ""),
	sep = '\t',
	quote=F
	)
plot1 <- VariableFeaturePlot(sobj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

jpeg(paste(outputDir, '/', 'ExpressionPlot.jpeg', sep = ""),
	unit = "in",
	res = 300,
	height = 6,
	width =  6)
	
plot2

dev.off()

#--- linear Scaling and plotting PCA
all.genes = rownames(sobj)
print("Go and drink coffee because it is going to take lots of time :-)")
print('')
print('')
print('')

#sobj <- ScaleData(sobj, vars.to.regress = "percent.mt", features = all.genes)
sobj <- ScaleData(sobj, features = all.genes)

sobj <- RunPCA(sobj, features = VariableFeatures(object = sobj))
sink(paste(outputDir, '/', 'PCAImportantGenes.txt', sep = ""))

print(sobj[['pca']])

sink()

print("Simulating to identify right number of the PCs")

sobj <- JackStraw(sobj, num.replicate = 100, dims = 50, verbose = T )
sobj <- ScoreJackStraw(sobj, dims = 1:50, verbose  = T )

jpeg(paste(outputDir, '/', 'ElbowPlot.jpeg', sep = ""),
	unit = "in",
	res = 300,
	height = 6,
	width =  6)
	
ElbowPlot(sobj, ndims = 50)

dev.off()

print(paste("Open ElbowPlot.jpeg in the your output folder and ",
"find number of eigenvectors required for the ",
"downstream analysis. The right number is on x-axis after which plateau reached.",
sep =""))


saveRDS(sobj, file = paste(outputDir, '/', 'sobj.rds', sep = ""))



